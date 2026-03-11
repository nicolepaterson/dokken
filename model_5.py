import os
import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from Bio.PDB import PDBParser, NeighborSearch
from Bio.Align import PairwiseAligner
from scipy.spatial import cKDTree
import numpy as np

# ------------------ CONFIG ------------------

PDB_DIR = "./pdbs"
DISTANCE_CUTOFF = 5.0
CLUSTER_DISTANCE = 40.0
BATCH_SIZE = 4
EPOCHS = 200
DEVICE = torch.device("cpu")

HEAD_START = 120
HEAD_END = 250

SIALIC_LIGANDS = {"SIA", "LSTa", "LSTb"}

REFERENCE_SEQ = (
    "MERIVIALAIISIVKGDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKEHNGKLCSIKGVRPLILKD"
    "CSVAGWLLGNPMCDEFLNVPEWSYIVEKDNPVNGLCYPGDFSDYEELKHLMSSTNHFEKIQIIPRSSWSN"
    "HDASSGVSSACPYNGRSSFFRNVVWLIKKNNAYPTIKRTYNNTNVEDLLIIWGIHHPNDAAEQTKLYQNS"
    "NTYVSVGTSTLNQRSIPEIATRPKVNGQSGRMEFFWTILRPNDAISFESNGNFIAPEYAYKIVKKGDSAI"
    "MRSELEYGNCDTKCQTPVGAINSSMPFHNVHPLTIGECPKYIKSDKLVLATGLRNVPQRETRGLFGAIAG"
    "FIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGITNKVNSIIDKMNTQFEAVGKEFNNLERRIE"
    "NLNRKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNEC"
    "MESVRNGTYDYPQYSEESRLNREEIDGVKLESMGTYQILSIYSTVASSLALAIMIAGLSFWMCSNGSLQC"
    "RICI"
)

# ------------------ AMINO ACID MAPPING ------------------

three_to_one = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
    'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TRP':'W','TYR':'Y'
}

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def one_hot_amino_acid(resname):
    vec = [0]*20
    if resname in AA_LIST:
        vec[AA_LIST.index(resname)] = 1
    return vec

# ------------------ ALIGNMENT ------------------

def pdb_to_seq_mapping(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)

    seq = []
    for chain in structure.get_chains():
        for res in chain:
            if res.get_id()[0] == " ":
                aa = three_to_one.get(res.get_resname())
                if aa:
                    seq.append((res.get_id(), aa))
    return seq


def map_to_reference(pdb_seq, reference_seq):
    pdb_string = "".join([aa for _, aa in pdb_seq])

    aligner = PairwiseAligner()
    aligner.mode = "global"  # important for HA numbering
    alignment = aligner.align(reference_seq, pdb_string)[0]

    ref_index = 0
    pdb_index = 0
    mapping = {}

    for r, p in zip(alignment.target, alignment.query):
        if p != "-":
            if r != "-":
                mapping[pdb_seq[pdb_index][0]] = ref_index
            pdb_index += 1
        if r != "-":
            ref_index += 1

    return mapping

# ------------------ INTERFACE DETECTION ------------------

def identify_interface_residues(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)

    interface_residues = set()

    for model in structure:
        ligand_atoms = [
            atom for chain in model
            for res in chain
            if res.get_id()[0] != " " and res.get_resname() in SIALIC_LIGANDS
            for atom in res
        ]

        receptor_atoms = [
            atom for chain in model
            for res in chain
            if res.get_id()[0] == " "
            for atom in res
        ]

        ns = NeighborSearch(receptor_atoms)

        for atom in ligand_atoms:
            neighbors = ns.search(atom.coord, DISTANCE_CUTOFF)
            for n in neighbors:
                interface_residues.add(n.get_parent().get_id())

    return interface_residues

# ------------------ GRAPH BUILDING ------------------

def build_graph(pdb_file, reference_seq):
    pdb_seq = pdb_to_seq_mapping(pdb_file)
    if len(pdb_seq) == 0:
        return None

    mapping = map_to_reference(pdb_seq, reference_seq)
    interface_residues = identify_interface_residues(pdb_file)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)

    valid_atoms = []
    x, y, ref_positions = [], [], []

    for atom in structure.get_atoms():

        if atom.get_id() != "CA":
            continue

        res = atom.get_parent()
        res_id = res.get_id()
        aa = three_to_one.get(res.get_resname())

        if aa is None or res_id not in mapping:
            continue

        ref_pos = mapping[res_id]

        aa_vec = one_hot_amino_acid(aa)
        ref_norm = ref_pos / len(reference_seq)
        aa_vec.append(ref_norm)

        x.append(aa_vec)
        y.append(1 if res_id in interface_residues else 0)
        ref_positions.append(ref_pos)
        valid_atoms.append(atom)

    if len(valid_atoms) == 0:
        return None

    x = torch.tensor(x, dtype=torch.float)
    y = torch.tensor(y, dtype=torch.long)
    ref_positions = torch.tensor(ref_positions, dtype=torch.long)

    # KDTree edge construction
    coords = np.array([atom.coord for atom in valid_atoms])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=DISTANCE_CUTOFF)

    edge_index = []
    for i, j in pairs:
        edge_index.append([i, j])
        edge_index.append([j, i])

    if edge_index:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)

    return Data(
        x=x,
        edge_index=edge_index,
        y=y,
        ref_positions=ref_positions,
        coords=torch.tensor(coords, dtype=torch.float)
    )

# ------------------ MODEL ------------------

class InterfaceGCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels=64):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, 2)

    def forward(self, data):
        x = F.relu(self.conv1(data.x, data.edge_index))
        x = self.conv2(x, data.edge_index)
        return x

# ------------------ SPATIAL CLUSTER FILTER ------------------

def spatial_cluster_filter(coords, probs, threshold=0.5, max_dist=40.0):
    coords = coords.cpu().numpy()
    pred_indices = np.where(probs > threshold)[0]

    if len(pred_indices) <= 1:
        return np.zeros_like(probs)

    keep = []
    for i in pred_indices:
        dists = np.linalg.norm(coords[pred_indices] - coords[i], axis=1)
        if np.sum((dists < max_dist) & (dists > 0)) > 0:
            keep.append(i)

    filtered = np.zeros_like(probs)
    filtered[keep] = probs[keep]
    return filtered

# ------------------ DATASET ------------------

def build_dataset():
    dataset = []
    pdb_files = [os.path.join(PDB_DIR, f) for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        data = build_graph(pdb_file, REFERENCE_SEQ)
        if data:
            dataset.append(data)

    return dataset

# ------------------ TRAIN + INFER ------------------
if __name__ == "__main__":

    dataset = build_dataset()
    split = int(0.8 * len(dataset))
    train_dataset = dataset[:split]
    test_dataset = dataset[split:]

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)

    model = InterfaceGCN(in_channels=21).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    # compute class weights
    total_pos, total_neg = 0, 0
    for data in train_dataset:
        mask = (data.ref_positions >= HEAD_START) & (data.ref_positions <= HEAD_END)
        labels = data.y[mask]
        total_pos += (labels == 1).sum().item()
        total_neg += (labels == 0).sum().item()

    pos_weight = min(total_neg / (total_pos + 1e-6), 30.0)
    class_weights = torch.tensor([1.0, pos_weight], dtype=torch.float).to(DEVICE)

    # TRAIN
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0

        for data in train_loader:
            data = data.to(DEVICE)
            optimizer.zero_grad()

            out = model(data)
            mask = (data.ref_positions >= HEAD_START) & (data.ref_positions <= HEAD_END)

            if mask.sum() == 0:
                continue

            ce_loss = F.cross_entropy(out[mask], data.y[mask], weight=class_weights)

            probs = F.softmax(out, dim=1)[:, 1]
            smooth_loss = ((probs[data.edge_index[0]] - probs[data.edge_index[1]])**2).mean()

            loss = ce_loss + 0.1 * smooth_loss
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        print(f"Epoch {epoch+1:03d} | Loss: {total_loss:.4f}")

    # MODEL
    torch.save(model.state_dict(), "interface_gcn_mapped.pth")
    print("Model saved to interface_gcn_mapped.pth")

    print("Training complete.")
