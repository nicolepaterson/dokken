import torch
import torch.nn.functional as F
import numpy as np
from Bio.PDB import PDBParser
from model_5 import InterfaceGCN, build_graph, REFERENCE_SEQ, DEVICE, pdb_to_seq_mapping
import os

pdb_dir = "/Users/nmp/Desktop/DOKKEN/naive_pdbs"
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]

if len(pdb_files) == 0:
    raise ValueError("No PDB files found in the directory.")

# Build a single graph to get feature size
tmp_data = build_graph(os.path.join(pdb_dir, pdb_files[0]), REFERENCE_SEQ)
in_channels = tmp_data.x.shape[1]  # number of features per node
print("in_channels inferred from graph:", in_channels)

# Now we can instantiate the model
model = InterfaceGCN(in_channels)
model.load_state_dict(torch.load("interface_gcn_mapped.pth", map_location=DEVICE))
model.to(DEVICE)
model.eval()

def predict_interface(model, pdb_file, reference_seq, threshold=0.5, top_k=25, device="cpu"):
    """
    Predict interface residues for an HA structure and compute geometric metrics safely.

    Returns:
        data: graph object with residue_ids and coords
        probs: per-residue probability tensor
        preds: thresholded predictions (True/False)
        metrics: dict with centroid, radius, compactness, mean_prob_RBS, top residues
    """
    # Build the graph for the PDB
    data = build_graph(pdb_file, reference_seq)
    data = data.to(device)

    # Store residue IDs
    res_mapping = pdb_to_seq_mapping(pdb_file)
    data.residue_ids = [res_id for res_id, aa in res_mapping] if res_mapping else []

    # Get CA coordinates
    coords = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("HA", pdb_file)
    for atom in structure.get_atoms():
        if atom.get_id() == "CA":
            coords.append(atom.coord)
    data.coords = torch.tensor(coords, dtype=torch.float, device=device) if coords else torch.empty((0,3), dtype=torch.float, device=device)

    # Run the model
    model.eval()
    with torch.no_grad():
        logits = model(data)
        probs = F.softmax(logits, dim=1)[:, 1].flatten() if logits.numel() > 0 else torch.tensor([], device=device)
        preds = probs > threshold if probs.numel() > 0 else torch.tensor([], dtype=torch.bool, device=device)

    # Adjust top_k if there are fewer residues
    top_k_residues = min(top_k, len(probs))
    top_indices = torch.topk(probs, top_k_residues).indices if top_k_residues > 0 else torch.tensor([], dtype=torch.long, device=device)
    rbs_coords = data.coords[top_indices] if top_indices.numel() > 0 else torch.empty((0,3), dtype=torch.float, device=device)

    # Compute metrics safely
    if rbs_coords.numel() == 0 or not torch.isfinite(rbs_coords).all():
        centroid = np.array([np.nan, np.nan, np.nan], dtype=float)
        radius = np.nan
        compactness = np.nan
        mean_prob_RBS = np.nan
    else:
        centroid = rbs_coords.mean(dim=0).cpu().numpy().astype(float)
        diffs = rbs_coords - torch.tensor(centroid, device=device)
        radius = torch.norm(diffs, dim=1).max().item()
        compactness = torch.norm(diffs, dim=1).mean().item()
        mean_prob_RBS = probs[top_indices].mean().item()

    top_residues = [data.residue_ids[i] if i < len(data.residue_ids) else None for i in top_indices.cpu().numpy()]

    metrics = {
        "centroid": centroid,
        "radius": radius,
        "compactness": compactness,
        "mean_prob_RBS": mean_prob_RBS,
        "top_indices": top_indices.cpu().numpy(),
        "top_residues": top_residues
    }

    return data, probs.cpu(), preds.cpu(), metrics


model = InterfaceGCN(in_channels)
model.load_state_dict(torch.load("interface_gcn_mapped.pth", map_location=DEVICE))
model.to(DEVICE)

for pdb_file in pdb_files:
    full_path = os.path.join(pdb_dir, pdb_file)
    data, probs, preds, metrics = predict_interface(
        model,
        full_path,
        REFERENCE_SEQ,
        threshold=0.6,   # slightly stricter than default
        top_k=25,
        device=DEVICE
    )

    print(f"\nPredicted interface residues for {pdb_file} (top 25):")
    print(metrics["top_residues"])
    print("Geometric metrics:")
    print("Centroid:", metrics["centroid"])
    print("Radius:", metrics["radius"])
    print("Compactness:", metrics["compactness"])
    print("Mean RBS probability:", metrics["mean_prob_RBS"])

