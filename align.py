#!/usr/bin/env python
"""
Influenza HA PDB Pipeline

1. Search RCSB for influenza PDBs
2. Download PDBs in parallel
3. Filter for sialic acid / LST ligands
4. Fetch sequences and compute alignment to reference
5. Select top matching PDB(s)
6. Copy final PDBs to ./pdbs
"""

import os
import shutil
import requests
from concurrent.futures import ThreadPoolExecutor
from Bio.PDB import PDBParser
from Bio import Align

# ----------------------------- CONFIG --------------------------------
SIALIC_LIGANDS = {"SIA", "LSTa", "LSTb"}
TMP_PDB_DIR = "./tmp_pdbs"
FINAL_PDB_DIR = "./pdbs"
MAX_WORKERS = 8

os.makedirs(TMP_PDB_DIR, exist_ok=True)
os.makedirs(FINAL_PDB_DIR, exist_ok=True)

QUERY_SEQUENCE = """MERIVIALAIISIVKGDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKEHNGKLCSIKGVRPLILKD
CSVAGWLLGNPMCDEFLNVPEWSYIVEKDNPVNGLCYPGDFSDYEELKHLMSSTNHFEKIQIIPRSSWSN
HDASSGVSSACPYNGRSSFFRNVVWLIKKNNAYPTIKRTYNNTNVEDLLIIWGIHHPNDAAEQTKLYQNS
NTYVSVGTSTLNQRSIPEIATRPKVNGQSGRMEFFWTILRPNDAISFESNGNFIAPEYAYKIVKKGDSAI
MRSELEYGNCDTKCQTPVGAINSSMPFHNVHPLTIGECPKYIKSDKLVLATGLRNVPQRETRGLFGAIAG
FIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGITNKVNSIIDKMNTQFEAVGKEFNNLERRIE
NLNRKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNEC
MESVRNGTYDYPQYSEESRLNREEIDGVKLESMGTYQILSIYSTVASSLALAIMIAGLSFWMCSNGSLQC
RICI""".replace("\n", "")

# ------------------------- FUNCTIONS --------------------------------

def search_influenza_pdbs(limit=5000):
    """Search RCSB for influenza PDB IDs."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                "operator": "contains_phrase",
                "value": "Influenza"
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    resp = requests.post(url, json=query)
    resp.raise_for_status()
    data = resp.json()
    pdb_ids = [r["identifier"] for r in data["result_set"]]
    return pdb_ids[:limit]


def download_pdb(pdb_id, folder=TMP_PDB_DIR):
    """Download a PDB file and return its path."""
    pdb_file = os.path.join(folder, f"{pdb_id}.pdb")
    if os.path.exists(pdb_file):
        return pdb_file

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url)
    if r.status_code == 200:
        with open(pdb_file, "w") as f:
            f.write(r.text)
        return pdb_file
    return None


def has_sialic_ligand(pdb_file, ligands=SIALIC_LIGANDS):
    """Check if a PDB structure contains sialic/LST ligands."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)

    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_id()[0] != " " and res.get_resname().upper() in ligands:
                    return True
    return False


def fetch_pdb_sequence(pdb_file):
    """Fetch PDB sequence from RCSB and return (pdb_id, seq)."""
    pdb_id = os.path.basename(pdb_file).replace(".pdb", "")
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    resp = requests.get(url)
    if resp.status_code != 200:
        return None
    fasta = resp.text.split("\n")
    seq = "".join(fasta[1:])
    return pdb_id, seq


def fetch_sequences_parallel(pdb_files, max_workers=MAX_WORKERS):
    """Fetch sequences in parallel and return dict {pdb_id: seq}."""
    sequences = {}

    def fetch_seq(pdb_file):
        result = fetch_pdb_sequence(pdb_file)
        if result is None:
            return None
        return result

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for result in executor.map(fetch_seq, pdb_files):
            if result is not None:
                pdb_id, seq = result
                sequences[pdb_id] = seq

    return sequences


def compute_alignment_score(query_seq, target_seq):
    """Compute normalized global alignment score."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    score = aligner.score(query_seq, target_seq)
    return score / max(len(query_seq), len(target_seq))


def filter_top_matches(query_seq, pdb_sequences, cutoff=0.7, top_n=1):
    """Return top PDB IDs matching the query sequence."""
    scores = {
        pdb_id: compute_alignment_score(query_seq, seq)
        for pdb_id, seq in pdb_sequences.items()
    }
    filtered = {pdb_id: s for pdb_id, s in scores.items() if s >= cutoff}
    ranked = sorted(filtered.items(), key=lambda x: x[1], reverse=True)
    return [pdb_id for pdb_id, _ in ranked[:top_n]]


def select_top_pdbs(query_seq, pdb_files, cutoff=0.7, top_n=1):
    """Return PDB file paths of top matches after sequence alignment."""
    id_to_file = {os.path.basename(p).replace(".pdb", ""): p for p in pdb_files}
    pdb_sequences = fetch_sequences_parallel(pdb_files)
    top_ids = filter_top_matches(query_seq, pdb_sequences, cutoff=cutoff, top_n=top_n)
    return [id_to_file[pdb_id] for pdb_id in top_ids]


def main():
    print("Searching for influenza PDBs...")
    influenza_pdbs = search_influenza_pdbs(limit=5000)
    print(f"Found {len(influenza_pdbs)} influenza PDBs.")

    # Download and filter for ligand
    pdb_files_with_ligand = []
    for pdb_id in influenza_pdbs:
        pdb_file = download_pdb(pdb_id)
        if pdb_file and has_sialic_ligand(pdb_file):
            pdb_files_with_ligand.append(pdb_file)
    print(f"{len(pdb_files_with_ligand)} PDBs have SIA/LST bound.")

    # Select top PDB(s) by alignment
    top_pdb_files = select_top_pdbs(
        QUERY_SEQUENCE,
        pdb_files_with_ligand,
        cutoff=0.7,
        top_n=1
    )

    print("Top matching PDBs:")
    for f in top_pdb_files:
        print(os.path.basename(f))

    # Copy to final directory
    for pdb_file in pdb_files_with_ligand:
        shutil.copy(
            pdb_file,
            os.path.join(FINAL_PDB_DIR, os.path.basename(pdb_file))
        )

    print(f"Final PDB(s) saved in {FINAL_PDB_DIR}.")


if __name__ == "__main__":
    main()
