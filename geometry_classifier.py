import os
import traceback
import numpy as np
import pandas as pd
import argparse
from Bio.PDB import PDBParser, NeighborSearch
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
import joblib
from xgboost import XGBClassifier
import re

LOGFILE = "ligand_classification_results.txt"

def log(msg):
    print(msg)
    with open(LOGFILE, "a") as f:
        f.write(msg + "\n")

# --- Ligand geometry ---
def read_ligand_geometry(pdb_file, pocket_centroid, sia_atom_name="C1"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ligand", pdb_file)
    
    sia_coord = None
    all_coords = []

    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res:
                    all_coords.append(atom.get_coord())
                    if atom.get_name() == sia_atom_name:
                        sia_coord = atom.get_coord()
    
    if sia_coord is None:
        raise ValueError(f"SIA atom '{sia_atom_name}' not found in {pdb_file}")
    
    all_coords = np.array(all_coords)
    dists = np.linalg.norm(all_coords - sia_coord, axis=1)
    terminal_coord = all_coords[np.argmax(dists)]

    vec_pocket_to_sia = sia_coord - np.array(pocket_centroid)
    vec_sia_to_terminal = terminal_coord - sia_coord
    distance_to_centroid = np.linalg.norm(vec_pocket_to_sia)
    angle_deg = np.degrees(
        np.arccos(
            np.clip(
                np.dot(vec_pocket_to_sia, vec_sia_to_terminal) /
                (np.linalg.norm(vec_pocket_to_sia) * np.linalg.norm(vec_sia_to_terminal)),
                -1.0, 1.0
            )
        )
    )
    vector_magnitude = np.linalg.norm(vec_sia_to_terminal)
    
    return distance_to_centroid, angle_deg, vector_magnitude, all_coords

# --- Extract Vina score ---
def extract_vina_score(pdbqt_file):
    with open(pdbqt_file) as f:
        for line in f:
            if line.strip().startswith("REMARK VINA RESULT:"):
                parts = line.strip().split()
                return float(parts[3])
    return None

def read_per_residue_probs(filename):
    """
    Reads per-residue probabilities from a file formatted like:
      1. RefPos  174 | Prob 0.923931
    Returns a dictionary: {resnum: probability}
    """
    probs = {}
    with open(filename) as f:
        for line in f:
            # Skip lines that don't contain "RefPos" and "Prob"
            if "RefPos" not in line or "Prob" not in line:
                continue
            # Extract numbers using regex
            m = re.search(r"RefPos\s+(\d+).*Prob\s+([0-9.]+)", line)
            if m:
                resnum = int(m.group(1))
                prob = float(m.group(2))
                probs[resnum] = prob
    return probs


# --- Compute interface features ---
def compute_interface_features(ligand_coords, receptor_pdb, per_res_probs, distance_cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", receptor_pdb)
    atoms = [atom for atom in structure.get_atoms()]
    ns = NeighborSearch(atoms)

    nearby_probs = []
    for atom_coord in ligand_coords:
        neighbors = ns.search(atom_coord, distance_cutoff, level="R")
        for res in neighbors:
            resnum = res.get_id()[1]
            if resnum in per_res_probs:
                nearby_probs.append(per_res_probs[resnum])
    
    if not nearby_probs:
        return 0.0, 0.0, 0.0
    return np.sum(nearby_probs), np.mean(nearby_probs), np.max(nearby_probs)

# --- Build combined dataset ---
def build_dataset(ligand_dir, pocket_centroid, receptor_pdb, per_res_probs, sia_atom_name="C1", label_map=None):
    rows = []
    ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith("_smina_out.pdbqt")]

    for f in ligand_files:
        pdb_path = os.path.join(ligand_dir, f)
        base_name = os.path.splitext(f)[0]
        output_prefix = f"{base_name}_dokken_out"
        
        try:
            dist, angle, mag, ligand_coords = read_ligand_geometry(pdb_path, pocket_centroid, sia_atom_name)
            score = extract_vina_score(pdb_path)
            interface_sum, interface_mean, interface_max = compute_interface_features(
                ligand_coords, receptor_pdb, per_res_probs
            )
            label = label_map.get(f, 1) if label_map else 1

            rows.append({
                "ligand": f,
                "distance_to_centroid": dist,
                "glycan_angle": angle,
                "terminal_vector_mag": mag,
                "vina_score": score,
                "interface_sum": interface_sum,
                "interface_mean": interface_mean,
                "interface_max": interface_max,
                "label": label
            })

            # Save per-ligand CSV
            pd.DataFrame([{
                "distance_to_centroid": dist,
                "glycan_angle": angle,
                "terminal_vector_mag": mag,
                "vina_score": score,
                "interface_sum": interface_sum,
                "interface_mean": interface_mean,
                "interface_max": interface_max,
                "label": label
            }]).to_csv(f"{output_prefix}.csv", index=False)

            log(f"[INFO] Processed {f}, output prefix: {output_prefix}")

        except Exception as e:
            log(f"[ERROR] Failed processing {f}: {e}")
            log(traceback.format_exc())

    return pd.DataFrame(rows)

# --- Main ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ligand_dir", required=False, default=os.getcwd(),help="Directory containing _smina_out.pdbqt files")
    parser.add_argument("--pocket_centroid", nargs=3, type=float, required=True, help="x y z coordinates")
    parser.add_argument("--receptor_pdb", required=False, help="Optional: prepped receptor PDB")
    parser.add_argument("--per_residue_probs", required=False, default="per_residue_probs.txt", help="Per-residue probabilities file")
    parser.add_argument("--sia_atom_name", default="C1", help="Atom name for glycan SIA")
    args = parser.parse_args()

    # Clear previous logfile
    if os.path.exists(LOGFILE):
        os.remove(LOGFILE)

    try:
        per_res_probs = read_per_residue_probs(args.per_residue_probs)

        # Default receptor PDB if not specified
        receptor_pdb = args.receptor_pdb if args.receptor_pdb else "receptor.pdb"

        # Build dataset from all ligands
        df = build_dataset(args.ligand_dir, args.pocket_centroid, receptor_pdb, per_res_probs, args.sia_atom_name)
        log(f"[INFO] Dataset built: {len(df)} ligands")

        # Train classifiers
        feature_cols = ["distance_to_centroid","glycan_angle","terminal_vector_mag","vina_score",
                        "interface_sum","interface_mean","interface_max"]
        X = df[feature_cols].values
        y = df["label"].values

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        classifiers = {
            "RandomForest": RandomForestClassifier(n_estimators=100, random_state=42),
            "XGBoost": XGBClassifier(n_estimators=100, use_label_encoder=False, eval_metric='logloss', random_state=42),
            "LogisticRegression": LogisticRegression(max_iter=500)
        }

        for name, clf in classifiers.items():
            try:
                clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)
                report = classification_report(y_test, y_pred)
                log(f"\n=== {name} ===\n{report}")
                joblib.dump(clf, f"ligand_classifier_{name}.pkl")
                log(f"[INFO] Saved model ligand_classifier_{name}.pkl")
            except Exception as e:
                log(f"[ERROR] Training {name} failed: {e}")
                log(traceback.format_exc())

    except Exception as e:
        log(f"[CRITICAL] Script failed: {e}")
        log(traceback.format_exc())
