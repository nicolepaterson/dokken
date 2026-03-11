#!/usr/bin/env python
"""
python smina.py \
 --pdb_dir pdbs \
 --output all_docks.csv \
--extra_ligands 2_3_smina.pdbqt 2_6_smina.pdbqt \
--pdb_dir ./pdb_pairs --output pdb_pairs_smina_scores.csv  \
--fasta_dir /Users/nmp/Desktop/DOKKEN/pdbs/fasta
"""

import os
import csv
import subprocess
import gc
from pathlib import Path
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

SIA_LIGANDS = {"SIA", "LSTa", "LSTb", "LSTc", "NANA", "NEU5AC", "3SLN", "3SLN-LN", "SIA1", "SA", "KDN"}
SUGAR_NAMES    = {"NAG", "MAN", "BMA", "FUC"}
RECEPTOR_PDBQT_DIR = Path("pdbqt_receptor")
LIGAND_PDBQT_DIR   = Path("pdbqt_ligand")

TIMEOUT_RECEPTOR_PREP = 4000
TIMEOUT_SMINA_SCORE   = 4000
TIMEOUT_SMINA_DOCK    = 8000  

def extract_chain(pdb_file, chain_id="A"):
    out_file = pdb_file.parent / f"{pdb_file.stem}_chainA_tmp.pdb"
    with open(pdb_file) as f, open(out_file, "w") as out:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain_id:
                continue
            resname = line[17:20].strip()
            if resname in SIA_LIGANDS or resname in SUGAR_NAMES:
                continue
            altloc = line[16]                        
            if altloc not in (" ", "A"):           
                continue
            out.write(line[:16] + " " + line[17:]) 
        out.write("END\n")
    return out_file

def strip_sugars(receptor_pdb):
    out_file = receptor_pdb.parent / f"{receptor_pdb.stem}_nosugar_tmp.pdb"
    with open(receptor_pdb) as f, open(out_file, "w") as out:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                if line[17:20].strip() in SUGAR_NAMES:
                    continue
            out.write(line)
    return out_file

def extract_sia_ligand(pdb_file):
    out_file = pdb_file.parent / f"{pdb_file.stem}_ligand_tmp.pdb"
    lines = [l for l in open(pdb_file) if l.startswith("HETATM") and l[17:20].strip() in SIA_LIGANDS]
    if not lines:
        return None
    with open(out_file, "w") as out:
        out.writelines(lines)
        out.write("END\n")
    return out_file

def prepare_smina_receptor(receptor_pdb, output_dir):
    output_dir.mkdir(exist_ok=True)
    out_stem = output_dir / receptor_pdb.stem
    cmd = [
    "mk_prepare_receptor.py",
    "-i", str(receptor_pdb),
    "-o", str(out_stem),
    "-p",
    "--default_altloc", "A",   
    "-a",                 
    ]
    print("[INFO] Preparing receptor PDBQT:", " ".join(cmd))
    subprocess.run(cmd, check=True, timeout=TIMEOUT_RECEPTOR_PREP)
    candidates = sorted(output_dir.glob(f"{receptor_pdb.stem}*.pdbqt"))
    if not candidates:
        raise FileNotFoundError(f"No receptor PDBQT produced in {output_dir} for {receptor_pdb.stem}")
    rigid = next((c for c in candidates if "_rigid" in c.name), candidates[0])
    print(f"[INFO] Using receptor PDBQT: {rigid}")
    return rigid

def prepare_flex_receptor(receptor_pdb, centroid, flex_residues, size=(30, 30, 30)):
    receptor_out_stem = receptor_pdb.parent / f"{receptor_pdb.stem}_flex"
    mk_cmd = [
        "mk_prepare_receptor.py",
        "-i", str(receptor_pdb),
        "-o", str(receptor_out_stem),
        "-p", "-g", "-v",
        "--box_size",   str(size[0]), str(size[1]), str(size[2]),
        "--box_center", str(centroid[0]), str(centroid[1]), str(centroid[2]),
        "-f", ",".join(flex_residues),
        "-a"
    ]
    print("[INFO] Preparing flex receptor:", " ".join(mk_cmd))
    subprocess.run(mk_cmd, check=True, timeout=TIMEOUT_RECEPTOR_PREP)
    rigid = receptor_out_stem.parent / f"{receptor_out_stem.name}_rigid.pdbqt"
    flex  = receptor_out_stem.parent / f"{receptor_out_stem.name}_flex.pdbqt"
    if not rigid.exists():
        candidates = sorted(receptor_out_stem.parent.glob(f"{receptor_out_stem.name}*.pdbqt"))
        if not candidates:
            raise FileNotFoundError(f"No receptor PDBQT found under {receptor_out_stem.parent}")
        rigid = candidates[0]
        flex  = candidates[1] if len(candidates) > 1 else None
    print(f"[INFO] Rigid PDBQT: {rigid}")
    print(f"[INFO] Flex  PDBQT: {flex}")
    return rigid, flex

def prep_meeko_pdbqt(pdb_path, output_dir):
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f"{pdb_path.stem}.pdbqt"
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit failed to read {pdb_path}")
    mol = Chem.AddHs(mol, addCoords=True)
    frags = rdmolops.GetMolFrags(mol, asMols=True)
    mol = max(frags, key=lambda m: m.GetNumAtoms())
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    if not setups:
        raise ValueError(f"Meeko prep failed for {pdb_path}")
    writer = PDBQTWriterLegacy()
    pdbqt_string = writer.write_string(setups[0])[0]
    if not isinstance(pdbqt_string, str):
        pdbqt_string = "\n".join(pdbqt_string)
    with open(output_path, "w") as f:
        f.write(pdbqt_string)
    return output_path

def calc_centroid(pdb_file):
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.mean(coords, axis=0) if coords else None

def assess_pose(coords, pocket_centroid, sia_idx=0):
    coords = np.asarray(coords)
    pocket_centroid = np.asarray(pocket_centroid)
    if coords.size == 0:
        return None, None
    sia      = coords[sia_idx]
    terminal = coords[np.argmax(np.linalg.norm(coords - sia, axis=1))]
    vec_c    = sia - pocket_centroid
    vec_t    = terminal - sia
    dist     = float(np.linalg.norm(vec_c))
    denom    = np.linalg.norm(vec_c) * np.linalg.norm(vec_t)
    angle_deg = 0.0
    if denom != 0:
        cos_theta = np.clip(np.dot(vec_c, vec_t) / denom, -1.0, 1.0)
        angle_deg = float(np.degrees(np.arccos(cos_theta)))
    return dist, angle_deg

def find_flexible_residues(receptor_pdb, centroid, cutoff=6.0):
    centroid = np.asarray(centroid, dtype=float)
    residues = set()
    with open(receptor_pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                chain = line[21]
                resnum = int(line[22:26])
                coord = np.array([x, y, z])
                if np.linalg.norm(coord - centroid) <= cutoff:
                    residues.add(f"{chain}:{resnum}")
            except Exception:
                continue
    residues = sorted(residues)
    print(f"[INFO] Geometric flexible residues (cutoff={cutoff}Å):")
    print(residues)
    return residues

def run_smina_score(receptor_pdbqt, ligand_pdbqt, center, size=(30, 30, 30)):
    v = Vina(sf_name="vina")
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=center.tolist(), box_size=list(size))
    energy = v.score()
    return energy[0]  

def run_flex_docking(receptor_rigid_pdbqt, receptor_flex_pdbqt,
                     centroid, ligand_path, size=(30, 30, 30)):
    dock_out = ligand_path.parent / f"{ligand_path.stem}_docked_to_{receptor_rigid_pdbqt.stem}.pdbqt"
    log_out  = ligand_path.parent / f"{ligand_path.stem}_docked_to_{receptor_rigid_pdbqt.stem}.txt"
    smina_cmd = [
        "smina",
        "--receptor",  str(receptor_rigid_pdbqt),
        "--ligand",    str(ligand_path),
        "--center_x",  str(centroid[0]),
        "--center_y",  str(centroid[1]),
        "--center_z",  str(centroid[2]),
        "--size_x",    str(size[0]),
        "--size_y",    str(size[1]),
        "--size_z",    str(size[2]),
        "--out",       str(dock_out),
        "--exhaustiveness", "12",
        "--num_modes", "5",
        "--log", str(log_out),
        "--seed","0"
    ]

    if receptor_flex_pdbqt and Path(receptor_flex_pdbqt).exists():
        smina_cmd += ["--flex", str(receptor_flex_pdbqt)]
    print("[INFO] Running Smina docking:", " ".join(smina_cmd))
    result = subprocess.run(smina_cmd, capture_output=True, text=True, timeout=TIMEOUT_SMINA_DOCK)
    with open(log_out, "w") as f:
        f.write("=== STDOUT ===\n"); f.write(result.stdout)
        f.write("\n=== STDERR ===\n"); f.write(result.stderr)
    if result.returncode != 0:
        print("[ERROR] Smina docking failed. Log:", log_out)
        return dock_out, None, []
    scores = []
    with open(dock_out) as f:
        for line in f:
            if "REMARK minimizedAffinity" in line:
                try:
                    scores.append(float(line.split()[2]))
                except:
                    pass
            elif "REMARK VINA RESULT:" in line:
                try:
                    scores.append(float(line.split()[3]))
                except:
                    pass
    best_score = min(scores) if scores else None
    print(f"[INFO] Docking done → {dock_out}")
    print(f"[INFO] Best score: {best_score}")
    print(f"[INFO] Pose count: {len(scores)}")
    return dock_out, best_score, scores

def read_sequence(fasta_dir, pdb_name):
    fasta_file = Path(fasta_dir) / f"{pdb_name}.fasta"
    if fasta_file.exists():
        with open(fasta_file) as f:
            return "".join(l.strip() for l in f if not l.startswith(">"))
    return None

def parse_smina_all_scores(pdbqt_file):
    scores = []
    try:
        with open(pdbqt_file) as f:
            for line in f:
                if "REMARK minimizedAffinity" in line:
                    try:
                        scores.append(float(line.split()[3]))
                    except:
                        pass
    except Exception:
        pass
    return scores

def boltzmann_free_energy(scores, temperature=298.15):
    if not scores:
        return None
    R = 0.001987  # kcal/mol/K
    RT = R * temperature
    scores = np.array(scores)
    G = -RT * np.log(np.sum(np.exp(-scores / RT)))
    return float(G)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", required=True)
    parser.add_argument("--extra_ligands", nargs="+", default=[
        "/Users/nmp/Desktop/DOKKEN/2_3_smina.pdbqt",
        "/Users/nmp/Desktop/DOKKEN/2_6_smina.pdbqt",
    ])
    parser.add_argument("--fasta_dir", default="/Users/nmp/Desktop/DOKKEN/pdbs/fasta")
    parser.add_argument("--output",    required=True)
    args = parser.parse_args()
    pdb_dir = "pdb_pairs" 
    #pdb_files = sorted(p for p in pdb_dir.glob("*.pdb") if "_tmp" not in p.name)
    pdb_dir   = Path(args.pdb_dir)
    fasta_dir = Path(args.fasta_dir)
    RECEPTOR_PDBQT_DIR.mkdir(exist_ok=True)
    LIGAND_PDBQT_DIR.mkdir(exist_ok=True)
    results = []
    pdb_files = sorted(p for p in pdb_dir.glob("*.pdb") if "_tmp" not in p.name)
    for pdb_file in pdb_files:
        print(f"\n{'='*60}\n[INFO] Processing {pdb_file.name}")
        receptor_pdb = None
        stripped_receptor_pdb = None
        ligand_pdb = None
        try:
            receptor_pdb = extract_chain(pdb_file)
            stripped_receptor_pdb = strip_sugars(receptor_pdb)
            ligand_pdb = extract_sia_ligand(pdb_file)
            if ligand_pdb is None:
                print("[WARN] No SIA ligand found — skipping")
                continue
            ligand_pdbqt = prep_meeko_pdbqt(ligand_pdb, LIGAND_PDBQT_DIR)
            centroid = calc_centroid(ligand_pdb)
            if centroid is None:
                print("[WARN] Could not compute centroid — skipping")
                continue
            receptor_pdbqt = prepare_smina_receptor(receptor_pdb, RECEPTOR_PDBQT_DIR)
            score_orig = run_smina_score(receptor_pdbqt, ligand_pdbqt, centroid)
            coords = np.array([
                [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in open(ligand_pdb) if l.startswith(("ATOM", "HETATM"))
            ])
            dist, angle_deg = assess_pose(coords, centroid)
            flex_residues = find_flexible_residues(receptor_pdb, centroid, cutoff=10.0)
            if not flex_residues:
                print("[WARN] No flexible residues found within cutoff")
            holo_rigid, holo_flex   = prepare_flex_receptor(receptor_pdb, centroid, flex_residues)
            empty_rigid, empty_flex = prepare_flex_receptor(stripped_receptor_pdb, centroid, flex_residues)
            extra_scores = {}
            for extra in args.extra_ligands:
                extra_pdbqt = Path(extra)
                if not extra_pdbqt.exists():
                    print(f"[WARN] Extra ligand not found: {extra_pdbqt}")
                    extra_scores[extra_pdbqt.stem + "_holo"]  = None
                    extra_scores[extra_pdbqt.stem + "_empty"] = None
                    continue
               #_, holo_score  = run_flex_docking(holo_rigid,  holo_flex,  centroid, extra_pdbqt)
                #_, empty_score = run_flex_docking(empty_rigid, empty_flex, centroid, extra_pdbqt)
                _, holo_score, holo_scores = run_flex_docking(holo_rigid,  holo_flex,  centroid, extra_pdbqt)
                _, empty_score, empty_scores = run_flex_docking(empty_rigid, empty_flex, centroid, extra_pdbqt)
                G_holo = boltzmann_free_energy(holo_scores)
                G_empty = boltzmann_free_energy(empty_scores)
                extra_scores[extra_pdbqt.stem + "_holo_best"]  = holo_score
                extra_scores[extra_pdbqt.stem + "_empty_best"] = empty_score
                extra_scores[extra_pdbqt.stem + "_holo_G"]  = G_holo
                extra_scores[extra_pdbqt.stem + "_empty_G"] = G_empty
            seq = read_sequence(fasta_dir, pdb_file.stem)
            row = {
                "pdb":        pdb_file.name,
                "vina_score": score_orig,
                "cx":         centroid[0],
                "cy":         centroid[1],
                "cz":         centroid[2],
                "angle_deg":  angle_deg,
                "sequence":   seq,
            }
            row.update(extra_scores)
            results.append(row)
        except subprocess.TimeoutExpired as e:
            print(f"[ERROR] Subprocess timed out for {pdb_file.name}: {e}")
        except Exception as e:
            print(f"[ERROR] {pdb_file.name}: {e}")
        finally:
            for tmp in [receptor_pdb, stripped_receptor_pdb, ligand_pdb]:
                if tmp and Path(tmp).exists():
                    Path(tmp).unlink(missing_ok=True)
            gc.collect()
    fieldnames = ["pdb","vina_score","cx","cy","cz","angle_deg","sequence"]
    for extra in args.extra_ligands:
        stem = Path(extra).stem
        fieldnames += [
            stem + "_holo_best",
            stem + "_empty_best",
            stem + "_holo_G",
            stem + "_empty_G"
        ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\n[DONE] Results written to {args.output}")
