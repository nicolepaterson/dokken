#!/usr/bin/env python
"""
mk_prepare_receptor.py

Prepare receptor PDB for flexible docking with AutoDock Vina and Smina.
- Compute flexible residues based on GCN-predicted interface residues and SASA.
- Generate receptor PDBQT files for rigid/flexible docking.
- Create docking config files.
- Submit docking jobs for Vina and Smina.

Usage:
python mk_prepare_receptor.py \
    --input receptor.pdb \
    --output out_receptor.pdb \
    --top_res top_residues.txt \
    --centroid 12.3,34.5,56.7 \
    --box_size 25,25,25 \
    --sasa_cutoff 0.3 \
    --vina_exhaustiveness 32 \
    --smina_exhaustiveness 16
"""

import argparse, os, subprocess
from Bio.PDB import PDBParser
import freesasa
import numpy as np
import subprocess
import os



def compute_flexible_residues(pdb_file, top_res, centroid, half_size=10.0, sasa_cutoff=0.3):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    flexible_residues = []

    # Compute SASA
    fs_structure = freesasa.Structure(pdb_file)
    fs_result = freesasa.calc(fs_structure)
    res_sasa_dict = fs_result.residueAreas()  # Returns dict: keys = (chain, resnum, resname), values = float

    for res in structure.get_residues():
        resnum = res.get_id()[1]
        chain_id = res.get_parent().get_id()
        if resnum not in top_res:
            continue
        if 'CA' not in res:
            continue
        ca_coord = np.array(res['CA'].coord)
        if np.all(np.abs(ca_coord - centroid) <= half_size):
            sasa = res_sasa_dict.get((chain_id, resnum, res.get_resname()), 0.0)
            if sasa >= sasa_cutoff:
                flexible_residues.append((resnum, chain_id))

    return flexible_residues

def write_resfile(resfile_path, flexible_residues):
    """Write flexible residues to a Rosetta-style resfile"""
    with open(resfile_path, "w") as f:
        f.write("NATRO\nstart\n")
        for resnum, chain_id in flexible_residues:
            f.write(f"{resnum} {chain_id} NATAA\n")

# -------------------- Docking Config --------------------
def write_vina_config(config_path, center, box_size):
    x, y, z = center
    sx, sy, sz = box_size
    with open(config_path, "w") as f:
        f.write(f"center_x = {x}\n")
        f.write(f"center_y = {y}\n")
        f.write(f"center_z = {z}\n")
        f.write(f"size_x = {sx}\n")
        f.write(f"size_y = {sy}\n")
        f.write(f"size_z = {sz}\n")
        f.write("exhaustiveness = 32\n")
        f.write("num_modes = 20\n")

def write_smina_config(config_path, center, box_size, flex_file=""):
    x, y, z = center
    sx, sy, sz = box_size
    with open(config_path, "w") as f:
        f.write(f"--center_x {x} --center_y {y} --center_z {z} \\\n")
        f.write(f"--size_x {sx} --size_y {sy} --size_z {sz} \\\n")
        if flex_file:
            f.write(f"--flex {flex_file} \\\n")
        f.write("--exhaustiveness 16 --num_modes 20\n")

# -------------------- Helper: Run docking --------------------
def run_vina(receptor, ligand, config, output):
    cmd = f"vina --receptor {receptor} --ligand {ligand} --config {config} --out {output}"
    print(f"Running Vina: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def run_smina(receptor, ligand, config, output):
    cmd = f"smina --receptor {receptor} --ligand {ligand} {config} --out {output}"
    print(f"Running Smina: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# -------------------- Main --------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare receptor for flexible docking.")
    parser.add_argument("--input", required=True, help="Input receptor PDB")
    parser.add_argument("--output", required=True, help="Output receptor PDB")
    parser.add_argument("--top_res", required=True, nargs='+', type=int,help="List of top residues for flexibility (space-separated)")
    parser.add_argument("--centroid", required=True, nargs=3, type=float,help="Centroid coordinates (x y z)")
    parser.add_argument("--box_size", nargs=3, type=float, default=[25,25,25],help="Bounding box size in Å (x y z)")
    parser.add_argument("--sasa_cutoff", type=float, default=0.3, help="SASA cutoff for flexibility")
    parser.add_argument("--vina_exhaustiveness", type=int, default=32, help="Exhaustiveness for Vina")
    parser.add_argument("--smina_exhaustiveness", type=int, default=16, help="Exhaustiveness for Smina")
    parser.add_argument("--ligands", required=True, nargs='+', help="Ligand pdbqt files")
    args = parser.parse_args()
    centroid = args.centroid
    top_res = args.top_res
    box_size = args.box_size 
    #box_size = [float(x) for x in args.box_size.split(",")]

    # Determine flexible residues
    flexible_residues = compute_flexible_residues(args.input, top_res, centroid, half_size=max(box_size)/2, sasa_cutoff=args.sasa_cutoff)
    print(f"Flexible residues: {flexible_residues}")

    resfile_path = args.output.replace(".pdb","_flex.resfile")
    write_resfile(resfile_path, flexible_residues)
    print(f"Resfile written: {resfile_path}")

    # Prepare docking config files
    vina_config = args.output.replace(".pdb","_vina_config.txt")
    write_vina_config(vina_config, centroid, box_size)

    smina_config = args.output.replace(".pdb","_smina_config.txt")
    write_smina_config(smina_config, centroid, box_size, flex_file=resfile_path)


# --- Convert PDB to PDBQT using Meeko ---
def prepare_pdbqt(pdb_file, flexible_residues, rigid_out, flex_out=None):
    # Load receptor
    receptor = PDBQTMolecule()
    receptor.read_pdb(pdb_file)

    # Write rigid receptor
    receptor.write(rigid_out)
    print(f"Rigid receptor written: {rigid_out}")

    # Write flexible receptor if flexible residues provided
    if flex_out and flexible_residues:
        receptor.write(flex_out, flexible_residues=flexible_residues)
        print(f"Flexible receptor written: {flex_out}")

# --- Docking helpers ---
def run_vina(receptor, ligand, config, output):
    cmd = f"vina --receptor {receptor} --ligand {ligand} --config {config} --out {output}"
    print(f"Running Vina: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def run_smina(receptor, ligand, config, output):
    cmd = f"smina --receptor {receptor} --ligand {ligand} --config {config} --out {output}"
    print(f"Running Smina: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# --- Usage inside main ---
rigid_pdbqt = args.output.replace(".pdb","_rigid.pdbqt")
flex_pdbqt = args.output.replace(".pdb","_flex.pdbqt")

prep = ReceptorPreparation()
prep.prepare(args.input)

# Write rigid
prep.write_pdbqt_file(rigid_pdbqt)

# Write flexible if residues exist
if flexible_residues:
    flex_list = [f"{chain}:{resnum}" for resnum, chain in flexible_residues]
    prep.write_pdbqt_file(flex_pdbqt, flexible_residues=flex_list)

# Run docking
for ligand in args.ligands:
    vina_out = ligand.replace(".pdbqt","_vina_out.pdbqt")
    run_vina(rigid_pdbqt, ligand, vina_config, vina_out)

    smina_out = ligand.replace(".pdbqt","_smina_out.pdbqt")
    run_smina(rigid_pdbqt, ligand, smina_config, smina_out)

    # Dock against ligands
    for ligand in args.ligands:
        vina_out = ligand.replace(".pdbqt","_vina_out.pdbqt")
        run_vina(rigid_pdbqt, ligand, vina_config, vina_out)

        smina_out = ligand.replace(".pdbqt","_smina_out.pdbqt")
        run_smina(rigid_pdbqt, ligand, smina_config, smina_out)

