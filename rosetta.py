#!usr/bin/env python 
# usage notes
# python rosetta.py \
#    --receptor receptor.pdb \
#    --flex_res flexible_residues.resfile \
#    --centroid x,y,z \
#    --log_file packing_log.txt

import argparse
import subprocess
import os

# -------------------- Helper: Run Rosetta --------------------
def run_rosetta_score(pdb_file, score_bin, log_file="packing_log.txt", label=""):
    cmd = [score_bin, "-s", pdb_file, "-out:file:score_only"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stdout
        with open(log_file, "a") as f:
            f.write(f"\n==== {label} {os.path.basename(pdb_file)} ====\n")
            f.write(output)
        print(f"Scored {pdb_file} ({label}). Results appended to {log_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error scoring {pdb_file}: {e.stderr}")

# -------------------- Helper: Run pocket_relax --------------------
def run_pocket_relax(pdb_file, relax_bin, flex_res_file, output_pdb):
    cmd = [
        relax_bin,
        "-s", pdb_file,
        "-relax:constrain_relax_to_start_coords",
        "-relax:coord_cst_stdev 0.5",
        "-packing:ex1",
        "-packing:ex2",
        "-resfile", flex_res_file,  # flexible residues
        "-out:file:o", output_pdb
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"Pocket relax completed: {output_pdb}")
    except subprocess.CalledProcessError as e:
        print(f"Error running pocket_relax: {e.stderr}")

# -------------------- Main --------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimize receptor sidechains and evaluate packing.")
    parser.add_argument("--receptor", required=True, help="Input receptor PDB")
    parser.add_argument("--flex_res", required=True, help="File with flexible residues (Rosetta resfile)")
    parser.add_argument("--centroid", required=True, help="Centroid coordinates (x,y,z)")
    parser.add_argument("--log_file", default="packing_log.txt", help="Log file for scores")
    args = parser.parse_args()

    # Paths to Rosetta executables
    score_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/score.static.macosclangrelease"
    pack_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/packstat.static.macosclangrelease"
    relax_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/pocket_relax.static.macosclangrelease"

    # Score before packing
    run_rosetta_score(args.receptor, score_bin, args.log_file, label="BEFORE score.static")
    run_rosetta_score(args.receptor, pack_bin, args.log_file, label="BEFORE packstat.static")

    # Output PDB after pocket_relax
    output_pdb = args.receptor.replace(".pdb", "_relaxed.pdb")
    run_pocket_relax(args.receptor, relax_bin, args.flex_res, output_pdb)

    # Score after packing
    run_rosetta_score(output_pdb, score_bin, args.log_file, label="AFTER score.static")
    run_rosetta_score(output_pdb, pack_bin, args.log_file, label="AFTER packstat.static")

