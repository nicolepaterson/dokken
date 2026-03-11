#!/usr/bin/env python
# usage:
# python rosetta.py \
#    --receptor receptor.pdb \
#    --flex_res flexible_residues.resfile \
#    --centroid x,y,z \
#    --log_file packing_log.txt

import argparse
import subprocess
import os
import logging

# -------------------- Setup Logging --------------------
def setup_logger(log_file):
    logger = logging.getLogger("pocket_relax")
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

# -------------------- Helper: Run Rosetta --------------------
def run_rosetta_score(pdb_file, score_bin, logger, label=""):
    cmd = [score_bin, "-s", pdb_file, "-out:file:score_only"]
    logger.info(f"Running {label} on {pdb_file}...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stdout
        # Save verbose output to log
        logger.debug(f"{label} output for {os.path.basename(pdb_file)}:\n{output}")

        # Parse the first numeric score in the output
        score = None
        for line in output.splitlines():
            if line.strip() and not line.startswith("#"):
                try:
                    score = float(line.strip().split()[1])
                    break
                except (IndexError, ValueError):
                    continue
        return score
    except subprocess.CalledProcessError as e:
        logger.error(f"Error scoring {pdb_file} with {label}: {e.stderr}")
        return None

# -------------------- Helper: Run pocket_relax --------------------
def run_pocket_relax(pdb_file, relax_bin, flex_res_file, output_pdb, logger):
    cmd = [
        relax_bin,
        "-s", pdb_file,
        "-relax:constrain_relax_to_start_coords",
        "-relax:coord_cst_stdev", "0.5",
        "-packing:ex1",
        "-packing:ex2",
        "-resfile", flex_res_file,
        "-out:file:o", output_pdb
    ]
    logger.info(f"Running pocket_relax on {pdb_file} -> {output_pdb}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Pocket relax completed: {output_pdb}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running pocket_relax on {pdb_file}: {e.stderr}")
        return False

# -------------------- Main --------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimize receptor sidechains and evaluate packing.")
    parser.add_argument("--receptor", required=True, help="Input receptor PDB")
    parser.add_argument("--flex_res", required=True, help="File with flexible residues (Rosetta resfile)")
    parser.add_argument("--centroid", required=True, help="Centroid coordinates (x,y,z)")
    parser.add_argument("--log_file", default="packing_log.txt", help="Log file for verbose output")
    args = parser.parse_args()

    logger = setup_logger(args.log_file)
    logger.info(f"Starting Pocket Relax for {args.receptor}")

    # Paths to Rosetta executables
    score_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/score.static.macosclangrelease"
    pack_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/packstat.static.macosclangrelease"
    relax_bin = "/Users/nmp/Desktop/rosetta/main/source/bin/pocket_relax.static.macosclangrelease"

    # Score before packing
    score_before = run_rosetta_score(args.receptor, score_bin, logger, label="BEFORE score.static")
    pack_before = run_rosetta_score(args.receptor, pack_bin, logger, label="BEFORE packstat.static")

    # Output PDB after pocket_relax
    output_pdb = args.receptor.replace(".pdb", "_relaxed.pdb")
    success = run_pocket_relax(args.receptor, relax_bin, args.flex_res, output_pdb, logger)
    if not success:
        logger.error("Pocket relax failed. Exiting.")
        exit(1)

    # Score after packing
    score_after = run_rosetta_score(output_pdb, score_bin, logger, label="AFTER score.static")
    pack_after = run_rosetta_score(output_pdb, pack_bin, logger, label="AFTER packstat.static")

    # Write results file
    results_file = f"pocket_relax_results_{os.path.basename(args.receptor).replace('.pdb','')}.txt"
    with open(results_file, "w") as f:
        f.write(f"Pocket Relax Results for {args.receptor}\n")
        f.write("-"*50 + "\n")
        f.write(f"BEFORE score.static: {score_before}\n")
        f.write(f"AFTER score.static: {score_after}\n")
        f.write(f"ΔScore (AFTER - BEFORE): {score_after - score_before if score_before is not None and score_after is not None else 'N/A'}\n\n")
        f.write(f"BEFORE packstat.static: {pack_before}\n")
        f.write(f"AFTER packstat.static: {pack_after}\n")
        f.write(f"ΔPackstat (AFTER - BEFORE): {pack_after - pack_before if pack_before is not None and pack_after is not None else 'N/A'}\n")
    logger.info(f"Results written to {results_file}")
    logger.info("Pocket Relax workflow completed successfully.")

