#!/usr/bin/env python

"""
python prepare_receptor2.py \
    --input rosetta_results/H3N6_Chicken_Nanchang_sample5_0001.pdb \
    --output H3N6_Chicken_Nanchang_sample5_prepped.pdb \
    --ligands 2_3_smina.pdbqt 2_6_smina.pdbqt \
    --top_res 124 241 200 243 123 132 153 125 199 133 152 209 122 208 198 168 205 \
    --centroid 10.31 -17.82 -44.84 \
    --box_size 25 25 25


Flexible docking pipeline using:
- freesasa (surface exposure filter)
- GCN predicted residues
- mk_prepare_receptor.py (Meeko CLI)
- Smina
"""
#!/usr/bin/env python

"""
prepare_receptor_smina.py

Flexible docking pipeline using:
- freesasa (surface exposure filter)
- GCN predicted residues
- mk_prepare_receptor.py (Meeko CLI)
- SMINA only (no Vina)

Outputs:
- rigid receptor PDBQT
- flexible receptor PDBQT (if residues selected)
- docked ligand files
"""

import argparse
import subprocess
from Bio.PDB import PDBParser
import freesasa
import numpy as np
import os


# ------------------------------------------------------------
# Flexible Residue Selection
# ------------------------------------------------------------
def compute_flexible_residues(
    pdb_file,
    #top_res,
    centroid,
    half_size=12.5,
    #rsa_cutoff=0.3,  # now using relative SASA
    output_prefix="receptor"
):
    """
    Identify flexible residues for docking based on:
    - Predicted top residues (top_res)
    - 3D box around centroid
    - Relative SASA cutoff (rsa_cutoff)

    Returns:
        flexible_residues: list of "Chain:Resnum" strings
    Also writes a SASA + RSA report: {output_prefix}_rsa_report.txt
    """

    # Max solvent accessible surface area (Å^2) per residue
    MAX_ASA = {
        "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193,
        "CYS": 167, "GLN": 223, "GLU": 225, "GLY": 104,
        "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
        "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
        "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", pdb_file)

    # Compute SASA
    fs_structure = freesasa.Structure(pdb_file)
    fs_result = freesasa.calc(fs_structure)
    res_sasa = fs_result.residueAreas()

    flexible_residues = []
    report_lines = []

    for res in structure.get_residues():
        chain = res.get_parent().get_id()
        resnum = res.get_id()[1]
        icode = res.get_id()[2].strip()
        resnum_str = f"{resnum}{icode}" if icode else str(resnum)
        resname = res.get_resname()

        # Skip residues not in MAX_ASA (unknown residues)
        if resname not in MAX_ASA:
            continue

        # Skip residues not in top list or missing CA
        if "CA" not in res:
        #resnum not in top_res or 
            in_box = False
            passes_rsa = False
            sasa_val = 0.0
            rsa_val = 0.0
        else:
            coord = np.array(res["CA"].coord)
            in_box = np.all(np.abs(coord - centroid) <= half_size)

            try:
                sasa_val = res_sasa[chain][resnum_str].total
            except KeyError:
                sasa_val = 0.0

            rsa_val = sasa_val / MAX_ASA[resname]  # Relative SASA
            passes_rsa = (rsa_val >= 0.1) and in_box

            if passes_rsa:
                flexible_residues.append(f"{chain}:{resnum}")

        report_lines.append(
            f"{chain}\t{resnum}\t{resname}\t{sasa_val:.4f}\t{rsa_val:.3f}\t{in_box}\t{passes_rsa}"
        )

    # Write RSA report
    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] RSA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues

    # Write SASA/RSA report
    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] RSA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues


    # Write RSA report
    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:l
    f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
    f.write("--------------------------------------------------------------\n")
    for line in report_lines:
        f.write(line + "\n")

        print(f"[INFO] RSA report written to {report_file}")
        print(f"[INFO] Flexible residues selected: {flexible_residues}")
    # Return Meeko-ready string
    return ",".join(flexible_residues)

    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] RSA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues


    # Write RSA report
    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] RSA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues

    # Write SASA/RSA report
    report_file = f"{output_prefix}_rsa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tRSA\tInBox\tPassesRSA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] RSA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues

    # Write SASA report
    report_file = f"{output_prefix}_sasa_report.txt"
    with open(report_file, "w") as f:
        f.write("Chain\tResnum\tResname\tSASA\tInBox\tPassesSASA\n")
        f.write("--------------------------------------------------------------\n")
        for line in report_lines:
            f.write(line + "\n")

    print(f"[INFO] SASA report written to {report_file}")
    print(f"[INFO] Flexible residues selected: {flexible_residues}")

    return flexible_residues

    # --------------------------------------------------------
    # Write SASA report
    # --------------------------------------------------------
    if output_prefix:

        sasa_file = f"{output_prefix}_sasa_report.txt"

        with open(sasa_file, "w") as f:
            f.write(
                "Chain\tResnum\tResname\tSASA\tInBox\tPassesSASA\n"
            )
            f.write(
                "--------------------------------------------------------------\n"
            )

            for row in report_rows:
                chain, resnum, resname, sasa, in_box, passes_sasa = row
                f.write(
                    f"{chain}\t{resnum}\t{resname}\t"
                    f"{sasa:.4f}\t{in_box}\t{passes_sasa}\n"
                )

        print(f"[INFO] SASA report written to {sasa_file}")

    return flexible

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument("--ligands", nargs="+", required=True)
    #parser.add_argument("--top_res", nargs="+", type=int, required=True)

    parser.add_argument("--centroid", nargs=3, type=float, required=True)
    parser.add_argument("--box_size", nargs=3, type=float, default=[25, 25, 25])

    parser.add_argument("--rsa_cutoff", type=float, default=0.1)
    parser.add_argument("--smina_exhaustiveness", type=int, default=16)
    parser.add_argument("--num_modes", type=int, default=20)

    args = parser.parse_args()

    centroid = np.array(args.centroid)
    half_size = max(args.box_size) / 2

    print("\nSelecting flexible residues...")

    flexible_residues = compute_flexible_residues(
    pdb_file=args.input,
    #top_res=args.top_res,
    centroid=centroid,
    half_size=half_size,
    #rsa_cutoff=args.rsa_cutoff,
    output_prefix=args.output.replace(".pdb","")
)
    print("Flexible residues:")
    print(flexible_residues)

    # --------------------------------------------------------
    # Prepare receptor with Meeko
    # --------------------------------------------------------
    
    base = args.output.replace(".pdb", "")
    rigid_pdbqt = f"{base}_rigid.pdbqt"
    flex_pdbqt = f"{base}_flex.pdbqt"

# Join flexible residues as a comma-separated string
    flex_str = ",".join(flexible_residues) if flexible_residues else ""

    print("\nRunning mk_prepare_receptor.py...")

# Build command
    cmd = [
    "mk_prepare_receptor.py",
    "-i", args.input,
    "-o", base,
    "-p",
    "-g",
    "-v",
    "--box_center",
    str(args.centroid[0]),
    str(args.centroid[1]),
    str(args.centroid[2]),
    "--box_size",
    str(args.box_size[0]),
    str(args.box_size[1]),
    str(args.box_size[2]),
    "-a"
]

# Add flexible residues if specified
    if flex_str:        
        cmd += ["-f", flex_str]

# Optionally print the command for debugging
        print("Command:", " ".join(cmd))

# Run the command
    import subprocess
    subprocess.run(cmd, check=True)


    # --------------------------------------------------------
    # Dock with SMINA
    # --------------------------------------------------------
    for ligand in args.ligands:
    # Output file
        ligand_base = os.path.splitext(os.path.basename(ligand))[0]
        output_file = f"{ligand_base}_smina_out.pdbqt"

        print(f"\nDocking {ligand} with Smina...")

    # Start building the docking command
        dock_cmd = [
        "smina",
        "--receptor", rigid_pdbqt,
        "--ligand", ligand,
        "--center_x", str(args.centroid[0]),
        "--center_y", str(args.centroid[1]),
        "--center_z", str(args.centroid[2]),
        "--size_x", str(args.box_size[0]),
        "--size_y", str(args.box_size[1]),
        "--size_z", str(args.box_size[2]),
        "--exhaustiveness", str(args.smina_exhaustiveness),
        "--num_modes", str(args.num_modes),
        "--out", output_file,
        ]

    # Only add flexible receptor if residues were selected
        if flexible_residues and os.path.exists(flex_pdbqt):
            dock_cmd += ["--flex", flex_pdbqt]

    # Print the final command for debugging
        print("Running:", " ".join(dock_cmd))
    
    # Run SMINA
        subprocess.run(dock_cmd, check=True)



    print("\nSMINA docking complete.\n")
