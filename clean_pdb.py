from Bio.PDB import PDBParser, PDBIO
import os
import sys
# Input/output files
input_pdb = sys.argv[1]
base = os.path.splitext(os.path.basename(input_pdb))[0]
output_pdb = f"{base}_cleaned.pdb"
# Parse structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("rec", input_pdb)

# Process residues
for model in structure:
    for chain in model:
        for res in list(chain):  # list() because we may detach residues
            resname = res.get_resname().strip()
            if resname == "HIS":
                res.resname = "HIE"          # protonation at pH 7
            elif resname == "UNK":
                chain.detach_child(res.id)   # remove unknown residues

# Save cleaned PDB
io = PDBIO()
io.set_structure(structure)
io.save(output_pdb)

print(f"[INFO] Cleaned PDB saved as {output_pdb}")

