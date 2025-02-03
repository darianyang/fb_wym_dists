"""
Convert a comma-separated list of PDB IDs to a newline-separated list.
"""

import sys

# string args
file_in = str(sys.argv[1])
file_out = str(sys.argv[2])

# Read the input file
with open(file_in, "r") as f:
    data = f.read()

# Split the text by commas and strip any whitespace
pdb_ids = [entry.strip() for entry in data.split(",")]

# Write to output file with each entry on a new line
with open(file_out, "w") as f:
    for pdb_id in pdb_ids:
        f.write(pdb_id + "\n")

print(f"{file_in} conversion complete! Check {file_out} for the formatted list.")
