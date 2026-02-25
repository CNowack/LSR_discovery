"""
Core logic for reconstructing attB and attP from MGE boundaries.

Given post-integration structure:
    B1 - D - P1 - [MGE with LSR] - P2 - D - B2

attB = B1 + D + B2   (50bp window around center)
attP = P2 + D + P1   (50bp window around center)
"""

import os
import sys
import pandas as pd

# 1. Grab inputs/outputs
boundaries_file = snakemake.input.boundaries
out_file = snakemake.output.att_sites

# 2. Check if the boundaries file is empty or only contains a header
try:
    df = pd.read_csv(boundaries_file, sep='\t')
    is_empty = df.empty
except pd.errors.EmptyDataError:
    is_empty = True

# 3. If no insertions were found, write an empty file and exit safely
if is_empty:
    with open(out_file, 'w') as f:
        # Write the headers your downstream scripts expect
        f.write("genome\tsite_type\tsequence\n")
    sys.exit(0)

def reconstruct_att_sites(boundary_record):
    b1 = boundary_record["left_flank"]
    b2 = boundary_record["right_flank"]
    d  = boundary_record["target_site_duplication"]
    p1 = boundary_record["left_mge_terminal"]
    p2 = boundary_record["right_mge_terminal"]

    attB_full = b1 + d + b2
    attP_full = p2 + d + p1

    # Extract 50bp around the center (dinucleotide core)
    center = len(attB_full) // 2
    attB = attB_full[center-25 : center+25]
    attP = attP_full[center-25 : center+25]

    return attB, attP