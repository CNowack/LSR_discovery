"""
Core logic for reconstructing attB and attP from MGE boundaries.

Given post-integration structure:
    B1 - D - P1 - [MGE with LSR] - P2 - D - B2

attB = B1 + D + B2   (50bp window around center)
attP = P2 + D + P1   (50bp window around center)
"""

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