"""Atom-type parameter tables.

Two independent parameter sets live here:

* ``AD_*``   -- AutoDock4-style force-field parameters (van der Waals radii and
                well depths, hydrogen-bond 12-10 parameters, atomic solvation
                volumes and parameters).  Published in Huey et al., J. Comput.
                Chem. 2007, 28, 1145.
* ``PLP_*``  -- piecewise-linear-potential breakpoints for the empirical
                scoring function, in the spirit of Gehlhaar et al., Chem. Biol.
                1995, 2, 317 and the ChemPLP steric term.

Plus a compact Amber-flavoured partial-charge table for standard protein
residues, used by the force-field scoring function.  It is deliberately
approximate: polar and charged atoms carry proper charges, apolar carbon is
near-neutral.  Good enough for pose ranking, not for free-energy work.
"""
from __future__ import annotations

# --------------------------------------------------------------------------- #
# AutoDock-style types
# --------------------------------------------------------------------------- #
# type -> (Rii [A], epsii [kcal/mol], vol [A^3], solpar)
AD_PARAMS: dict[str, tuple[float, float, float, float]] = {
    "C":  (4.00, 0.150, 33.5103, -0.00143),
    "A":  (4.00, 0.150, 33.5103, -0.00052),
    "N":  (3.50, 0.160, 22.4493, -0.00162),
    "NA": (3.50, 0.160, 22.4493, -0.00162),
    "O":  (3.20, 0.200, 17.1573, -0.00251),
    "OA": (3.20, 0.200, 17.1573, -0.00251),
    "S":  (4.00, 0.200, 33.5103, -0.00214),
    "SA": (4.00, 0.200, 33.5103, -0.00214),
    "H":  (2.00, 0.020,  0.0000,  0.00051),
    "HD": (2.00, 0.020,  0.0000,  0.00051),
    "F":  (3.09, 0.080, 15.4480, -0.00110),
    "Cl": (4.09, 0.276, 35.8235, -0.00110),
    "Br": (4.33, 0.389, 42.5661, -0.00110),
    "I":  (4.72, 0.550, 55.0585, -0.00110),
    "P":  (4.20, 0.200, 38.7924, -0.00110),
    "Mg": (1.30, 0.875,  1.5600, -0.00110),
    "Zn": (1.48, 0.550,  1.7000, -0.00110),
    "Ca": (1.98, 0.550,  2.7700, -0.00110),
    "Mn": (1.30, 0.875,  2.1400, -0.00110),
    "Fe": (1.30, 0.010,  1.8400, -0.00110),
}

# Directional hydrogen bond, 12-10 form.
AD_HB_R = 1.90       # A, optimal HD...acceptor separation
AD_HB_EPS = 5.0      # kcal/mol well depth

# AutoDock 4.2 free-energy coefficients.
AD_W_VDW = 0.1662
AD_W_HBOND = 0.1209
AD_W_ELEC = 0.1406
AD_W_DESOLV = 0.1322
AD_W_TORS = 0.2983

# Mehler-Solmajer distance-dependent dielectric.
MS_A = -8.5525
MS_LAMBDA = 0.003627
MS_K = 7.7839
MS_B = 78.4 - MS_A
MS_SIGMA = 0.003627

# Smoothing half-width applied to the vdW/H-bond minima (AutoDock uses 0.25 A).
AD_SMOOTH = 0.25

# --------------------------------------------------------------------------- #
# Empirical (ChemScore-like) coefficients
# --------------------------------------------------------------------------- #
# Eldridge et al., J. Comput. Aided Mol. Des. 1997, 11, 425.  kJ/mol.
CS_DG0 = -5.48
CS_DG_HBOND = -3.34
CS_DG_METAL = -6.03
CS_DG_LIPO = -0.117
CS_DG_ROT = 2.56

# Block-function breakpoints.
CS_HB_R_IDEAL = 1.85     # A, H...acceptor
CS_HB_DR1 = 0.25         # full weight below this deviation
CS_HB_DR2 = 0.65         # zero weight above this deviation
CS_HB_DA1 = 30.0         # deg, deviation of D-H...A from 180
CS_HB_DA2 = 80.0
CS_HB_ACC_A1 = 80.0      # deg, deviation of H...A-AA from the ideal lone-pair
CS_HB_ACC_A2 = 100.0

CS_METAL_R1 = 2.60       # A
CS_METAL_R2 = 3.00

CS_LIPO_R1 = 4.10        # A, heavy-atom separation
CS_LIPO_R2 = 7.10

# Steric repulsion.  Contact distances depend on the polarity of both partners:
# two polar atoms at 2.6 A are hydrogen bonded, two apolar atoms are clashing.
CLASH_R_NP_NP = 3.45     # A, apolar vs apolar
CLASH_R_NP_P = 3.25      # A, apolar vs polar
CLASH_R_P_P = 2.60       # A, polar vs polar
CS_CLASH_WEIGHT = 20.0   # kJ/mol per A^2 of overlap
CS_INTERNAL_WEIGHT = 10.0
CS_OUT_OF_BOX = 25.0     # kJ/mol per atom outside the grid

# Intramolecular repulsion onset, as a fraction of the summed vdW radii.
INTERNAL_SCALE = 0.82

# --------------------------------------------------------------------------- #
# Piecewise linear potential (steric and polar channels)
# --------------------------------------------------------------------------- #
# Each channel is (A, B, C, D, E, F): energy is E at r<A, rises linearly to
# -F... see scoring/plp.py for the exact piecewise form.
#            A     B     C     D      well    wall
PLP_STERIC = (3.4, 3.6, 4.5, 5.5, -0.4, 20.0)
PLP_POLAR = (2.3, 2.6, 3.1, 3.4, -2.0, 20.0)
PLP_REPULSIVE = (3.2, 3.2, 3.2, 3.2, 0.1, 20.0)
PLP_BURIED = (3.4, 3.6, 4.5, 5.5, -0.1, 20.0)

# --------------------------------------------------------------------------- #
# Protein partial charges (approximate, Amber-flavoured)
# --------------------------------------------------------------------------- #
BACKBONE_CHARGES = {
    "N": -0.42, "H": 0.27, "CA": 0.03, "HA": 0.08,
    "C": 0.60, "O": -0.57, "OXT": -0.80, "CB": -0.10,
}

SIDECHAIN_CHARGES: dict[str, dict[str, float]] = {
    "ASP": {"CG": 0.80, "OD1": -0.80, "OD2": -0.80},
    "GLU": {"CD": 0.80, "OE1": -0.80, "OE2": -0.80},
    "LYS": {"NZ": -0.30, "HZ1": 0.33, "HZ2": 0.33, "HZ3": 0.33, "CE": 0.20},
    "ARG": {"CZ": 0.64, "NE": -0.53, "HE": 0.33, "NH1": -0.80, "NH2": -0.80,
            "HH11": 0.46, "HH12": 0.46, "HH21": 0.46, "HH22": 0.46},
    "HIS": {"ND1": -0.55, "HD1": 0.33, "NE2": -0.55, "HE2": 0.33,
            "CE1": 0.20, "CD2": 0.13},
    "SER": {"OG": -0.65, "HG": 0.40},
    "THR": {"OG1": -0.65, "HG1": 0.40, "CB": 0.20},
    "CYS": {"SG": -0.31, "HG": 0.19},
    "TYR": {"OH": -0.56, "HH": 0.40, "CZ": 0.32},
    "ASN": {"CG": 0.60, "OD1": -0.57, "ND2": -0.62, "HD21": 0.33, "HD22": 0.33},
    "GLN": {"CD": 0.60, "OE1": -0.57, "NE2": -0.62, "HE21": 0.33, "HE22": 0.33},
    "TRP": {"NE1": -0.34, "HE1": 0.34},
    "MET": {"SD": -0.27},
}

# Residue/atom names that are rotatable polar hydrogens in the receptor.
# GOLD rotates these during docking; we do the same.
ROTATABLE_OH: dict[str, tuple[str, str, str]] = {
    # residue -> (pivot heavy atom, its parent, the hydrogen)
    "SER": ("OG", "CB", "HG"),
    "THR": ("OG1", "CB", "HG1"),
    "TYR": ("OH", "CZ", "HH"),
    "CYS": ("SG", "CB", "HG"),
}

METALS = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI", "CU", "CD", "NA", "K"}

ELEMENT_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "Se": 1.90,
}
