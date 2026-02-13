#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example usage and test cases for plotlyMol3D.

This module demonstrates the various ways to use the plotlyMol3D visualization
functions, including SMILES input, XYZ files, cube files with orbitals, and
combined visualizations.

Created on Thu Jun 27 11:41:46 2024
@author: benjaminlear
"""

from pathlib import Path

from plotlyMol3D import draw_3D_mol, draw_3D_rep, cubefile_to_xyzblock
from rdkit import Chem

# Get the directory containing this file for relative paths to sample data
PACKAGE_DIR = Path(__file__).parent.resolve()

# %% smiles test
# Draw a molecule directly from a SMILES string - no file paths needed
moldraw = draw_3D_rep(smiles="CCNCOCSC", mode="ball", ambient=0.1)

# %% xyztest
# Draw from an XYZ file - uses sample data included with the package
# Note: cube.xyz is included in the package; for your own files, provide the path
moldraw = draw_3D_rep(
    xyzfile=str(PACKAGE_DIR / "cube.xyz"), mode="ball+stick", ambient=0.1
)

# %% Cubetest
# Draw a molecule with orbital visualization from cube file
# Uses sample cube and mol files included with the package
moldraw = draw_3D_rep(
    cubefile=str(PACKAGE_DIR / "anto_occ_1-min2.cube"),
    molfile=str(PACKAGE_DIR / "cube.mol"),
    mode="ball+stick",
    ambient=0.1,
    cubedraw="orbitals",
    orbital_opacity=0.25,
    orbital_colors=["darkorange", "darkblue"],
)

# %% multidraw test
# Combine multiple input types in a single visualization
# Note: For xyz file, using the sample file included with the package
moldraw = draw_3D_rep(
    smiles="CCNCOCSC",
    xyzfile=str(PACKAGE_DIR / "cube.xyz"),
    mode="ball+stick",
    ambient=0.1,
)

# %% Load mol file directly with RDKit
m = Chem.MolFromMolFile(str(PACKAGE_DIR / "cube.mol"))

# %% Extract xyz coordinates from a cube file
testblock = cubefile_to_xyzblock(str(PACKAGE_DIR / "anto_occ_1-min2.cube"))


# =============================================================================
# Notes and archived code snippets
# =============================================================================

# The following code demonstrates XYZ processing approaches that were explored
# but had issues (kept for reference):

# def process_xyz_coords(xyz, charge=0):
#     """Process XYZ coordinates using RDKit bond determination.
#
#     Note: rdDetermineBonds.DetermineConnectivity can hang on import in some
#     environments. This approach may need alternative implementations.
#     """
#     raw_mol = Chem.MolFromXYZBlock(xyz)
#     from rdkit.Chem import rdDetermineBonds  # <-- can hang on import
#     conn_mol = Chem.Mol(raw_mol)
#
#     rdDetermineBonds.DetermineConnectivity(conn_mol)
#     rdDetermineBonds.DetermineBondOrders(conn_mol, charge=charge)
#
#     atoms = conn_mol.GetAtoms()
#     bonds = conn_mol.GetBonds()
#
#     return atoms, bonds

# Example XYZ block format for reference:
# testblock = '''9
# 	Energy:       4.3719968
# C         -5.60141        3.84674       -0.07080
# C         -4.10301        3.83876       -0.08074
# H         -3.72298        4.87682        0.04590
# H         -3.75288        3.44948       -1.06129
# N         -3.59011        2.98389        0.98445
# H         -3.85276        3.38964        1.91098
# H         -2.54705        2.97307        0.92691
# O         -6.21202        4.25745        0.90282
# H         -6.15086        3.49213       -0.93685
# '''
