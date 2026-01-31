#!/usr/bin/env python3
"""
Interactive demo script for plotlyMol3D.

Run this script to test the molecule visualization capabilities locally.
Each section demonstrates different features - uncomment the ones you want to test.

Usage:
    python demo_visualizations.py

Requirements:
    - pip install -e .  (install package in development mode)
    - Or: pip install plotly numpy rdkit
"""
from pathlib import Path

from plotly.subplots import make_subplots

from plotlymol3d import draw_3D_mol, draw_3D_rep, smiles_to_rdkitmol

# Path to sample data files included with the package
PACKAGE_DIR = Path(__file__).parent / "plotlymol3d"


def demo_smiles_visualization():
    """Demo 1: Visualize molecules from SMILES strings.

    SMILES (Simplified Molecular Input Line Entry System) is a notation
    that represents molecular structures as text strings.
    """
    print("\n=== Demo 1: SMILES Visualization ===")

    # Common molecules to try:
    molecules = {
        "ethanol": "CCO",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "benzene": "c1ccccc1",
        "glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "water": "O",
        "methane": "C",
        "alanine": "CC(N)C(=O)O",
    }

    # Change this to try different molecules
    mol_name = "caffeine"
    smiles = molecules[mol_name]

    print(f"Visualizing {mol_name}: {smiles}")
    fig = draw_3D_rep(smiles=smiles, mode="ball+stick", ambient=0.1)
    return fig


def demo_visualization_modes():
    """Demo 2: Compare different visualization modes.

    Available modes:
    - "ball+stick": Atoms as small spheres, bonds as cylinders (default)
    - "ball": Atoms only as small spheres (scaled VDW radii)
    - "stick": Bonds only
    - "vdw": Van der Waals space-filling model
    """
    print("\n=== Demo 2: Visualization Modes ===")

    smiles = "c1ccccc1"  # Benzene - good for showing modes

    # Try each mode (uncomment the one you want)
    # mode = "ball+stick"  # Default - shows both atoms and bonds
    # mode = "ball"        # Just atoms
    # mode = "stick"       # Just bonds (minimal atoms)
    mode = "vdw"  # Space-filling model

    print(f"Mode: {mode}")
    fig = draw_3D_rep(smiles=smiles, mode=mode, ambient=0.2)
    return fig


def demo_xyz_file():
    """Demo 3: Load molecule from XYZ file.

    XYZ files contain atomic coordinates - commonly output from
    quantum chemistry calculations.
    """
    print("\n=== Demo 3: XYZ File Visualization ===")

    xyz_file = PACKAGE_DIR / "cube.xyz"
    if xyz_file.exists():
        print(f"Loading: {xyz_file}")
        fig = draw_3D_rep(xyzfile=str(xyz_file), mode="ball+stick", ambient=0.1)
        return fig
    else:
        print(f"XYZ file not found: {xyz_file}")
        return None


def demo_mol_file():
    """Demo 4: Load molecule from MOL file.

    MOL files contain both coordinates AND bonding information,
    so they're more reliable than XYZ files.
    """
    print("\n=== Demo 4: MOL File Visualization ===")

    mol_file = PACKAGE_DIR / "cube.mol"
    if mol_file.exists():
        print(f"Loading: {mol_file}")
        fig = draw_3D_rep(molfile=str(mol_file), mode="ball+stick", ambient=0.1)
        return fig
    else:
        print(f"MOL file not found: {mol_file}")
        return None


def demo_orbital_visualization():
    """Demo 5: Visualize molecular orbitals from cube file.

    Cube files contain volumetric data (electron density, orbitals)
    from quantum chemistry calculations. This demo shows how to
    render orbital isosurfaces.
    """
    print("\n=== Demo 5: Orbital Visualization ===")

    cube_file = PACKAGE_DIR / "anto_occ_1-min2.cube"
    mol_file = PACKAGE_DIR / "cube.mol"

    if cube_file.exists() and mol_file.exists():
        print(f"Loading cube: {cube_file}")
        fig = draw_3D_rep(
            cubefile=str(cube_file),
            molfile=str(mol_file),  # Use MOL file for better bond detection
            mode="ball+stick",
            ambient=0.1,
            cubedraw="orbitals",  # Draw orbital isosurfaces
            orbital_opacity=0.25,
            orbital_colors=["darkorange", "darkblue"],  # +/- lobes
        )
        return fig
    else:
        print("Cube/MOL files not found")
        return None


def demo_lighting_options():
    """Demo 6: Experiment with lighting parameters.

    Adjust these to get publication-quality renders:
    - ambient: Background lighting (0-1)
    - diffuse: Scattered light (0-1)
    - specular: Shiny highlights (0-1)
    - roughness: Surface texture (0-1)
    - fresnel: Edge glow effect (0-1)
    """
    print("\n=== Demo 6: Lighting Options ===")

    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin

    fig = draw_3D_rep(
        smiles=smiles,
        mode="vdw",  # Space-filling shows lighting best
        ambient=0.3,
        diffuse=0.8,
        specular=0.5,
        roughness=0.3,
        fresnel=0.2,
        lightx=500,
        lighty=500,
        lightz=1000,
    )
    return fig


def demo_custom_workflow():
    """Demo 7: Custom workflow with lower-level functions.

    For more control, you can use the lower-level functions directly:
    - smiles_to_rdkitmol() - Convert SMILES to RDKit molecule
    - draw_3D_mol() - Draw on existing figure
    """
    print("\n=== Demo 7: Custom Workflow ===")

    # Create figure manually
    fig = make_subplots()

    # Convert SMILES to RDKit molecule
    mol = smiles_to_rdkitmol("CCCC")  # Butane

    # Draw with custom settings
    fig = draw_3D_mol(fig, mol, mode="ball+stick", resolution=64)

    # You can add multiple molecules, annotations, etc.

    fig.show("browser")
    return fig


# =============================================================================
# Run a demo - uncomment the one you want to try
# =============================================================================

if __name__ == "__main__":
    print("plotlyMol3D Demo")
    print("================")
    print("Uncomment a demo function below to test different features.\n")

    # Basic demos
    demo_smiles_visualization()  # Most common use case
    # demo_visualization_modes()       # Compare ball+stick, vdw, etc.

    # File-based demos
    # demo_xyz_file()                  # Load from XYZ coordinates
    # demo_mol_file()                  # Load from MOL file
    # demo_orbital_visualization()     # Quantum chemistry orbitals

    # Advanced demos
    # demo_lighting_options()          # Publication-quality rendering
    # demo_custom_workflow()           # Lower-level API access
