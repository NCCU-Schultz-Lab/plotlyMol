#!/usr/bin/env python3
"""
Streamlit GUI for plotlyMol3D - Visual Testing & Demo App.

A simple web interface for testing and demonstrating molecule visualizations.
Useful for visual confirmation during development and as a demo for users.

Run with:
    streamlit run gui_app.py

Requirements:
    pip install streamlit
"""
import streamlit as st
from pathlib import Path

# Must be first Streamlit command
st.set_page_config(
    page_title="plotlyMol3D Viewer",
    page_icon="🧪",
    layout="wide",
)

# Import visualization functions
from plotlymol3d import (
    draw_3D_mol,
    smiles_to_rdkitmol,
    xyzfile_to_xyzblock,
    xyzblock_to_rdkitmol,
    cubefile_to_xyzblock,
    rdkitmol_to_atoms_bonds_lists,
    draw_atoms,
    draw_bonds,
    format_figure,
    format_lighting,
)
from plotlymol3d.cube import draw_cube_orbitals
from plotly.subplots import make_subplots
from rdkit import Chem
import tempfile
import os

# Sample data directory
PACKAGE_DIR = Path(__file__).parent / "plotlymol3d"

# =============================================================================
# Sidebar Controls
# =============================================================================

st.sidebar.title("🧪 plotlyMol3D")
st.sidebar.markdown("Interactive 3D Molecular Visualization")

# Input method selection
input_method = st.sidebar.radio(
    "Input Method",
    ["SMILES", "MOL File", "XYZ File", "Cube File", "Sample Molecules"],
    index=0,
)

# Visualization mode
mode = st.sidebar.selectbox(
    "Visualization Mode",
    ["ball+stick", "ball", "vdw", "stick"],
    index=0,
    help="ball+stick: atoms and bonds | ball: atoms only | vdw: space-filling | stick: thin atoms",
)

# Lighting settings in expander
with st.sidebar.expander("⚡ Lighting Settings", expanded=False):
    ambient = st.slider("Ambient", 0.0, 1.0, 0.2, 0.05)
    diffuse = st.slider("Diffuse", 0.0, 1.0, 0.8, 0.05)
    specular = st.slider("Specular", 0.0, 1.0, 0.3, 0.05)
    roughness = st.slider("Roughness", 0.0, 1.0, 0.5, 0.05)
    fresnel = st.slider("Fresnel", 0.0, 1.0, 0.1, 0.05)

# Resolution setting
resolution = st.sidebar.slider("Resolution", 8, 64, 32, 8, help="Higher = smoother spheres")

# =============================================================================
# Helper Functions
# =============================================================================

def create_figure_from_mol(rdkitmol, mode, resolution, ambient, diffuse, specular, roughness, fresnel):
    """Create a Plotly figure from an RDKit molecule."""
    fig = make_subplots()
    fig = format_figure(fig)
    fig = draw_3D_mol(fig, rdkitmol, mode=mode, resolution=resolution)
    fig = format_lighting(
        fig,
        ambient=ambient,
        diffuse=diffuse,
        specular=specular,
        roughness=roughness,
        fresnel=fresnel,
    )
    # Update layout for Streamlit
    fig.update_layout(
        height=600,
        margin=dict(l=0, r=0, t=30, b=0),
    )
    return fig


def display_molecule_info(rdkitmol):
    """Display molecule information in sidebar."""
    st.sidebar.markdown("---")
    st.sidebar.markdown("### Molecule Info")
    st.sidebar.write(f"**Atoms:** {rdkitmol.GetNumAtoms()}")
    st.sidebar.write(f"**Bonds:** {rdkitmol.GetNumBonds()}")
    
    # Get atom counts by element
    atom_counts = {}
    for atom in rdkitmol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    formula = "".join(f"{sym}{cnt if cnt > 1 else ''}" for sym, cnt in sorted(atom_counts.items()))
    st.sidebar.write(f"**Formula:** {formula}")


# =============================================================================
# Main Content Area
# =============================================================================

st.title("Molecule Viewer")

# Initialize molecule variable
rdkitmol = None
show_orbitals = False
cube_path = None

# --- SMILES Input ---
if input_method == "SMILES":
    st.markdown("### Enter SMILES String")
    
    # Initialize session state for SMILES using the widget key directly
    if "smiles_input" not in st.session_state:
        st.session_state.smiles_input = "c1ccccc1"  # Benzene as default
    
    def set_random_smiles():
        import random
        examples = [
            ("CCO", "Ethanol"),
            ("c1ccccc1", "Benzene"),
            ("CC(=O)O", "Acetic acid"),
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine"),
            ("CC(N)C(=O)O", "Alanine"),
            ("C1CCCCC1", "Cyclohexane"),
            ("c1ccc2ccccc2c1", "Naphthalene"),
            ("CCCCCCCC", "Octane"),
            ("C=C", "Ethene"),
            ("C#C", "Ethyne"),
            ("C1=CC=C(C=C1)C=O", "Benzaldehyde"),
            ("CC(=O)OC1=CC=CC=C1C(=O)O", "Aspirin"),
        ]
        choice = random.choice(examples)
        st.session_state.smiles_input = choice[0]
        st.session_state.random_molecule_name = choice[1]
    
    col1, col2 = st.columns([3, 1])
    with col1:
        smiles_input = st.text_input(
            "SMILES",
            placeholder="Enter a SMILES string (e.g., CCO for ethanol)",
            label_visibility="collapsed",
            key="smiles_input",
        )
    with col2:
        if st.button("🎲 Random", help="Try a random molecule"):
            set_random_smiles()
            st.rerun()
    
    # Show toast after rerun if we just selected a random molecule
    if "random_molecule_name" in st.session_state:
        st.toast(f"🎲 Selected: {st.session_state.random_molecule_name}")
        del st.session_state.random_molecule_name
    
    if st.session_state.smiles_input:
        try:
            rdkitmol = smiles_to_rdkitmol(st.session_state.smiles_input)
            st.success(f"✅ Parsed: {Chem.MolToSmiles(rdkitmol)}")
        except Exception as e:
            st.error(f"❌ Invalid SMILES: {e}")

# --- MOL File Upload ---
elif input_method == "MOL File":
    st.markdown("### Upload MOL File")
    
    uploaded_file = st.file_uploader("Choose a .mol file", type=["mol", "sdf"])
    
    if uploaded_file is not None:
        try:
            mol_content = uploaded_file.read().decode("utf-8")
            rdkitmol = Chem.MolFromMolBlock(mol_content)
            if rdkitmol is None:
                st.error("❌ Could not parse MOL file")
            else:
                st.success(f"✅ Loaded: {uploaded_file.name}")
        except Exception as e:
            st.error(f"❌ Error reading file: {e}")

# --- XYZ File Upload ---
elif input_method == "XYZ File":
    st.markdown("### Upload XYZ File")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        uploaded_file = st.file_uploader("Choose a .xyz file", type=["xyz"])
    with col2:
        charge = st.number_input("Molecular Charge", value=0, min_value=-5, max_value=5)
    
    if uploaded_file is not None:
        try:
            xyz_content = uploaded_file.read().decode("utf-8")
            rdkitmol = xyzblock_to_rdkitmol(xyz_content, charge=charge)
            st.success(f"✅ Loaded: {uploaded_file.name}")
        except Exception as e:
            st.error(f"❌ Error: {e}")
            st.info("💡 Tip: XYZ bond detection can be tricky. Try specifying the correct charge, or use a MOL file instead.")

# --- Cube File Upload ---
elif input_method == "Cube File":
    st.markdown("### Upload Cube File (Orbital Visualization)")
    
    uploaded_file = st.file_uploader("Choose a .cube file", type=["cube", "cub"])
    
    col1, col2 = st.columns(2)
    with col1:
        show_molecule = st.checkbox("Show Molecule", value=True)
    with col2:
        show_orbitals = st.checkbox("Show Orbitals", value=True)
    
    if show_orbitals:
        col1, col2 = st.columns(2)
        with col1:
            orbital_opacity = st.slider("Orbital Opacity", 0.1, 1.0, 0.3, 0.05)
        with col2:
            pos_color = st.color_picker("Positive Lobe", "#FF8C00")
            neg_color = st.color_picker("Negative Lobe", "#1E90FF")
    
    if uploaded_file is not None:
        try:
            # Save to temp file (cube processing needs file path)
            with tempfile.NamedTemporaryFile(delete=False, suffix=".cube") as tmp:
                tmp.write(uploaded_file.read())
                cube_path = tmp.name
            
            if show_molecule:
                xyzblock, cube_charge = cubefile_to_xyzblock(cube_path)
                rdkitmol = xyzblock_to_rdkitmol(xyzblock, charge=cube_charge)
            
            st.success(f"✅ Loaded: {uploaded_file.name}")
        except Exception as e:
            st.error(f"❌ Error: {e}")

# --- Sample Molecules ---
elif input_method == "Sample Molecules":
    st.markdown("### Select a Sample Molecule")
    
    samples = {
        "Ethanol": "CCO",
        "Benzene": "c1ccccc1",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Alanine": "CC(N)C(=O)O",
        "Cyclohexane": "C1CCCCC1",
        "Naphthalene": "c1ccc2ccccc2c1",
        "Methane": "C",
        "Water": "O",
    }
    
    selected = st.selectbox("Choose molecule", list(samples.keys()))
    smiles = samples[selected]
    st.code(smiles, language=None)
    
    try:
        rdkitmol = smiles_to_rdkitmol(smiles)
    except Exception as e:
        st.error(f"❌ Error: {e}")

# =============================================================================
# Render Visualization
# =============================================================================

if rdkitmol is not None:
    # Show molecule info
    display_molecule_info(rdkitmol)
    
    # Create and display figure
    fig = create_figure_from_mol(
        rdkitmol, mode, resolution,
        ambient, diffuse, specular, roughness, fresnel
    )
    
    # Add orbitals if cube file
    if show_orbitals and cube_path is not None:
        try:
            draw_cube_orbitals(fig, cube_path, orbital_opacity, [pos_color, neg_color])
        except Exception as e:
            st.warning(f"⚠️ Could not render orbitals: {e}")
    
    # Display the plot
    st.plotly_chart(fig, use_container_width=True)
    
    # Clean up temp file
    if cube_path and os.path.exists(cube_path):
        os.unlink(cube_path)

elif input_method not in ["Sample Molecules"]:
    st.info("👆 Enter a molecule above to visualize it")

# =============================================================================
# Footer
# =============================================================================

st.sidebar.markdown("---")
st.sidebar.markdown(
    """
    **Quick Examples:**
    - `CCO` - Ethanol
    - `c1ccccc1` - Benzene  
    - `CC(=O)O` - Acetic acid
    - `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` - Caffeine
    """
)

st.sidebar.markdown("---")
st.sidebar.caption("plotlyMol3D v0.1.0 | [GitHub](https://github.com/jonathanschultzNU/plotlyMol)")
