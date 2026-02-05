"""
Streamlit GUI for plotlyMol3D - Visual Testing & Demo App.

Run with:
    streamlit run examples/gui_app.py
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
import json

import streamlit as st
from plotly.subplots import make_subplots
import plotly.io as pio
from rdkit import Chem

from plotlymol3d import (
    cubefile_to_xyzblock,
    draw_3D_mol,
    format_figure,
    format_lighting,
    smiles_to_rdkitmol,
    xyzblock_to_rdkitmol,
    parse_vibrations,
    add_vibrations_to_figure,
    create_vibration_animation,
)
from plotlymol3d.cube import draw_cube_orbitals

CONFIG_PATH = Path(__file__).resolve().parents[2] / ".plotlymol3d_config.json"


@st.cache_resource(show_spinner=False)
def cached_smiles_to_mol(smiles: str):
    return smiles_to_rdkitmol(smiles)


@st.cache_resource(show_spinner=False)
def cached_xyzblock_to_mol(xyzblock: str, charge: int = 0):
    return xyzblock_to_rdkitmol(xyzblock, charge=charge)


@st.cache_resource(show_spinner=False)
def cached_molblock_to_mol(molblock: str):
    return Chem.MolFromMolBlock(molblock)


@st.cache_data(show_spinner=False)
def cached_cube_bytes_to_xyzblock(cube_bytes: bytes):
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".cube") as tmp:
            tmp.write(cube_bytes)
            temp_path = tmp.name
        return cubefile_to_xyzblock(temp_path)
    finally:
        if temp_path and os.path.exists(temp_path):
            os.unlink(temp_path)


@st.cache_resource(show_spinner=False)
def cached_parse_vibrations(vib_bytes: bytes, filename: str):
    """Parse vibration file from bytes."""
    temp_path = None
    try:
        # Determine file extension
        suffix = Path(filename).suffix
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(vib_bytes)
            temp_path = tmp.name
        return parse_vibrations(temp_path)
    finally:
        if temp_path and os.path.exists(temp_path):
            os.unlink(temp_path)


@st.cache_data(show_spinner=False)
def cached_figure_from_mol_pickle(
    mol_pkl: bytes,
    mode: str,
    resolution: int,
    ambient: float,
    diffuse: float,
    specular: float,
    roughness: float,
    fresnel: float,
):
    """Cache figures using pickled mol to preserve hydrogens."""
    import pickle

    mol = pickle.loads(mol_pkl)
    return create_figure_from_mol(
        mol,
        mode,
        resolution,
        ambient,
        diffuse,
        specular,
        roughness,
        fresnel,
    )


@st.cache_data(show_spinner=False)
def cached_image_bytes(
    mol_pkl: bytes,
    mode: str,
    resolution: int,
    ambient: float,
    diffuse: float,
    specular: float,
    roughness: float,
    fresnel: float,
    width: int,
    height: int,
    fmt: str,
):
    fig = cached_figure_from_mol_pickle(
        mol_pkl,
        mode,
        resolution,
        ambient,
        diffuse,
        specular,
        roughness,
        fresnel,
    )
    fig.update_layout(width=width, height=height)
    return pio.to_image(fig, format=fmt)


def create_figure_from_mol(
    rdkitmol,
    mode,
    resolution,
    ambient,
    diffuse,
    specular,
    roughness,
    fresnel,
):
    """Create a Plotly figure from an RDKit molecule."""
    fig = make_subplots()
    fig = format_figure(fig)  # Transparent background to match theme
    fig = draw_3D_mol(fig, rdkitmol, mode=mode, resolution=resolution)
    fig = format_lighting(
        fig,
        ambient=ambient,
        diffuse=diffuse,
        specular=specular,
        roughness=roughness,
        fresnel=fresnel,
    )
    fig.update_layout(
        height=600,
        margin={"l": 0, "r": 0, "t": 30, "b": 0},
    )
    return fig


def display_molecule_info(rdkitmol):
    """Display molecule information in sidebar."""
    st.sidebar.markdown("---")
    st.sidebar.markdown("### Molecule Info")
    st.sidebar.write(f"**Atoms:** {rdkitmol.GetNumAtoms()}")
    st.sidebar.write(f"**Bonds:** {rdkitmol.GetNumBonds()}")

    atom_counts = {}
    for atom in rdkitmol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    formula = "".join(
        f"{sym}{cnt if cnt > 1 else ''}" for sym, cnt in sorted(atom_counts.items())
    )
    st.sidebar.write(f"**Formula:** {formula}")


def main():
    """Run the Streamlit app."""
    st.set_page_config(
        page_title="plotlyMol3D Viewer",
        page_icon="⚛",
        layout="wide",
    )

    st.sidebar.title("plotlyMol3D")
    st.sidebar.markdown("**Interactive 3D Molecular Visualization**")

    input_method = st.sidebar.radio(
        "Input Method",
        ["SMILES", "MOL File", "XYZ File", "Cube File", "Sample Molecules"],
        index=0,
    )

    mode = st.sidebar.selectbox(
        "Visualization Mode",
        ["ball+stick", "ball", "vdw", "stick"],
        index=0,
        help="ball+stick: atoms and bonds | ball: atoms only | vdw: space-filling | stick: thin atoms",
    )

    with st.sidebar.expander("Lighting Settings", expanded=False):
        # Initialize defaults only if not present (fixes Session State warning)
        if "ambient" not in st.session_state:
            st.session_state["ambient"] = 0.2
        if "diffuse" not in st.session_state:
            st.session_state["diffuse"] = 0.8
        if "specular" not in st.session_state:
            st.session_state["specular"] = 0.3
        if "roughness" not in st.session_state:
            st.session_state["roughness"] = 0.5
        if "fresnel" not in st.session_state:
            st.session_state["fresnel"] = 0.1

        presets = {}
        if CONFIG_PATH.exists():
            try:
                presets = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                presets = {}

        lighting_presets = presets.get("lighting_presets", {})
        preset_names = sorted(lighting_presets.keys())

        if preset_names:
            selected_preset = st.selectbox(
                "Load preset",
                preset_names,
                key="lighting_preset",
            )

            def apply_preset():
                data = lighting_presets.get(selected_preset, {})
                if data:
                    st.session_state["ambient"] = data.get("ambient", 0.2)
                    st.session_state["diffuse"] = data.get("diffuse", 0.8)
                    st.session_state["specular"] = data.get("specular", 0.3)
                    st.session_state["roughness"] = data.get("roughness", 0.5)
                    st.session_state["fresnel"] = data.get("fresnel", 0.1)

            st.button("Load preset", on_click=apply_preset)
        else:
            st.caption("No saved presets yet.")

        # Use sliders without default value since we're using session state
        ambient = st.slider("Ambient", 0.0, 1.0, key="ambient", step=0.05)
        diffuse = st.slider("Diffuse", 0.0, 1.0, key="diffuse", step=0.05)
        specular = st.slider("Specular", 0.0, 1.0, key="specular", step=0.05)
        roughness = st.slider("Roughness", 0.0, 1.0, key="roughness", step=0.05)
        fresnel = st.slider("Fresnel", 0.0, 1.0, key="fresnel", step=0.05)

        st.markdown("---")
        preset_name = st.text_input(
            "Preset name",
            value="default",
            help="Save the current lighting settings under this name",
        )
        if st.button("Save lighting preset"):
            presets = {}
            if CONFIG_PATH.exists():
                try:
                    presets = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    presets = {}

            lighting_presets = presets.get("lighting_presets", {})
            lighting_presets[preset_name] = {
                "ambient": ambient,
                "diffuse": diffuse,
                "specular": specular,
                "roughness": roughness,
                "fresnel": fresnel,
            }
            presets["lighting_presets"] = lighting_presets
            CONFIG_PATH.write_text(json.dumps(presets, indent=2), encoding="utf-8")
            st.success(f"Saved preset: {preset_name}")

    with st.sidebar.expander("Settings", expanded=False):
        perf_mode = st.selectbox(
            "Mode",
            ["Balanced", "Performance"],
            index=0,
            help="Balanced keeps full quality; Performance reduces render cost.",
        )

    # Vibration Settings Section
    vib_data = None
    vib_mode_number = None
    vib_display_type = None
    vib_amplitude = 1.0
    vib_arrow_color = "#FF0000"
    vib_arrow_scale = 1.0
    vib_heatmap_colorscale = "Reds"
    vib_n_frames = 20

    with st.sidebar.expander("Vibration Settings", expanded=False):
        st.markdown("Upload a vibrational frequency calculation file")

        vib_file = st.file_uploader(
            "Vibration File",
            type=["log", "out", "molden"],
            help="Gaussian .log, ORCA .out, or Molden .molden files",
        )

        if vib_file is not None:
            try:
                with st.spinner("Parsing vibration file..."):
                    vib_bytes = vib_file.read()
                    vib_data = cached_parse_vibrations(vib_bytes, vib_file.name)

                st.success(
                    f"Loaded {len(vib_data.modes)} modes from {vib_data.program.upper()} file"
                )

                # Mode selection dropdown
                mode_options = []
                for mode in vib_data.modes:
                    freq_str = f"{mode.frequency:.1f} cm⁻¹"
                    if mode.is_imaginary:
                        freq_str += " (imaginary)"

                    if mode.ir_intensity is not None:
                        mode_label = f"Mode {mode.mode_number}: {freq_str} (IR: {mode.ir_intensity:.1f})"
                    else:
                        mode_label = f"Mode {mode.mode_number}: {freq_str}"

                    mode_options.append(mode_label)

                selected_mode = st.selectbox(
                    "Select Vibrational Mode",
                    options=range(len(mode_options)),
                    format_func=lambda i: mode_options[i],
                    index=0,
                )
                vib_mode_number = vib_data.modes[selected_mode].mode_number

                # Display type selection
                vib_display_type = st.radio(
                    "Display Type",
                    ["Static arrows", "Animation", "Heatmap", "Arrows + Heatmap"],
                    index=0,
                    help="How to visualize the vibrational mode",
                )

                st.markdown("---")
                st.markdown("**Visualization Parameters**")

                # Common parameters
                vib_amplitude = st.slider(
                    "Amplitude",
                    0.1,
                    5.0,
                    1.0,
                    0.1,
                    help="Displacement amplitude multiplier",
                )

                # Arrow-specific parameters
                if vib_display_type in ["Static arrows", "Arrows + Heatmap"]:
                    vib_arrow_color = st.color_picker(
                        "Arrow Color",
                        value="#FF0000",
                        help="Color for displacement arrows",
                    )
                    vib_arrow_scale = st.slider(
                        "Arrow Size",
                        0.1,
                        3.0,
                        1.0,
                        0.1,
                        help="Visual scale for arrow size",
                    )

                # Heatmap-specific parameters
                if vib_display_type in ["Heatmap", "Arrows + Heatmap"]:
                    vib_heatmap_colorscale = st.selectbox(
                        "Heatmap Colorscale",
                        ["Reds", "Blues", "Viridis", "Plasma", "Hot", "YlOrRd"],
                        index=0,
                        help="Color scheme for displacement magnitude",
                    )

                # Animation-specific parameters
                if vib_display_type == "Animation":
                    vib_n_frames = st.slider(
                        "Animation Frames",
                        5,
                        50,
                        20,
                        5,
                        help="Number of frames (more = smoother but slower)",
                    )

            except Exception as e:
                st.error(f"Error parsing vibration file: {e}")
                vib_data = None

    resolution = st.sidebar.slider(
        "Resolution", 8, 64, 32, 8, help="Higher = smoother spheres"
    )
    resolution_used = 16 if perf_mode == "Performance" else resolution
    if perf_mode == "Performance":
        st.sidebar.caption("Performance mode uses lower resolution for faster UI.")

    st.title("Molecule Viewer")

    rdkitmol = None
    show_orbitals = False
    cube_path = None

    if input_method == "SMILES":
        st.markdown("### Enter SMILES String")

        if "smiles_input" not in st.session_state:
            st.session_state.smiles_input = "c1ccccc1"

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
            st.session_state["smiles_input"] = choice[0]
            st.session_state["random_molecule_name"] = choice[1]

        col1, col2 = st.columns([3, 1])
        with col1:
            smiles_input = st.text_input(
                "SMILES",
                placeholder="Enter a SMILES string (e.g., CCO for ethanol)",
                label_visibility="collapsed",
                key="smiles_input",
            )
        with col2:
            st.button(
                "Random",
                help="Try a random molecule",
                on_click=set_random_smiles,
            )

        if "random_molecule_name" in st.session_state:
            st.toast(f"Selected: {st.session_state.random_molecule_name}")
            del st.session_state.random_molecule_name

        if st.session_state.smiles_input:
            try:
                with st.spinner("Parsing SMILES..."):
                    rdkitmol = cached_smiles_to_mol(st.session_state.smiles_input)
                st.success(f"Parsed: {Chem.MolToSmiles(rdkitmol)}")
            except Exception as e:
                st.error(f"Invalid SMILES: {e}")

    elif input_method == "MOL File":
        st.markdown("### Upload MOL File")

        uploaded_file = st.file_uploader("Choose a .mol file", type=["mol", "sdf"])

        if uploaded_file is not None:
            try:
                with st.spinner("Loading MOL file..."):
                    mol_content = uploaded_file.read().decode("utf-8")
                    rdkitmol = cached_molblock_to_mol(mol_content)
                if rdkitmol is None:
                    st.error("Could not parse MOL file")
                else:
                    st.success(f"Loaded: {uploaded_file.name}")
            except Exception as e:
                st.error(f"Error reading file: {e}")

    elif input_method == "XYZ File":
        st.markdown("### Upload XYZ File")

        col1, col2 = st.columns([3, 1])
        with col1:
            uploaded_file = st.file_uploader("Choose a .xyz file", type=["xyz"])
        with col2:
            charge = st.number_input(
                "Molecular Charge", value=0, min_value=-5, max_value=5
            )

        if uploaded_file is not None:
            try:
                with st.spinner("Parsing XYZ file..."):
                    xyz_content = uploaded_file.read().decode("utf-8")
                    rdkitmol = cached_xyzblock_to_mol(xyz_content, charge=charge)
                st.success(f"Loaded: {uploaded_file.name}")
            except Exception as e:
                st.error(f"Error: {e}")
                st.info(
                    "Tip: XYZ bond detection can be tricky. Try specifying the correct charge, or use a MOL file instead."
                )

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
                with st.spinner("Processing cube file..."):
                    cube_bytes = uploaded_file.read()
                    with tempfile.NamedTemporaryFile(
                        delete=False, suffix=".cube"
                    ) as tmp:
                        tmp.write(cube_bytes)
                        cube_path = tmp.name

                    if show_molecule:
                        xyzblock, cube_charge = cached_cube_bytes_to_xyzblock(
                            cube_bytes
                        )
                        rdkitmol = cached_xyzblock_to_mol(xyzblock, charge=cube_charge)

                st.success(f"Loaded: {uploaded_file.name}")
            except Exception as e:
                st.error(f"Error: {e}")

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
            with st.spinner("Generating 3D structure..."):
                rdkitmol = cached_smiles_to_mol(smiles)
        except Exception as e:
            st.error(f"Error: {e}")

    if rdkitmol is not None:
        display_molecule_info(rdkitmol)

        # Handle vibration animation separately (creates new figure)
        if vib_data is not None and vib_display_type == "Animation":
            try:
                with st.spinner("Creating vibration animation..."):
                    fig = create_vibration_animation(
                        vib_data=vib_data,
                        mode_number=vib_mode_number,
                        mol=rdkitmol,
                        amplitude=vib_amplitude,
                        n_frames=vib_n_frames,
                        mode=mode,
                        resolution=resolution_used,
                    )
                    # Apply lighting to animation frames (optional)
                    fig = format_lighting(
                        fig,
                        ambient=ambient,
                        diffuse=diffuse,
                        specular=specular,
                        roughness=roughness,
                        fresnel=fresnel,
                    )
            except Exception as e:
                st.error(f"Error creating animation: {e}")
                # Fall back to regular figure
                with st.spinner("Rendering 3D visualization..."):
                    import pickle

                    mol_pkl = pickle.dumps(rdkitmol)
                    fig = cached_figure_from_mol_pickle(
                        mol_pkl,
                        mode,
                        resolution_used,
                        ambient,
                        diffuse,
                        specular,
                        roughness,
                        fresnel,
                    )
        else:
            # Regular figure with optional vibration overlays
            with st.spinner("Rendering 3D visualization..."):
                import pickle

                mol_pkl = pickle.dumps(rdkitmol)
                fig = cached_figure_from_mol_pickle(
                    mol_pkl,
                    mode,
                    resolution_used,
                    ambient,
                    diffuse,
                    specular,
                    roughness,
                    fresnel,
                )

            # Apply vibrations if loaded (non-animation modes)
            if vib_data is not None and vib_mode_number is not None:
                try:
                    with st.spinner("Adding vibration visualization..."):
                        # Map display type to internal format
                        display_map = {
                            "Static arrows": "arrows",
                            "Heatmap": "heatmap",
                            "Arrows + Heatmap": "both",
                        }
                        display_type_internal = display_map.get(vib_display_type)

                        if display_type_internal:
                            fig = add_vibrations_to_figure(
                                fig=fig,
                                vib_data=vib_data,
                                mode_number=vib_mode_number,
                                display_type=display_type_internal,
                                amplitude=vib_amplitude,
                                arrow_scale=vib_arrow_scale,
                                arrow_color=vib_arrow_color,
                                heatmap_colorscale=vib_heatmap_colorscale,
                                show_colorbar=True,
                            )
                except Exception as e:
                    st.warning(f"Could not add vibrations: {e}")

        if show_orbitals and cube_path is not None:
            try:
                with st.spinner("Rendering orbitals..."):
                    draw_cube_orbitals(
                        fig, cube_path, orbital_opacity, [pos_color, neg_color]
                    )
            except Exception as e:
                st.warning(f"Could not render orbitals: {e}")

        st.plotly_chart(fig, use_container_width=True)

        with st.sidebar.expander("Save Image", expanded=False):
            preset = st.selectbox(
                "Preset",
                ["Small (800x600)", "HD (1280x720)", "Large (1920x1080)"],
                index=1,
            )
            fmt = st.selectbox("Format", ["png", "svg"], index=0)

            preset_sizes = {
                "Small (800x600)": (800, 600),
                "HD (1280x720)": (1280, 720),
                "Large (1920x1080)": (1920, 1080),
            }
            width, height = preset_sizes[preset]

            file_name = st.text_input(
                "File name",
                value=f"plotlymol3d_{width}x{height}.{fmt}",
            )

            try:
                with st.spinner("Preparing image..."):
                    import pickle

                    mol_pkl = pickle.dumps(rdkitmol)
                    image_bytes = cached_image_bytes(
                        mol_pkl,
                        mode,
                        resolution_used,
                        ambient,
                        diffuse,
                        specular,
                        roughness,
                        fresnel,
                        width,
                        height,
                        fmt,
                    )

                st.download_button(
                    "Download image",
                    data=image_bytes,
                    file_name=file_name,
                    mime=f"image/{fmt}",
                )
            except Exception as e:
                st.error(
                    f"Image export failed: {e}. If this is a missing dependency, install kaleido."
                )

        if cube_path and os.path.exists(cube_path):
            os.unlink(cube_path)

    elif input_method not in ["Sample Molecules"]:
        st.info("Enter a molecule above to visualize it")

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
    st.sidebar.caption(
        "plotlyMol3D v0.1.0 | [GitHub](https://github.com/jonathanschultzNU/plotlyMol)"
    )


if __name__ == "__main__":
    main()
