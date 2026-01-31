# plotlyMol

[![Tests](https://github.com/NCCU-Schultz-Lab/plotlyMol/actions/workflows/test.yml/badge.svg)](https://github.com/NCCU-Schultz-Lab/plotlyMol/actions/workflows/test.yml)
[![Lint](https://github.com/NCCU-Schultz-Lab/plotlyMol/actions/workflows/lint.yml/badge.svg)](https://github.com/NCCU-Schultz-Lab/plotlyMol/actions/workflows/lint.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Interactive molecular visualizations with Plotly. Supports SMILES, XYZ, MOL/PDB, and cube orbitals.

## Features

- 3D ball-and-stick, stick, and VDW representations
- SMILES-to-3D embedding via RDKit
- XYZ, MOL/SDF (single), and PDB input support
- Cube file orbital isosurfaces
- Streamlit GUI for interactive exploration

## Installation

### From source (recommended for now)

```bash
git clone https://github.com/NCCU-Schultz-Lab/plotlyMol.git
cd plotlyMol
pip install -e .
```

### Development dependencies

```bash
pip install -r requirements.txt
```

## Quick start

```python
from plotlymol3d import draw_3D_rep

# Draw a molecule from SMILES
fig = draw_3D_rep(smiles="CCNCOCSC", mode="ball+stick", ambient=0.1)
fig.show()

# Draw from XYZ file
fig = draw_3D_rep(xyzfile="path/to/file.xyz", mode="ball+stick", ambient=0.1)
fig.show()
```

### Orbitals from cube files

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
	cubefile="path/to/file.cube",
	molfile="path/to/file.mol",
	mode="ball+stick",
	ambient=0.1,
	cubedraw="orbitals",
	orbital_opacity=0.25,
	orbital_colors=["darkorange", "darkblue"],
)
fig.show()
```

## GUI

Launch the Streamlit app for interactive controls:

```bash
streamlit run examples/gui_app.py
```

## Examples

- Demo script: `python examples/demo_visualizations.py`
- Package data includes sample XYZ/MOL/CUBE files under src/plotlymol3d/

## Repository layout

```
plotlyMol/
├─ src/
│  └─ plotlymol3d/        # Library package code + sample data files
├─ examples/              # Demo scripts and GUI app
├─ tests/                 # Pytest suite
├─ docs/                  # Roadmap and documentation assets
├─ pyproject.toml          # Packaging and tooling configuration
├─ requirements.txt        # Consolidated dependencies
└─ README.md
```

## Roadmap

See the current roadmap in [docs/ROADMAP.md](docs/ROADMAP.md).

