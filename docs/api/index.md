# API Reference

Complete API documentation for plotlyMol.

## Overview

plotlyMol's API is organized into several modules:

- **[Core Functions](core.md)** - Main visualization functions (`draw_3D_rep`, `draw_3D_mol`, etc.)
- **[Atom Properties](properties.md)** - Atomic data (colors, radii, symbols)
- **[Orbital Rendering](orbitals.md)** - Molecular orbital visualization from cube files

## Quick Reference

### Main Function

The primary entry point is `draw_3D_rep()`:

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
    smiles=None,              # SMILES string
    xyzfile=None,             # Path to XYZ file
    xyzblock=None,            # XYZ text block
    molfile=None,             # Path to MOL/SDF file
    pdbfile=None,             # Path to PDB file
    cubefile=None,            # Path to Gaussian cube file
    mol=None,                 # RDKit Mol object
    mode="ball+stick",        # Visualization mode
    resolution=32,            # Sphere tessellation
    ambient=0.5,              # Ambient lighting
    diffuse=0.8,              # Diffuse reflection
    specular=0.5,             # Specular highlights
    roughness=0.5,            # Surface roughness
    bgcolor="white",          # Background color
    title=None,               # Figure title
    cubedraw=None,            # "orbitals" for orbital viz
    orbital_isovalue=0.02,    # Isosurface threshold
    orbital_colors=["darkorange", "darkblue"],
    orbital_opacity=0.25      # Orbital transparency
)
```

### Input Methods

plotlyMol accepts multiple input formats:

=== "SMILES"
    ```python
    fig = draw_3D_rep(smiles="CCO")
    ```

=== "XYZ File"
    ```python
    fig = draw_3D_rep(xyzfile="molecule.xyz")
    ```

=== "MOL/SDF"
    ```python
    fig = draw_3D_rep(molfile="molecule.mol")
    ```

=== "PDB"
    ```python
    fig = draw_3D_rep(pdbfile="protein.pdb")
    ```

=== "RDKit Mol"
    ```python
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")
    fig = draw_3D_rep(mol=mol)
    ```

### Visualization Modes

Three rendering modes are available:

| Mode | Description | Use Case |
|------|-------------|----------|
| `"ball+stick"` | Full-size atoms + bonds | General visualization |
| `"stick"` | Small atoms + bonds | Complex molecules |
| `"vdw"` | Van der Waals spheres only | Space-filling models |

### Helper Functions

#### Formatting

```python
from plotlymol3d import format_figure, format_lighting

# Apply scene formatting
format_figure(fig, bgcolor="black", title="My Molecule")

# Adjust lighting
format_lighting(fig, ambient=0.2, diffuse=0.9)
```

#### Input Processing

```python
from plotlymol3d import (
    smiles_to_rdkitmol,
    xyzfile_to_xyzblock,
    xyzblock_to_rdkitmol,
    molfile_to_rdkitmol,
    pdbfile_to_rdkitmol
)

# Convert SMILES to RDKit Mol
mol = smiles_to_rdkitmol("CCO")

# Read XYZ file
xyz_text = xyzfile_to_xyzblock("molecule.xyz")
mol = xyzblock_to_rdkitmol(xyz_text, charge=0)
```

#### Mesh Generation

```python
from plotlymol3d import (
    fibonacci_sphere,
    cylinder_mesh,
    make_atom_mesh_trace,
    make_bond_mesh_trace
)

# Generate sphere vertices
vertices = fibonacci_sphere(samples=32, radius=1.0)

# Create cylinder mesh
vertices, faces = cylinder_mesh(
    start=[0, 0, 0],
    end=[1, 1, 1],
    radius=0.15,
    resolution=16
)
```

## Data Structures

### Atom Dataclass

```python
from dataclasses import dataclass
from typing import List

@dataclass
class Atom:
    atom_id: int          # 0-indexed atom identifier
    atom_number: int      # Atomic number (6=C, 8=O, etc.)
    atom_symbol: str      # Element symbol ("C", "O", etc.)
    atom_xyz: List[float] # [x, y, z] coordinates
    atom_vdw: float       # Van der Waals radius
```

### Bond Dataclass

```python
@dataclass
class Bond:
    a1_id: int            # First atom ID
    a2_id: int            # Second atom ID
    a1_number: int        # First atom number
    a2_number: int        # Second atom number
    a1_xyz: List[float]   # First atom [x, y, z]
    a2_xyz: List[float]   # Second atom [x, y, z]
    a1_vdw: float         # First atom VDW radius
    a2_vdw: float         # Second atom VDW radius
    bond_order: float     # 1.0, 1.5, 2.0, 3.0
```

## Constants

```python
from plotlymol3d import DEFAULT_RESOLUTION, DEFAULT_BOND_RADIUS

print(DEFAULT_RESOLUTION)   # 32
print(DEFAULT_BOND_RADIUS)  # 0.15
```

## Type Hints

plotlyMol uses comprehensive type hints:

```python
from typing import Optional, List, Tuple
import plotly.graph_objects as go
from rdkit.Chem import Mol

def draw_3D_rep(
    smiles: Optional[str] = None,
    mode: str = "ball+stick",
    resolution: int = 32,
    # ...
) -> go.Figure:
    """Returns a Plotly Figure object"""
```

## Error Handling

```python
from plotlymol3d import draw_3D_rep

try:
    fig = draw_3D_rep(smiles="invalid_smiles")
except ValueError as e:
    print(f"Invalid input: {e}")
```

Common exceptions:
- `ValueError` - Invalid input or parameters
- `FileNotFoundError` - File not found
- RDKit exceptions - Molecule parsing errors

## Next Steps

- [Core Functions Reference](core.md) - Detailed function documentation
- [Atom Properties](properties.md) - Atomic data reference
- [Orbital Rendering](orbitals.md) - Cube file visualization
