# Basic Usage

Comprehensive guide to using plotlyMol for molecular visualization.

## Installation

See the [Installation Guide](../installation.md) for setup instructions.

## Quick Start

The simplest way to use plotlyMol:

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(smiles="CCO")
fig.show()
```

This creates an interactive 3D visualization of ethanol that you can rotate, zoom, and pan.

## The `draw_3D_rep` Function

The main entry point is `draw_3D_rep()`, which handles all input formats and visualization options.

### Function Signature

```python
def draw_3D_rep(
    smiles: Optional[str] = None,
    xyzfile: Optional[str] = None,
    xyzblock: Optional[str] = None,
    molfile: Optional[str] = None,
    pdbfile: Optional[str] = None,
    cubefile: Optional[str] = None,
    mol: Optional[Chem.Mol] = None,
    mode: str = "ball+stick",
    resolution: int = 32,
    ambient: float = 0.5,
    diffuse: float = 0.8,
    specular: float = 0.5,
    roughness: float = 0.5,
    bgcolor: str = "white",
    title: Optional[str] = None,
    cubedraw: Optional[str] = None,
    orbital_isovalue: float = 0.02,
    orbital_colors: List[str] = ["darkorange", "darkblue"],
    orbital_opacity: float = 0.25
) -> go.Figure:
```

### Parameters

#### Input Sources (Choose One)

| Parameter | Type | Description |
|-----------|------|-------------|
| `smiles` | str | SMILES string (e.g., "CCO") |
| `xyzfile` | str | Path to XYZ file |
| `xyzblock` | str | XYZ text block |
| `molfile` | str | Path to MOL/SDF file |
| `pdbfile` | str | Path to PDB file |
| `cubefile` | str | Path to Gaussian cube file |
| `mol` | Chem.Mol | RDKit Mol object |

#### Visualization Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | str | "ball+stick" | "ball+stick", "stick", or "vdw" |
| `resolution` | int | 32 | Sphere tessellation (16-64) |
| `ambient` | float | 0.5 | Ambient lighting (0.0-1.0) |
| `diffuse` | float | 0.8 | Diffuse reflection (0.0-1.0) |
| `specular` | float | 0.5 | Specular highlights (0.0-1.0) |
| `roughness` | float | 0.5 | Surface roughness (0.0-1.0) |
| `bgcolor` | str | "white" | Background color |
| `title` | str | None | Figure title |

#### Orbital Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubedraw` | str | None | Set to "orbitals" to draw |
| `orbital_isovalue` | float | 0.02 | Isosurface threshold |
| `orbital_colors` | List[str] | ["darkorange", "darkblue"] | [positive, negative] |
| `orbital_opacity` | float | 0.25 | Transparency (0.0-1.0) |

### Return Value

Returns a `plotly.graph_objects.Figure` object that can be:

- Displayed with `.show()`
- Saved as HTML with `.write_html()`
- Exported as image with `.write_image()`
- Further customized with Plotly methods

## Input Methods

### From SMILES

SMILES (Simplified Molecular Input Line Entry System) is the most convenient input:

```python
# Simple molecules
fig = draw_3D_rep(smiles="CCO")  # Ethanol
fig = draw_3D_rep(smiles="c1ccccc1")  # Benzene

# Complex molecules
fig = draw_3D_rep(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
```

**Advantages:**
- Compact text representation
- No file needed
- Automatic 3D generation

**Limitations:**
- RDKit must generate valid 3D coordinates
- May not preserve specific conformations
- Some complex structures may fail

### From XYZ Files

XYZ files contain atomic coordinates:

```python
fig = draw_3D_rep(xyzfile="path/to/molecule.xyz")
```

**XYZ Format:**
```
3
Ethanol
C   0.000   0.000   0.000
C   1.520   0.000   0.000
O   2.020   1.320   0.000
```

**Advantages:**
- Exact coordinates
- Simple format

**Limitations:**
- RDKit must infer bonds
- May fail for charged molecules

### From MOL/SDF Files

MOL (MDL Molfile) and SDF (Structure Data File) include bond information:

```python
fig = draw_3D_rep(molfile="molecule.mol")
```

**Advantages:**
- Explicit bond orders
- Preserves stereochemistry
- Most reliable for complex molecules

**Recommended when:**
- XYZ bond perception fails
- Stereochemistry matters
- Formal charges present

### From PDB Files

Protein Data Bank format:

```python
fig = draw_3D_rep(pdbfile="protein.pdb")
```

Useful for:
- Protein/peptide structures
- Biomolecular systems
- Crystallographic data

### From RDKit Mol Objects

For integration with RDKit workflows:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Create RDKit Mol
mol = Chem.MolFromSmiles("CCO")
AllChem.EmbedMolecule(mol)

# Visualize
fig = draw_3D_rep(mol=mol)
```

**Use cases:**
- Custom molecule manipulation
- Property calculations
- Conformer generation
- Advanced cheminformatics

## Working with Figures

### Display

```python
fig = draw_3D_rep(smiles="CCO")
fig.show()  # Opens in default browser
```

### Save as HTML

```python
fig = draw_3D_rep(smiles="CCO")
fig.write_html("molecule.html")
```

The HTML file is:
- Fully interactive
- Self-contained
- Shareable
- No server needed

### Export as Image

```python
fig = draw_3D_rep(smiles="CCO")

# PNG
fig.write_image("molecule.png", width=800, height=600)

# High-resolution for publication
fig.write_image("molecule.png", width=1600, height=1200, scale=2)

# Other formats
fig.write_image("molecule.pdf")
fig.write_image("molecule.svg")
```

Requires `kaleido` package (installed by default).

### Customize with Plotly

The returned Figure can be modified using Plotly's API:

```python
fig = draw_3D_rep(smiles="CCO")

# Update layout
fig.update_layout(
    title="Ethanol Molecule",
    font=dict(family="Arial", size=14),
    scene=dict(
        bgcolor="lightgray",
        xaxis_title="X (Å)",
        yaxis_title="Y (Å)",
        zaxis_title="Z (Å)"
    )
)

# Update camera
fig.update_layout(
    scene_camera=dict(
        eye=dict(x=1.5, y=1.5, z=1.5),
        center=dict(x=0, y=0, z=0)
    )
)

# Hide axes
fig.update_layout(
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False)
    )
)

fig.show()
```

## Error Handling

### Invalid Input

```python
try:
    fig = draw_3D_rep(smiles="invalid_smiles")
except ValueError as e:
    print(f"Error: {e}")
```

### File Not Found

```python
try:
    fig = draw_3D_rep(xyzfile="nonexistent.xyz")
except FileNotFoundError as e:
    print(f"File not found: {e}")
```

### RDKit Errors

```python
from rdkit import Chem

try:
    mol = Chem.MolFromSmiles("invalid")
    if mol is None:
        print("Failed to parse SMILES")
except Exception as e:
    print(f"RDKit error: {e}")
```

## Best Practices

### Choosing Input Format

1. **Use SMILES for:**
   - Simple molecules
   - Quick visualizations
   - When coordinates don't matter

2. **Use XYZ for:**
   - Optimized geometries
   - Specific conformations
   - Simple molecules

3. **Use MOL/SDF for:**
   - Complex molecules
   - Charged species
   - When XYZ fails

4. **Use PDB for:**
   - Proteins/peptides
   - Biomolecules
   - Crystallographic structures

### Performance Tips

1. **Lower resolution for preview:**
   ```python
   fig = draw_3D_rep(smiles="large_molecule", resolution=16)
   ```

2. **Higher resolution for publication:**
   ```python
   fig = draw_3D_rep(smiles="molecule", resolution=64)
   ```

3. **Use stick mode for large molecules:**
   ```python
   fig = draw_3D_rep(smiles="protein", mode="stick")
   ```

### Quality Settings

**For presentations:**
```python
fig = draw_3D_rep(
    smiles="molecule",
    resolution=48,
    ambient=0.15,
    bgcolor="white"
)
```

**For publication:**
```python
fig = draw_3D_rep(
    smiles="molecule",
    resolution=64,
    ambient=0.1,
    diffuse=0.9,
    specular=0.6,
    roughness=0.3
)
fig.write_image("figure.png", width=1600, height=1200, scale=2)
```

## Next Steps

- [Input Formats](input-formats.md) - Detailed format information
- [Visualization Modes](visualization-modes.md) - Rendering options
- [Customization](customization.md) - Advanced styling
- [API Reference](../api/core.md) - Complete function documentation
