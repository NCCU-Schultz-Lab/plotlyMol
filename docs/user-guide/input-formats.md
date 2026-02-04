# Input Formats

plotlyMol supports multiple molecular file formats and input methods.

## SMILES

**Simplified Molecular Input Line Entry System**

### Overview

SMILES is the most convenient input format - a text-based representation of molecular structure.

### Usage

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(smiles="CCO")
```

### Common SMILES Examples

| Molecule | SMILES | Description |
|----------|--------|-------------|
| Methane | `C` | Single carbon |
| Ethane | `CC` | Two carbons |
| Ethanol | `CCO` | Carbon chain with OH |
| Benzene | `c1ccccc1` | Aromatic ring |
| Acetic acid | `CC(=O)O` | Double bond to oxygen |
| Aspirin | `CC(=O)OC1=CC=CC=C1C(=O)O` | Complex molecule |

### SMILES Syntax

- **Atoms**: Element symbols (C, N, O, etc.)
- **Bonds**: Single (-), double (=), triple (#)
- **Aromatics**: Lowercase letters (c for aromatic carbon)
- **Rings**: Numbers indicate ring closures
- **Branches**: Parentheses for branching

### 3D Generation

plotlyMol automatically:
1. Parses SMILES with RDKit
2. Adds explicit hydrogens
3. Generates 3D coordinates (UFF force field)
4. Optimizes geometry

### Limitations

- May not preserve specific conformations
- Some complex structures may fail 3D generation
- Stereochemistry may not be perfect

For details on SMILES syntax, see [DaylightSMILES](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html).

## XYZ Files

**Cartesian Coordinate Format**

### Format Specification

```
N_atoms
Comment line
Element1 X1 Y1 Z1
Element2 X2 Y2 Z2
...
```

### Example

```
3
Water molecule
O   0.000   0.000   0.000
H   0.757   0.586   0.000
H  -0.757   0.586   0.000
```

### Usage

```python
fig = draw_3D_rep(xyzfile="molecule.xyz")
```

### Advantages

- Exact atomic coordinates
- Simple, human-readable format
- Widely supported

### Limitations

- **No bond information** - RDKit must infer bonds
- May fail for:
  - Charged molecules
  - Unusual bonding
  - Transition metal complexes

### Troubleshooting

If XYZ fails:

```python
# Try specifying charge
from plotlymol3d import xyzblock_to_rdkitmol

with open("molecule.xyz") as f:
    xyz_text = f.read()

mol = xyzblock_to_rdkitmol(xyz_text, charge=-1)
```

Or use MOL/SDF format instead.

## MOL/SDF Files

**MDL Molfile / Structure Data File**

### Overview

MOL files include explicit bond information, making them more reliable than XYZ.

### Usage

```python
fig = draw_3D_rep(molfile="molecule.mol")
```

### Advantages

- **Explicit bonds** - No inference needed
- Preserves bond orders
- Handles charged species
- Includes stereochemistry

### When to Use

Prefer MOL/SDF over XYZ when:
- XYZ bond perception fails
- Molecule has formal charges
- Stereochemistry matters
- Complex functional groups present

### SDF Format

SDF (Structure Data File) can contain multiple molecules:

```python
# Currently plotlyMol handles single molecules
# Multi-molecule support coming in Phase 5
fig = draw_3D_rep(molfile="single_molecule.sdf")
```

## PDB Files

**Protein Data Bank Format**

### Overview

PDB format is designed for biomolecules from X-ray crystallography or NMR.

### Usage

```python
fig = draw_3D_rep(pdbfile="protein.pdb")
```

### Use Cases

- Protein structures
- Peptides
- Nucleic acids
- Protein-ligand complexes
- Crystallographic structures

### Limitations

- Large files may render slowly
- Use `mode="stick"` for better performance

## Cube Files

**Gaussian Cube Format**

### Overview

Cube files contain volumetric data from quantum chemistry calculations.

### Usage

```python
fig = draw_3D_rep(
    molfile="molecule.mol",  # Structure
    cubefile="orbital.cube",  # Volumetric data
    cubedraw="orbitals"
)
```

See [Orbital Visualization Guide](orbitals.md) for details.

## RDKit Mol Objects

### Overview

For integration with RDKit workflows.

### Usage

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Create Mol object
mol = Chem.MolFromSmiles("CCO")
AllChem.EmbedMolecule(mol)

# Visualize
fig = draw_3D_rep(mol=mol)
```

### Use Cases

- Custom molecule manipulation
- Conformer generation
- Property calculations
- Integration with existing code

## Format Comparison

| Format | Bonds | Charges | Stereo | 3D | Best For |
|--------|-------|---------|--------|----|----|
| SMILES | Implicit | Yes | Yes | Auto-generated | Quick visualization |
| XYZ | Inferred | Limited | No | Yes | Optimized geometries |
| MOL/SDF | Explicit | Yes | Yes | Yes | Complex molecules |
| PDB | Explicit | Yes | Yes | Yes | Biomolecules |
| Cube | N/A | N/A | N/A | Yes | Orbitals + structure |

## File Reading Functions

### Low-Level API

For advanced users:

```python
from plotlymol3d import (
    smiles_to_rdkitmol,
    xyzfile_to_xyzblock,
    xyzblock_to_rdkitmol,
    molfile_to_rdkitmol,
    pdbfile_to_rdkitmol
)

# Read SMILES
mol = smiles_to_rdkitmol("CCO")

# Read XYZ
xyz_text = xyzfile_to_xyzblock("molecule.xyz")
mol = xyzblock_to_rdkitmol(xyz_text, charge=0)

# Read MOL
mol = molfile_to_rdkitmol("molecule.mol")

# Read PDB
mol = pdbfile_to_rdkitmol("protein.pdb")
```

## Best Practices

### Choosing Format

1. **Start with SMILES** if possible
2. **Use XYZ** for optimized geometries
3. **Switch to MOL** if XYZ fails
4. **Use PDB** for biomolecules

### Troubleshooting

**"Failed to perceive bonds":**
- Switch from XYZ to MOL format
- Verify atomic charges
- Check geometry quality

**"Invalid SMILES":**
- Check syntax
- Verify element symbols
- Try online SMILES validator

**"RDKit error":**
- Molecule may be too complex
- Try different input format
- Check for unusual bonding

## Next Steps

- [Visualization Modes](visualization-modes.md)
- [Customization](customization.md)
- [API Reference](../api/core.md)
