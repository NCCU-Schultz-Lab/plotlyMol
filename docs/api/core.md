# Core Functions

Auto-generated API documentation for plotlyMol's core visualization functions.

## Main Visualization Function

::: plotlymol3d.plotlyMol3D.draw_3D_rep
    options:
      show_source: true
      heading_level: 3

## Drawing from RDKit Mol

::: plotlymol3d.plotlyMol3D.draw_3D_mol
    options:
      show_source: true
      heading_level: 3

## Input Processing Functions

### SMILES Processing

::: plotlymol3d.plotlyMol3D.smiles_to_rdkitmol
    options:
      show_source: true
      heading_level: 4

### XYZ File Processing

::: plotlymol3d.plotlyMol3D.xyzfile_to_xyzblock
    options:
      show_source: true
      heading_level: 4

::: plotlymol3d.plotlyMol3D.xyzblock_to_rdkitmol
    options:
      show_source: true
      heading_level: 4

### Cube File Processing

::: plotlymol3d.plotlyMol3D.cubefile_to_xyzblock
    options:
      show_source: true
      heading_level: 4

## Molecule Processing

::: plotlymol3d.plotlyMol3D.rdkitmol_to_atoms_bonds_lists
    options:
      show_source: true
      heading_level: 3

## Mesh Generation Functions

### Sphere Tessellation

::: plotlymol3d.plotlyMol3D.make_fibonacci_sphere
    options:
      show_source: true
      heading_level: 4

### Cylinder Meshes

::: plotlymol3d.plotlyMol3D.generate_cylinder_mesh_rectangles
    options:
      show_source: true
      heading_level: 4

## Trace Creation Functions

### Atom Traces

::: plotlymol3d.plotlyMol3D.make_atom_mesh_trace
    options:
      show_source: true
      heading_level: 4

::: plotlymol3d.plotlyMol3D.draw_atoms
    options:
      show_source: true
      heading_level: 4

### Bond Traces

::: plotlymol3d.plotlyMol3D.make_bond_mesh_trace
    options:
      show_source: true
      heading_level: 4

::: plotlymol3d.plotlyMol3D.draw_bonds
    options:
      show_source: true
      heading_level: 4

## Formatting Functions

### Figure Formatting

::: plotlymol3d.plotlyMol3D.format_figure
    options:
      show_source: true
      heading_level: 4

### Lighting Control

::: plotlymol3d.plotlyMol3D.format_lighting
    options:
      show_source: true
      heading_level: 4

## Data Classes

### Atom

```python
@dataclass
class Atom:
    """Represents an atom with position and properties.

    Attributes:
        atom_id: 0-indexed unique identifier
        atom_number: Atomic number (e.g., 6 for carbon)
        atom_symbol: Element symbol (e.g., "C")
        atom_xyz: [x, y, z] coordinates in Angstroms
        atom_vdw: Van der Waals radius in Angstroms
    """
    atom_id: int
    atom_number: int
    atom_symbol: str
    atom_xyz: List[float]
    atom_vdw: float
```

### Bond

```python
@dataclass
class Bond:
    """Represents a bond between two atoms.

    Attributes:
        a1_id, a2_id: Atom IDs
        a1_number, a2_number: Atomic numbers
        a1_xyz, a2_xyz: Atomic coordinates
        a1_vdw, a2_vdw: Van der Waals radii
        bond_order: 1.0 (single), 1.5 (aromatic), 2.0 (double), 3.0 (triple)
    """
    a1_id: int
    a2_id: int
    a1_number: int
    a2_number: int
    a1_xyz: List[float]
    a2_xyz: List[float]
    a1_vdw: float
    a2_vdw: float
    bond_order: float
```
