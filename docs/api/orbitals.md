# Orbital Visualization

Documentation for molecular orbital visualization from Gaussian cube files.

## Overview

plotlyMol can visualize molecular orbitals and electron density from quantum chemistry calculations stored in Gaussian cube file format.

## Module Reference

::: plotlymol3d.cube.draw_cube_orbitals
    options:
      show_source: true
      heading_level: 3

## Cube File Format

Gaussian cube files contain volumetric data on a 3D grid:

```
Comment line 1
Comment line 2
N_atoms X_origin Y_origin Z_origin
N_x dX_x dX_y dX_z
N_y dY_x dY_y dY_z
N_z dZ_x dZ_y dZ_z
Atom1_number Atom1_charge X1 Y1 Z1
Atom2_number Atom2_charge X2 Y2 Z2
...
Volumetric data (N_x * N_y * N_z values)
```

## Basic Usage

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
    molfile="molecule.mol",
    cubefile="orbital.cube",
    mode="ball+stick",
    cubedraw="orbitals",
    orbital_isovalue=0.02,
    orbital_colors=["red", "blue"],
    orbital_opacity=0.3
)
fig.show()
```

## Parameters

### Isovalue

The `orbital_isovalue` parameter controls the isosurface threshold:

- **Lower values** (0.01-0.015): Larger, more diffuse orbitals
- **Medium values** (0.02-0.03): Balanced visualization (recommended)
- **Higher values** (0.04+): Smaller, more compact orbitals

```python
# Diffuse orbital
fig = draw_3D_rep(
    cubefile="orbital.cube",
    molfile="mol.mol",
    cubedraw="orbitals",
    orbital_isovalue=0.01
)

# Compact orbital
fig = draw_3D_rep(
    cubefile="orbital.cube",
    molfile="mol.mol",
    cubedraw="orbitals",
    orbital_isovalue=0.05
)
```

### Colors

Orbitals have two phases (positive and negative), colored separately:

```python
# Classic red/blue
orbital_colors=["red", "blue"]

# Orange/dark blue (default)
orbital_colors=["darkorange", "darkblue"]

# Custom colors
orbital_colors=["#FF6B6B", "#4ECDC4"]
```

### Opacity

Control orbital transparency:

```python
# More transparent (better for seeing molecule)
orbital_opacity=0.2

# More opaque (emphasize orbital)
orbital_opacity=0.5
```

## Marching Cubes Algorithm

plotlyMol uses the marching cubes algorithm to generate isosurfaces from volumetric data:

1. **Grid Traversal**: Scan through 3D grid cells
2. **Surface Detection**: Identify cells intersecting the isovalue
3. **Vertex Interpolation**: Calculate exact intersection points
4. **Triangulation**: Generate triangular mesh representing surface
5. **Normal Calculation**: Compute surface normals for lighting

The implementation handles:
- Arbitrary grid orientations
- Non-cubic voxels
- Both positive and negative isosurfaces
- Mesh smoothing

## Performance Considerations

### Grid Size

Cube file size affects performance:

| Grid Size | Voxels | Build Time | Quality |
|-----------|--------|------------|---------|
| 50³ | 125K | ~1s | Low |
| 100³ | 1M | ~5s | Medium |
| 150³ | 3.4M | ~15s | High |
| 200³ | 8M | ~40s | Very High |

### Optimization Tips

1. **Use appropriate isovalue**: Higher values = fewer triangles = faster

2. **Reduce opacity**: Lower opacity can hide mesh artifacts

3. **Simplify molecular structure**: Use stick mode for complex molecules

4. **Cache results**: If viewing same orbital multiple times

## Complete Example

```python
from plotlymol3d import draw_3D_rep

# HOMO visualization
fig = draw_3D_rep(
    molfile="benzene.mol",
    cubefile="benzene_HOMO.cube",
    mode="stick",
    resolution=32,
    ambient=0.1,
    bgcolor="white",
    cubedraw="orbitals",
    orbital_isovalue=0.02,
    orbital_colors=["#FF6B6B", "#4ECDC4"],
    orbital_opacity=0.35
)

fig.update_layout(
    title="Benzene HOMO",
    scene=dict(
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
    )
)

fig.write_html("benzene_HOMO.html")
fig.show()
```

## Generating Cube Files

Cube files can be generated from various quantum chemistry packages:

### Gaussian

```
%chk=molecule.chk
#P B3LYP/6-31G(d) Pop=Full

...

(After calculation)
formchk molecule.chk molecule.fchk
cubegen 0 MO=HOMO molecule.fchk homo.cube 80 h
```

### ORCA

```
! B3LYP 6-31G(d)

%output
  Print[P_MOs] 1
end

...

(Use orca_plot to generate cube)
```

### Psi4

```python
import psi4

# Run calculation
energy, wfn = psi4.energy('b3lyp/6-31g(d)', return_wfn=True)

# Generate cube file
psi4.cubeprop(wfn)
```

### Q-Chem

```
$rem
  ...
  make_cube_files true
$end
```

## Troubleshooting

### Orbital Not Visible

- Check isovalue (try lower values like 0.01)
- Verify cube file contains data
- Check orbital opacity

### Blocky Appearance

- Cube file resolution too low
- Need finer grid in QM calculation

### Slow Rendering

- Reduce cube grid size
- Increase isovalue
- Use stick mode for molecule

### Wrong Orbital Shown

- Check cube file contains desired orbital
- Verify file path is correct
- Some QM packages number orbitals differently

## Additional Resources

- [Gaussian Cube Format Specification](http://paulbourke.net/dataformats/cube/)
- [Marching Cubes Algorithm](https://en.wikipedia.org/wiki/Marching_cubes)
- [RDKit Cube File Handling](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDraw2D.html)
