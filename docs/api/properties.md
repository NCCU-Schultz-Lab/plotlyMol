# Atom Properties

Reference for atomic properties used in molecular visualization.

## Overview

The `atomProperties` module provides atomic data used for rendering molecules:

- **CPK Colors**: Standard element coloring scheme
- **VDW Radii**: Van der Waals radii in Angstroms
- **Element Symbols**: Atomic number to symbol mapping

## Module Reference

::: plotlymol3d.atomProperties
    options:
      show_source: true
      members:
        - CPK_COLORS
        - VDW_RADII
        - SYMBOLS

## CPK Color Scheme

Standard CPK (Corey-Pauling-Koltun) colors for elements:

| Element | Symbol | Color | RGB |
|---------|--------|-------|-----|
| Hydrogen | H | White | (255, 255, 255) |
| Carbon | C | Gray | (144, 144, 144) |
| Nitrogen | N | Blue | (48, 80, 248) |
| Oxygen | O | Red | (255, 13, 13) |
| Fluorine | F | Green | (144, 224, 80) |
| Phosphorus | P | Orange | (255, 128, 0) |
| Sulfur | S | Yellow | (255, 255, 48) |
| Chlorine | Cl | Green | (31, 240, 31) |

**Usage:**
```python
from plotlymol3d.atomProperties import CPK_COLORS

# Get color for carbon (atomic number 6)
carbon_color = CPK_COLORS[6]  # "rgb(144, 144, 144)"
```

## Van der Waals Radii

Atomic radii in Angstroms:

| Element | Symbol | VDW Radius (Å) |
|---------|--------|----------------|
| H | Hydrogen | 1.20 |
| C | Carbon | 1.70 |
| N | Nitrogen | 1.55 |
| O | Oxygen | 1.52 |
| F | Fluorine | 1.47 |
| P | Phosphorus | 1.80 |
| S | Sulfur | 1.80 |
| Cl | Chlorine | 1.75 |

**Usage:**
```python
from plotlymol3d.atomProperties import VDW_RADII

# Get VDW radius for oxygen (atomic number 8)
oxygen_radius = VDW_RADII[8]  # 1.52
```

## Element Symbols

Mapping of atomic numbers to element symbols:

**Usage:**
```python
from plotlymol3d.atomProperties import SYMBOLS

# Get symbol for atomic number 6
symbol = SYMBOLS[6]  # "C"
```

## Customization

You can override default properties when needed:

```python
from plotlymol3d import draw_3D_rep
from plotlymol3d.atomProperties import CPK_COLORS

# Modify colors (not recommended - creates side effects)
# Better to use Plotly's figure update methods

fig = draw_3D_rep(smiles="CCO")

# Update trace colors after creation
for trace in fig.data:
    if hasattr(trace, 'color'):
        # Customize individual trace colors
        pass
```

## Complete Element List

The module includes data for elements 1-118 (Hydrogen through Oganesson).

**Common Elements in Organic Chemistry:**

| Z | Symbol | Name | VDW (Å) | Color |
|---|--------|------|---------|-------|
| 1 | H | Hydrogen | 1.20 | White |
| 6 | C | Carbon | 1.70 | Gray |
| 7 | N | Nitrogen | 1.55 | Blue |
| 8 | O | Oxygen | 1.52 | Red |
| 9 | F | Fluorine | 1.47 | Green |
| 15 | P | Phosphorus | 1.80 | Orange |
| 16 | S | Sulfur | 1.80 | Yellow |
| 17 | Cl | Chlorine | 1.75 | Green |
| 35 | Br | Bromine | 1.85 | Brown |
| 53 | I | Iodine | 1.98 | Purple |

**Transition Metals:**

| Z | Symbol | Name | VDW (Å) | Color |
|---|--------|------|---------|-------|
| 26 | Fe | Iron | 2.00 | Orange |
| 29 | Cu | Copper | 1.40 | Orange |
| 30 | Zn | Zinc | 1.39 | Gray |

## Default Values

If an element is not found in the dictionaries, default values are used:

- **Color**: `"rgb(255, 20, 147)"` (deep pink) - indicates missing data
- **VDW Radius**: 1.70 Å (carbon-like default)
- **Symbol**: "X" (unknown)

## References

- CPK colors: Corey, R. B.; Pauling, L. (1953)
- VDW radii: Bondi, A. (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68 (3): 441–451.
