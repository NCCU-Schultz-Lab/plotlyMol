# plotlyMol - Context for Claude

**Last Updated:** 2026-02-03
**Version:** 0.1.0
**Python Package:** plotlymol
**Primary Module:** plotlymol3d

## Project Overview

plotlyMol is a Python package for creating **interactive 3D molecular visualizations** using Plotly. It enables chemists, researchers, and students to visualize molecular structures with various rendering modes and includes support for orbital visualization from quantum chemistry calculations.

### Core Capabilities

- **3D Molecular Visualization**: Ball-and-stick, stick-only, and van der Waals (VDW) representations
- **Multiple Input Formats**: SMILES strings, XYZ files, MOL/SDF files, PDB files, and Gaussian cube files
- **SMILES-to-3D**: Automatic 3D coordinate generation from SMILES via RDKit
- **Orbital Visualization**: Isosurface rendering from cube files using marching cubes algorithm
- **Vibrational Mode Visualization**: 🆕 Visualize molecular vibrations from quantum chemistry calculations
  - Support for Gaussian (.log), ORCA (.out), and Molden (.molden) formats
  - Three visualization modes: static displacement arrows, animated vibrations, heatmap coloring
  - Interactive controls for mode selection and parameters
- **Bond Order Display**: Visual differentiation of single, double, triple, and aromatic bonds
- **Interactive GUI**: Streamlit-based web interface for visual exploration
- **Export Options**: Interactive HTML, static PNG (via Kaleido)

### Key Features

- Uses Plotly's `Mesh3d` traces for performant 3D rendering
- Fibonacci sphere algorithm for smooth atom spheres
- Cylinder mesh generation for bonds
- Customizable lighting (ambient, diffuse, specular, roughness)
- Aromatic bonds displayed with dashed rendering
- Hover tooltips showing atom/bond information

---

## Repository Structure

```
plotlyMol/
├── src/
│   └── plotlymol3d/               # Main package (note: package name vs module name)
│       ├── __init__.py            # Exports: draw_3D_rep, draw_3D_mol, vibration functions, etc.
│       ├── plotlyMol3D.py         # Core visualization module (~800 lines)
│       ├── atomProperties.py      # Atomic data (CPK colors, VDW radii, symbols, symbol_to_number mapping)
│       ├── cube.py                # Marching cubes for orbital isosurfaces
│       ├── vibrations.py          # 🆕 Vibrational mode visualization (~1000 lines)
│       ├── app.py                 # Streamlit GUI application (with vibration settings)
│       ├── test.py                # Legacy test script (use pytest instead)
│       ├── Cube_to_Blender v3.py  # Legacy Blender export (excluded from linting)
│       └── *.{xyz,mol,pdb,cube}   # Sample molecular data files
│
├── examples/
│   ├── demo_visualizations.py     # Demo script showing various features
│   └── gui_app.py                 # Streamlit GUI (can be run directly)
│
├── tests/
│   ├── __init__.py
│   ├── conftest.py                # pytest fixtures (sample molecules + vibration files)
│   ├── test_input_processing.py   # Tests for file parsing and conversion
│   ├── test_visualization.py      # Tests for rendering functions
│   ├── test_vibrations.py         # 🆕 Tests for vibration parsers and visualization
│   └── fixtures/                  # 🆕 Test data directory
│       ├── water_gaussian.log     # Sample Gaussian frequency calculation
│       ├── water_orca.out         # Sample ORCA frequency calculation
│       └── water.molden           # Sample Molden format file
│
├── docs/
│   └── ROADMAP.md                 # Detailed development roadmap
│
├── .github/
│   └── workflows/
│       ├── test.yml               # CI: pytest on 3 OSes, 4 Python versions
│       └── lint.yml               # CI: black, ruff, mypy
│
├── pyproject.toml                 # Modern Python packaging config
├── requirements.txt               # All dependencies (runtime + dev)
├── README.md                      # User-facing documentation (with vibration examples)
├── CHANGELOG.md                   # Version history
├── LICENSE                        # MIT License
├── .gitignore                     # Comprehensive Python gitignore
├── .pre-commit-config.yaml        # Local pre-commit hooks
├── launch_app.{bat,vbs}          # Windows scripts to launch Streamlit GUI
└── stop_app.{bat,vbs}            # Windows scripts to stop Streamlit GUI
```

### Important Note on Naming

- **Package name (PyPI/pip):** `plotlymol` (no "3d")
- **Module name (import):** `plotlymol3d` (with "3d")
- Usage: `from plotlymol3d import draw_3D_rep`

---

## Core Modules

### 1. plotlyMol3D.py - Main Visualization Engine

**Location:** [src/plotlymol3d/plotlyMol3D.py](src/plotlymol3d/plotlyMol3D.py)

This is the heart of the package (~800 lines). Key components:

#### Data Classes
```python
@dataclass
class Atom:
    atom_id: int
    atom_number: int
    atom_symbol: str
    atom_xyz: List[float]
    atom_vdw: float

@dataclass
class Bond:
    a1_id: int
    a2_id: int
    a1_number: int
    a2_number: int
    a1_xyz: List[float]
    a2_xyz: List[float]
    a1_vdw: float
    a2_vdw: float
    bond_order: float  # 1.0=single, 2.0=double, 3.0=triple, 1.5=aromatic
```

#### Input Processing Functions
- `smiles_to_rdkitmol(smiles)` - Parse SMILES, add H, generate 3D coords
- `xyzfile_to_xyzblock(xyzfile)` - Read XYZ file to text block
- `xyzblock_to_rdkitmol(xyzblock, charge)` - Convert XYZ to RDKit Mol (uses rdDetermineBonds)
- `cubefile_to_xyzblock(cubefile)` - Extract coords from Gaussian cube file
- `molfile_to_rdkitmol(molfile)` - Read MOL/SDF file
- `pdbfile_to_rdkitmol(pdbfile)` - Read PDB file

#### Molecule Processing
- `rdkitmol_to_atoms_bonds_lists(mol)` - Extract Atom and Bond dataclasses from RDKit Mol

#### 3D Mesh Generation
- `fibonacci_sphere(samples, radius)` - Generate sphere vertices using Fibonacci spiral
- `cylinder_mesh(start, end, radius, resolution)` - Generate cylinder for bond
- `dashed_cylinder_mesh(...)` - Generate dashed cylinder for aromatic bonds

#### Rendering Functions
- `make_atom_mesh_trace(atom, resolution, mode)` - Create Mesh3d for atom
- `make_bond_mesh_trace(bond, resolution, mode)` - Create Mesh3d for bond(s)
- `draw_3D_mol(mol, mode, resolution)` - Generate all traces from RDKit Mol
- `draw_3D_rep(...)` - **Main entry point** - handles all input types

#### Formatting
- `format_figure(fig, bgcolor, title)` - Apply scene/layout settings
- `format_lighting(fig, ambient, diffuse, specular, roughness)` - Customize lighting

#### Visualization Modes
- `"ball+stick"` - Full-size atoms + bonds
- `"stick"` - Small atoms + bonds
- `"vdw"` - Van der Waals spheres only (no bonds)

### 2. atomProperties.py - Atomic Data

**Location:** [src/plotlymol3d/atomProperties.py](src/plotlymol3d/atomProperties.py)

Contains dictionaries mapping atomic numbers to:
- CPK color schemes (standard element colors)
- Van der Waals radii
- Element symbols

### 3. cube.py - Orbital Visualization

**Location:** [src/plotlymol3d/cube.py](src/plotlymol3d/cube.py)

Implements marching cubes algorithm for generating isosurfaces from volumetric data (electron density, molecular orbitals).

**Key Function:**
- `draw_cube_orbitals(cubefile, isovalue, colors, opacity)` - Generate orbital mesh

**Note:** This file has relaxed linting rules due to its mathematical complexity and legacy code style.

### 4. app.py - Streamlit GUI

**Location:** [src/plotlymol3d/app.py](src/plotlymol3d/app.py)

Interactive web-based GUI for visual testing and demonstrations. Features:
- Multiple input methods (SMILES, file upload, random molecules)
- Real-time parameter adjustment (lighting, mode, resolution)
- Orbital visualization from cube files
- 🆕 **Vibration Settings** expandable section with file upload and visualization controls
- Configuration persistence to `.plotlymol3d_config.json`
- Caching for performance (`@st.cache_resource`)

**Run with:** `streamlit run src/plotlymol3d/app.py` or `streamlit run examples/gui_app.py`

### 5. vibrations.py - Vibrational Mode Visualization 🆕

**Location:** [src/plotlymol3d/vibrations.py](src/plotlymol3d/vibrations.py)

Comprehensive module for visualizing molecular vibrations from quantum chemistry calculations (~1000 lines).

#### Data Classes
```python
@dataclass
class VibrationalMode:
    mode_number: int              # 1-based index
    frequency: float              # cm⁻¹ (negative if imaginary)
    ir_intensity: Optional[float] # km/mol (if available)
    displacement_vectors: np.ndarray  # (n_atoms, 3) Cartesian displacements
    is_imaginary: bool            # True for imaginary frequencies

@dataclass
class VibrationalData:
    coordinates: np.ndarray       # (n_atoms, 3) in Angstroms
    atomic_numbers: List[int]     # Atomic numbers
    modes: List[VibrationalMode]  # All vibrational modes
    source_file: str              # Original filename
    program: str                  # "gaussian", "orca", or "molden"

    def get_mode(self, mode_number: int) -> Optional[VibrationalMode]
    def get_displacement_magnitudes(self, mode_number: int) -> np.ndarray
```

#### Parser Functions
- `parse_gaussian_vibrations(filepath)` - Parse Gaussian .log files
  - Extracts coordinates from last "Standard orientation"
  - Parses "Harmonic frequencies" section (3 modes per block)
  - Extracts frequencies, IR intensities, and displacement vectors
- `parse_orca_vibrations(filepath)` - Parse ORCA .out files
  - Extracts "CARTESIAN COORDINATES (ANGSTROEM)"
  - Parses "VIBRATIONAL FREQUENCIES" section
  - Filters first 6 translation/rotation modes
  - Extracts displacement vectors from "NORMAL MODES" (6 modes per block)
- `parse_molden_vibrations(filepath)` - Parse Molden .molden files
  - Parses [Atoms], [FREQ], [INT], [FR-NORM-COORD] sections
  - Handles Angs/AU unit conversion
  - Well-structured format with clear section markers
- `parse_vibrations(filepath)` - Auto-detect format and route to appropriate parser

#### Visualization Functions
- `create_displacement_arrows(vib_data, mode_number, amplitude, ...)` - Generate 3D arrow traces
  - Uses Plotly Cone traces for vector field visualization
  - Filters small displacements by threshold
  - Custom hover info with mode and displacement data
- `create_vibration_animation(vib_data, mode_number, mol, n_frames, ...)` - Create animated vibration
  - Sinusoidal motion: coords(t) = coords_eq + A·sin(2πt)·displacement
  - Plotly frames with play/pause controls and slider
  - Typical usage: 20-30 frames for smooth motion
- `create_heatmap_colored_figure(fig, vib_data, mode_number, colorscale, ...)` - Color by displacement
  - Modifies existing atom traces with displacement magnitude coloring
  - Normalizes magnitudes to 0-1 range
  - Applies Plotly colorscale (Reds, Blues, Viridis, etc.)
- `add_vibrations_to_figure(fig, vib_data, mode_number, display_type, ...)` - Main integration function
  - Entry point following `draw_cube_orbitals` pattern
  - Handles "arrows", "heatmap", or "both" display types
  - Updates figure title with mode info

#### Key Features
- **Format Support**: Gaussian, ORCA, Molden with auto-detection
- **Three Visualization Modes**: Arrows, Animation, Heatmap
- **Performance**: Caching with `@st.cache_resource` in Streamlit
- **Error Handling**: Informative messages for parsing failures
- **Testing**: Comprehensive test suite with 21 tests covering all parsers and visualization modes

---

## Architecture & Design Decisions

### 1. RDKit as Core Dependency

RDKit provides:
- SMILES parsing (`Chem.MolFromSmiles`)
- 3D coordinate generation (`AllChem.EmbedMolecule`, `UFFOptimizeMolecule`)
- Bond perception from XYZ (`rdDetermineBonds.DetermineBonds`)
- File format readers (MOL, PDB)
- Molecular graph operations

**Known Limitation:** XYZ bond perception can fail for:
- Charged molecules without proper charge specification
- Complex functional groups (NITRO groups, some transition metals)
- Molecules with unusual bonding

**Workaround:** Use MOL/SDF files instead when RDKit bond perception fails.

### 2. Plotly Mesh3d Traces

Instead of using scatter plots with markers, the package uses `Mesh3d` traces:
- **Performance:** More efficient for complex molecules
- **Visual Quality:** Smooth surfaces with proper lighting
- **Interactivity:** Native Plotly zoom, rotate, pan

Each atom and bond is a separate `Mesh3d` trace, enabling:
- Individual coloring (half-bonds colored by atom)
- Custom hover info per atom/bond
- Show/hide capabilities (future feature)

### 3. Fibonacci Sphere Algorithm

Atoms are rendered using Fibonacci sphere tessellation:
- Distributes points evenly on sphere surface
- Creates triangular mesh with good uniformity
- Adjustable resolution (default: 32 points)

### 4. Bond Order Visualization

- **Single bonds:** One cylinder
- **Double bonds:** Two parallel cylinders with offset
- **Triple bonds:** Three cylinders arranged triangularly
- **Aromatic bonds:** One solid + one dashed cylinder

### 5. Half-Bond Coloring

Bonds are split at the midpoint and colored by the atom at each end:
- Visual clarity for heteroatom bonds (C-O, C-N)
- Each half is a separate mesh

**Future Enhancement:** Weight split point by VDW radii instead of 50/50.

---

## Development Setup

### Installation from Source

```bash
# Clone repository
git clone https://github.com/NCCU-Schultz-Lab/plotlyMol.git
cd plotlyMol

# Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in editable mode with all dependencies
pip install -r requirements.txt
pip install -e .
```

### Dependencies

**Runtime (required):**
- `plotly>=5.0.0` - Interactive plotting
- `numpy>=1.20.0` - Numerical arrays
- `rdkit>=2022.3.1` - Chemistry toolkit
- `kaleido>=0.2.1` - Static image export

**Development (optional):**
- `pytest>=7.0.0` - Testing framework
- `pytest-cov>=4.0.0` - Coverage reporting
- `black>=23.0.0` - Code formatter
- `ruff>=0.1.0` - Fast linter
- `flake8>=6.0.0` - Additional linting
- `mypy>=1.0.0` - Type checker
- `pre-commit>=3.0.0` - Git hooks
- `streamlit>=1.30.0` - GUI framework

### Python Version Support

- **Minimum:** Python 3.8
- **Tested:** Python 3.9, 3.10, 3.11, 3.12
- **CI Matrix:** Tests on Ubuntu, macOS, Windows

---

## Testing Strategy

### Test Suite Structure

**Location:** [tests/](tests/)

- `conftest.py` - Fixtures providing sample molecules and vibration files
- `test_input_processing.py` - Input format parsers
- `test_visualization.py` - Rendering functions
- `test_vibrations.py` - 🆕 Vibration parsers and visualization (21 tests)

### Current Coverage

- **Overall:** Significantly improved with vibration tests
- **plotlyMol3D.py:** ~73%
- **vibrations.py:** 🆕 ~95% (comprehensive test coverage)
- **cube.py:** Low (marching cubes needs more tests)
- **Target:** >60% achieved with vibration module

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=plotlymol3d --cov-report=term-missing

# Run specific test file
pytest tests/test_input_processing.py

# Run specific test
pytest tests/test_visualization.py::test_make_atom_mesh_trace
```

### Test Status

**Passing:** 47/47 tests ✅

Key tests:
- ✅ SMILES parsing
- ✅ XYZ file reading
- ✅ XYZ to molecule conversion (with known charge limitation)
- ✅ Cube file parsing
- ✅ Atom/bond extraction
- ✅ Mesh generation (sphere, cylinder)
- ✅ Bond order handling
- 🆕 ✅ Gaussian vibration parsing (3 tests)
- 🆕 ✅ ORCA vibration parsing (2 tests)
- 🆕 ✅ Molden vibration parsing (2 tests)
- 🆕 ✅ Displacement arrow generation (3 tests)
- 🆕 ✅ Vibration animation (2 tests)
- 🆕 ✅ Heatmap coloring (2 tests)
- 🆕 ✅ Dataclass methods (3 tests)
- 🆕 ✅ Integration with draw_3D_rep (4 tests)

---

## CI/CD Pipeline

### GitHub Actions Workflows

**Location:** [.github/workflows/](.github/workflows/)

#### 1. test.yml - Automated Testing

**Triggers:** Push/PR to `main` or `develop` branches

**Matrix:**
- **OS:** Ubuntu, macOS, Windows
- **Python:** 3.9, 3.10, 3.11, 3.12
- **Total:** 10 combinations (some excluded to reduce CI time)

**Steps:**
1. Checkout code
2. Set up Python with pip caching
3. Install dependencies from `requirements.txt`
4. Install package in editable mode
5. Run pytest with coverage
6. Upload coverage to Codecov (Ubuntu 3.11 only)

#### 2. lint.yml - Code Quality

**Triggers:** Push/PR to any branch

**Checks:**
1. **Black** - Code formatting (line length: 88)
2. **Ruff** - Fast linting (replaces flake8)
3. **Mypy** - Static type checking

**Files Excluded from Linting:**
- `cube.py` - Mathematical code with relaxed rules
- `Cube_to_Blender*.py` - Legacy code
- `test.py` - Old test file

### Pre-commit Hooks

**Location:** [.pre-commit-config.yaml](.pre-commit-config.yaml)

Local hooks run before each commit:
- Trailing whitespace removal
- End-of-file fixer
- YAML validation
- Black formatting
- Ruff linting

**Setup:**
```bash
pre-commit install
```

---

## Known Issues & Limitations

### Current Issues

1. **XYZ Bond Perception Failures**
   - **Issue:** RDKit's `rdDetermineBonds` fails on some charged molecules
   - **Example:** NITRO groups, charged transition metal complexes
   - **Workaround:** Use MOL/SDF files with explicit bond information
   - **Status:** Phase 5 roadmap item

2. **Half-Bond Positioning**
   - **Issue:** Bonds split at 50/50 midpoint, not weighted by atom size
   - **Impact:** Visual accuracy could be improved
   - **Solution:** Implement VDW-weighted midpoint calculation
   - **Status:** Phase 5 roadmap item

3. **Single Molecule Limitation**
   - **Issue:** Can only visualize one molecule per call
   - **Requested:** Handle lists of SMILES or structures
   - **Use Case:** Comparing multiple conformers or related structures
   - **Status:** Phase 5 roadmap item

4. **Orbital Integration**
   - **Issue:** Orbital drawing exists but needs better integration
   - **Status:** Basic cube file support works, needs polish

5. **Test Coverage**
   - **Issue:** Only ~27% overall coverage (cube.py has minimal tests)
   - **Target:** >60%
   - **Status:** Ongoing improvement

### Platform-Specific Notes

**Windows:**
- Uses `.bat` and `.vbs` scripts for Streamlit GUI launching
- File paths use backslashes (handled by pathlib)

**macOS/Linux:**
- Standard Unix paths work correctly
- Shell scripts would be more natural (could add in future)

---

## Future Roadmap

**Full details in:** [docs/ROADMAP.md](docs/ROADMAP.md)

### Phase Status
- ✅ **Phase 1:** Project Foundation (Complete)
- ✅ **Phase 2:** Code Quality (Complete)
- ✅ **Phase 3:** Testing & CI/CD (Complete)
- ✅ **Phase 4:** Documentation (Complete)
- 🔄 **Phase 5:** Feature Development (In Progress)
  - ✅ Vibrational Mode Visualization (Complete)
- ⏳ **Phase 6:** Advanced Features (Pending)
- ⏳ **Phase 7:** Community & Distribution (Pending)

### Recent Milestones (Phase 4-5) ✅

1. **Documentation** ✅
   - Comprehensive README with vibration examples
   - CHANGELOG updates
   - CLAUDE.md context updates

2. **Vibrational Mode Visualization** ✅ COMPLETED
   - Three file format parsers (Gaussian, ORCA, Molden)
   - Three visualization modes (arrows, animation, heatmap)
   - Streamlit UI integration with interactive controls
   - Comprehensive test suite (21 tests, 95% coverage)

### Near-Term Priorities (Phase 5-6)

1. **Vibration Feature Enhancements**
   - 🎯 Test with real quantum chemistry data
   - 🎯 Create example Jupyter notebooks
   - 🎯 IR spectrum viewer with clickable peaks
   - 🎯 Export animations as GIF/MP4
   - 🎯 Additional formats (ADF, Q-Chem, NWChem)
   - 🎯 Raman intensity support
   - 🎯 Transition state reaction coordinate visualization

2. **Documentation Expansion**
   - API documentation (Sphinx or MkDocs)
   - Tutorial notebooks showcasing vibration features
   - Example gallery with quantum chemistry workflows
   - CONTRIBUTING.md

3. **Core Feature Enhancements**
   - VDW-weighted bond splitting
   - Partial charge coloring (Gasteiger)
   - Enhanced hover tooltips
   - SDF multi-molecule support
   - 2D structure rendering (ChemDraw-like)

4. **RDKit Integration Expansion**
   - More input formats (MOL2, InChI)
   - Substructure highlighting
   - Conformer generation/viewing
   - Property calculations

### Long-Term Goals (Phase 6-7)

- Molecular dynamics trajectory visualization (4D animations)
- Interactive measurement tools (distances, angles, dihedrals)
- Performance optimization for large molecules (>1000 atoms)
- Side-by-side mode comparison for vibrational analysis
- PyPI publication
- conda-forge package
- Community building (GitHub Discussions, website)

---

## Code Style & Standards

### Python Style

**Formatter:** Black (line length: 88)
**Linter:** Ruff (replaces flake8)
**Type Checker:** Mypy

### Docstring Style

**Format:** Google Style

Example:
```python
def draw_3D_rep(smiles: Optional[str] = None, ...) -> go.Figure:
    """Create interactive 3D molecular visualization.

    Args:
        smiles: SMILES string for molecule.
        mode: Visualization mode ("ball+stick", "stick", "vdw").
        resolution: Sphere tessellation resolution (default: 32).

    Returns:
        Plotly Figure object with 3D molecular visualization.

    Raises:
        ValueError: If no input is provided or input is invalid.

    Example:
        >>> fig = draw_3D_rep(smiles="CCO", mode="ball+stick")
        >>> fig.show()
    """
```

### Type Hints

**Required:** All functions should have type hints
**Style:** Use `typing` module for complex types

```python
from typing import List, Optional, Tuple, Union

def example_func(
    param1: str,
    param2: Optional[int] = None,
    param3: List[float] = None
) -> Tuple[int, str]:
    ...
```

### Import Organization

**Order (enforced by ruff):**
1. Standard library
2. Third-party packages
3. Local modules

```python
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from rdkit import Chem

from .atomProperties import *
from .cube import *
```

---

## Key Files Reference

### Configuration Files

| File | Purpose |
|------|---------|
| `pyproject.toml` | Modern Python packaging, tool configs (pytest, black, ruff, mypy, coverage) |
| `requirements.txt` | All dependencies (runtime + dev consolidated) |
| `.pre-commit-config.yaml` | Pre-commit hook configurations |
| `.gitignore` | Comprehensive Python exclusions |

### Package Files

| File | Lines | Purpose |
|------|-------|---------|
| `plotlyMol3D.py` | ~800 | Main visualization engine |
| `atomProperties.py` | ~290 | Atomic data (colors, radii, symbols, symbol_to_number mapping) |
| `cube.py` | ~400 | Marching cubes for orbitals |
| `vibrations.py` | 🆕 ~1030 | Vibrational mode visualization (parsers + visualization) |
| `app.py` | ~527 | Streamlit GUI application (with vibration settings) |
| `__init__.py` | ~14 | Package exports (includes vibration functions) |

### Documentation Files

| File | Purpose |
|------|---------|
| `README.md` | User-facing documentation, quick start |
| `ROADMAP.md` | Detailed development plan (25+ phases) |
| `CHANGELOG.md` | Version history |
| `LICENSE` | MIT License |
| `CLAUDE.md` | This file - AI assistant context |

### Sample Data Files

Located in `src/plotlymol3d/`:
- `cube.xyz` - Sample XYZ molecular coordinates
- `cube.mol` - Sample MOL file
- `cube.pdb` - Sample PDB file
- `anto_occ_1-min2.cube` - Sample Gaussian cube file with orbital data

---

## Git Workflow

### Branch Strategy

- **main** - Stable releases
- **develop** - Development branch (if used)
- Feature branches for new work

### Current Git Status

```
M  pyproject.toml         # Modified
M  requirements.txt        # Modified
M  src/plotlymol3d/app.py  # Modified
```

### Commit Message Style

Based on repository history, commits follow conventional format:
- `Add <feature>` - New functionality
- `Update <file/feature>` - Modifications
- `Fix <issue>` - Bug fixes
- `Merge pull request #N` - PR merges

Example:
```
Add Streamlit app module and Windows scripts

- Created app.py for interactive GUI
- Added launch/stop batch files for Windows
- Configured persistent settings

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

### Repository URLs

- **GitHub:** https://github.com/NCCU-Schultz-Lab/plotlyMol
- **Issues:** https://github.com/NCCU-Schultz-Lab/plotlyMol/issues
- **CI Badges:** Tests and Lint workflows visible in README

---

## Quick Start for AI Assistants

### Common Tasks

**1. Adding a new visualization feature:**
- Modify `plotlyMol3D.py`
- Add corresponding test in `tests/test_visualization.py`
- Update docstrings
- Run `black` and `ruff --fix`
- Run `pytest`

**2. Supporting a new input format:**
- Add parser function in `plotlyMol3D.py` (section: Input Processing Functions)
- Add test in `tests/test_input_processing.py`
- Update `draw_3D_rep()` to accept new parameter
- Update README with example

**3. Fixing a bug:**
- Write a failing test first
- Implement fix
- Verify test passes
- Check coverage didn't decrease

**4. Improving documentation:**
- Update docstrings in code
- Update README.md for user-facing changes
- Update ROADMAP.md for completed items
- Update CHANGELOG.md for version history

### Testing Commands

```bash
# Format code
black src/plotlymol3d tests

# Lint code
ruff check src/plotlymol3d tests
ruff check --fix src/plotlymol3d tests  # Auto-fix

# Type check
mypy src/plotlymol3d

# Run tests
pytest -v

# Run tests with coverage
pytest --cov=plotlymol3d --cov-report=term-missing

# Run GUI for visual testing
streamlit run src/plotlymol3d/app.py
```

---

## Design Patterns & Best Practices

### 1. Separation of Concerns

- **Input Processing** → `*_to_rdkitmol()` functions
- **Data Extraction** → `rdkitmol_to_atoms_bonds_lists()`
- **Mesh Generation** → `fibonacci_sphere()`, `cylinder_mesh()`
- **Rendering** → `make_*_trace()` functions
- **Formatting** → `format_figure()`, `format_lighting()`

### 2. Dataclasses for Clarity

Use `@dataclass` instead of dicts or tuples for atom/bond data:
- Type safety
- Named fields
- Default values
- Self-documenting

### 3. Optional Parameters with Sensible Defaults

Most functions use optional parameters with defaults:
```python
def draw_3D_rep(
    smiles: Optional[str] = None,
    mode: str = "ball+stick",
    resolution: int = DEFAULT_RESOLUTION,
    ambient: float = 0.5,
    ...
)
```

### 4. Error Handling Philosophy

- Use `ValueError` for invalid inputs
- Let RDKit exceptions propagate (better error messages)
- Document potential failures in docstrings

### 5. Caching in GUI

Streamlit app uses `@st.cache_resource` for expensive operations:
- SMILES parsing
- XYZ conversion
- File reading

Prevents recomputation on every widget interaction.

---

## Performance Considerations

### 1. Mesh Resolution

- **Default:** 32 samples per sphere
- **Trade-off:** Higher resolution = smoother but slower
- **Range:** 16-64 typically sufficient

### 2. Large Molecules

- Each atom/bond is a separate trace
- Molecules with 100+ atoms may render slowly
- Future optimization: Combine meshes, WebGL improvements

### 3. Marching Cubes (Orbitals)

- Most computationally expensive operation
- Large cube files (>100³ grid) take several seconds
- Potential improvements: Numba JIT, parallel processing

---

## Debugging Tips

### 1. RDKit Issues

If bond perception fails:
```python
# Check RDKit Mol object
mol = xyzblock_to_rdkitmol(xyzblock, charge=0)
if mol is None:
    print("RDKit failed to create molecule")
else:
    print(f"Atoms: {mol.GetNumAtoms()}, Bonds: {mol.GetNumBonds()}")
```

### 2. Visualization Problems

If rendering looks wrong:
- Check `mode` parameter ("ball+stick", "stick", "vdw")
- Adjust `resolution` (try 16 for fast preview)
- Check lighting parameters (too high ambient = washed out)

### 3. Import Errors

If `from plotlymol3d import *` fails:
- Verify package installed: `pip list | grep plotlymol`
- Check working directory
- Try explicit import: `from plotlymol3d.plotlyMol3D import draw_3D_rep`

---

## Authors & License

**Authors:**
- Jonathan Schultz (jonathanschultzNU@users.noreply.github.com)
- Benjamin Lear

**License:** MIT (see [LICENSE](LICENSE) file)

**Repository:** https://github.com/NCCU-Schultz-Lab/plotlyMol

---

## Additional Resources

- **RDKit Documentation:** https://www.rdkit.org/docs/
- **Plotly Python Documentation:** https://plotly.com/python/
- **Plotly 3D Mesh:** https://plotly.com/python/3d-mesh/
- **SMILES Tutorial:** https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
- **Gaussian Cube Format:** http://paulbourke.net/dataformats/cube/
- **Marching Cubes:** https://en.wikipedia.org/wiki/Marching_cubes

---

**End of Context Document**

This document should be updated when:
- Major architectural changes occur
- New modules are added
- Development phase changes (see ROADMAP.md)
- API significantly changes
- New dependencies are added
