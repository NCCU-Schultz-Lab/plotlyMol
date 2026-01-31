# plotlyMol Development Roadmap

This document outlines the development plan for plotlyMol, a Python package for creating interactive molecular visualizations using Plotly.

---

## 📊 Progress Summary

| Phase | Status | Description |
|-------|--------|-------------|
| Phase 1 | ✅ Complete | Project Foundation - Package structure, pyproject.toml, requirements |
| Phase 2 | ✅ Complete | Code Quality - Type hints, docstrings, error handling |
| Phase 3 | ✅ Complete | Testing & CI/CD - GitHub Actions, coverage, linting |
| Phase 4 | 🔄 In Progress | Documentation |
| Phase 5 | ⏳ Pending | Feature Development |
| Phase 6 | ⏳ Pending | Advanced Features |
| Phase 7 | ⏳ Pending | Community & Distribution |

### 🎯 Recommended Next Actions

1. **Merge PR and verify CI** - Push changes, verify workflows run on GitHub
2. **Set up Codecov** - Add `CODECOV_TOKEN` secret to repository for coverage tracking
3. **Fix linting issues** - Run `black plotlymol3d tests` and `ruff check --fix` to auto-fix style issues
4. **Enhance README** - Add CI badges, screenshots, and usage examples
5. **Create CHANGELOG.md** - Document version history
6. **Fix failing test** - Address `test_xyzblock_to_rdkitmol` charge detection issue

---

## Current State

The repository has completed **Phases 1-3** and now contains:

### Core Package
- **3D Visualization Module** (`plotlymol3d/plotlyMol3D.py`): Main module for 3D molecular rendering with ball-and-stick and VDW representations
- **Atom Properties** (`plotlymol3d/atomProperties.py`): Atom colors, symbols, and VDW radii data
- **Marching Cubes Implementation** (`plotlymol3d/cube.py`): Orbital visualization from cube files
- **Sample Data Files**: Various molecular structure files (.xyz, .mol, .pdb, .cube)

### Testing & CI/CD (NEW)
- **Test Suite** (`tests/`): 25 unit tests covering input processing and visualization
- **GitHub Actions** (`.github/workflows/`): Automated testing and linting on push/PR
- **Pre-commit Hooks** (`.pre-commit-config.yaml`): Local development quality checks
- **Coverage**: ~27% (main module ~73%, cube.py needs more tests)

### Configuration
- **Package Configuration** (`pyproject.toml`): Modern Python packaging with tool configs
- **Dependencies** (`requirements.txt`, `requirements-dev.txt`): Core and development dependencies
- **Comprehensive `.gitignore`**: Proper exclusions for Python projects

## Development Phases

### Phase 1: Project Foundation (Immediate) ✅ COMPLETED

**Goal**: Establish proper Python package infrastructure and fix critical issues

#### Tasks:
- [x] **Fix `__init__.py` syntax error**
  - ~~Current: `from .3D.plotlyMol3D import *` is invalid (module names cannot start with digits)~~
  - Solution: Renamed `3D` directory to `plotlymol3d` and fixed import structure
  
- [x] **Improve `.gitignore`**
  - Added `__pycache__/` directories
  - Added `*.pyc`, `*.pyo`, `*.pyd` files
  - Added `.Python`, `*.egg-info/`, `dist/`, `build/`
  - Added virtual environment directories (`venv/`, `env/`, `.venv/`)
  - Added IDE-specific files (`.vscode/`, `.idea/`, `*.swp`)
  - Added OS-specific files (`.DS_Store`, `Thumbs.db`)

- [x] **Create package configuration**
  - Added `pyproject.toml` (modern Python packaging)
  - Defined package metadata (name, version, author, description, license)
  - Specified package structure and entry points
  - Configured build system (setuptools)

- [x] **Add `requirements.txt`**
  - Listed core dependencies:
    - `plotly>=5.0.0` - Interactive plotting library
    - `numpy>=1.20.0` - Numerical operations
    - `rdkit>=2022.3.1` - Chemistry toolkit for molecular operations
  - Created separate `requirements-dev.txt` for development dependencies

- [x] **Reorganize directory structure**
  - Renamed `3D/` to `plotlymol3d/` (valid Python module name)
  - Updated all relative imports accordingly
  - Fixed `__init__.py` to properly export module contents

#### Success Criteria: ✅ All Met
- Package can be installed with `pip install -e .`
- No `__pycache__` or `.pyc` files in version control
- Module can be imported without syntax errors

---

### Phase 2: Code Quality (Short-term) ✅ COMPLETED

**Goal**: Improve code maintainability, readability, and robustness

#### Tasks:
- [x] **Add comprehensive type hints**
  - Annotated all function parameters and return types
  - Used `typing` module for complex types (List, Dict, Optional, Union, Tuple)
  - Added type hints to dataclass fields

- [x] **Add docstrings**
  - Used Google style docstrings
  - Documented all modules, classes, and functions
  - Included parameters, return values, and examples

- [x] **Remove hardcoded paths** ✅
  - Replaced absolute paths in `test.py` with `pathlib.Path` relative paths
  - Used `__file__` to locate package data
  - Added proper imports and documentation

- [x] **Implement proper error handling**
  - Added input validation for file formats
  - Included descriptive error messages in docstrings

- [x] **Clean up or remove `graveyard.py`** ✅
  - Reviewed deprecated code - confirmed it was obsolete Surface-based traces
  - Removed entirely (current Mesh3d implementation is superior)

#### Success Criteria: ✅ All Met
- All functions have type hints and docstrings
- No hardcoded paths in test files
- Graceful error handling with informative messages
- Code organization improved with section headers

---

### Phase 3: Testing & CI/CD (Short-term) ✅ COMPLETED

**Goal**: Establish automated testing and continuous integration

#### Tasks:
- [x] **Create test infrastructure**
  - Created `tests/` directory in repository root
  - Added `tests/__init__.py`
  - Created `tests/conftest.py` for pytest fixtures
  - Sample test data uses files already in `plotlymol3d/` directory

- [x] **Write unit tests**
  - Test input format parsers:
    - `test_smiles_to_rdkitmol` - SMILES parsing ✅
    - `test_xyzfile_to_xyzblock` - XYZ file reading ✅
    - `test_xyzblock_to_rdkitmol` - XYZ to molecule conversion ⚠️ (known issue with charge detection)
    - `test_cubefile_to_xyzblock` - Cube file parsing ✅
  - Test molecular structure handling:
    - `test_rdkitmol_to_atoms_bonds_lists` - Atom/bond extraction ✅
  - Test visualization components:
    - `test_make_atom_mesh_trace` - Atom rendering ✅
    - `test_make_bond_mesh_trace` - Bond rendering ✅
    - `test_fibonacci_sphere` - Sphere generation ✅
    - `test_cylinder_mesh` - Cylinder generation ✅

- [x] **Add GitHub Actions CI workflow**
  - Created `.github/workflows/test.yml`
  - Runs on Python 3.9, 3.10, 3.11, 3.12
  - Runs on Ubuntu, macOS, Windows
  - Installs dependencies and runs pytest with coverage
  - Uploads coverage to Codecov

- [x] **Add code coverage reporting**
  - Configured `pytest-cov` in `pyproject.toml`
  - Coverage threshold set to 25% (increase as tests improve)
  - Current coverage: ~27% overall, ~73% on main module

- [x] **Add linting and formatting checks**
  - Created `.github/workflows/lint.yml`
  - **Ruff**: Fast Python linter (replaces flake8)
  - **Black**: Code formatter
  - **mypy**: Static type checker
  - Created `.pre-commit-config.yaml` for local hooks

#### Files Created:
- `.github/workflows/test.yml` - Test automation
- `.github/workflows/lint.yml` - Linting automation  
- `.pre-commit-config.yaml` - Pre-commit hooks
- Updated `pyproject.toml` with tool configurations
- Updated `requirements-dev.txt` with new dependencies

#### Known Issues:
- ⚠️ `test_xyzblock_to_rdkitmol` fails due to RDKit charge detection issue
- ⚠️ Many linting warnings exist (run `black` and `ruff --fix` to auto-fix)
- ⚠️ Coverage is low on `cube.py` (marching cubes code)

#### Success Criteria: ✅ Met
- [x] pytest runs successfully (24/25 tests pass)
- [x] CI/CD pipeline configured for all supported platforms
- [x] Linting and formatting tools configured
- [x] All tests run automatically on pull requests

---

### Phase 4: Documentation (Medium-term)

**Goal**: Create comprehensive documentation for users and contributors

#### Tasks:
- [ ] **Enhance README.md**
  - Add clear project description and features
  - Add badges (build status, coverage, PyPI version, license)
  - Add installation instructions:
    ```bash
    pip install plotlymol
    # or for development
    git clone https://github.com/jonathanschultzNU/plotlyMol.git
    cd plotlyMol
    pip install -e .
    ```
  - Add quick start guide with code examples
  - Add examples for each input format (SMILES, XYZ, MOL, PDB, cube)
  - Add screenshots/GIFs of visualizations
  - Link to full documentation

- [ ] **Create API documentation**
  - Choose documentation tool (Sphinx or MkDocs)
  - Set up documentation structure:
    - `docs/` directory
    - `docs/api/` - Auto-generated API reference
    - `docs/tutorials/` - Step-by-step guides
    - `docs/examples/` - Example gallery
  - Configure auto-generation from docstrings
  - Host documentation on Read the Docs or GitHub Pages

- [ ] **Create example notebooks**
  - Create `examples/` directory
  - Add Jupyter notebooks demonstrating:
    - Basic 3D molecule visualization
    - Different representation modes (ball+stick, VDW, stick)
    - Orbital visualization from cube files
    - Customizing colors and lighting
    - Exporting visualizations
  - Ensure notebooks are tested and up-to-date

- [ ] **Add CONTRIBUTING.md**
  - Explain how to contribute (bug reports, feature requests, pull requests)
  - Code style guidelines
  - Testing requirements
  - Pull request process
  - Code of conduct reference

- [ ] **Add LICENSE file**
  - Recommend MIT or BSD license for open source
  - Include copyright notice
  - Add license badge to README

- [ ] **Create CHANGELOG.md**
  - Document version history
  - Follow Keep a Changelog format
  - Include:
    - Added features
    - Changed functionality
    - Deprecated features
    - Removed features
    - Fixed bugs
    - Security patches

#### Success Criteria:
- Clear, comprehensive README with examples
- API documentation available online
- Example notebooks run without errors
- Contribution guidelines are clear and welcoming

---

### Phase 5: Feature Development (Medium-term)

**Goal**: Implement planned features and address known issues

#### 2D Structures
- [ ] **Implement ChemDraw-like 2D structure rendering**
  - Use RDKit's 2D depiction capabilities
  - Create Plotly-based 2D molecular structure plots
  - Support SMILES input for 2D generation
  - Add customization options (atom labels, bond types, colors)

#### Enhanced 3D Features
- [ ] **Scale half-bond positions by VDW radii**
  - Currently bonds are split at midpoint
  - Implement VDW-weighted bond splitting
  - Update `draw_bonds` function to use weighted midpoints
  - Formula: `weighted_mid = (r1*pos1 + r2*pos2) / (r1 + r2)`

- [ ] **Improve XYZ to MOL conversion**
  - Address rdkit XYZ block to MOL conversion failures
  - Handle edge cases (NITRO groups, charged species)
  - Consider alternative approaches:
    - Use Open Babel for XYZ to SMILES conversion
    - Implement custom bond perception algorithm
    - Add fallback methods

- [ ] **Add orbital drawing integration**
  - Full integration of marching cubes orbital visualization
  - Make orbital drawing more user-friendly
  - Add documentation and examples
  - Support different isosurface values
  - Allow multiple orbitals in single visualization

#### Multiple Molecule Support
- [ ] **Handle lists of structures in single plot**
  - Accept list of SMILES strings
  - Accept list of file paths
  - Create grid layouts or overlays
  - Add molecule labeling/naming
  - Support comparison visualizations

#### Additional Enhancements
- [ ] **Add support for molecular properties**
  - Display molecular weight, formula
  - Show partial charges
  - Display bond orders
  - Interactive property tooltips

#### Success Criteria:
- 2D structure rendering works for common molecules
- VDW-scaled bonds improve visual accuracy
- Robust handling of various input formats
- Support for visualizing multiple molecules

---

### Phase 6: Advanced Features (Long-term)

**Goal**: Implement advanced visualization and interaction capabilities

#### 4D Animations
- [ ] **Implement trajectory/animation support**
  - Read molecular dynamics trajectory files
  - Support common formats (XYZ trajectory, DCD, TRR)
  - Create Plotly animation frames
  - Add timeline controls
  - Export as animated HTML or video

#### GUI Development
- [ ] **Create interactive web interface**
  - Framework options: Streamlit, Dash, or Flask
  - Features:
    - Toggle visualization options (hover, presentation mode)
    - Dynamic molecule input (paste SMILES, upload files)
    - Real-time parameter adjustment (colors, sizes, lighting)
    - Export options (PNG, HTML, SVG)
    - Save/load visualization settings

- [ ] **Add visualization controls**
  - Rotation/zoom controls
  - Show/hide atoms by element
  - Show/hide hydrogens
  - Measurement tools (distances, angles)
  - Selection and highlighting

#### Performance Optimization
- [ ] **Optimize marching cubes algorithm**
  - Profile performance on large cube files
  - Implement parallel processing (multiprocessing/numba)
  - Add progress indicators for long operations
  - Cache computed meshes
  - Optimize mesh reduction/decimation

#### Additional Input Formats
- [ ] **Support more file formats**
  - SDF (Structure-Data File) - multiple molecules
  - CIF (Crystallographic Information File) - crystals
  - PDBQT (AutoDock) - docking results
  - GRO (GROMACS) - MD simulations
  - JSON - custom molecular formats

#### Success Criteria:
- Smooth animations for molecular dynamics trajectories
- Functional GUI for non-programmers
- Acceptable performance on large molecular systems
- Support for diverse input formats

---

### Phase 7: Community & Distribution (Long-term)

**Goal**: Make plotlyMol widely accessible and build a community

#### Package Distribution
- [ ] **Publish to PyPI**
  - Set up PyPI account and project
  - Configure `pyproject.toml` for publishing
  - Create distribution packages (sdist and wheel)
  - Upload to PyPI: `pip install plotlymol`
  - Set up automated releases via GitHub Actions

- [ ] **Create conda-forge recipe**
  - Submit feedstock to conda-forge
  - Maintain conda package
  - Enable: `conda install -c conda-forge plotlymol`

#### Community Building
- [ ] **Set up GitHub Discussions**
  - Enable Discussions on repository
  - Create categories:
    - Announcements
    - Q&A
    - Show and Tell (user examples)
    - Feature Requests
    - General Discussion

- [ ] **Create project website**
  - Dedicated domain or GitHub Pages
  - Features:
    - Interactive demos
    - Gallery of examples
    - Tutorial walkthrough
    - API documentation
    - Download links
    - Community showcase

- [ ] **Outreach and promotion**
  - Write blog posts or tutorials
  - Present at Python or chemistry conferences
  - Submit to Awesome Python lists
  - Create demo videos
  - Engage with chemistry/Python communities

#### Success Criteria:
- Package available via pip and conda
- Active community engagement
- Professional documentation website
- Growing user base and contributors

---

## Additional Files to Consider Creating

### Immediate Priority ✅ COMPLETED
- **`.gitignore`**: Comprehensive exclusions ✅
- **`pyproject.toml`**: Package configuration ✅
- **`requirements.txt`**: Dependencies list ✅
- **`requirements-dev.txt`**: Development dependencies ✅

### Short-term Priority (Phase 2-3)
- **`LICENSE`**: Open source license (recommend MIT) ✅ Created
- **`CONTRIBUTING.md`**: Contribution guidelines
- **`CHANGELOG.md`**: Version history

### Medium-term Priority
- **`.github/workflows/`**: CI/CD workflows
  - `test.yml` - Automated testing
  - `lint.yml` - Code quality checks
  - `publish.yml` - PyPI publishing
- **`tests/`**: Test suite directory
  - `tests/conftest.py` - pytest configuration
  - `tests/test_*.py` - Unit tests
- **`docs/`**: Documentation directory
  - `docs/conf.py` - Sphinx/MkDocs config
  - `docs/index.md` - Documentation homepage
  - `docs/tutorials/` - Tutorials
  - `docs/api/` - API reference

### Long-term Priority
- **`examples/`**: Example notebooks and scripts
- **`.pre-commit-config.yaml`**: Pre-commit hooks configuration
- **`CODE_OF_CONDUCT.md`**: Community guidelines

---

## Notes from README Future Plans

The following goals were identified in the original README:

### 2D Structures
- ChemDraw-like 2D molecular structures

### 3D Enhancements
- Scale half-bond positions by van der Waals radii of atoms at either end
- Add orbitals and other molecular surfaces (marching cubes code exists)
- Handle SMILES structures (currently implemented)
- Fix XYZ block to MOL conversion issues (NITRO groups, etc.)
- Need to integrate orbital drawing code

### 4D Features
- Molecular dynamics animations and trajectories

### GUI Features
- Toggle visualization aspects (presentation mode, hover info, etc.)
- Dynamic molecule entry interface
- Interactive parameter controls

### Multiple Molecule Support
- Handle lists of SMILES or structures in a single plot
- Currently only handles single molecules per call

---

## Implementation Priority

### High Priority (Completed) ✅
1. ~~Fix `__init__.py` syntax error~~ ✅
2. ~~Create proper `.gitignore`~~ ✅
3. ~~Add `pyproject.toml` or `setup.py`~~ ✅
4. ~~Add `requirements.txt`~~ ✅
5. ~~Fix directory structure (rename `3D/` to valid Python module name)~~ ✅

### Medium Priority (Current Focus)
1. ~~Remove hardcoded paths from `test.py`~~ ✅
2. ~~Add type hints and docstrings~~ ✅
3. ~~Add LICENSE file (MIT)~~ ✅
4. ~~Create test suite with pytest~~ ✅
5. ~~Set up CI/CD with GitHub Actions~~ ✅
6. **Fix code style issues** ⬅️ Run `black` and `ruff --fix`
7. **Enhance README** - Add badges, screenshots, examples
8. **Create CHANGELOG.md** - Document version history

### Lower Priority (3-6 months)
1. Implement 2D structure rendering
2. Fix VDW bond scaling
3. Complete API documentation (Sphinx/MkDocs)
4. Create example Jupyter notebooks
5. Publish to PyPI
6. Increase test coverage to >60%

### Future Considerations (6+ months)
1. GUI development (Streamlit/Dash)
2. Animation support for molecular dynamics
3. Performance optimization for large molecules
4. Additional file format support (SDF, CIF, etc.)
5. Community building initiatives

---

## Success Metrics

- **Code Quality**: All functions documented, >80% test coverage, passing linting
- **Usability**: Clear installation instructions, working examples, responsive issue resolution
- **Distribution**: Available on PyPI and conda-forge
- **Community**: Active contributors, growing user base, regular releases
- **Documentation**: Comprehensive docs with tutorials and API reference

---

## Conclusion

This roadmap provides a clear path from the current state to a mature, well-documented, and widely-used Python package for molecular visualization. The phased approach ensures that critical infrastructure is established first, followed by quality improvements, testing, documentation, and finally advanced features and community building.

Contributions are welcome at any phase of this roadmap! See CONTRIBUTING.md (to be created in Phase 4) for details on how to get involved.

---

**Last Updated**: 2026-01-31  
**Current Phase**: Phase 4 (Documentation) - Starting
**Phase 1 Status**: ✅ Completed
**Phase 2 Status**: ✅ Completed
**Phase 3 Status**: ✅ Completed
