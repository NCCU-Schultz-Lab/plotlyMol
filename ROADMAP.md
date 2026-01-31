# plotlyMol Development Roadmap

This document outlines the development plan for plotlyMol, a Python package for creating interactive molecular visualizations using Plotly.

---

## 📊 Progress Summary

| Phase | Status | Description |
|-------|--------|-------------|
| Phase 1 | ✅ Complete | Project Foundation - Package structure, pyproject.toml, requirements |
| Phase 2 | ✅ Complete | Code Quality - Type hints, docstrings, error handling |
| Phase 3 | 🔄 In Progress | Testing & CI/CD |
| Phase 4 | ⏳ Pending | Documentation |
| Phase 5 | ⏳ Pending | Feature Development |
| Phase 6 | ⏳ Pending | Advanced Features |
| Phase 7 | ⏳ Pending | Community & Distribution |

### 🎯 Recommended Next Actions

1. **Set up GitHub Actions CI** - Add automated testing workflow
2. **Add code coverage** - Configure pytest-cov and Codecov
3. **Add linting checks** - Configure flake8/ruff and black
4. **Enhance README** - Add badges, screenshots, and more examples
5. **Create CHANGELOG.md** - Document version history

---

## Current State

The repository has completed **Phase 1** and now contains:
- **3D Visualization Module** (`plotlymol3d/plotlyMol3D.py`): Main module for 3D molecular rendering with ball-and-stick and VDW representations
- **Atom Properties** (`plotlymol3d/atomProperties.py`): Atom colors, symbols, and VDW radii data
- **Marching Cubes Implementation** (`plotlymol3d/cube.py`): Orbital visualization from cube files
- **Test Examples** (`plotlymol3d/test.py`): Basic usage examples (still has hardcoded paths - to fix in Phase 2)
- **Sample Data Files**: Various molecular structure files (.xyz, .mol, .pdb, .cube)
- **Deprecated Code** (`graveyard.py`): Old/unused code (to clean up in Phase 2)
- **Documentation** (`README.md`): Updated with installation and usage instructions
- **Package Configuration** (`pyproject.toml`): Modern Python packaging with metadata
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

### Phase 3: Testing & CI/CD (Short-term) 🔄 IN PROGRESS

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
    - `test_xyzblock_to_rdkitmol` - XYZ to molecule conversion ✅
    - `test_cubefile_to_xyzblock` - Cube file parsing ✅
  - Test molecular structure handling:
    - `test_rdkitmol_to_atoms_bonds_lists` - Atom/bond extraction ✅
  - Test visualization components:
    - `test_make_atom_mesh_trace` - Atom rendering ✅
    - `test_make_bond_mesh_trace` - Bond rendering ✅
    - `test_fibonacci_sphere` - Sphere generation ✅
    - `test_cylinder_mesh` - Cylinder generation ✅

- [ ] **Add GitHub Actions CI workflow**
  - Create `.github/workflows/test.yml` ✅
  - Run tests on:
    - Multiple Python versions (3.9, 3.10, 3.11, 3.12) ✅
    - Multiple operating systems (Ubuntu, macOS, Windows) ✅
  - Install dependencies and run pytest ✅
  - Upload test results as artifacts

- [x] **Add code coverage reporting**
  - Install `pytest-cov` ✅
  - Configure coverage in `pyproject.toml` ✅
  - Add coverage reporting to CI workflow ✅
  - Upload to Codecov ✅
  - Add coverage badge to README

- [x] **Add linting and formatting checks**
  - **Linting**: Add `ruff` for style checking ✅
  - **Formatting**: Add `black` for code formatting ✅
  - **Type checking**: Add `mypy` for static type analysis ✅
  - Create `.github/workflows/lint.yml` for automated checks ✅
  - Add pre-commit hooks ✅

#### Success Criteria:
- `pytest` runs successfully with >80% code coverage
- CI/CD pipeline passes on all supported platforms
- Linting and formatting checks pass
- All tests run automatically on pull requests

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

### High Priority (Start Immediately)
1. ~~Fix `__init__.py` syntax error~~ ✅
2. ~~Create proper `.gitignore`~~ ✅
3. ~~Add `pyproject.toml` or `setup.py`~~ ✅
4. ~~Add `requirements.txt`~~ ✅
5. ~~Fix directory structure (rename `3D/` to valid Python module name)~~ ✅

### Medium Priority (Next Steps - 1-3 months)
1. ~~Remove hardcoded paths from `test.py`~~ ✅
2. ~~Add type hints and docstrings~~ ✅
3. ~~Add LICENSE file (MIT)~~ ✅
4. ~~Create test suite with pytest~~ ✅
5. Set up CI/CD with GitHub Actions ⬅️ **START HERE**
6. Enhance README with examples and badges

### Lower Priority (3-6 months)
1. Implement 2D structure rendering
2. Fix VDW bond scaling
3. Complete API documentation
4. Create example notebooks
5. Publish to PyPI

### Future Considerations (6+ months)
1. GUI development
2. Animation support
3. Performance optimization
4. Additional file format support
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
**Current Phase**: Phase 3 (Testing & CI/CD) - In Progress
**Phase 1 Status**: ✅ Completed
**Phase 2 Status**: ✅ Completed
