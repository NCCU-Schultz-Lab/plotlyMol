# plotlyMol Development Roadmap

This document outlines the development plan for plotlyMol, a Python package for creating interactive molecular visualizations using Plotly.

## Current State

The repository currently contains:
- **3D Visualization Module** (`3D/plotlyMol3D.py`): Main module for 3D molecular rendering with ball-and-stick and VDW representations
- **Atom Properties** (`3D/atomProperties.py`): Atom colors, symbols, and VDW radii data
- **Marching Cubes Implementation** (`3D/cube.py`): Orbital visualization from cube files
- **Test Examples** (`3D/test.py`): Basic usage examples (with hardcoded paths)
- **Sample Data Files**: Various molecular structure files (.xyz, .mol, .pdb, .cube)
- **Deprecated Code** (`graveyard.py`): Old/unused code
- **Basic Documentation** (`README.md`): Current functionality and future plans

## Development Phases

### Phase 1: Project Foundation (Immediate)

**Goal**: Establish proper Python package infrastructure and fix critical issues

#### Tasks:
- [ ] **Fix `__init__.py` syntax error**
  - Current: `from .3D.plotlyMol3D import *` is invalid (module names cannot start with digits)
  - Solution: Either rename `3D` directory to `plotlymol3d` or fix import structure
  
- [ ] **Improve `.gitignore`**
  - Add `__pycache__/` directories
  - Add `*.pyc`, `*.pyo`, `*.pyd` files
  - Add `.Python`, `*.egg-info/`, `dist/`, `build/`
  - Add virtual environment directories (`venv/`, `env/`, `.venv/`)
  - Add IDE-specific files (`.vscode/`, `.idea/`, `*.swp`)
  - Add OS-specific files (`.DS_Store`, `Thumbs.db`)

- [ ] **Create package configuration**
  - Add `pyproject.toml` (modern Python packaging) OR `setup.py` (traditional)
  - Define package metadata (name, version, author, description, license)
  - Specify package structure and entry points
  - Configure build system (setuptools, hatchling, or poetry)

- [ ] **Add `requirements.txt`**
  - List core dependencies:
    - `plotly` - Interactive plotting library
    - `numpy` - Numerical operations
    - `rdkit` - Chemistry toolkit for molecular operations
  - Consider pinning versions for reproducibility
  - Separate development dependencies (`requirements-dev.txt`)

- [ ] **Reorganize directory structure**
  - Option A: Rename `3D/` to `plotlymol3d/` (valid Python module name)
  - Option B: Create `plotlymol/` directory and move `3D/` inside
  - Ensure consistent naming conventions throughout
  - Update all relative imports accordingly

#### Success Criteria:
- Package can be installed with `pip install -e .`
- No `__pycache__` or `.pyc` files in version control
- Module can be imported without syntax errors

---

### Phase 2: Code Quality (Short-term)

**Goal**: Improve code maintainability, readability, and robustness

#### Tasks:
- [ ] **Add comprehensive type hints**
  - Annotate all function parameters and return types
  - Use `typing` module for complex types (List, Dict, Optional, Union)
  - Add type hints to dataclass fields
  - Example: `def make_atom_mesh_trace(atom: Atom, radius: Union[float, str] = DEFAULT_RADIUS, ...) -> go.Mesh3d:`

- [ ] **Add docstrings**
  - Choose documentation style (Google or NumPy style recommended)
  - Document all modules, classes, and functions
  - Include:
    - Brief description of functionality
    - Parameters with types and descriptions
    - Return values with types
    - Raises (exceptions)
    - Examples where appropriate
  - Example format (Google style):
    ```python
    def smiles_to_rdkitmol(smiles: str) -> Chem.Mol:
        """Convert SMILES string to RDKit molecule with 3D coordinates.
        
        Args:
            smiles: SMILES representation of the molecule
            
        Returns:
            RDKit molecule object with optimized 3D coordinates
            
        Raises:
            ValueError: If SMILES string is invalid
        """
    ```

- [ ] **Remove hardcoded paths**
  - Replace absolute paths in `test.py` with relative paths
  - Use `pathlib.Path` for cross-platform compatibility
  - Create fixtures or configuration for test data locations
  - Use `__file__` or `pkg_resources` to locate package data

- [ ] **Implement proper error handling**
  - Add input validation for all file formats
  - Raise descriptive exceptions with helpful error messages
  - Handle edge cases (empty files, invalid formats, missing data)
  - Add logging for debugging (using Python `logging` module)
  - Example:
    ```python
    if not os.path.exists(xyzfile):
        raise FileNotFoundError(f"XYZ file not found: {xyzfile}")
    ```

- [ ] **Clean up or remove `graveyard.py`**
  - Review deprecated code for useful components
  - Either remove entirely or move to `docs/archived/` with explanation
  - Document any breaking changes if removing functionality

#### Success Criteria:
- All functions have type hints and docstrings
- No hardcoded paths in test files
- Graceful error handling with informative messages
- Code passes linting checks (flake8/ruff)

---

### Phase 3: Testing & CI/CD (Short-term)

**Goal**: Establish automated testing and continuous integration

#### Tasks:
- [ ] **Create test infrastructure**
  - Create `tests/` directory in repository root
  - Add `tests/__init__.py`
  - Create `tests/conftest.py` for pytest fixtures
  - Add sample test data in `tests/data/` directory

- [ ] **Write unit tests**
  - Test input format parsers:
    - `test_smiles_to_rdkitmol` - SMILES parsing
    - `test_xyzfile_to_xyzblock` - XYZ file reading
    - `test_xyzblock_to_rdkitmol` - XYZ to molecule conversion
    - `test_cubefile_to_xyzblock` - Cube file parsing
  - Test molecular structure handling:
    - `test_rdkitmol_to_atoms_bonds_lists` - Atom/bond extraction
  - Test visualization components:
    - `test_make_atom_mesh_trace` - Atom rendering
    - `test_make_bond_mesh_trace` - Bond rendering
  - Test main functions:
    - `test_draw_3D_mol` - 3D molecule drawing
    - `test_draw_3D_rep` - Full representation with various inputs

- [ ] **Add GitHub Actions CI workflow**
  - Create `.github/workflows/test.yml`
  - Run tests on:
    - Multiple Python versions (3.8, 3.9, 3.10, 3.11, 3.12)
    - Multiple operating systems (Ubuntu, macOS, Windows)
  - Install dependencies and run pytest
  - Upload test results as artifacts

- [ ] **Add code coverage reporting**
  - Install `pytest-cov`
  - Configure coverage in `pyproject.toml` or `.coveragerc`
  - Add coverage reporting to CI workflow
  - Upload to Codecov or Coveralls
  - Add coverage badge to README

- [ ] **Add linting and formatting checks**
  - **Linting**: Add `flake8` or `ruff` for style checking
  - **Formatting**: Add `black` for code formatting
  - **Type checking**: Add `mypy` for static type analysis
  - Create `.github/workflows/lint.yml` for automated checks
  - Add pre-commit hooks (optional but recommended)

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

### Immediate Priority
- **`.gitignore`**: Comprehensive exclusions (already exists, needs improvement)
- **`pyproject.toml` or `setup.py`**: Package configuration
- **`requirements.txt`**: Dependencies list

### Short-term Priority
- **`CONTRIBUTING.md`**: Contribution guidelines
- **`LICENSE`**: Open source license (recommend MIT or BSD)
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
1. Fix `__init__.py` syntax error
2. Create proper `.gitignore`
3. Add `pyproject.toml` or `setup.py`
4. Add `requirements.txt`
5. Fix directory structure (rename `3D/` to valid Python module name)

### Medium Priority (1-3 months)
1. Add type hints and docstrings
2. Remove hardcoded paths
3. Create test suite with pytest
4. Set up CI/CD with GitHub Actions
5. Enhance README with examples
6. Add LICENSE file

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

**Last Updated**: 2026-01-09  
**Current Phase**: Phase 1 (Project Foundation)
