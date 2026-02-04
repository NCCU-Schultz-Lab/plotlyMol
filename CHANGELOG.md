# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Vibrational Mode Visualization** - Complete system for visualizing molecular vibrations from quantum chemistry calculations
  - Three file format parsers: Gaussian (.log), ORCA (.out), Molden (.molden) with auto-detection
  - Three visualization modes:
    - Static displacement arrows using Plotly Cone traces
    - Animated vibrations with interactive controls (play/pause, frame slider)
    - Heatmap coloring by displacement magnitude
  - New `vibrations.py` module (~1000 lines) with comprehensive dataclasses and functions
  - Streamlit "📊 Vibration Settings" section with file upload and interactive controls
  - 21 new tests achieving ~95% coverage of vibration module
  - Test fixtures for all three file formats (water molecule examples)
- Exported vibration functions in `__init__.py` for public API access
- Comprehensive vibration documentation in README with code examples
- Symbol-to-atomic-number mapping (`symbol_to_number`) in `atomProperties.py`

### Changed
- Enhanced Streamlit app with vibration file uploader and parameter controls
- Updated README with vibration visualization examples and available parsers
- Expanded test suite from 26 to 47 tests
- Updated CLAUDE.md with vibration module documentation and architecture details

### Previous Changes
- Expanded README with badges, features list, quick start, orbital example, and GUI instructions.
- Added this changelog.
- Applied Black formatting and Ruff auto-fixes across plotlymol3d and tests.
- Switched to a src/ layout and moved demos into examples/.
- Consolidated dev and runtime dependencies into a single requirements.txt.

## [0.1.0] - 2026-01-31

### Added
- Core 3D molecular visualization with Plotly and RDKit integration.
- Input support for SMILES, XYZ, MOL/PDB, and cube files.
- Streamlit GUI for interactive visualization.
- Test suite and CI workflows.
