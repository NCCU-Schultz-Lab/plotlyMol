"""
Pytest configuration and fixtures for plotlyMol tests.
"""

from pathlib import Path

import pytest


@pytest.fixture
def sample_data_dir():
    """Return the path to the sample data directory in the package."""
    return Path(__file__).parent.parent / "src" / "plotlymol3d"


@pytest.fixture
def sample_cube_file(sample_data_dir):
    """Return the path to the sample cube file."""
    return str(sample_data_dir / "anto_occ_1-min2.cube")


@pytest.fixture
def sample_mol_file(sample_data_dir):
    """Return the path to the sample MOL file."""
    return str(sample_data_dir / "cube.mol")


@pytest.fixture
def sample_xyz_file(sample_data_dir):
    """Return the path to the sample XYZ file."""
    return str(sample_data_dir / "cube.xyz")


@pytest.fixture
def sample_smiles():
    """Return a sample SMILES string for testing."""
    return "CCO"  # Ethanol


@pytest.fixture
def sample_smiles_complex():
    """Return a more complex SMILES string for testing."""
    return "CCNCOCSC"  # The example used in test.py


@pytest.fixture
def fixtures_dir():
    """Return the path to the test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def water_gaussian_log(fixtures_dir):
    """Return the path to the sample Gaussian log file for water."""
    return str(fixtures_dir / "water_gaussian.log")


@pytest.fixture
def water_orca_out(fixtures_dir):
    """Return the path to the sample ORCA output file for water."""
    return str(fixtures_dir / "water_orca.out")


@pytest.fixture
def water_molden(fixtures_dir):
    """Return the path to the sample Molden file for water."""
    return str(fixtures_dir / "water.molden")
