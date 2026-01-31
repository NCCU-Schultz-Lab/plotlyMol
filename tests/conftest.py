"""
Pytest configuration and fixtures for plotlyMol tests.
"""
import pytest
from pathlib import Path


@pytest.fixture
def sample_data_dir():
    """Return the path to the sample data directory in the package."""
    return Path(__file__).parent.parent / "plotlymol3d"


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
