"""
Unit tests for input processing functions in plotlyMol3D.

Tests cover:
- SMILES to RDKit molecule conversion
- XYZ file/block processing
- Cube file parsing
- RDKit molecule to atom/bond list conversion
"""
import pytest
from pathlib import Path

from plotlymol3d.plotlyMol3D import (
    smiles_to_rdkitmol,
    xyzfile_to_xyzblock,
    xyzblock_to_rdkitmol,
    cubefile_to_xyzblock,
    rdkitmol_to_atoms_bonds_lists,
    Atom,
    Bond,
)


class TestSmilesToRdkitmol:
    """Tests for smiles_to_rdkitmol function."""
    
    def test_simple_smiles(self, sample_smiles):
        """Test conversion of a simple SMILES string (ethanol)."""
        mol = smiles_to_rdkitmol(sample_smiles)
        
        assert mol is not None
        # Ethanol: C2H5OH = 9 atoms (2C + 6H + 1O)
        assert mol.GetNumAtoms() == 9
        # Should have a conformer with 3D coordinates
        assert mol.GetNumConformers() == 1
    
    def test_complex_smiles(self, sample_smiles_complex):
        """Test conversion of a more complex SMILES string."""
        mol = smiles_to_rdkitmol(sample_smiles_complex)
        
        assert mol is not None
        assert mol.GetNumAtoms() > 0
        assert mol.GetNumConformers() == 1
    
    def test_deterministic_coordinates(self, sample_smiles):
        """Test that coordinates are deterministic (same random seed)."""
        mol1 = smiles_to_rdkitmol(sample_smiles)
        mol2 = smiles_to_rdkitmol(sample_smiles)
        
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        # Coordinates should be identical due to fixed random seed
        for i in range(mol1.GetNumAtoms()):
            pos1 = conf1.GetAtomPosition(i)
            pos2 = conf2.GetAtomPosition(i)
            assert abs(pos1.x - pos2.x) < 1e-6
            assert abs(pos1.y - pos2.y) < 1e-6
            assert abs(pos1.z - pos2.z) < 1e-6


class TestXyzProcessing:
    """Tests for XYZ file processing functions."""
    
    def test_xyzfile_to_xyzblock(self, sample_xyz_file):
        """Test reading XYZ file to block string."""
        if not Path(sample_xyz_file).exists():
            pytest.skip("Sample XYZ file not found")
        
        xyzblock = xyzfile_to_xyzblock(sample_xyz_file)
        
        assert xyzblock is not None
        assert len(xyzblock) > 0
        # XYZ files start with atom count
        lines = xyzblock.strip().split('\n')
        assert len(lines) >= 3  # At least: count, comment, one atom
    
    def test_xyzblock_to_rdkitmol(self, sample_xyz_file):
        """Test converting XYZ block to RDKit molecule."""
        if not Path(sample_xyz_file).exists():
            pytest.skip("Sample XYZ file not found")
        
        xyzblock = xyzfile_to_xyzblock(sample_xyz_file)
        mol = xyzblock_to_rdkitmol(xyzblock, charge=0)
        
        assert mol is not None
        assert mol.GetNumAtoms() > 0


class TestCubefileProcessing:
    """Tests for cube file processing."""
    
    def test_cubefile_to_xyzblock(self, sample_cube_file):
        """Test extracting XYZ coordinates from cube file."""
        if not Path(sample_cube_file).exists():
            pytest.skip("Sample cube file not found")
        
        xyzblock, charge = cubefile_to_xyzblock(sample_cube_file)
        
        assert xyzblock is not None
        assert len(xyzblock) > 0
        assert isinstance(charge, int)


class TestRdkitmolToAtomsLists:
    """Tests for rdkitmol_to_atoms_bonds_lists function."""
    
    def test_atom_list_structure(self, sample_smiles):
        """Test that atom list contains proper Atom objects."""
        mol = smiles_to_rdkitmol(sample_smiles)
        atomList, bondList = rdkitmol_to_atoms_bonds_lists(mol)
        
        assert len(atomList) == mol.GetNumAtoms()
        
        for atom in atomList:
            assert isinstance(atom, Atom)
            assert atom.atom_id >= 0
            assert atom.atom_number > 0
            assert len(atom.atom_symbol) > 0
            assert len(atom.atom_xyz) == 3
            assert atom.atom_vdw > 0
    
    def test_bond_list_structure(self, sample_smiles):
        """Test that bond list contains proper Bond objects."""
        mol = smiles_to_rdkitmol(sample_smiles)
        atomList, bondList = rdkitmol_to_atoms_bonds_lists(mol)
        
        assert len(bondList) == mol.GetNumBonds()
        
        for bond in bondList:
            assert isinstance(bond, Bond)
            assert bond.a1_id >= 0
            assert bond.a2_id >= 0
            assert bond.a1_id != bond.a2_id
            assert len(bond.a1_xyz) == 3
            assert len(bond.a2_xyz) == 3
    
    def test_ethanol_structure(self, sample_smiles):
        """Test ethanol has expected atoms and bonds."""
        mol = smiles_to_rdkitmol(sample_smiles)
        atomList, bondList = rdkitmol_to_atoms_bonds_lists(mol)
        
        # Count atom types
        symbols = [a.atom_symbol for a in atomList]
        assert symbols.count('C') == 2
        assert symbols.count('O') == 1
        assert symbols.count('H') == 6
