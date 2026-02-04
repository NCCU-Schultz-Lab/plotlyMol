"""Tests for molecular vibration visualization module."""

import os

import numpy as np
import pytest

from plotlymol3d.vibrations import (
    VibrationalData,
    VibrationalMode,
    parse_gaussian_vibrations,
    parse_orca_vibrations,
    parse_molden_vibrations,
    parse_vibrations,
    create_displacement_arrows,
    create_vibration_animation,
    create_heatmap_colored_figure,
    add_vibrations_to_figure,
)
from plotlymol3d import draw_3D_rep


@pytest.fixture
def water_gaussian_log():
    """Path to sample Gaussian log file for water."""
    return os.path.join(
        os.path.dirname(__file__), "fixtures", "water_gaussian.log"
    )


@pytest.fixture
def sample_vib_data():
    """Create minimal VibrationalData for testing."""
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    atomic_numbers = [8, 1, 1]

    mode1 = VibrationalMode(
        mode_number=1,
        frequency=1595.0,
        ir_intensity=75.0,
        displacement_vectors=np.array(
            [[0.0, 0.07, 0.0], [0.0, -0.56, -0.43], [0.0, -0.56, 0.43]]
        ),
        is_imaginary=False,
    )

    mode2 = VibrationalMode(
        mode_number=2,
        frequency=3657.0,
        ir_intensity=10.0,
        displacement_vectors=np.array(
            [[0.0, 0.0, 0.06], [0.0, -0.59, 0.43], [0.0, 0.59, 0.43]]
        ),
        is_imaginary=False,
    )

    mode3 = VibrationalMode(
        mode_number=3,
        frequency=3757.0,
        ir_intensity=57.0,
        displacement_vectors=np.array(
            [[0.0, -0.06, 0.0], [0.0, 0.47, -0.60], [0.0, 0.47, 0.60]]
        ),
        is_imaginary=False,
    )

    return VibrationalData(
        coordinates=coords,
        atomic_numbers=atomic_numbers,
        modes=[mode1, mode2, mode3],
        source_file="test.log",
        program="test",
    )


class TestGaussianParser:
    """Test Gaussian file parser."""

    def test_parse_gaussian_vibrations(self, water_gaussian_log):
        """Test Gaussian log file parsing."""
        vib_data = parse_gaussian_vibrations(water_gaussian_log)

        assert vib_data.program == "gaussian"
        assert len(vib_data.atomic_numbers) == 3  # H2O
        assert vib_data.atomic_numbers == [8, 1, 1]  # O, H, H
        assert len(vib_data.modes) == 3  # 3N-6 for water

        # Check coordinates shape
        assert vib_data.coordinates.shape == (3, 3)

        # Check first mode
        mode1 = vib_data.modes[0]
        assert mode1.mode_number == 1
        assert mode1.frequency == pytest.approx(1595.0338, rel=1e-4)
        assert mode1.ir_intensity == pytest.approx(75.0326, rel=1e-4)
        assert not mode1.is_imaginary
        assert mode1.displacement_vectors.shape == (3, 3)

        # Check second mode
        mode2 = vib_data.modes[1]
        assert mode2.mode_number == 2
        assert mode2.frequency == pytest.approx(3656.9382, rel=1e-4)
        assert mode2.ir_intensity == pytest.approx(10.2941, rel=1e-4)

        # Check third mode
        mode3 = vib_data.modes[2]
        assert mode3.mode_number == 3
        assert mode3.frequency == pytest.approx(3756.8229, rel=1e-4)
        assert mode3.ir_intensity == pytest.approx(57.3682, rel=1e-4)

    def test_auto_detect_format(self, water_gaussian_log):
        """Test automatic format detection."""
        vib_data = parse_vibrations(water_gaussian_log)
        assert vib_data.program == "gaussian"
        assert len(vib_data.modes) == 3

    def test_missing_file(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_gaussian_vibrations("nonexistent.log")


class TestOrcaParser:
    """Test ORCA file parser."""

    def test_parse_orca_vibrations(self, water_orca_out):
        """Test ORCA output file parsing."""
        vib_data = parse_orca_vibrations(water_orca_out)

        assert vib_data.program == "orca"
        assert len(vib_data.atomic_numbers) == 3  # H2O
        assert vib_data.atomic_numbers == [8, 1, 1]  # O, H, H
        assert len(vib_data.modes) == 3  # 3N-6 for water

        # Check coordinates shape
        assert vib_data.coordinates.shape == (3, 3)

        # Check first mode
        mode1 = vib_data.modes[0]
        assert mode1.mode_number == 1
        assert mode1.frequency == pytest.approx(1595.03, rel=1e-4)
        assert mode1.ir_intensity is None  # ORCA doesn't include IR in std output
        assert not mode1.is_imaginary
        assert mode1.displacement_vectors.shape == (3, 3)

        # Check second mode
        mode2 = vib_data.modes[1]
        assert mode2.mode_number == 2
        assert mode2.frequency == pytest.approx(3656.94, rel=1e-4)

        # Check third mode
        mode3 = vib_data.modes[2]
        assert mode3.mode_number == 3
        assert mode3.frequency == pytest.approx(3756.82, rel=1e-4)

    def test_auto_detect_orca(self, water_orca_out):
        """Test automatic format detection for ORCA."""
        vib_data = parse_vibrations(water_orca_out)
        assert vib_data.program == "orca"
        assert len(vib_data.modes) == 3


class TestMoldenParser:
    """Test Molden file parser."""

    def test_parse_molden_vibrations(self, water_molden):
        """Test Molden file parsing."""
        vib_data = parse_molden_vibrations(water_molden)

        assert vib_data.program == "molden"
        assert len(vib_data.atomic_numbers) == 3  # H2O
        assert vib_data.atomic_numbers == [8, 1, 1]  # O, H, H
        assert len(vib_data.modes) == 3  # 3N-6 for water

        # Check coordinates shape
        assert vib_data.coordinates.shape == (3, 3)

        # Check first mode
        mode1 = vib_data.modes[0]
        assert mode1.mode_number == 1
        assert mode1.frequency == pytest.approx(1595.0338, rel=1e-4)
        assert mode1.ir_intensity == pytest.approx(75.0326, rel=1e-4)
        assert not mode1.is_imaginary
        assert mode1.displacement_vectors.shape == (3, 3)

        # Check second mode
        mode2 = vib_data.modes[1]
        assert mode2.mode_number == 2
        assert mode2.frequency == pytest.approx(3656.9382, rel=1e-4)
        assert mode2.ir_intensity == pytest.approx(10.2941, rel=1e-4)

        # Check third mode
        mode3 = vib_data.modes[2]
        assert mode3.mode_number == 3
        assert mode3.frequency == pytest.approx(3756.8229, rel=1e-4)
        assert mode3.ir_intensity == pytest.approx(57.3682, rel=1e-4)

    def test_auto_detect_molden(self, water_molden):
        """Test automatic format detection for Molden."""
        vib_data = parse_vibrations(water_molden)
        assert vib_data.program == "molden"
        assert len(vib_data.modes) == 3


class TestVibrationalData:
    """Test VibrationalData dataclass methods."""

    def test_get_mode(self, sample_vib_data):
        """Test mode retrieval by number."""
        mode1 = sample_vib_data.get_mode(1)
        assert mode1 is not None
        assert mode1.mode_number == 1
        assert mode1.frequency == pytest.approx(1595.0)

        mode2 = sample_vib_data.get_mode(2)
        assert mode2 is not None
        assert mode2.mode_number == 2

        # Test non-existent mode
        mode_none = sample_vib_data.get_mode(999)
        assert mode_none is None

    def test_get_displacement_magnitudes(self, sample_vib_data):
        """Test displacement magnitude calculation."""
        magnitudes = sample_vib_data.get_displacement_magnitudes(1)

        assert len(magnitudes) == 3  # 3 atoms
        assert magnitudes[0] == pytest.approx(0.07)  # O atom
        assert magnitudes[1] > 0.5  # H atom has larger displacement
        assert magnitudes[2] > 0.5  # H atom has larger displacement

    def test_get_displacement_magnitudes_invalid_mode(self, sample_vib_data):
        """Test displacement magnitudes with invalid mode."""
        magnitudes = sample_vib_data.get_displacement_magnitudes(999)
        assert len(magnitudes) == 0


class TestVibrationalMode:
    """Test VibrationalMode dataclass."""

    def test_imaginary_frequency(self):
        """Test handling of imaginary frequencies."""
        mode = VibrationalMode(
            mode_number=1,
            frequency=-150.0,
            ir_intensity=None,
            displacement_vectors=np.zeros((3, 3)),
            is_imaginary=True,
        )

        assert mode.is_imaginary
        assert mode.frequency < 0

    def test_missing_ir_intensity(self):
        """Test modes without IR intensity data."""
        mode = VibrationalMode(
            mode_number=1,
            frequency=1500.0,
            ir_intensity=None,
            displacement_vectors=np.zeros((3, 3)),
            is_imaginary=False,
        )

        assert mode.ir_intensity is None
        assert mode.frequency > 0


class TestVisualization:
    """Test visualization functions."""

    def test_create_displacement_arrows(self, sample_vib_data):
        """Test arrow trace generation."""
        traces = create_displacement_arrows(
            vib_data=sample_vib_data, mode_number=1, amplitude=1.0
        )

        assert len(traces) > 0
        assert traces[0].type == "cone"
        assert len(traces[0].x) == 3  # 3 atoms

    def test_create_displacement_arrows_with_threshold(self, sample_vib_data):
        """Test arrow filtering by displacement threshold."""
        traces = create_displacement_arrows(
            vib_data=sample_vib_data,
            mode_number=1,
            amplitude=1.0,
            show_small_displacements=False,
            displacement_threshold=0.4,
        )

        # Should filter out O atom (displacement 0.07)
        assert len(traces[0].x) < 3

    def test_create_displacement_arrows_invalid_mode(self, sample_vib_data):
        """Test error handling for invalid mode number."""
        with pytest.raises(ValueError, match="Mode .* not found"):
            create_displacement_arrows(vib_data=sample_vib_data, mode_number=999)

    def test_add_vibrations_to_figure(self, sample_vib_data):
        """Test main integration function."""
        fig = draw_3D_rep(smiles="O", mode="ball+stick")
        original_trace_count = len(fig.data)

        fig_with_vib = add_vibrations_to_figure(
            fig=fig, vib_data=sample_vib_data, mode_number=1, display_type="arrows"
        )

        # Should have added arrow traces
        assert len(fig_with_vib.data) > original_trace_count

        # Title should be updated
        assert "Mode 1" in fig_with_vib.layout.title.text

    def test_draw_3D_rep_with_vibrations(self, water_gaussian_log):
        """Test end-to-end vibration visualization via draw_3D_rep."""
        fig = draw_3D_rep(
            smiles="O",
            vibration_file=water_gaussian_log,
            vibration_mode=1,
            vibration_display="arrows",
            vibration_amplitude=1.5,
        )

        # Should have molecular traces + vibration arrows
        assert len(fig.data) > 3  # More than just water molecule
        # Should have mode info in title
        assert "Mode 1" in fig.layout.title.text
        assert "cm" in fig.layout.title.text  # Frequency unit

    def test_create_heatmap_colored_figure(self, sample_vib_data, water_gaussian_log):
        """Test heatmap coloring of atoms by displacement magnitude."""
        # Use actual vibration data to ensure coordinates match
        vib_data = parse_gaussian_vibrations(water_gaussian_log)

        # Create figure from water molecule
        fig = draw_3D_rep(smiles="O", mode="ball+stick")
        original_trace_count = len(fig.data)

        fig_with_heatmap = create_heatmap_colored_figure(
            fig=fig,
            vib_data=vib_data,
            mode_number=1,
            colorscale="Reds",
            show_colorbar=True,
        )

        # Should have same number of traces (just modified colors)
        assert len(fig_with_heatmap.data) == original_trace_count

        # Check that at least one trace has intensity values (heatmap applied)
        has_intensity = any(
            hasattr(trace, "intensity") and trace.intensity is not None
            for trace in fig_with_heatmap.data
        )
        assert has_intensity, "No traces have intensity values (heatmap not applied)"

    def test_create_vibration_animation(self, sample_vib_data):
        """Test animation generation for vibrational mode."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Create a simple water molecule for testing
        mol = Chem.MolFromSmiles("O")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        # Create animation with few frames for speed
        fig = create_vibration_animation(
            vib_data=sample_vib_data,
            mode_number=1,
            mol=mol,
            amplitude=0.5,
            n_frames=5,
            mode="stick",
            resolution=16,
        )

        # Should have frames
        assert hasattr(fig, "frames"), "Figure should have frames"
        assert len(fig.frames) == 5, "Should have 5 frames"

        # Should have animation controls (play/pause buttons)
        assert hasattr(fig.layout, "updatemenus"), "Should have animation controls"
        assert len(fig.layout.updatemenus) > 0, "Should have at least one update menu"

        # Should have slider
        assert hasattr(fig.layout, "sliders"), "Should have slider"
        assert len(fig.layout.sliders) > 0, "Should have at least one slider"

        # Title should contain mode info
        assert "Mode 1" in fig.layout.title.text

    def test_create_heatmap_invalid_mode(self, sample_vib_data):
        """Test heatmap with invalid mode number."""
        fig = draw_3D_rep(smiles="O", mode="ball+stick")

        with pytest.raises(ValueError, match="Mode .* not found"):
            create_heatmap_colored_figure(
                fig=fig, vib_data=sample_vib_data, mode_number=999
            )

    def test_create_animation_invalid_mode(self, sample_vib_data):
        """Test animation with invalid mode number."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles("O")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        with pytest.raises(ValueError, match="Mode .* not found"):
            create_vibration_animation(
                vib_data=sample_vib_data, mode_number=999, mol=mol
            )
