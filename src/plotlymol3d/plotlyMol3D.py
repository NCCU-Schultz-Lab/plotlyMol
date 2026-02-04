"""
plotlyMol3D - Interactive 3D molecular visualization using Plotly.

This module provides functions to create interactive 3D visualizations of molecules
using Plotly's Mesh3d traces. It supports multiple input formats including SMILES
strings, XYZ files, MOL files, PDB files, and cube files for orbital visualization.

Example:
    >>> from plotlymol3d import draw_3D_rep
    >>> fig = draw_3D_rep(smiles="CCO", mode="ball+stick")
"""

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem
from .atomProperties import *
from .cube import *

DEFAULT_RESOLUTION = 32
DEFAULT_RADIUS = 0.1


@dataclass
class Atom:
    """Represents an atom with its properties for 3D visualization.

    Attributes:
        atom_id: Unique identifier for the atom within the molecule.
        atom_number: Atomic number (e.g., 6 for Carbon, 8 for Oxygen).
        atom_symbol: Element symbol (e.g., "C", "O", "N").
        atom_xyz: 3D coordinates [x, y, z] in Angstroms.
        atom_vdw: Van der Waals radius in Angstroms.
    """

    atom_id: int
    atom_number: int = field(default=0)
    atom_symbol: str = field(default="unknown")
    atom_xyz: List[float] = field(default_factory=list)
    atom_vdw: float = field(default=1.70)


@dataclass
class Bond:
    """Represents a bond between two atoms for 3D visualization.

    Attributes:
        a1_id: Atom ID of the first atom.
        a2_id: Atom ID of the second atom.
        a1_number: Atomic number of the first atom.
        a2_number: Atomic number of the second atom.
        a1_xyz: 3D coordinates [x, y, z] of the first atom.
        a2_xyz: 3D coordinates [x, y, z] of the second atom.
        a1_vdw: Van der Waals radius of the first atom.
        a2_vdw: Van der Waals radius of the second atom.
        bond_order: Bond order (1=single, 2=double, 3=triple, 1.5=aromatic).
    """

    a1_id: int
    a2_id: int
    a1_number: int
    a2_number: int
    a1_xyz: List[float] = field(default_factory=list)
    a2_xyz: List[float] = field(default_factory=list)
    a1_vdw: float = field(default=1.70)
    a2_vdw: float = field(default=1.70)
    bond_order: float = field(default=1.0)


# =============================================================================
# Input Processing Functions
# =============================================================================


def cubefile_to_xyzblock(cubefile: str) -> Tuple[str, int]:
    """Extract atomic coordinates from a Gaussian cube file.

    Parses a cube file and extracts the molecular geometry in XYZ format.
    Cube files contain volumetric data (e.g., electron density, orbitals)
    along with atomic coordinates.

    Args:
        cubefile: Path to the cube file.

    Returns:
        A tuple containing:
            - xyzblock: String in XYZ format with atom count, blank line,
              and atomic coordinates.
            - total_charge: Sum of nuclear charges (integer).

    Example:
        >>> xyzblock, charge = cubefile_to_xyzblock("orbital.cube")
        >>> print(xyzblock[:50])
        12

        C          0.00000       0.00000       0.00000
    """
    total_charge: float = 0.0
    xyzblock = ""
    with open(cubefile, "r") as cf:
        for i, line in enumerate(cf):
            if i == 2:
                num_atoms = int(line.strip().split()[0])
                xyzblock = xyzblock + str(num_atoms) + "\n \n"
                stopat = 2 + 3 + num_atoms

            elif i > 5 and i <= stopat:
                if i <= stopat:
                    parts = line.strip().split()
                    atom_symbol = atom_symbols[int(parts[0])]
                    total_charge = total_charge + float(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])

                    xyzblock = (
                        xyzblock
                        + f"{atom_symbol:<3} {x:>14.5f} {y:>14.5f} {z:>14.5f} \n"
                    )
                else:
                    break
    xyzblock = xyzblock + "\n"
    print(xyzblock)
    print(f"total charge = {total_charge}")
    return xyzblock, int(total_charge)


def rdkitmol_to_atoms_bonds_lists(mol: Chem.Mol) -> Tuple[List[Atom], List[Bond]]:
    """Convert an RDKit molecule to lists of Atom and Bond objects.

    Extracts atom and bond information from an RDKit molecule object,
    including 3D coordinates, element types, and van der Waals radii.

    Args:
        mol: RDKit molecule object with 3D coordinates (must have a conformer).

    Returns:
        A tuple containing:
            - atomList: List of Atom dataclass objects.
            - bondList: List of Bond dataclass objects.

    Raises:
        ValueError: If the molecule has no conformer (3D coordinates).

    Example:
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> AllChem.EmbedMolecule(mol)
        >>> atoms, bonds = rdkitmol_to_atoms_bonds_lists(mol)
        >>> len(atoms)
        9
    """
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    conf = mol.GetConformer()

    atomList = []
    for a in atoms:
        atomList.append(
            Atom(
                atom_id=a.GetIdx(),
                atom_number=a.GetAtomicNum(),
                atom_symbol=a.GetSymbol(),
                atom_xyz=[
                    conf.GetAtomPosition(a.GetIdx()).x,
                    conf.GetAtomPosition(a.GetIdx()).y,
                    conf.GetAtomPosition(a.GetIdx()).z,
                ],
                atom_vdw=vdw_radii[a.GetAtomicNum()],
            )
        )

    bondList = []
    for b in bonds:
        # Get bond order: SINGLE=1, DOUBLE=2, TRIPLE=3, AROMATIC=1.5
        bond_type = b.GetBondType()
        if bond_type == Chem.BondType.SINGLE:
            bond_order = 1.0
        elif bond_type == Chem.BondType.DOUBLE:
            bond_order = 2.0
        elif bond_type == Chem.BondType.TRIPLE:
            bond_order = 3.0
        elif bond_type == Chem.BondType.AROMATIC:
            bond_order = 1.5
        else:
            bond_order = 1.0  # Default to single for other types

        bondList.append(
            Bond(
                a1_id=b.GetBeginAtomIdx(),
                a2_id=b.GetEndAtomIdx(),
                a1_number=b.GetBeginAtom().GetAtomicNum(),
                a2_number=b.GetEndAtom().GetAtomicNum(),
                a1_xyz=atomList[b.GetBeginAtomIdx()].atom_xyz,
                a2_xyz=atomList[b.GetEndAtomIdx()].atom_xyz,
                a1_vdw=atomList[b.GetBeginAtomIdx()].atom_vdw,
                a2_vdw=atomList[b.GetEndAtomIdx()].atom_vdw,
                bond_order=bond_order,
            )
        )

    return atomList, bondList


def smiles_to_rdkitmol(smiles: str) -> Chem.Mol:
    """Convert a SMILES string to an RDKit molecule with 3D coordinates.

    Parses the SMILES string, adds hydrogens, embeds in 3D space,
    and optimizes the geometry using the Universal Force Field (UFF).

    Args:
        smiles: SMILES representation of the molecule.

    Returns:
        RDKit molecule object with optimized 3D coordinates.

    Raises:
        ValueError: If the SMILES string is invalid.

    Example:
        >>> mol = smiles_to_rdkitmol("CCO")
        >>> mol.GetNumAtoms()
        9
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def xyzfile_to_xyzblock(file: str) -> str:
    """Read an XYZ file and return its contents as a string.

    Args:
        file: Path to the XYZ file.

    Returns:
        Contents of the XYZ file as a string block.

    Example:
        >>> xyzblock = xyzfile_to_xyzblock("molecule.xyz")
    """
    xyzblock = ""
    with open(file, "r") as f:
        for line in f:
            xyzblock = xyzblock + line
    return xyzblock


from rdkit.Chem import rdDetermineBonds


def xyzblock_to_rdkitmol(xyzblock: str, charge: int = 0) -> Chem.Mol:
    """Convert an XYZ coordinate block to an RDKit molecule with bonds.

    Parses XYZ coordinates and uses RDKit's bond perception algorithms
    to determine connectivity and bond orders.

    Args:
        xyzblock: String containing XYZ format coordinates.
        charge: Total molecular charge (used for bond order determination).

    Returns:
        RDKit molecule object with perceived bonds.

    Note:
        Bond perception may fail for certain functional groups (e.g., nitro groups).
        For problematic molecules, consider using a MOL file instead.

    Example:
        >>> xyzblock = "3\\n\\nO  0.0 0.0 0.0\\nH  0.96 0.0 0.0\\nH -0.24 0.93 0.0"
        >>> mol = xyzblock_to_rdkitmol(xyzblock)
    """
    raw_mol = Chem.MolFromXYZBlock(xyzblock)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    rdDetermineBonds.DetermineBondOrders(conn_mol, charge=charge)
    return conn_mol


# =============================================================================
# Atom Drawing Functions
# =============================================================================

DEFAULT_RADIUS = 0.1
DEFAULT_RESOLUTION = 32
a_res_scale = 10


def make_fibonacci_sphere(
    center: Sequence[float],
    radius: float = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate points on a sphere using the Fibonacci lattice method.

    Creates evenly distributed points on a sphere surface, which produces
    better visual results than latitude/longitude grids.

    Args:
        center: Center point [x, y, z] of the sphere.
        radius: Radius of the sphere.
        resolution: Number of points to generate on the sphere.

    Returns:
        Tuple of (x, y, z) numpy arrays containing point coordinates.
    """
    num_points = resolution
    indices = np.arange(0, num_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / num_points)
    theta = np.pi * (1 + 5**0.5) * indices

    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    return x, y, z


def make_atom_mesh_trace(
    atom: Atom,
    radius: Union[float, str] = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
    color: str = "grey",
) -> go.Mesh3d:
    """Create a Plotly Mesh3d trace for a single atom.

    Generates a spherical mesh representing an atom at its 3D position.

    Args:
        atom: Atom object containing position and element information.
        radius: Sphere radius. Can be a float, "vdw" for van der Waals radius,
            or "ball" for scaled VDW radius (0.2x).
        resolution: Number of points for sphere generation.
        color: Fallback color (actual color is determined by element).

    Returns:
        Plotly Mesh3d trace object for the atom.
    """
    if radius == "vdw":
        radius_value = atom.atom_vdw
    elif radius == "ball":
        radius_value = atom.atom_vdw * 0.2
    else:
        radius_value = float(radius)

    x, y, z = make_fibonacci_sphere(
        atom.atom_xyz, radius=radius_value, resolution=resolution * a_res_scale
    )

    atom_trace = go.Mesh3d(
        x=x,
        y=y,
        z=z,
        color=atom_colors[atom.atom_number],
        opacity=1,
        alphahull=0,
        name=f"{atom.atom_symbol}{atom.atom_id}",
        hoverinfo="name",
    )
    return atom_trace


def draw_atoms(
    fig: go.Figure,
    atomList: List[Atom],
    resolution: int = DEFAULT_RESOLUTION,
    radius: Union[float, str] = DEFAULT_RADIUS,
) -> go.Figure:
    """Add atom traces to a Plotly figure.

    Args:
        fig: Plotly figure to add atoms to.
        atomList: List of Atom objects to draw.
        resolution: Sphere resolution for each atom.
        radius: Atom sphere radius (float, "vdw", or "ball").

    Returns:
        The figure with atom traces added.
    """
    for a in atomList:
        a_trace = make_atom_mesh_trace(a, resolution=resolution, radius=radius)
        fig.add_trace(a_trace)
    return fig


# =============================================================================
# Bond Drawing Functions
# =============================================================================


def generate_cylinder_mesh_rectangles(
    point1: Union[List[float], np.ndarray],
    point2: Union[List[float], np.ndarray],
    radius: float = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate mesh vertices for a cylinder between two points.

    Creates vertices for the top and bottom circles of a cylinder,
    which can be used to draw bonds between atoms.

    Args:
        point1: Starting point [x, y, z] of the cylinder.
        point2: Ending point [x, y, z] of the cylinder.
        radius: Radius of the cylinder.
        resolution: Number of vertices around the circular cross-section.

    Returns:
        Tuple of (x, y, z) numpy arrays containing vertex coordinates.
    """
    point1 = np.array(point1)
    point2 = np.array(point2)

    v = point2 - point1
    height = np.linalg.norm(v)
    v = v / height  # Normalize the vector

    # Find two vectors orthogonal to the axis of the cylinder
    if np.allclose(v, np.array([0, 0, 1])) or np.allclose(v, np.array([0, 0, -1])):
        not_v = np.array([1, 0, 0])
    else:
        not_v = np.array([0, 0, 1])

    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)

    # Generate the angles for the circular cross-section
    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    circle = np.array([np.cos(theta), np.sin(theta)])

    # Generate the points for the bottom and top circles of the cylinder
    bottom_circle = point1[:, None] + radius * (
        n1[:, None] * circle[0] + n2[:, None] * circle[1]
    )
    top_circle = point2[:, None] + radius * (
        n1[:, None] * circle[0] + n2[:, None] * circle[1]
    )

    x = np.concatenate([bottom_circle[0], top_circle[0]])
    y = np.concatenate([bottom_circle[1], top_circle[1]])
    z = np.concatenate([bottom_circle[2], top_circle[2]])

    return x, y, z


def make_bond_mesh_trace(
    point1: Union[List[float], np.ndarray],
    point2: Union[List[float], np.ndarray],
    radius: float = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
    color: str = "grey",
) -> go.Mesh3d:
    """Create a Plotly Mesh3d trace for a bond (cylinder).

    Args:
        point1: Starting point [x, y, z] of the bond.
        point2: Ending point [x, y, z] of the bond.
        radius: Radius of the bond cylinder.
        resolution: Number of vertices around the cylinder.
        color: Color of the bond.

    Returns:
        Plotly Mesh3d trace object for the bond segment.
    """
    x, y, z = generate_cylinder_mesh_rectangles(point1, point2, radius, resolution)

    # Create the faces for the cylinder using rectangles
    i, j, k, l = [], [], [], []
    num_vertices = resolution
    for n in range(num_vertices):
        next_n = (n + 1) % num_vertices
        i.extend([n, next_n, next_n, n])
        j.extend([n, n, n + num_vertices, n + num_vertices])
        k.extend(
            [
                n + num_vertices,
                n + num_vertices,
                next_n + num_vertices,
                next_n + num_vertices,
            ]
        )
        l.extend([next_n, next_n + num_vertices, next_n + num_vertices, next_n])

    bond_trace = go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        color=color,
        opacity=1,
        hoverinfo="skip",
    )
    return bond_trace


def draw_bonds(
    fig: go.Figure,
    bondList: List[Bond],
    resolution: int = DEFAULT_RESOLUTION,
    radius: Union[float, str] = DEFAULT_RADIUS,
) -> go.Figure:
    """Add bond traces to a Plotly figure.

    Draws bonds as two half-cylinders colored by each atom's element.
    Double and triple bonds are shown as multiple parallel cylinders.
    Aromatic bonds are shown as 1.5 bonds (one full + one thinner).

    Args:
        fig: Plotly figure to add bonds to.
        bondList: List of Bond objects to draw.
        resolution: Cylinder resolution for each bond.
        radius: Bond cylinder radius. Can be float or "ball" for ball+stick mode.

    Returns:
        The figure with bond traces added.
    """
    # Convert string radius to numeric value
    if isinstance(radius, str):
        if radius == "ball":
            radius = DEFAULT_RADIUS  # Use default for ball+stick mode
        else:
            radius = DEFAULT_RADIUS

    for bond in bondList:
        # Calculate bond vector and midpoint
        a1 = np.array(bond.a1_xyz)
        a2 = np.array(bond.a2_xyz)
        bond_vec = a2 - a1
        midpoint = (a1 + a2) / 2

        # Get perpendicular offset vector for multiple bonds
        # Find a vector perpendicular to the bond
        if np.allclose(bond_vec / np.linalg.norm(bond_vec), [0, 0, 1]) or np.allclose(
            bond_vec / np.linalg.norm(bond_vec), [0, 0, -1]
        ):
            perp = np.array([1, 0, 0])
        else:
            perp = np.cross(bond_vec, [0, 0, 1])
        perp = perp / np.linalg.norm(perp)

        # Determine bond offsets based on bond order
        bond_order = bond.bond_order
        offset_distance = radius * 1.8  # Spacing between parallel bonds

        # Initialize dashed flags (default: all solid)
        is_dashed = None

        if bond_order == 1.0:
            # Single bond: one cylinder at center
            offsets = [np.zeros(3)]
            radii = [radius]
        elif bond_order == 2.0:
            # Double bond: two parallel cylinders
            offsets = [perp * offset_distance * 0.5, -perp * offset_distance * 0.5]
            radii = [radius * 0.7, radius * 0.7]
        elif bond_order == 3.0:
            # Triple bond: three parallel cylinders
            offsets = [
                np.zeros(3),
                perp * offset_distance * 0.7,
                -perp * offset_distance * 0.7,
            ]
            radii = [radius * 0.6, radius * 0.6, radius * 0.6]
        elif bond_order == 1.5:
            # Aromatic bond: one solid + one dashed (indicating resonance)
            offsets = [perp * offset_distance * 0.3, -perp * offset_distance * 0.3]
            radii = [radius * 0.7, radius * 0.5]
            is_dashed = [False, True]  # Second bond is dashed for aromatic
        else:
            # Default to single
            offsets = [np.zeros(3)]
            radii = [radius]

        # If no dashed flags set, default to all solid
        if is_dashed is None:
            is_dashed = [False] * len(offsets)

        # Draw each sub-bond
        for idx, (offset, r) in enumerate(zip(offsets, radii)):
            p1 = a1 + offset
            p2 = a2 + offset
            mid = midpoint + offset

            if is_dashed[idx]:
                # Dashed bond: draw segments with gaps
                num_dashes = 5  # Number of dash segments per half-bond

                # First half of bond (atom 1 color) - dashed
                for dash_idx in range(num_dashes):
                    t_start = dash_idx / num_dashes
                    t_end = (dash_idx + 0.6) / num_dashes  # 60% dash, 40% gap
                    dash_start = p1 + (mid - p1) * t_start
                    dash_end = p1 + (mid - p1) * t_end
                    bond_trace = make_bond_mesh_trace(
                        dash_start.tolist(),
                        dash_end.tolist(),
                        color=atom_colors[bond.a1_number],
                        resolution=resolution,
                        radius=r,
                    )
                    fig.add_trace(bond_trace)

                # Second half of bond (atom 2 color) - dashed
                for dash_idx in range(num_dashes):
                    t_start = dash_idx / num_dashes
                    t_end = (dash_idx + 0.6) / num_dashes
                    dash_start = mid + (p2 - mid) * t_start
                    dash_end = mid + (p2 - mid) * t_end
                    bond_trace = make_bond_mesh_trace(
                        dash_start.tolist(),
                        dash_end.tolist(),
                        color=atom_colors[bond.a2_number],
                        resolution=resolution,
                        radius=r,
                    )
                    fig.add_trace(bond_trace)
            else:
                # Solid bond: single cylinder per half
                # First half of bond (atom 1 color)
                bond_trace = make_bond_mesh_trace(
                    p1.tolist(),
                    mid.tolist(),
                    color=atom_colors[bond.a1_number],
                    resolution=resolution,
                    radius=r,
                )
                fig.add_trace(bond_trace)

                # Second half of bond (atom 2 color)
                bond_trace = make_bond_mesh_trace(
                    mid.tolist(),
                    p2.tolist(),
                    color=atom_colors[bond.a2_number],
                    resolution=resolution,
                    radius=r,
                )
                fig.add_trace(bond_trace)

    return fig


# =============================================================================
# Figure Formatting Functions
# =============================================================================


def format_lighting(
    fig: go.Figure,
    ambient: float = 0,
    diffuse: float = 1,
    specular: float = 0,
    roughness: float = 1,
    fresnel: float = 0,
    lightx: float = 1000,
    lighty: float = 1000,
    lightz: float = 1000,
) -> go.Figure:
    """Configure lighting for 3D mesh traces.

    Args:
        fig: Plotly figure to configure.
        ambient: Ambient light intensity (0-1).
        diffuse: Diffuse light intensity (0-1).
        specular: Specular highlight intensity (0-1).
        roughness: Surface roughness (0-1).
        fresnel: Fresnel effect intensity (0-1).
        lightx: Light position x-coordinate.
        lighty: Light position y-coordinate.
        lightz: Light position z-coordinate.

    Returns:
        The figure with updated lighting settings.
    """
    fig.update_traces(
        lighting=dict(
            ambient=ambient,
            diffuse=diffuse,
            specular=specular,
            roughness=roughness,
            fresnel=fresnel,
        ),
        lightposition=dict(x=lightx, y=lighty, z=lightz),
    )

    return fig


def format_figure(fig: go.Figure) -> go.Figure:
    """Apply default formatting to a molecular visualization figure.

    Hides axes and grid lines for a clean molecular visualization.
    Sets aspectmode='data' to ensure spheres aren't distorted.

    Args:
        fig: Plotly figure to format.

    Returns:
        The formatted figure.
    """
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                visible=False, showbackground=False, showgrid=False, zeroline=False
            ),
            yaxis=dict(
                visible=False, showbackground=False, showgrid=False, zeroline=False
            ),
            zaxis=dict(
                visible=False, showbackground=False, showgrid=False, zeroline=False
            ),
            aspectmode="data",  # Ensure equal scaling on all axes
        ),
        margin=dict(l=0, r=0, t=0, b=0),
    )

    return fig


# =============================================================================
# Main Drawing Functions
# =============================================================================


def draw_3D_mol(
    fig: go.Figure,
    rdkitmol: Chem.Mol,
    resolution: int = DEFAULT_RESOLUTION,
    radius: Union[float, str] = DEFAULT_RADIUS,
    mode: str = "ball+stick",
) -> go.Figure:
    """Draw a 3D molecule representation on a Plotly figure.

    Args:
        fig: Plotly figure to draw on.
        rdkitmol: RDKit molecule object with 3D coordinates.
        resolution: Mesh resolution for atoms and bonds.
        radius: Atom/bond radius (float, "vdw", or "ball").
        mode: Visualization mode - "ball+stick", "ball", "stick", or "vdw".

    Returns:
        The figure with the molecule drawn.

    Example:
        >>> fig = make_subplots()
        >>> mol = smiles_to_rdkitmol("CCO")
        >>> fig = draw_3D_mol(fig, mol, mode="ball+stick")
    """
    atomList, bondList = rdkitmol_to_atoms_bonds_lists(rdkitmol)

    if "ball" in mode:
        fig = draw_atoms(fig, atomList, resolution=resolution, radius="ball")
        if "stick" in mode:
            fig = draw_bonds(fig, bondList, resolution=resolution, radius="ball")
    elif "stick" == mode:
        fig = draw_atoms(fig, atomList, resolution=resolution, radius=radius)
    elif "vdw" == mode:
        fig = draw_atoms(fig, atomList, resolution=resolution * 4, radius="vdw")

    return fig


def draw_3D_rep(
    smiles: Optional[str] = None,
    xyzfile: Optional[str] = None,
    charge: int = 0,
    cubefile: Optional[str] = None,
    molfile: Optional[str] = None,
    pbdfile: Optional[str] = None,
    resolution: int = DEFAULT_RESOLUTION,
    radius: Union[float, str] = DEFAULT_RADIUS,
    mode: str = "ball+stick",
    orbital_opacity: float = 0.25,
    orbital_colors: Optional[List[str]] = None,
    cubedraw: str = "orbitals",
    vibration_file: Optional[str] = None,
    vibration_mode: Optional[int] = None,
    vibration_display: str = "arrows",
    vibration_amplitude: float = 1.0,
    vibration_arrow_scale: float = 1.0,
    vibration_arrow_color: str = "red",
    ambient: float = 0,
    diffuse: float = 1,
    specular: float = 0,
    roughness: float = 1,
    fresnel: float = 0,
    lightx: float = 1000,
    lighty: float = 1000,
    lightz: float = 1000,
) -> go.Figure:
    """Create a complete 3D molecular visualization from various input formats.

    This is the main entry point for creating molecular visualizations.
    Accepts multiple input formats and combines them into a single figure.

    Args:
        smiles: SMILES string for the molecule.
        xyzfile: Path to an XYZ coordinate file.
        charge: Molecular charge (used for XYZ bond perception).
        cubefile: Path to a Gaussian cube file (for orbitals).
        molfile: Path to a MOL file.
        pbdfile: Path to a PDB file (not yet implemented).
        resolution: Mesh resolution for rendering.
        radius: Atom/bond radius (float, "vdw", or "ball").
        mode: Visualization mode - "ball+stick", "ball", "stick", or "vdw".
        orbital_opacity: Opacity for orbital isosurfaces (0-1).
        orbital_colors: Colors for positive/negative orbital lobes.
            Defaults to ["darkorange", "skyblue"].
        cubedraw: What to draw from cube file - "orbitals", "molecule", or both.
        vibration_file: Path to vibration file (.log, .out, .molden).
        vibration_mode: Mode number to visualize (1-based).
        vibration_display: "arrows", "heatmap", or "both".
        vibration_amplitude: Displacement amplitude scale.
        vibration_arrow_scale: Visual scale for arrows.
        vibration_arrow_color: Color for displacement arrows.
        ambient: Ambient light intensity (0-1).
        diffuse: Diffuse light intensity (0-1).
        specular: Specular highlight intensity (0-1).
        roughness: Surface roughness (0-1).
        fresnel: Fresnel effect intensity (0-1).
        lightx: Light position x-coordinate.
        lighty: Light position y-coordinate.
        lightz: Light position z-coordinate.

    Returns:
        Plotly figure with the molecular visualization.

    Example:
        >>> fig = draw_3D_rep(smiles="CCO", mode="ball+stick", ambient=0.1)
        >>> fig = draw_3D_rep(cubefile="orbital.cube", cubedraw="orbitals")
        >>> fig = draw_3D_rep(smiles="O", vibration_file="water.log", vibration_mode=1)
    """
    if orbital_colors is None:
        orbital_colors = ["darkorange", "skyblue"]

    fig = make_subplots()
    fig = format_figure(fig)

    if smiles is not None:
        rdkitmol = smiles_to_rdkitmol(smiles)
        draw_3D_mol(fig, rdkitmol)
    if xyzfile is not None:
        xyzblock = xyzfile_to_xyzblock(xyzfile)
        rdkitmol = xyzblock_to_rdkitmol(xyzblock, charge=0)
        draw_3D_mol(fig, rdkitmol)
    if molfile is not None:
        rdkitmol = Chem.MolFromMolFile(molfile)
        draw_3D_mol(fig, rdkitmol)
    if cubefile is not None:
        if "molecule" in cubedraw:
            xyzblock, cubecharge = cubefile_to_xyzblock(cubefile)
            print(cubecharge)
            rdkitmol = xyzblock_to_rdkitmol(xyzblock, charge=cubecharge)
            draw_3D_mol(fig, rdkitmol)
        if "orbitals" in cubedraw:
            draw_cube_orbitals(fig, cubefile, orbital_opacity, orbital_colors)

    # Add vibration visualization if requested
    if vibration_file is not None and vibration_mode is not None:
        from .vibrations import parse_vibrations, add_vibrations_to_figure

        vib_data = parse_vibrations(vibration_file)
        fig = add_vibrations_to_figure(
            fig=fig,
            vib_data=vib_data,
            mode_number=vibration_mode,
            display_type=vibration_display,
            amplitude=vibration_amplitude,
            arrow_scale=vibration_arrow_scale,
            arrow_color=vibration_arrow_color,
        )

    format_lighting(
        fig,
        ambient=ambient,
        diffuse=diffuse,
        specular=specular,
        roughness=roughness,
        fresnel=fresnel,
        lightx=lightx,
        lighty=lighty,
        lightz=lightz,
    )
    fig.show("browser")

    return fig
