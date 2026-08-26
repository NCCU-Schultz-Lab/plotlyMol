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
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple, Union, Sequence

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

# Aromatic (bond_order == 1.5) bonds are drawn as one solid cylinder plus
# one dashed cylinder to indicate resonance.
#
# AROMATIC_NUM_DASHES: dash segments per half of the dashed line.
# AROMATIC_DASH_OFFSET_FACTOR (multiplied by offset_distance = radius * 1.8):
#   how far the dashed line sits from the solid one -- raised above the
#   multi-bond default of 0.7 so the dash clears the solid bond instead of
#   hugging it.
# AROMATIC_DASH_DUTY_CYCLE: fraction of each dash+gap unit that is solid
#   material, the rest is gap.
# AROMATIC_DASH_RESOLUTION: cylinder sides for dash segments. These are
#   thin decorative rods, not load-bearing geometry, so they don't need
#   the full bond `resolution` -- capped low regardless of it.
AROMATIC_NUM_DASHES = 3
AROMATIC_DASH_OFFSET_FACTOR = 1.3
AROMATIC_DASH_DUTY_CYCLE = 0.6
AROMATIC_DASH_RESOLUTION = 10


# -----------------------------------------------------------------------
# Precomputed unit primitives
#
# Atom spheres are built from a subdivided icosahedron with explicit
# triangle faces, rather than an unstructured point cloud. Supplying
# faces up front means the browser never has to run a convex-hull
# triangulation (Plotly's `alphahull`) per atom at render time, and the
# vertex count no longer needs to be inflated to give Qhull enough points
# to work with.
# -----------------------------------------------------------------------


# `functools.cache` would collide with cube.py's wildcard-imported `cache` dict.
@lru_cache(maxsize=None)  # noqa: UP033
def _unit_icosphere(subdivisions: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """Build a unit-radius icosphere centered at the origin.

    Args:
        subdivisions: Number of times to subdivide each triangular face.
            0 gives the base icosahedron (12 verts / 20 faces), 1 gives 42
            verts / 80 faces, 2 gives 162 verts / 320 faces, 3 gives 642
            verts / 1280 faces.

    Returns:
        Tuple of (vertices, faces): vertices is an (N, 3) array on the
        unit sphere, faces is an (M, 3) array of vertex indices.
    """
    t = (1 + 5**0.5) / 2
    verts = [
        np.array(v, dtype=float)
        for v in [
            (-1, t, 0),
            (1, t, 0),
            (-1, -t, 0),
            (1, -t, 0),
            (0, -1, t),
            (0, 1, t),
            (0, -1, -t),
            (0, 1, -t),
            (t, 0, -1),
            (t, 0, 1),
            (-t, 0, -1),
            (-t, 0, 1),
        ]
    ]
    verts = [v / np.linalg.norm(v) for v in verts]
    faces = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),
        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),
        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),
        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ]

    midpoint_cache: Dict[Tuple[int, int], int] = {}

    def midpoint(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        if key in midpoint_cache:
            return midpoint_cache[key]
        m = verts[a] + verts[b]
        m = m / np.linalg.norm(m)
        verts.append(m)
        idx = len(verts) - 1
        midpoint_cache[key] = idx
        return idx

    for _ in range(subdivisions):
        new_faces = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new_faces

    return np.array(verts), np.array(faces, dtype=int)


def _subdivisions_for_resolution(resolution: int) -> int:
    """Map a legacy sphere `resolution` value onto an icosphere subdivision level."""
    if resolution <= 8:
        return 0
    elif resolution <= 20:
        return 1
    elif resolution <= 48:
        return 2
    else:
        return 3


def _icosphere_for_resolution(resolution: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return cached (vertices, faces) for a unit icosphere at the given resolution."""
    return _unit_icosphere(_subdivisions_for_resolution(resolution))


class _ColorMeshGroup:
    """Accumulates mesh geometry keyed by render color into merged Mesh3d traces.

    Rather than emitting one Mesh3d trace per atom/bond primitive (each of
    which is a separate WebGL draw call), geometry sharing a color is
    concatenated into shared vertex/face buffers and emitted as a single
    trace per color. This keeps the atom/bond count from driving the trace
    count -- a small molecule and a large one differ only in vertex count,
    not in the number of draw calls the browser has to make.
    """

    def __init__(self) -> None:
        self._verts: Dict[str, List[np.ndarray]] = {}
        self._faces: Dict[str, List[np.ndarray]] = {}
        self._offsets: Dict[str, int] = {}
        self._tags: Dict[str, List[np.ndarray]] = {}
        self._meta: Dict[str, List[Any]] = {}

    def add(
        self,
        verts: np.ndarray,
        faces: np.ndarray,
        color: str,
        tag_meta: Optional[Any] = None,
    ) -> None:
        """Add one piece of geometry to the group sharing `color`.

        `tag_meta` is an optional small metadata payload (e.g. an atom's
        world position) identifying this piece of geometry. Rather than
        repeating it on every vertex, it is appended once to a per-color
        lookup table (the trace's `meta`), and every vertex added in this
        call is tagged with a compact integer index into that table (the
        trace's `customdata`) -- so downstream code can recover which atom
        a vertex came from after merging, without bloating the payload.
        Geometry added without `tag_meta` (e.g. bonds) leaves `customdata`
        unset for that color's trace.
        """
        offset = self._offsets.get(color, 0)
        self._verts.setdefault(color, []).append(verts)
        self._faces.setdefault(color, []).append(faces + offset)
        self._offsets[color] = offset + len(verts)
        if tag_meta is not None:
            meta_list = self._meta.setdefault(color, [])
            idx = len(meta_list)
            meta_list.append(tag_meta)
            self._tags.setdefault(color, []).append(
                np.full(len(verts), idx, dtype=np.int32)
            )

    def add_traces(self, fig: go.Figure, **mesh_kwargs) -> go.Figure:
        for color in self._verts:
            V = np.vstack(self._verts[color])
            F = np.vstack(self._faces[color])
            tags = self._tags.get(color)
            customdata = np.concatenate(tags) if tags else None
            meta = self._meta.get(color)
            fig.add_trace(
                go.Mesh3d(
                    x=np.round(V[:, 0], 4),
                    y=np.round(V[:, 1], 4),
                    z=np.round(V[:, 2], 4),
                    i=F[:, 0],
                    j=F[:, 1],
                    k=F[:, 2],
                    color=color,
                    opacity=1,
                    customdata=customdata,
                    meta=meta,
                    **mesh_kwargs,
                )
            )
        return fig


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

    V, F = _icosphere_for_resolution(resolution)
    verts = V * radius_value + np.asarray(atom.atom_xyz)

    atom_trace = go.Mesh3d(
        x=verts[:, 0],
        y=verts[:, 1],
        z=verts[:, 2],
        i=F[:, 0],
        j=F[:, 1],
        k=F[:, 2],
        color=atom_colors[atom.atom_number],
        opacity=1,
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
    V0, F0 = _icosphere_for_resolution(resolution)
    group = _ColorMeshGroup()
    hover_xyz = []
    hover_text = []

    for a in atomList:
        if radius == "vdw":
            r = a.atom_vdw
        elif radius == "ball":
            r = a.atom_vdw * 0.2
        else:
            r = float(radius)

        center = np.asarray(a.atom_xyz)
        # Tag this atom's sphere with its own center, so downstream code
        # (e.g. vibration heatmap coloring) can recover which atom a
        # vertex belongs to even after merging by color.
        group.add(
            V0 * r + center,
            F0,
            atom_colors[a.atom_number],
            tag_meta=tuple(center.tolist()),
        )
        hover_xyz.append(center)
        hover_text.append(f"{a.atom_symbol}{a.atom_id}")

    group.add_traces(fig, hoverinfo="skip")

    if hover_xyz:
        # Per-atom hover names are lost once atoms are merged into shared
        # per-color meshes, so a single cheap Scatter3d trace of atom
        # centers restores hover-to-identify without adding a trace per atom.
        pts = np.array(hover_xyz)
        fig.add_trace(
            go.Scatter3d(
                x=pts[:, 0],
                y=pts[:, 1],
                z=pts[:, 2],
                mode="markers",
                marker=dict(size=1, opacity=0.01, color="black"),
                text=hover_text,
                hoverinfo="text",
                showlegend=False,
            )
        )

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


def _cylinder_mesh(
    point1: Union[List[float], np.ndarray],
    point2: Union[List[float], np.ndarray],
    radius: float = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
    add_caps: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build vertex/face arrays for a single bond cylinder.

    Shared by `make_bond_mesh_trace` (a standalone trace per call) and
    `draw_bonds` (which accumulates many of these into merged per-color
    traces instead of adding one trace per call).
    """
    p1 = np.asarray(point1)
    p2 = np.asarray(point2)
    x, y, z = generate_cylinder_mesh_rectangles(p1, p2, radius, resolution)
    V = np.column_stack([x, y, z])

    res = resolution
    faces = []
    for n in range(res):
        nxt = (n + 1) % res
        faces.append((n, n + res, nxt + res))
        faces.append((n, nxt + res, nxt))

    if add_caps:
        c_bottom = len(V)
        c_top = len(V) + 1
        V = np.vstack([V, p1, p2])
        for n in range(res):
            nxt = (n + 1) % res
            faces.append((c_bottom, nxt, n))
            faces.append((c_top, n + res, nxt + res))

    return V, np.array(faces, dtype=int)


def make_bond_mesh_trace(
    point1: Union[List[float], np.ndarray],
    point2: Union[List[float], np.ndarray],
    radius: float = DEFAULT_RADIUS,
    resolution: int = DEFAULT_RESOLUTION,
    color: str = "grey",
    add_caps: bool = True,
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
    V, F = _cylinder_mesh(point1, point2, radius, resolution, add_caps)
    return go.Mesh3d(
        x=V[:, 0],
        y=V[:, 1],
        z=V[:, 2],
        i=F[:, 0],
        j=F[:, 1],
        k=F[:, 2],
        color=color,
        opacity=1,
        hoverinfo="skip",
    )


def _oval_cap_mesh(
    center: np.ndarray,
    bond_dir: np.ndarray,
    perp_major: np.ndarray,
    semi_a: float,
    semi_b: float,
    resolution: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build vertex/face arrays for a flat elliptical end cap (multi-bond termini)."""
    perp_minor = np.cross(bond_dir, perp_major)
    norm = np.linalg.norm(perp_minor)
    if norm > 0:
        perp_minor /= norm

    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    rim = (
        center[:, None]
        + semi_a * perp_major[:, None] * np.cos(theta)
        + semi_b * perp_minor[:, None] * np.sin(theta)
    )

    V = np.vstack([rim.T, center])
    c_idx = resolution
    faces = [(c_idx, n, (n + 1) % resolution) for n in range(resolution)]
    return V, np.array(faces, dtype=int)


def _make_oval_cap(
    center: np.ndarray,
    bond_dir: np.ndarray,
    perp_major: np.ndarray,
    semi_a: float,
    semi_b: float,
    resolution: int,
    color: str,
) -> go.Mesh3d:
    """Flat elliptical end cap for multi-bond termini, as a standalone trace."""
    V, F = _oval_cap_mesh(center, bond_dir, perp_major, semi_a, semi_b, resolution)
    return go.Mesh3d(
        x=V[:, 0],
        y=V[:, 1],
        z=V[:, 2],
        i=F[:, 0],
        j=F[:, 1],
        k=F[:, 2],
        color=color,
        opacity=1,
        hoverinfo="skip",
    )


def _shrink_toward(
    anchor: np.ndarray, other: np.ndarray, trim_dist: float
) -> np.ndarray:
    """Move `anchor` toward `other` by `trim_dist`, without passing it.

    Used to pull a dashed aromatic bond's atom-adjacent endpoint back to
    the atom's sphere surface, so it doesn't spend geometry on a segment
    that's fully hidden inside the sphere anyway.
    """
    if trim_dist <= 0:
        return anchor
    vec = other - anchor
    length = np.linalg.norm(vec)
    if length <= 1e-9:
        return anchor
    trim_dist = min(trim_dist, length * 0.9)  # never collapse the segment
    return anchor + vec / length * trim_dist


def _sphere_line_entry_distance(atom_radius: float, offset_mag: float) -> float:
    """Distance along a dash line, from its closest point to an atom center,
    at which the line enters the atom's sphere.

    The dashed line runs parallel to the true bond axis, offset sideways
    from it by `offset_mag` (perpendicular to the bond, hence also
    perpendicular to the dash's own direction). That makes the dash's
    starting point already the closest point on the line to the atom
    center, so this is an exact sphere-line intersection with no need to
    search: 0 if the line passes outside the sphere already.
    """
    return float(np.sqrt(max(atom_radius**2 - offset_mag**2, 0.0)))


def draw_bonds(
    fig: go.Figure,
    bondList: List[Bond],
    atomList: Optional[List[Atom]] = None,
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
        atomList: List of Atom objects (used for calculating ring centers).
        resolution: Cylinder resolution for each bond.
        radius: Bond cylinder radius. Can be float or "ball" for ball+stick mode.

    Returns:
        The figure with bond traces added.
    """
    # Convert string radius to numeric value. In "ball" mode atoms are
    # drawn at atom_vdw * 0.2 (see draw_atoms), independent of the stick
    # radius; in "stick" mode atoms are drawn at this same numeric radius.
    # is_ball_mode is kept so dashed aromatic bonds can trim back to
    # whichever atom radius is actually in effect.
    is_ball_mode = isinstance(radius, str) and radius == "ball"
    if isinstance(radius, str):
        radius = DEFAULT_RADIUS

    group = _ColorMeshGroup()

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
            # Calculate ring center to determine correct offset direction
            ring_center_direction = perp  # Default direction

            if atomList is not None and bondList is not None:
                # Find atoms connected to both bond atoms via AROMATIC bonds only
                # This excludes substituents and gives a better ring center
                connected_atoms = []
                for other_bond in bondList:
                    # Only consider other aromatic bonds (same ring)
                    if other_bond.bond_order != 1.5:
                        continue

                    # Find atoms connected to a1 (excluding a2)
                    if (
                        other_bond.a1_id == bond.a1_id
                        and other_bond.a2_id != bond.a2_id
                    ):
                        connected_atoms.append(np.array(other_bond.a2_xyz))
                    elif (
                        other_bond.a2_id == bond.a1_id
                        and other_bond.a1_id != bond.a2_id
                    ):
                        connected_atoms.append(np.array(other_bond.a1_xyz))
                    # Find atoms connected to a2 (excluding a1)
                    if (
                        other_bond.a1_id == bond.a2_id
                        and other_bond.a2_id != bond.a1_id
                    ):
                        connected_atoms.append(np.array(other_bond.a2_xyz))
                    elif (
                        other_bond.a2_id == bond.a2_id
                        and other_bond.a1_id != bond.a1_id
                    ):
                        connected_atoms.append(np.array(other_bond.a1_xyz))

                # Calculate average position of connected atoms (ring center approximation)
                if len(connected_atoms) >= 2:
                    ring_center_approx = np.mean(connected_atoms, axis=0)
                    to_ring_center = ring_center_approx - midpoint

                    # Determine which perpendicular direction points toward ring center
                    # Use dot product to see if perp points toward or away from ring center
                    if np.dot(perp, to_ring_center) < 0:
                        ring_center_direction = -perp  # Flip direction
                    else:
                        ring_center_direction = perp

            # Place solid at center and dashed offset inward toward ring center.
            # AROMATIC_DASH_OFFSET_FACTOR pushes the dashed line further from
            # the solid one so the two are clearly separated rather than
            # hugging each other.
            offsets = [
                np.zeros(3),
                ring_center_direction * offset_distance * AROMATIC_DASH_OFFSET_FACTOR,
            ]
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
                # Dashed bond: draw segments with gaps, trimmed back from
                # each atom so no geometry is spent on the part of the
                # dash that would render fully hidden inside the sphere.
                num_dashes = AROMATIC_NUM_DASHES
                dash_resolution = min(resolution, AROMATIC_DASH_RESOLUTION)
                offset_mag = float(np.linalg.norm(offset))
                atom_r1 = bond.a1_vdw * 0.2 if is_ball_mode else radius
                atom_r2 = bond.a2_vdw * 0.2 if is_ball_mode else radius
                trim1 = _sphere_line_entry_distance(atom_r1, offset_mag)
                trim2 = _sphere_line_entry_distance(atom_r2, offset_mag)
                dash_p1 = _shrink_toward(p1, mid, trim1)
                dash_p2 = _shrink_toward(p2, mid, trim2)

                # First half of bond (atom 1 color) - dashed
                for dash_idx in range(num_dashes):
                    t_start = dash_idx / num_dashes
                    t_end = (dash_idx + AROMATIC_DASH_DUTY_CYCLE) / num_dashes
                    dash_start = dash_p1 + (mid - dash_p1) * t_start
                    dash_end = dash_p1 + (mid - dash_p1) * t_end
                    V, F = _cylinder_mesh(
                        dash_start, dash_end, r, dash_resolution, add_caps=True
                    )
                    group.add(V, F, atom_colors[bond.a1_number])

                # Second half of bond (atom 2 color) - dashed
                for dash_idx in range(num_dashes):
                    t_start = dash_idx / num_dashes
                    t_end = (dash_idx + AROMATIC_DASH_DUTY_CYCLE) / num_dashes
                    dash_start = mid + (dash_p2 - mid) * t_start
                    dash_end = mid + (dash_p2 - mid) * t_end
                    V, F = _cylinder_mesh(
                        dash_start, dash_end, r, dash_resolution, add_caps=True
                    )
                    group.add(V, F, atom_colors[bond.a2_number])
            else:
                # Solid bond: single cylinder per half
                use_oval_caps = bond_order in (2.0, 3.0)
                # First half of bond (atom 1 color)
                V, F = _cylinder_mesh(
                    p1, mid, r, resolution, add_caps=not use_oval_caps
                )
                group.add(V, F, atom_colors[bond.a1_number])

                # Second half of bond (atom 2 color)
                V, F = _cylinder_mesh(
                    mid, p2, r, resolution, add_caps=not use_oval_caps
                )
                group.add(V, F, atom_colors[bond.a2_number])

        # Oval end caps for double and triple bonds
        if bond_order in (2.0, 3.0):
            bond_dir = bond_vec / np.linalg.norm(bond_vec)
            max_offset = max(np.linalg.norm(o) for o in offsets)
            r0 = radii[0]
            semi_a = max_offset + r0
            semi_b = r0
            for center, color_num in [(a1, bond.a1_number), (a2, bond.a2_number)]:
                V, F = _oval_cap_mesh(
                    center, bond_dir, perp, semi_a, semi_b, resolution
                )
                group.add(V, F, atom_colors[color_num])

    group.add_traces(fig, hoverinfo="skip")
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
        selector=dict(type="mesh3d"),
    )

    return fig


def format_figure(fig: go.Figure, bgcolor: str = "rgba(0,0,0,0)") -> go.Figure:
    """Apply default formatting to a molecular visualization figure.

    Hides axes and grid lines for a clean molecular visualization.
    Sets aspectmode='data' to ensure spheres aren't distorted.

    Args:
        fig: Plotly figure to format.
        bgcolor: Background color for the 3D scene (default: transparent).

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
            bgcolor=bgcolor,  # Transparent background to match theme
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
            fig = draw_bonds(
                fig, bondList, atomList, resolution=resolution, radius="ball"
            )
    elif "stick" == mode:
        fig = draw_atoms(fig, atomList, resolution=resolution, radius=radius)
        fig = draw_bonds(fig, bondList, atomList, resolution=resolution, radius=radius)
    elif "vdw" == mode:
        fig = draw_atoms(fig, atomList, resolution=resolution * 4, radius="vdw")

    return fig


def create_trajectory_animation(
    xyzblocks: List[str],
    energies_hartree: Optional[List[float]] = None,
    charge: int = 0,
    mode: str = "ball+stick",
    resolution: int = 16,
    title: str = "Geometry Optimization Trajectory",
) -> go.Figure:
    """Create an animated Plotly figure stepping through geometry optimization frames.

    Each XYZ block is one step in a BFGS or similar optimization trajectory.
    The first frame is the starting geometry; the last is the optimized structure.

    Args:
        xyzblocks: List of XYZ-format strings (count line + title line + coord lines).
            Must contain at least 2 entries.
        energies_hartree: SCF energy in Hartrees for each frame (same length as
            xyzblocks). Used to annotate frame labels. Optional.
        charge: Molecular charge for bond-order perception.
        mode: Visualization mode - "ball+stick", "stick", or "vdw".
        resolution: Sphere/cylinder mesh resolution (lower = faster).
        title: Figure title.

    Returns:
        Plotly Figure with animation frames and play/step controls.

    Raises:
        ValueError: If fewer than 2 xyzblocks are provided.

    Example:
        >>> blocks = [xyz_block_1, xyz_block_2, xyz_block_3]
        >>> energies = [-75.0, -75.5, -75.6]
        >>> fig = create_trajectory_animation(blocks, energies)
        >>> fig.show()
    """
    if len(xyzblocks) < 2:
        raise ValueError(
            f"create_trajectory_animation requires at least 2 frames, "
            f"got {len(xyzblocks)}"
        )

    n_frames = len(xyzblocks)

    # Parse the reference mol (first frame) for bond connectivity.
    ref_mol = xyzblock_to_rdkitmol(xyzblocks[0], charge=charge)

    frames = []
    first_frame_data = None

    for i, xyzblock in enumerate(xyzblocks):
        # Update atom positions from this frame's XYZ block.
        raw = Chem.MolFromXYZBlock(xyzblock)
        frame_mol = Chem.RWMol(ref_mol)
        conf = frame_mol.GetConformer()
        raw_conf = raw.GetConformer()
        for atom_idx in range(frame_mol.GetNumAtoms()):
            pos = raw_conf.GetAtomPosition(atom_idx)
            conf.SetAtomPosition(atom_idx, (pos.x, pos.y, pos.z))

        # Build traces for this frame.
        empty_fig = go.Figure()
        fig_frame = draw_3D_mol(
            empty_fig, frame_mol.GetMol(), mode=mode, resolution=resolution
        )
        frame_traces = list(fig_frame.data)

        # Build frame label.
        if energies_hartree is not None and i < len(energies_hartree):
            e_label = f"Step {i}: E = {energies_hartree[i]:.6f} Hₐ"
        else:
            e_label = f"Step {i}"

        frames.append(
            go.Frame(
                data=frame_traces,
                name=f"frame_{i}",
                layout=go.Layout(title_text=f"{title} — {e_label}"),
            )
        )
        if i == 0:
            first_frame_data = frame_traces

    fig = go.Figure(data=first_frame_data, frames=frames)

    # Animation controls: play/pause button + step slider.
    fig.update_layout(
        title=title,
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {
                        "label": "▶ Play",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "frame": {"duration": 300, "redraw": True},
                                "fromcurrent": False,
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "⏸ Pause",
                        "method": "animate",
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": False},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                ],
                "x": 0.1,
                "y": 0.0,
                "xanchor": "left",
                "yanchor": "bottom",
            }
        ],
        sliders=[
            {
                "active": 0,
                "steps": [
                    {
                        "args": [
                            [f"frame_{k}"],
                            {
                                "frame": {"duration": 0, "redraw": True},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                        "label": str(k),
                        "method": "animate",
                    }
                    for k in range(n_frames)
                ],
                "x": 0.1,
                "len": 0.85,
                "xanchor": "left",
                "y": 0.0,
                "yanchor": "top",
                "pad": {"b": 10, "t": 50},
                "currentvalue": {
                    "visible": True,
                    "prefix": "Step: ",
                    "xanchor": "right",
                    "font": {"size": 14},
                },
                "transition": {"duration": 0},
            }
        ],
        scene={
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "aspectmode": "data",
        },
    )

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
        draw_3D_mol(fig, rdkitmol, resolution=resolution, radius=radius, mode=mode)
    if xyzfile is not None:
        xyzblock = xyzfile_to_xyzblock(xyzfile)
        rdkitmol = xyzblock_to_rdkitmol(xyzblock, charge=0)
        draw_3D_mol(fig, rdkitmol, resolution=resolution, radius=radius, mode=mode)
    if molfile is not None:
        rdkitmol = Chem.MolFromMolFile(molfile)
        draw_3D_mol(fig, rdkitmol, resolution=resolution, radius=radius, mode=mode)
    if cubefile is not None:
        if "molecule" in cubedraw:
            xyzblock, cubecharge = cubefile_to_xyzblock(cubefile)
            print(cubecharge)
            rdkitmol = xyzblock_to_rdkitmol(xyzblock, charge=cubecharge)
            draw_3D_mol(fig, rdkitmol, resolution=resolution, radius=radius, mode=mode)
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

    return fig
