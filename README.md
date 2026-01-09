# plotlyMol

A package to create interactive molecular visualizations using Plotly.

## Installation

To install in development mode:

```bash
git clone https://github.com/jonathanschultzNU/plotlyMol.git
cd plotlyMol
pip install -e .
```

## Usage

```python
from plotlymol3d import draw_3D_rep

# Draw a molecule from SMILES
moldraw = draw_3D_rep(smiles="CCNCOCSC", mode="ball+stick", ambient=0.1)

# Draw from XYZ file
moldraw = draw_3D_rep(xyzfile="path/to/file.xyz", mode="ball+stick", ambient=0.1)

# Draw from cube file with orbitals
moldraw = draw_3D_rep(
    cubefile="path/to/file.cube",
    molfile="path/to/file.mol",
    mode="ball+stick",
    ambient=0.1,
    cubedraw="orbitals",
    orbital_opacity=0.25,
    orbital_colors=["darkorange", "darkblue"]
)
```

## Future Plans

consider the ability to be able to pass multiple inputs of the same time... [smiles, smiles, etc] <- need to handle single inputs and lists of inputs 

For the future:

- 2D: give chemdraw like structures
- 3D: 
	- draw 3D structures 
		- add in the ability to scale the position of the half-bond by the van der waals radii of the atoms at either end. 
	- add orbitals and other surfaces
	- for smiles: need to be able to handle add structures
	- PROBLEM: the rdkit xyzblock --> mol object can fail, for things like NITRO.  Need to find a new way to handle xyz (maybe bable xyz --> smiles???)
	- Need to add in the ability to draw orbitals.  Code works, just needs to be implemented. 
- 4D: Animations
- GUI
	- ability to toggle aspects of the presentation/hover/etc
	- ability to dynamically enter mol information

