# viewer3d
Display `PackedPose` objects, `Pose` objects, or `.pdb` files within a Jupyter notebook and Google Colab.

# Description
The `viewer3d` macromolecular viewer quickly renders `.pdb` files, dynamically instantiating `Pose` objects if required for certain visualization modules (matching the name `viewer3d.set*`). So when adding visualization modules to the viewer or using presets, passing `Pose` or `PackedPose` objects to the viewer is suggested for quicker rendering. If a `Pose` object or `list`, `tuple`, or `set` of `Pose` objects are provided to the viewer, the `Pose`(s) pointer location(s) in memory remain fixed, and so the viewer can dynamically update upon `Pose` conformational changes by calling the following methods:

`add_pose(pose, index=None, update_viewer=True)`
`add_pdbstring(pdbstring, index=None, update_viewer=True)`
`remove_pose(index=None, model=None, update_viewer=True)`
`remove_pdbstring(index=None, model=None, update_viewer=True)`
`update_pose(pose, index=None, model=None, update_viewer=True)`
`update_pdbstring(pdbstring, index=None, model=None, update_viewer=True)`
`update_poses(poses, index=None, update_viewer=True)`
`update_pdbstrings(pdbstrings, index=None, update_viewer=True)`

See docstrings for more information. The viewer applies visualization modules in the same order they are added (from left to right), so layering different styles (and `ResidueSelector`s) on top of one another becomes possible. The user must have already initialized PyRosetta providing `.params` files for any ligands and non-canonical residues in the input molecule(s), otherwise `pyrosetta.distributed` automatically initializes PyRosetta with default command line options.

# Contributing
Please open a pull request with an updated unit test from a custom branch.

# Unit Tests
`python -m unittest`

# PyPI Releases
1. Update version number in `setup.py`
2. `pip install --user --upgrade setuptools wheel twine`
3. `rm -rf dist/*`
4. `python setup.py sdist bdist_wheel`
5. `python -m twine upload dist/*`


# Usage Examples:
```
import viewer3d
```

Example Jupyter notebook commands:

```
v = viewer3d.init("path/to/pdbfile.pdb")
v.show()
```

```
import logging
logging.basicConfig(level=logging.WARNING)
import pyrosetta
pyrosetta.init("-mute all")

pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("5BVL")

view = viewer3d.init(pose, window_size=(800, 600), backend="nglview")
view() # Equivalent to view.show()
```

```
poses = [pyrosetta.toolbox.rcsb.pose_from_rcsb(pdbid) for pdbid in ["5BVL", "6MSR", "1QCQ"]]

view = viewer3d.init(poses) \
+ viewer3d.setStyle(colorscheme="lightgreyCarbon") \
+ viewer3d.setHydrogenBonds()
view()
```

```
import pyrosetta.distributed.io as io
packed_pose = io.to_packed(pyrosetta.toolbox.pose_from_rcsb("2FD7"))
polar_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(
    pyrosetta.rosetta.core.chemical.ResidueProperty(52)
)

view = viewer3d.init(packed_pose)
view.add(viewer3d.setStyle(radius=0.1))
view.add(viewer3d.setStyle(residue_selector=polar_residue_selector, colorscheme="whiteCarbon", radius=0.25, label=False))
view.add(viewer3d.setHydrogens(color="white", polar_only=True, radius=0.1))
view.add(viewer3d.setHydrogenBonds(color="black"))
view.add(viewer3d.setDisulfides(radius=0.1))
view()
```

```
view = sum(
    [
        viewer3d.init(packed_pose),
        viewer3d.setStyle(cartoon=False, style="sphere", radius=1.5, colorscheme="darkgreyCarbon"),
        viewer3d.setZoom(factor=1.5)
    ]
)
view.show()
```

```
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("6MSR")
chA = ChainSelector("A")
chB = ChainSelector("B")

view = sum(
    [
        viewer3d.init(pose, backend=1, gui=True),
        viewer3d.setStyle(cartoon_color="lightgrey", radius=0.25),
        viewer3d.setSurface(residue_selector=chA, colorscheme="greenCarbon", opacity=0.65, surface_type="VDW"),
        viewer3d.setSurface(residue_selector=chB, color="blue", opacity=1.0, surface_type="SAS"),
        viewer3d.setDisulfides(radius=0.25),
        viewer3d.setZoom(factor=1.5)
    ]
)
view()
```

```
import pyrosetta.distributed.io as io
import viewer3d
from pyrosetta.rosetta.core.select.residue_selector import (
    LayerSelector,
    ChainSelector,
)

poses = [io.pose_from_sequence("TEST" * i) for i in range(1, 6)]
core_selector = LayerSelector()
core_selector.set_layers(True, False, False)
boundary_selector = LayerSelector()
boundary_selector.set_layers(False, True, False)
surface_selector = LayerSelector()
surface_selector.set_layers(False, False, True)
view = (
    viewer3d.init(poses, (1200, 800.01), delay=0, backend=1) \
    + viewer3d.setStyle()  #residue_selector=surface_selector, style="stick", colorscheme="element", label=False) #+ viewer3d.setSurface(opacity=0.5)
    + viewer3d.setDisulfides() \
    + viewer3d.setBackgroundColor(color="lightgray")
    + viewer3d.setHydrogenBonds()
    + viewer3d.setHydrogens(polar_only=True, color="white")
    + viewer3d.setZoomTo(ChainSelector("A")) \
    + viewer3d.setZoom(-1) \
)
view.show()
```

```
from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector
helix_selector = SecondaryStructureSelector("H")
sheet_selector = SecondaryStructureSelector("E")
loop_selector = SecondaryStructureSelector("L")

modules = [
    viewer3d.setBackgroundColor(color="grey"),
    viewer3d.setStyle(residue_selector=helix_selector, cartoon_color="blue", label=False, radius=0),
    viewer3d.setStyle(residue_selector=sheet_selector, cartoon_color="red", label=False, radius=0),
    viewer3d.setStyle(residue_selector=loop_selector, cartoon_color="white", label=False, radius=0)
]

view = viewer3d.init(poses, window_size=(1200, 600), modules=modules, continuous_update=True)
view()
```

```
view.clear() # Subtract all visualization modules previously added to the Viewer
view()
```

View a live trajectory:
```
pose = pyrosetta.toolbox.pose_from_rcsb("2FD7")
v = viewer3d.init(pose, delay=0.1) + viewer3d.setStyle(radius=0.1) + viewer3d.setDisulfides(radius=0.1)
backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()

for _ in range(100):
    backrub.apply(pose)
    minimize.apply(pose)
    v.update_pose(pose)
```

Display preset custom viewers for routine visualizations:
```
viewer3d.presets.coreBoundarySurface(poses, window_size=(800, 600), continuous_update=True)
```

```
view = viewer3d.presets.makeBundle(backend=1)
view.show()
```

Visualize the residue 1 <i>psi</i>-space of a 20-residue polyvaline extended peptide:
```
from bokeh.palettes import Viridis
n = 11 # 256
v = viewer3d.init(backend=1)
pose = pyrosetta.io.pose_from_sequence("V" * 20)
for i, hex_str in enumerate(Viridis[n]):
    v.set_modules(viewer3d.setStyle(cartoon_color=hex_str, radius=0))
    pose.set_psi(1, i * 360 / n)
    v.add_pose(pose)
```

Overlay multiple poses:
```
pose1 = pyrosetta.toolbox.rcsb.pose_from_rcsb("6V67")
pose2 = pyrosetta.toolbox.rcsb.pose_from_rcsb("1ATP")
v = viewer3d.init()
v += viewer3d.setStyle(colorscheme='redCarbon')
v += viewer3d.setSurface()
v.add_pose(pose1) # Automatically updates the viewer with the currently set modules
v.reset() # Subtract all modules
v += viewer3d.setStyle(colorscheme='blueCarbon')
v.add_pose(pose2)
```

View different sets of different overlaid poses:
```
pose = pyrosetta.toolbox.pose_from_rcsb("2FD7")
v = viewer3d.init(pose, delay=1, backend=1) #+ viewer3d.setStyle(radius=0.1) + viewer3d.setDisulfides(radius=0.1)
backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
v.set_modules([viewer3d.setStyle(), viewer3d.setDisulfides()])

for h in range(5):
    for i in range(5):
        for j in range(5):
            backrub.apply(pose)
            minimize.apply(pose)
        v.add_pose(pose.clone(), index=h, update_viewer=False)
v.show()
```
