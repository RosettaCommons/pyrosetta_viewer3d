"""
Display PackedPose or Pose objects, or .pdb files, in py3Dmol within a Jupyter notebook.

Usage:
import viewer3d

Example Jupyter notebook commands:
--------------------------------------------------------------------------------

view = viewer3d.init("path/to/pdbfile.pdb")
view.show()

--------------------------------------------------------------------------------

import logging
logging.basicConfig(level=logging.WARNING)
import pyrosetta
pyrosetta.init("-mute all")

pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("5BVL")

view = viewer3d.init(pose, window_size=(800, 600))
view() # Equivalent to view.show()

--------------------------------------------------------------------------------

poses = [pyrosetta.toolbox.rcsb.pose_from_rcsb(id) for id in ["5BVL", "6MSR", "1QCQ"]]

view = viewer3d.init(poses) \
+ viewer3d.setStyle(colorscheme="lightgreyCarbon") \
+ viewer3d.setHydrogenBonds()
view()

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

view = sum(
    [
        viewer3d.init(packed_pose),
        viewer3d.setStyle(cartoon=False, style="sphere", radius=1.5, colorscheme="darkgreyCarbon"),
        viewer3d.setZoom(factor=1.5)
    ]
)
view.show()

--------------------------------------------------------------------------------

pose = pyrosetta.toolbox.rcsb.pose_from_rcsb("6MSR")
chA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
chB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")

view = sum(
    [
        viewer3d.init(pose),
        viewer3d.setStyle(cartoon_color="lightgrey", radius=0.25),
        viewer3d.setSurface(residue_selector=chA, colorscheme="greenCarbon", opacity=0.65, surface_type="VDW"),
        viewer3d.setSurface(residue_selector=chB, color="blue", opacity=1.0, surface_type="SAS"),
        viewer3d.setDisulfides(radius=0.25),
        viewer3d.setZoom(factor=1.5)
    ]
)
view()

--------------------------------------------------------------------------------

helix_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("H")
sheet_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("E")
loop_selector = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector("L")

modules = [
    viewer3d.setBackgroundColor(color="grey"),
    viewer3d.setStyle(residue_selector=helix_selector, cartoon_color="blue", label=False, radius=0),
    viewer3d.setStyle(residue_selector=sheet_selector, cartoon_color="red", label=False, radius=0),
    viewer3d.setStyle(residue_selector=loop_selector, cartoon_color="white", label=False, radius=0)
]

view = viewer3d.init(poses, window_size=(1200, 600), modules=modules)
view()

--------------------------------------------------------------------------------

view.reinit() # Subtract all visualization modules previously added to the Viewer
view()

--------------------------------------------------------------------------------

# View live trajectory:

pose = pyrosetta.toolbox.pose_from_rcsb("2FD7")
view = viewer3d.init(pose, delay=0.15) + viewer3d.setStyle(radius=0.1) + viewer3d.setDisulfides(radius=0.1)
backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()

for _ in range(100):
    backrub.apply(pose)
    minimize.apply(pose)
    view.show()

--------------------------------------------------------------------------------

# Display preset custom viewers for routine visualizations:

viewer3d.presets.coreBoundarySurface(poses, window_size=(800, 600), continuous_update=True)

--------------------------------------------------------------------------------
"""

import warnings

from viewer3d.core import expand_notebook, init
from viewer3d.modules import (
    setBackgroundColor,
    setDisulfides,
    setHydrogenBonds,
    setHydrogens,
    setStyle,
    setSurface,
    setZoom,
    setZoomTo,
)

__all__ = [
    "expand_notebook",
    "init",
    "presets",
    "setBackgroundColor",
    "setDisulfides",
    "setHydrogenBonds",
    "setHydrogens",
    "setStyle",
    "setSurface",
    "setZoom",
    "setZoomTo",
]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        get_ipython().Completer.limit_to__all__ = True
    except:
        pass

expand_notebook()
