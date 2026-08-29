"""
Interactively visualize `PackedPose` objects, `Pose` objects, PDB files and PDB strings
within Jupyter Notebook, JupyterLab, and Google Colab with py3Dmol, NGLview, and PyMOL
backends with a unified application programming interface (API).

## Description

The `viewer3d` bio-macromolecular viewer quickly renders PDB files and PDB strings,
dynamically instantiating `Pose` objects if required for certain visualization modules
(matching the object name `viewer3d.set*`). Therefore, when adding visualization modules to
the `Viewer` object or using the provided visualization presets, passing `Pose` or
`PackedPose` objects to a `Viewer` object is recommended for quicker rendering. If a `Pose`
or `PackedPose` object or iterable of `Pose` or `PackedPose` objects are provided to
`viewer3d`, the underlying `Pose`(s) pointer location(s) in memory remain fixed, and
therefore the visualization can dynamically update upon `Pose` conformational changes by
calling the methods on the instantiated `Py3DmolViewer`, `NGLViewViewer`, or `PyMOLViewer`
object.

## Examples

### Instantiating a `Viewer` object

```
v = viewer3d.init(poses)
```

The `Viewer` object can manage multiple `Pose` objects and PDB strings to allow overlaying
decoys for structural comparisons. Each `index` corresponds to a separate visualization
controlled by the interactive "Decoys" slider, and one or more `model` numbers may be
overlaid within each `index`. Models are internally cached according to the following scheme:

```
{
    index_0: [model_0, model_1, model_2, ...], # First visualization
    index_1: [model_0, model_1, model_2, ...], # Second visualization
    index_2: [model_0, model_1, model_2, ...], # Third visualization
    ...
}
```

### Modifying `Pose` objects

Add a `Pose` object, and optionally update the `Viewer` object in-place. If the `index` is
`None`, the currently visualized `index` is used:

```
v.add_pose(pose, index=None, update_viewer=True)
```

Remove a `Pose` object, and optionally update the `Viewer` object in-place. If the `index`
is `None`, the currently visualized `index` is used. If the `model` is `None` or out of
range, the last model at the `index` is removed:

```
v.remove_pose(index=None, model=None, update_viewer=True)
```

Update a `Pose` object, and optionally update the `Viewer` object in-place. If the `index`
is `None`, the currently visualized `index` is used. If the `model` is `None` or out of
range, the first model at the `index` is updated:

```
v.update_pose(pose, index=None, model=None, update_viewer=True)
```

Update all `Pose` objects, and optionally update the `Viewer` object in-place. If the
`index` is `None`, the currently visualized `index` is used:

```
v.update_poses(poses, index=None, update_viewer=True)
```

### Modifying PDB strings

Add a PDB string, and optionally update the `Viewer` object in-place. If the `index` is
`None`, the currently visualized `index` is used:

```
v.add_pdbstring(pdbstring, index=None, update_viewer=True)
```

Remove a PDB string, and optionally update the `Viewer` object in-place. If the `index` is
`None`, the currently visualized `index` is used. If the `model` is `None` or out of range,
the last model at the `index` is removed:

```
v.remove_pdbstring(index=None, model=None, update_viewer=True)
```

Update a PDB string, and optionally update the `Viewer` object in-place. If the `index` is
`None`, the currently visualized `index` is used. If the `model` is `None` or out of range,
the first model at the `index` is updated:

```
v.update_pdbstring(pdbstring, index=None, model=None, update_viewer=True)
```

Update all PDB strings, and optionally update the `Viewer` object in-place. If the `index`
is `None`, the currently visualized `index` is used:

```
v.update_pdbstrings(pdbstrings, index=None, update_viewer=True)
```

### Modifying visualization modules

To programmatically add visualization modules, simply add  (`+`) them to the instantiated
`Viewer` object. To programatically set visualization modules, either initialize the
`Viewer` object with `viewer3d.init(modules=[...])` syntax, or call the
`Viewer.set_modules(...)` method to overwrite the current `list` of visualization modules.
Otherwise, call the `Viewer.clear_modules()` method, then add the new visualization
modules. After `Pose` conformational changes and setting the visualization modules, call
the `Viewer.update_viewer()` method if the `Pose` object's pointer location in memory
remains fixed, otherwise call the
`Viewer.update_pose(pose, index=..., model=..., update_viewer=...)` method to update the
`Pose` object in the `Viewer` and optionally update the visualization.

The `Viewer` object applies visualization modules in the same order in which they are added
(from left to right), so layering different styles (and `ResidueSelector` objects) on top
of one another becomes possible. The user must have already initialized PyRosetta providing
Rosetta topology files and/or Rosetta patch files for any ligands and/or non-canonical
residues in the input molecule(s), otherwise the `viewer3d` framework automatically
initializes PyRosetta with default command line options.
"""

__author__ = "Jason C. Klima"

try:
    import pyrosetta
except ImportError:
    print(
        "To use `viewer3d`, please install the 'pyrosetta' package into your virtual environment. "
        "A PyRosetta license is required in order to download and use PyRosetta. "
        "Licensing is free for academic and non-profit institutions and is available "
        "to commercial users for a fee. Academic and commercial licensing of PyRosetta "
        "is handled with the standard Rosetta license through RosettaCommons.\n"
        "For more information, please visit:\n"
        "    https://www.pyrosetta.org/home/licensing-pyrosetta\n"
        "For Jupyter Notebook and Google Colab installation instructions, please visit:\n"
        "    https://github.com/RosettaCommons/PyRosetta.notebooks#chapter-10-how-to-get-started\n"
    )
    raise

import warnings

from viewer3d.base3d import expand_notebook
from viewer3d.converters import (
    atom_array_to_pdbstring,
    atom_array_to_pose,
    get_matplotlib_cmap,
    pdbstring_from_alphafold_id,
    pdbstring_from_pdb_id,
    pose_from_alphafold_id,
    pose_from_pdb_id,
)
from viewer3d.core import init
from viewer3d.modules import (
    setBackgroundColor,
    setDisulfides,
    setHydrogenBonds,
    setHydrogens,
    setPerResidueRealMetric,
    setStyle,
    setSurface,
    setZoom,
    setZoomTo,
)
from viewer3d.presets import (
    alphaFoldPLDDT,
    coreBoundarySurface,
    ligandsAndMetals,
    makeBundle,
    perResidueClashMetric,
    perResidueEnergyMetric,
    perResidueSasaMetric,
    rosettaViewer,
    unsatSelector,
)

__all__ = [
    "alphaFoldPLDDT",
    "atom_array_to_pdbstring",
    "atom_array_to_pose",
    "coreBoundarySurface",
    "expand_notebook",
    "init",
    "ligandsAndMetals",
    "makeBundle",
    "pdbstring_from_alphafold_id",
    "pdbstring_from_pdb_id",
    "perResidueClashMetric",
    "perResidueEnergyMetric",
    "perResidueSasaMetric",
    "pose_from_alphafold_id",
    "pose_from_pdb_id",
    "presets",
    "rosettaViewer",
    "setBackgroundColor",
    "setDisulfides",
    "setHydrogenBonds",
    "setHydrogens",
    "setPerResidueRealMetric",
    "setStyle",
    "setSurface",
    "setZoom",
    "setZoomTo",
    "unsatSelector",
]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        get_ipython().Completer.limit_to__all__ = True
    except:
        pass

expand_notebook()
