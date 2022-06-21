"""
Display PackedPose or Pose objects, or .pdb files, in py3Dmol
and NGLView within a Jupyter notebook, JupyterLab, and Google Colab.
"""
try:
    import pyrosetta
except ImportError:
    print(
        "To use `viewer3d`, please install `pyrosetta` into the python environment."
        + "A Rosetta license is required in order to download and use PyRosetta. "
        + "Licensing is free for academic and non-profit institutions and is available "
        + "to commercial users for a fee. Academic and commercial licensing of PyRosetta "
        + "is handled with the standard Rosetta license through Rosetta Commons.\n"
        + "For more information, please visit:\n"
        + "https://www.pyrosetta.org/home/licensing-pyrosetta \n"
        + "For Jupyter notebook and Google Colab installation instructions, please visit:\n"
        + "https://github.com/RosettaCommons/PyRosetta.notebooks#chapter-10-how-to-get-started \n"
    )
    raise

import warnings

from viewer3d.base import expand_notebook
from viewer3d.core import init
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
import viewer3d.presets

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
