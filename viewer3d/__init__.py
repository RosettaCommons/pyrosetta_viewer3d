"""
Display PackedPose or Pose objects, or .pdb files, in py3Dmol
and NGLView within a Jupyter notebook, JupyterLab, and Google Colab.
"""
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
