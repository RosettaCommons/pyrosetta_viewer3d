"""
Display simple preset custom viewers for routine visualizations.
"""

__author__ = "Jason C. Klima"

import logging

from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Any,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.presets.base")


@requires_init
def templatePreset(*args: Any, **kwargs: Any) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Add a description of the preset `Viewer` object here.
    """
    __author__ = ""

    view = init(*args, **kwargs)

    # Add custom Viewer commands here

    return view
