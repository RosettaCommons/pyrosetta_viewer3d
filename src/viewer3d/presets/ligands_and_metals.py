__author__ = "Jason C. Klima"

import inspect
import logging

from ipywidgets.widgets import Output
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    OrResidueSelector,
    ResiduePropertySelector,
)

from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.modules import (
    setBackgroundColor,
    setDisulfides,
    setHydrogenBonds,
    setHydrogens,
    setStyle,
    setSurface,
    setZoomTo,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Any,
    Union,
)

out = Output()
_logger: logging.Logger = logging.getLogger("viewer3d.presets.ligands_and_metals")


@requires_init
def ligandsAndMetals(*args: Any, **kwargs: Any) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Display residues with `ResidueProperty.LIGAND` as `"brownCarbon"` sticks with opaque
    surface, and `ResidueProperty.METAL` as `"chainHetatm"` spheres, with `"spectrum"` cartoon
    representation, disulfide bonds, polar hydrogens, and dashed hydrogen bonds.

    Args:
        *args: Variable length argument list passed to `viewer3d.init`.

        **kwargs: Arbitrary keyword arguments passed to `viewer3d.init`.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    with out:
        metals_selector = ResiduePropertySelector(ResidueProperty.METAL)
        ligands_selector = ResiduePropertySelector(ResidueProperty.LIGAND)
        ligands_or_metals_selector = OrResidueSelector(metals_selector, ligands_selector)

    def get_backend(func, *args, **kwargs):
        sig = inspect.signature(func)
        bound = sig.bind_partial(*args, **kwargs)
        bound.apply_defaults()
        return bound.arguments.get("backend")

    backend = get_backend(init, *args, **kwargs)

    view = (
        init(*args, **kwargs)
        + setBackgroundColor()
        + setStyle(
            style="stick",
            colorscheme="grey70" if backend == "pymol" else "lightgreyCarbon",
            radius=0.15,
        )
        + setStyle(
            residue_selector=ligands_selector,
            style="stick",
            colorscheme="brownCarbon",
            radius=0.5,
            label=True,
        )
        + setStyle(
            residue_selector=metals_selector,
            style="sphere",
            colorscheme=(
                "chainHetatm"
                if backend == "py3Dmol"
                else "moleculetype"
                if backend == "nglview"
                else "gray50"
            ),
            radius=1.5,
            label=True,
        )
        + setHydrogenBonds()
        + setDisulfides(radius=0.15)
        + setHydrogens(
            color="white",
            radius=0.033 if backend != "pymol" else 0.2,
            polar_only=True,
        )
        + setSurface(
            residue_selector=ligands_selector,
            surface_type="VDW",
            opacity=0.5,
            color="brown" if backend != "py3Dmol" else None,
            colorscheme="brownCarbon" if backend == "py3Dmol" else None,
        )
        + setZoomTo(residue_selector=ligands_or_metals_selector)
    )

    return view
