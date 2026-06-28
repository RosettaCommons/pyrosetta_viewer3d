__author__ = "Jason C. Klima"

import logging
import pyrosetta

from ipywidgets.widgets import Output
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector,
    NotResidueSelector,
    OrResidueSelector,
)
from pyrosetta.rosetta.protocols.hbnet import UnsatSelector
from viewer3d.converters import _to_backend
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
    setZoomTo,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Optional,
    Union,
)

out = Output()
_logger: logging.Logger = logging.getLogger("viewer3d.presets.unsat_selector")


@requires_init
def unsatSelector(
    pose: Pose,
    scorefxn: Optional[ScoreFunction] = None,
    hbond_energy_cutoff: float = -0.5,
    backend: Union[int, str] = 0,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Visualize residues with unsatisfied backbone amine and backbone carbonyl hydrogen bonds.

    Residues with an unsatisfied backbone nitrogen are colored blue, residues with an unsatisfied
    backbone oxygen are colored red, and residues with both an unsatisfied backbone nitrogen and
    unsatisfied backbone oxygen are colored yellow. Residues with satisfied backbone hydrogen
    bonds are colored black. Cartoon representation is shown and hydrogen bonds are shown in
    black dashed lines.

    Args:
        pose: A required `Pose` object to display.

        scorefxn: An optional `ScoreFunction` object to use for scoring. If `None`, then use
            the default `"ref2015"` scorefunction.

        hbond_energy_cutoff: An optional energy cutoff for selecting unsatisfied hydrogen
            bonds.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")
    unsat_amine_selector = UnsatSelector()
    unsat_amine_selector.set_scorefxn(scorefxn)
    unsat_amine_selector.set_consider_mainchain_only(False)
    unsat_amine_selector.set_hbond_energy_cutoff(hbond_energy_cutoff)
    unsat_amine_selector.set_legacy(False)
    unsat_amine_selector.set_mode(False)
    with out:
        unsat_amine_selector.apply(pose)

    unsat_carbonyl_selector = UnsatSelector()
    unsat_carbonyl_selector.set_scorefxn(scorefxn)
    unsat_carbonyl_selector.set_consider_mainchain_only(False)
    unsat_carbonyl_selector.set_hbond_energy_cutoff(hbond_energy_cutoff)
    unsat_carbonyl_selector.set_legacy(False)
    unsat_carbonyl_selector.set_mode(True)
    with out:
        unsat_carbonyl_selector.apply(pose)

    unsat_carbonyl_and_amine_selector = AndResidueSelector(
        unsat_amine_selector, unsat_carbonyl_selector
    )
    unsat_selector = OrResidueSelector(unsat_amine_selector, unsat_carbonyl_selector)
    not_unsat_selector = NotResidueSelector(unsat_selector)

    backend = _to_backend(backend)
    h_radius = 0.1 if backend == "pymol" else 0.033
    view = (
        init(pose, backend=backend)
        + setBackgroundColor("white")
        + setStyle(
            residue_selector=not_unsat_selector,
            style="stick",
            colorscheme="blackCarbon",
            radius=0.15,
            cartoon=backend != "nglview",
            cartoon_color="black",
            label=False,
        )
        + setHydrogens(
            residue_selector=not_unsat_selector,
            color="gray",
            radius=h_radius,
            polar_only=True,
        )
        + setStyle(
            residue_selector=unsat_amine_selector,
            style="stick",
            colorscheme="blueCarbon",
            radius=0.15,
            cartoon=backend != "nglview",
            cartoon_color="blue",
            label=True,
            label_fontsize=16,
        )
        + setHydrogens(
            residue_selector=unsat_amine_selector,
            color="white",
            radius=h_radius,
            polar_only=True,
        )
        + setStyle(
            residue_selector=unsat_carbonyl_selector,
            style="stick",
            colorscheme="redCarbon",
            radius=0.15,
            cartoon=backend != "nglview",
            cartoon_color="red",
            label=True,
            label_fontsize=16,
        )
        + setStyle(
            residue_selector=unsat_carbonyl_and_amine_selector,
            style="stick",
            colorscheme="yellowCarbon",
            radius=0.15,
            cartoon=backend != "nglview",
            cartoon_color="yellow",
            label=True,
            label_fontsize=16,
        )
        + setHydrogens(
            residue_selector=unsat_carbonyl_selector,
            color="lightgray" if backend != "pymol" else "gray70",
            radius=h_radius,
            polar_only=True,
        )
        + setHydrogenBonds()
        + setDisulfides(radius=0.15)
        + setZoomTo(residue_selector=unsat_selector)
    )
    if backend == "nglview":
        view += setStyle(
            cartoon=True,
            cartoon_color="black",
            cartoon_radius=0.1,
            cartoon_opacity=0.5,
            radius=0.0,
            label=False,
        )

    return view
