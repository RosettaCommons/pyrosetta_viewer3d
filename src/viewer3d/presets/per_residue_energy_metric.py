__author__ = "Jason C. Klima"

import bokeh.palettes
import logging
import pyrosetta

from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueEnergyMetric
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
    setPerResidueRealMetric,
)
from viewer3d.pose import apply_metric_to_poses
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.presets.per_residue_energy_metric")


@requires_init
def perResidueEnergyMetric(
    poses: Union[Pose, Iterable[Pose]],
    scorefxn: Optional[ScoreFunction] = None,
    vmin: Union[int, float] = -5,
    vmax: Union[int, float] = 5,
    log: Optional[Union[int, float]] = None,
    palette: Iterable[Union[int, str]] = list(bokeh.palettes.Greens256)
    + list(reversed(bokeh.palettes.Reds256)),
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    backend: Union[int, str] = 1,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Score the input pose(s) with `PerResidueEnergyMetric` and color sidechains by
    per-residue total score, with cartoon backbone, polar hydrogens, hydrogen bonds,
    and disulfide bonds also shown.

    Args:
        poses: A `Pose` object or iterable of `Pose` objects to score and display.

        scorefxn: An optional `ScoreFunction` object to use for scoring. If `None`, then use
            the default `"ref2015"` scorefunction.

        vmin: A `float` or `int` object representing the minimum energy score value for color
            map. If `None`, set `vmin` to the minimum scoretype value.

        vmax: A `float` or `int` object representing the maximum energy score value for color
            map. If `None`, set `vmin` to the maximum scoretype value.

        log: `None` to map colors spaced evenly on a linear scale between `vmin` to `vmax`. If
            an `int` or `float` object is provided, map colors spaced evenly on a log scale
            with the base provided.

        palette: an iterable of `str` (or `int`) objects representing a color map.

        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")
    e = PerResidueEnergyMetric()
    e.set_scorefunction(scorefxn)
    e.set_output_as_pdb_nums(output_as_pdb_nums=True)
    apply_metric_to_poses(e, poses)
    backend = _to_backend(backend)
    v = init(poses, backend=backend)
    v += setBackgroundColor(color="white")
    scorefxn_name = scorefxn.get_name()
    weights_ext = ".wts"
    if scorefxn_name.endswith(weights_ext):
        scorefxn_name = scorefxn_name[: -len(weights_ext)]
    v += setPerResidueRealMetric(
        scoretype="res_energy",
        colorbar_label=f"Per-Residue Energy ({scorefxn_name})",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(True, True),
    )
    v += setHydrogens(polar_only=True, color="lightgray" if backend != "pymol" else "gray70")
    v += setHydrogenBonds()
    v += setDisulfides()

    return v
