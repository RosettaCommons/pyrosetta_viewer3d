__author__ = "Jason C. Klima"

import bokeh.palettes
import logging

from pyrosetta import Pose

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueClashMetric
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

_logger: logging.Logger = logging.getLogger("viewer3d.presets.per_residue_clash_metric")


@requires_init
def perResidueClashMetric(
    poses: Union[Pose, Iterable[Pose]],
    vmin: Union[int, float] = 0,
    vmax: Union[int, float] = 10,
    log: Optional[Union[int, float]] = None,
    palette: Iterable[Union[int, str]] = list(reversed(bokeh.palettes.Reds256)),
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    backend: Union[int, str] = 1,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Score the input `Pose` object(s) with `PerResidueClashMetric` and color sidechains by
    per-residue clash score, with cartoon backbone, polar hydrogens, hydrogen bonds, and
    disulfide bonds also shown.

    Args:
        poses: A `Pose` object or iterable of `Pose` objects to score and display.

        vmin: A `float` or `int` object representing the minimum clash score value for color
            map. If `None`, set `vmin` to the minimum scoretype value.

        vmax: A `float` or `int` object representing the maximum clash score value for color
            map. If `None`, set `vmin` to the maximum scoretype value.

        log: `None` to map colors spaced evenly on a linear scale between `vmin` to `vmax`. If
            an `int` or `float` object is provided, map colors spaced evenly on a log scale
            with the base provided.

        palette: An iterable of `str` (or `int`) objects representing a color map.

        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    c = PerResidueClashMetric()
    c.set_output_as_pdb_nums(output_as_pdb_nums=True)
    c.set_residue_selector(TrueResidueSelector())
    c.set_secondary_residue_selector(TrueResidueSelector())
    c.set_soft_dampening(dampening=0.33)
    c.set_use_hydrogens(use_hydrogens=True)
    c.set_use_soft_clash(soft_clash_check=True)
    apply_metric_to_poses(c, poses)
    backend = _to_backend(backend)
    v = init(poses, window_size=window_size, backend=backend)
    v += setBackgroundColor(color="white")
    v += setPerResidueRealMetric(
        scoretype="atomic_clashes",
        colorbar_label="Per-Residue Clashes",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(False, True),
        colorbar_discrete_ticks=True,
        colorbar_nticks=len(range(vmin, vmax + 1)),
    )
    v += setHydrogens(polar_only=True, color="lightgray" if backend != "pymol" else "gray70")
    v += setHydrogenBonds()
    v += setDisulfides()

    return v
