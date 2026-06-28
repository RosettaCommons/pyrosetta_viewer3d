__author__ = "Jason C. Klima"

import bokeh.palettes
import logging

from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring.sasa import SasaMethodHPMode
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueSasaMetric
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

_logger: logging.Logger = logging.getLogger("viewer3d.presets.per_residue_sasa_metric")


@requires_init
def perResidueSasaMetric(
    poses: Union[Pose, Iterable[Pose]],
    mode: int = 0,
    vmin: Optional[Union[int, float]] = None,
    vmax: Optional[Union[int, float]] = None,
    log: Optional[Union[int, float]] = None,
    palette: Iterable[Union[int, str]] = list(bokeh.palettes.Viridis256),
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    backend: Union[int, str] = 1,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Score the input `Pose` object(s) with `PerResidueSasaMetric` and color sidechains by
    per-residue SASA score, with cartoon backbone, polar hydrogens, hydrogen bonds, and
    disulfide bonds also shown.

    Args:
        poses: A `Pose` object or iterable of `Pose` objects to score and display.

        mode: An `int` object to set the SASA mode: `0`: all SASA; `1`: hydrophobic only; `2`:
            polar only.

        vmin: A `float` or `int` object representing the minimum SASA value for color map. If
            `None`, set `vmin` to the minimum scoretype value.

        vmax: A `float` or `int` object representing the maximum SASA value for color map.
            If `None`, set `vmin` to the maximum scoretype value.

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

    if mode == 0:
        sasa_mode = SasaMethodHPMode.ALL_SASA
    elif mode == 1:
        sasa_mode = SasaMethodHPMode.HYDROPHOBIC_SASA
    elif mode == 2:
        sasa_mode = SasaMethodHPMode.POLAR_SASA
    else:
        raise ValueError(f"The 'mode' argument value must be an `int` object in {list(range(3))}.")
    s = PerResidueSasaMetric()
    s.set_mode(sasa_mode)
    s.set_output_as_pdb_nums(True)
    s.set_residue_selector(TrueResidueSelector())
    apply_metric_to_poses(s, poses)
    backend = _to_backend(backend)
    v = init(poses, backend=backend)
    v += setBackgroundColor(color="white")
    v += setPerResidueRealMetric(
        scoretype="res_sasa",
        colorbar_label="Per-Residue SASA (Å$^{2}$)",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(False, True),
    )
    v += setHydrogens(polar_only=True, color="lightgray" if backend != "pymol" else "gray70")
    v += setHydrogenBonds()
    v += setDisulfides()

    return v
