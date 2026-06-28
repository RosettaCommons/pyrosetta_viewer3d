__author__ = "Jason C. Klima"

import attr
import logging
import matplotlib.colors as mcolors
import pyrosetta

from pyrosetta import Pose

from viewer3d.converters import _to_backend
from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.modules import (
    setBackgroundColor,
    setHydrogenBonds,
    setHydrogens,
    setPerResidueRealMetric,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.presets.alphafold_plddt")


@attr.define(kw_only=True, slots=True, frozen=True)
class _PerResidueBfactorMetric:
    def apply(self, pose):
        has_cache = hasattr(pose, "cache")
        for res in range(1, pose.size() + 1):
            ca_plddt = pose.pdb_info().temperature(
                res=res,
                atom_index=pose.residue(res).atom_index("CA"),
            )
            if has_cache:
                pose.cache[f"Bfact_{res}"] = ca_plddt
            else:
                pose.scores[f"Bfact_{res}"] = ca_plddt


@requires_init
def alphaFoldPLDDT(
    poses: Union[Pose, Iterable[Pose]],
    cartoon_color: Union[int, str] = "black",
    rescale: bool = False,
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    continuous_update: bool = True,
    backend: Union[int, str] = 0,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Visualize C-alpha pLDDT with AlphaFold coloring.

    Args:
        poses: A `Pose` object or an iterable of `Pose` objects.

        cartoon_color: A `str` or `int` representing the cartoon color.

        rescale: A `bool` object. If `True`, plot pLDDT values in the range [0, 1]. If
            `False`, plot pLDDT values in the range [0, 100].

        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.

        continuous_update: a `bool` object. When using the interactive slider widget,
            `False` restricts rendering to mouse button release events.

        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    try:
        xml_obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
            """<SIMPLE_METRICS><PerResidueBfactorMetric name="b" atom_type="CA"/></SIMPLE_METRICS>"""
        ).get_simple_metric("b")
    except RuntimeError:
        xml_obj = _PerResidueBfactorMetric()

    _msg = "The 'poses' argument value must be a `Pose` object or an iterable of `Pose` objects. Received: {0}"
    if isinstance(poses, Pose):
        xml_obj.apply(poses)
    elif isinstance(poses, Iterable):
        for pose in poses:
            if isinstance(pose, Pose):
                xml_obj.apply(pose)
            else:
                raise ValueError(_msg.format(type(pose)))
    else:
        raise ValueError(_msg.format(type(poses)))

    palette = list(
        map(
            mcolors.to_hex,
            [
                [0.992, 0.490, 0.302],
                [0.996, 0.851, 0.212],
                [0.996, 0.851, 0.212],
                [0.416, 0.796, 0.945],
                [0.416, 0.796, 0.945],
                [0.416, 0.796, 0.945],
                [0.051, 0.341, 0.827],
            ],
        )
    )

    backend = _to_backend(backend)
    v = init(
        poses,
        window_size=window_size,
        continuous_update=continuous_update,
        backend=backend,
    )
    v += setBackgroundColor()
    v += setPerResidueRealMetric(
        scoretype="Bfact",
        vmin=0.4 if rescale else 40.0,
        vmax=1.0 if rescale else 100.0,
        palette=palette,
        log=None,
        style="stick",
        radius=0.4,
        show_hydrogens=False,
        bonds="symmetric",
        cartoon=True,
        cartoon_color=cartoon_color,
        cartoon_radius=0.25,
        cartoon_opacity=0.9,
        colorbar=True,
        colorbar_extremes=(True, False),
        colorbar_label="pLDDT",
        colorbar_fontsize=22,
        colorbar_nticks=7,
    )
    v += setHydrogenBonds(
        color="grey",
        dashed=True,
        radius=None,
    )
    v += setHydrogens(
        color="lightgray" if backend != "pymol" else "gray70",
        polar_only=True,
        radius=0.1,
        residue_selector=None,
    )

    return v
