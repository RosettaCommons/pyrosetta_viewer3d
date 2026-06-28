__author__ = "Jason C. Klima"

import itertools
import logging
import sys

from IPython.display import display
from ipywidgets.widgets import (
    Dropdown,
    Image,
    Label,
    interactive,
)
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose

from viewer3d.config import BACKENDS, COLORBAR_ATTR
from viewer3d.converters import (
    _to_backend,
    _to_poses_pdbstrings,
)
from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.presets.ligands_and_metals import ligandsAndMetals
from viewer3d.presets.per_residue_energy_metric import perResidueEnergyMetric
from viewer3d.presets.per_residue_clash_metric import perResidueClashMetric
from viewer3d.presets.per_residue_sasa_metric import perResidueSasaMetric
from viewer3d.presets.unsat_selector import unsatSelector
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.presets.rosetta_viewer")


@requires_init
def rosettaViewer(
    packed_and_poses_and_pdbs: Union[PackedPose, Pose, str, Iterable[Union[PackedPose, Pose, str]]],
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    continuous_update: bool = True,
    backend: Union[int, str] = 1,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Interactively visualize the following `viewer3d` presets:

        (0) ligandsAndMetals
        (1) unsatSelector
        (2) perResidueEnergyMetric
        (3) perResidueClashMetric
        (4) perResidueSasaMetric

    *Warning*: this `viewer3d` preset resets and sets scoring data cached in input `Pose`
    and `PackedPose` objects. It is recommended to pass cloned `Pose` and `PackedPose` objects
    to prevent clearing any important scoring data.

    Args:
        packed_and_poses_and_pdbs: A `PackedPose`, `Pose`, or `str` of PDB string or a valid
            filesystem path to a ".pdb" file, or an iterable of these objects.

        window_size: An optional `list` or `tuple` of `int` or `float` values for the
            `(width, height)` dimensions of the displayed window screen size.

        continuous_update: A `bool` object. When using the interactive slider widget, `False`
            restricts rendering to mouse button release events.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    presets = (
        ligandsAndMetals,
        unsatSelector,
        perResidueEnergyMetric,
        perResidueClashMetric,
        perResidueSasaMetric,
    )  # Presets to display using an IntSlider widget
    _preset_scoretypes = ("res_energy", "atomic_clashes", "res_sasa")
    backend = _to_backend(backend)
    view = init(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs,
        window_size=window_size,
        continuous_update=continuous_update,
        backend=backend,
    )
    _has_cache = hasattr(Pose(), "cache")

    def on_preset(i):
        preset = presets[i]
        _poses, _pdbstrings = _to_poses_pdbstrings(packed_and_poses_and_pdbs)
        # Set up poses
        poses = list(itertools.chain(*_poses.values()))
        # Clear scores
        for pose in poses:
            if _has_cache:
                pose.cache.clear()
            else:
                pose.scores.clear()

        # Score
        v = preset(
            poses,
            backend=backend,
        )
        view.poses = v.poses
        view.pdbstrings = v.pdbstrings
        view.set_modules(v.get_modules())
        view._invalidate_viewer_state()
        view.update_decoy(index=view.get_decoy_widget_index())

        if preset.__name__ not in (
            "perResidueEnergyMetric",
            "perResidueClashMetric",
            "perResidueSasaMetric",
        ):
            if hasattr(view, "viewer"):
                if hasattr(view.viewer, COLORBAR_ATTR):
                    if (backend != BACKENDS[2]) or (
                        (backend == BACKENDS[2])
                        and not isinstance(
                            getattr(view.viewer, COLORBAR_ATTR),
                            sys.modules["xmlrpc.client"]._Method,
                        )
                    ):
                        delattr(view.viewer, COLORBAR_ATTR)
        if hasattr(view, "viewer"):
            if hasattr(view.viewer, COLORBAR_ATTR):
                _value = getattr(view.viewer, COLORBAR_ATTR)
                if (backend != BACKENDS[2]) or (
                    (backend == BACKENDS[2])
                    and not isinstance(_value, sys.modules["xmlrpc.client"]._Method)
                ):
                    display(Image(value=_value))
                else:
                    display(Label(value=""))
            else:
                display(Label(value=""))

    dropdown_options = [(preset.__name__, i) for (i, preset) in enumerate(presets, start=0)]
    dropdown = Dropdown(
        options=dropdown_options,
        value=0,
        description="Preset",
    )
    preset_widget = interactive(
        on_preset,
        i=dropdown,
    )
    view.set_widgets(preset_widget)

    return view
