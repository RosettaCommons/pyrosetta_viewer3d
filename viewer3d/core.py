import attr
import collections
import logging
import os
import pyrosetta.distributed.io as io
import time

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Iterable, Tuple, Union

from viewer3d.base import ViewerBase
from viewer3d.config import BACKENDS
from viewer3d.converters import _to_float, _to_poses_pdbstrings
from viewer3d.modules import ModuleBase
from viewer3d.validators import _validate_int_float, _validate_window_size


_logger = logging.getLogger("viewer3d.core")


@attr.s(kw_only=True, slots=False, frozen=False)
class Py3DmolViewer(ViewerBase):
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
    n_decoys = attr.ib(type=int)
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)
    _was_show_called = attr.ib(type=bool, default=False, init=False)

    def __attrs_post_init__(self):
        self._toggle_window(self.window_size)
        self.py3Dmol = self._maybe_import_backend()
        self.viewer = self.py3Dmol.view(
            width=self.window_size[0],
            height=self.window_size[1],
        )
        self.surface_types_dict = {
            "VDW": self.py3Dmol.VDW,
            "MS": self.py3Dmol.MS,
            "SES": self.py3Dmol.SES,
            "SAS": self.py3Dmol.SAS,
        }

    def update(self, _pose, _pdbstring):
        """Update Py3DmolViewer in Jupyter notebook."""
        time.sleep(self.delay)
        self.viewer.removeAllLabels()
        self.viewer.removeAllModels()
        self.viewer.removeAllShapes()
        self.viewer.removeAllSurfaces()

        if _pose is not None:
            self.viewer.addModels(io.to_pdbstring(_pose), "pdb")
        else:
            self.viewer.addModels(_pdbstring, "pdb")

        for module in self.modules:
            self.viewer = module.apply_py3Dmol(
                self.viewer,
                _pose,
                _pdbstring,
                surface_types_dict=self.surface_types_dict,
            )

        if self._was_show_called:
            self.viewer.update()

    def show(self):
        """Display Py3DmolViewer in Jupyter notebook."""
        self.setup()
        self.viewer.show()
        self._was_show_called = True


@attr.s(kw_only=True, slots=False, frozen=False)
class NGLviewViewer(ViewerBase):
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
    n_decoys = attr.ib(type=int)
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)

    def __attrs_post_init__(self):
        self.nglview = self._maybe_import_backend()
        self.viewer = self.nglview.widget.NGLWidget()

    def update(self, _pose, _pdbstring):
        """Update NGLviewViewer in Jupyter notebook."""
        time.sleep(self.delay)

        for component_id in self.viewer._ngl_component_ids:
            self.viewer.remove_component(component_id)

        if _pose is not None:
            structure = self.nglview.adaptor.RosettaStructure(_pose)
        else:
            structure = self.nglview.adaptor.TextStructure(_pdbstring, ext="pdb")
        self.viewer.add_structure(structure)

        for module in self.modules:
            self.viewer = module.apply_nglview(
                self.viewer,
                _pose,
                _pdbstring,
            )

    def show(self):
        """Display NGLviewViewer in Jupyter notebook."""
        self.setup()
        self.viewer.display(gui=True, style="ngl")
        self.viewer._ipython_display_()


@attr.s(kw_only=True, slots=False, frozen=False)
class PyMOLViewer(ViewerBase):
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
    n_decoys = attr.ib(type=int)
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)

    def __attrs_post_init__(self):
        self.pymol = self._maybe_import_backend()

        raise NotImplementedError(
            f"{self.__class__.__name__} is not currently supported."
        )

    def update(self, _pose, _pdbstring):
        """Update PyMOLViewer."""
        pass

    def show(self):
        """Display PyMOLViewer."""
        _viewer = None  # TODO
        for module in self.modules:
            _viewer = module.apply_pymol(
                _viewer,
                _pose,
                _pdbstring,
            )


@attr.s(kw_only=True, slots=False, frozen=False)
class SetupViewer:
    packed_and_poses_and_pdbs = attr.ib(
        type=Union[PackedPose, Pose, Iterable[Union[PackedPose, Pose]], None],
        default=None,
    )
    poses = attr.ib(type=Iterable[Pose], init=False)
    pdbstrings = attr.ib(type=Iterable[str], init=False)
    window_size = attr.ib(
        type=Tuple[Union[int, float], Union[int, float]],
        default=None,
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of((int, float)),
                iterable_validator=attr.validators.instance_of(
                    collections.abc.Iterable
                ),
            ),
            _validate_window_size,
        ],
        converter=attr.converters.default_if_none(default=(1200, 800)),
    )
    modules = attr.ib(
        type=list,
        default=None,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(ModuleBase),
            iterable_validator=attr.validators.instance_of(list),
        ),
        converter=attr.converters.default_if_none(default=[]),
    )
    delay = attr.ib(
        type=float,
        default=None,
        validator=_validate_int_float,
        converter=attr.converters.pipe(
            attr.converters.default_if_none(0.25), _to_float
        ),
    )
    continuous_update = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    backend = attr.ib(
        type=str,
        default=None,
        validator=[attr.validators.instance_of(str), attr.validators.in_(BACKENDS)],
        converter=attr.converters.default_if_none(default=BACKENDS[0]),
    )

    def __attrs_post_init__(self):
        self.poses, self.pdbstrings = _to_poses_pdbstrings(
            self.packed_and_poses_and_pdbs
        )
        self.n_decoys = len(self.pdbstrings)
        self.viewer_kwargs = dict(
            poses=self.poses,
            pdbstrings=self.pdbstrings,
            window_size=self.window_size,
            modules=self.modules.copy(),
            delay=self.delay,
            continuous_update=self.continuous_update,
            backend=self.backend,
            n_decoys=self.n_decoys,
        )

    def initialize_viewer(self):
        if self.backend == BACKENDS[0]:
            viewer = Py3DmolViewer(**self.viewer_kwargs)
        elif self.backend == BACKENDS[1]:
            viewer = NGLviewViewer(**self.viewer_kwargs)
        elif self.backend == BACKENDS[2]:
            viewer = PyMOLViewer(**self.viewer_kwargs)

        return viewer


def init(
    packed_and_poses_and_pdbs=None,
    window_size=None,
    modules=None,
    delay=None,
    continuous_update=None,
    backend=None,
):
    """
    Initialize the Viewer object.

    Parameters
    ----------
    first : required
        `packed_and_poses_and_pdbs`

        `PackedPose`, `Pose`, or `str` of a valid path to a .pdb file, or a `list`, `set`, or `tuple` of these objects.

    second : optional
        `window_size`

        `list` or `tuple` of `int` or `float` values for the (width, height) dimensions of the displayed window screen size.
        Default: (1200, 800)

    third : optional
        `modules`

        `list` of instantiated visualization modules to run upon changing amongst `packed_and_poses_and_pdbs` objects
        with the slider, matching the namespace `viewer3d.set*`
        Default: []

    fourth : optional
        `delay`

        `float` or `int` time delay in seconds before rendering the Viewer in a Jupyter notebook, which is useful to prevent
        overburdening the Jupyter notebook client if `for` looping over quick modifications to a `Pose`, and should be >= 0.
        Default: 0.25

    fifth : optional
        `continuous_update`

        `True` or `False`. When using the interactive slider widget, `False` restricts rendering to mouse release events.
        Default: False

    sixth : optional
        `backend`

        The viewer backend to for the visualization. Supported backends are 'py3Dmol', 'nglview', and 'pymol'.
        Default: 'py3Dmol'


    Returns
    -------
    A Viewer instance.
    """
    viewer = SetupViewer(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs,
        window_size=window_size,
        modules=modules,
        delay=delay,
        continuous_update=continuous_update,
        backend=backend,
    ).initialize_viewer()

    return viewer
