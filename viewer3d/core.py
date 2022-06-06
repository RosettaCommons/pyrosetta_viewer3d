import attr
import collections
import logging
import os
import pyrosetta.distributed.io as io
import time


from ipywidgets import interact, IntSlider
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
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)

    def __attrs_post_init__(self):
        self._toggle_window(self.window_size)
        self.py3Dmol = self._maybe_import_backend()
        self.surface_types_dict = {
            "VDW": self.py3Dmol.VDW,
            "MS": self.py3Dmol.MS,
            "SES": self.py3Dmol.SES,
            "SAS": self.py3Dmol.SAS,
        }

    def show(self):
        """Display Py3DmolViewer in Jupyter notebook."""

        def view(i=0):
            _viewer = self.py3Dmol.view(
                width=self.window_size[0],
                height=self.window_size[1],
            )
            _pose = self.poses[i]
            _pdbstring = self.pdbstrings[i]

            if _pose is not None:
                _viewer.addModels(io.to_pdbstring(_pose), "pdb")
            else:
                _viewer.addModels(_pdbstring, "pdb")
            _viewer.zoomTo()

            for module in self.modules:
                _viewer = module.apply_py3Dmol(
                    _viewer,
                    _pose,
                    _pdbstring,
                    surface_types_dict=self.surface_types_dict,
                )

            self._clear_output()

            if _pose is not None and _pose.pdb_info() and _pose.pdb_info().name():
                _logger.debug("Decoy {0}: {1}".format(i, _pose.pdb_info().name()))

            return _viewer.show()

        time.sleep(self.delay)

        num_decoys = len(self.pdbstrings)
        if num_decoys > 1:
            s_widget = IntSlider(
                min=0,
                max=num_decoys - 1,
                description="Decoys",
                continuous_update=self.continuous_update,
            )
            widget = interact(view, i=s_widget)
        else:
            widget = view()

        return widget


@attr.s(kw_only=True, slots=False, frozen=False)
class NGLviewViewer(ViewerBase):
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)

    def __attrs_post_init__(self):
        # self._toggle_window(self.window_size)
        self.nglview = self._maybe_import_backend()

    def show(self):
        """Display NGLviewViewer in Jupyter notebook."""

        def view(i=0):
            _pose = self.poses[i]
            _pdbstring = self.pdbstrings[i]

            if _pose is not None:
                _viewer = self.nglview.show_rosetta(_pose)
            else:
                raise NotImplementedError(
                    f"PDB strings are currently not supported using the `{backend}` backend."
                )

            for module in self.modules:
                _viewer = module.apply_nglview(
                    _viewer,
                    _pose,
                    _pdbstring,
                )

            return _viewer.display(gui=True)

        time.sleep(self.delay)

        num_decoys = len(self.pdbstrings)
        if num_decoys > 1:
            s_widget = IntSlider(
                min=0,
                max=num_decoys - 1,
                description="Decoys",
                continuous_update=self.continuous_update,
            )
            widget = interact(view, i=s_widget)
        else:
            widget = view()

        return widget


@attr.s(kw_only=True, slots=False, frozen=False)
class PyMOLViewer(ViewerBase):
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
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
        self.viewer_kwargs = dict(
            poses=self.poses,
            pdbstrings=self.pdbstrings,
            window_size=self.window_size,
            modules=self.modules.copy(),
            delay=self.delay,
            continuous_update=self.continuous_update,
            backend=self.backend,
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
