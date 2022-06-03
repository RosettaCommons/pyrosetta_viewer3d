import attr
import collections
import logging
import os
import pyrosetta.distributed.io as io
import sys
import time

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Iterable, Tuple, TypeVar, Union

from viewer3d.converters import _to_float, _to_poses_pdbstrings
from viewer3d.validators import _validate_int_float, _validate_window_size
from viewer3d.base import BACKENDS
from viewer3d.modules import ModuleBase


_logger = logging.getLogger("viewer3d.core")
ModuleType = TypeVar("M", bound=ModuleBase)
import py3Dmol

try:
    import numpy
    from ipywidgets import interact, IntSlider
except ImportError:
    print(
        "Importing 'viewer3d' requires the third-party packages "
        + "'numpy', and 'ipywidgets' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/numpy/\n"
        + "https://ipywidgets.readthedocs.io/en/latest/user_install.html"
    )
    raise

from viewer3d.base import ViewerBase


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
        self._toggle_scrolling()
        if self.backend not in sys.modules:
            try:
                import py3Dmol
            except ImportError:
                print(
                    f"Using the '{self.backend}' backend requires the third-party package `{self.backend}`.\n"
                    + "Please install the package into your python environment. "
                    + "For installation instructions, visit:\n"
                    + "https://pypi.org/project/py3Dmol/\n"
                )
                raise

    def show(self):
        """Display Viewer in Jupyter notebook."""

        def view(i=0):

            _viewer = py3Dmol.view(
                width=self.window_size[0],
                height=self.window_size[1],
            )
            _pose = self.poses[i]
            _pdbstring = self.pdbstrings[i]

            if _pose:
                _viewer.addModels(io.to_pdbstring(_pose), "pdb")
            else:
                _viewer.addModels(_pdbstring, "pdb")
            _viewer.zoomTo()

            for module in self.modules:
                _viewer = module.apply(viewer=_viewer, pose=_pose, pdbstring=_pdbstring)

            self._clear_output()

            if _pose and _pose.pdb_info() and _pose.pdb_info().name():
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
        if self.backend not in sys.modules:
            try:
                import nglview
            except ImportError:
                print(
                    f"Using the '{self.backend}' backend requires the third-party package `{self.backend}`.\n"
                    + "Please install the package into your python environment. "
                    + "For installation instructions, visit:\n"
                    + "https://pypi.org/project/nglview/\n"
                )
                raise
        raise NotImplementedError(
            f"{self.__class__.__name__} is not currently supported."
        )


@attr.s(kw_only=True, slots=False, frozen=False)
class PyMOLViewer:
    poses = attr.ib(type=Pose)
    pdbstrings = attr.ib(type=PackedPose)
    window_size = attr.ib(type=Tuple[Union[int, float]])
    modules = attr.ib(type=list)
    delay = attr.ib(type=float)
    continuous_update = attr.ib(type=bool)
    backend = attr.ib(type=str)

    def __attrs_post_init__(self):
        if self.backend not in sys.modules:
            try:
                import pymol
            except ImportError:
                print(
                    f"Using the '{self.backend}' backend requires the third-party package `{self.backend}`.\n"
                    + "Please install the package into your python environment. "
                    + "For installation instructions, visit:\n"
                    + "https://anaconda.org/schrodinger/pymol\n"
                )
                raise
        raise NotImplementedError(
            f"{self.__class__.__name__} is not currently supported."
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
        type=Tuple[int, int],
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
            member_validator=attr.validators.instance_of(ModuleType),
            iterable_validator=attr.validators.instance_of(collections.abc.Iterable),
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
        converter=attr.converters.default_if_none(default="py3Dmol"),
    )

    def __attrs_post_init__(self):
        self.poses, self.pdbstrings = _to_poses_pdbstrings(
            self.packed_and_poses_and_pdbs
        )

    def initialize_viewer(self):
        viewer_kwargs = dict(
            poses=self.poses,
            pdbstrings=self.pdbstrings,
            window_size=self.window_size,
            modules=self.modules,
            delay=self.delay,
            continuous_update=self.continuous_update,
            backend=self.backend,
        )
        if self.backend == BACKENDS[0]:
            return Py3DmolViewer(**viewer_kwargs)
        elif self.backend == BACKENDS[1]:
            return NGLviewViewer(**viewer_kwargs)
        elif self.backend == BACKENDS[2]:
            return PyMOLViewer(**viewer_kwargs)


def init(
    packed_and_poses_and_pdbs=None,
    window_size=None,
    modules=None,
    delay=None,
    continuous_update=None,
    backend=None,
    *args,
    **kwargs,
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

    TODO: backend

    Returns
    -------
    A Viewer instance.
    """
    return SetupViewer(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs,
        window_size=window_size,
        modules=modules,
        delay=delay,
        continuous_update=continuous_update,
        backend=backend,
    ).initialize_viewer()
