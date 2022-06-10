import attr
import collections
import logging
import os
import pyrosetta.distributed.io as io

from ipywidgets.widgets import Widget
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Iterable, List, Optional, Tuple, Union

from viewer3d.base import ViewerBase, SetupBase
from viewer3d.config import BACKENDS
from viewer3d.converters import _to_poses_pdbstrings


_logger = logging.getLogger("viewer3d.core")


@attr.s(kw_only=True, slots=False)
class Py3DmolViewer(ViewerBase):
    _displayed = attr.ib(type=bool, default=False, init=False)

    def setup(self):
        self.py3Dmol = self._maybe_import_backend()
        self.viewer = self.py3Dmol.view(
            width=self.window_size[0],
            height=self.window_size[1],
        )

    def add_objects(self, _poses, _pdbstrings):
        for _model in range(len(_poses)):
            _pose = _poses[_model]
            if _pose is not None:
                self.viewer.addModels(io.to_pdbstring(_pose), "pdb")
            else:
                _pdbstring = _pdbstrings[_model]
                self.viewer.addModels(_pdbstring, "pdb")

    def remove_objects(self):
        self.viewer.removeAllLabels()
        self.viewer.removeAllModels()
        self.viewer.removeAllShapes()
        self.viewer.removeAllSurfaces()

    def update_viewer(self, _poses, _pdbstrings):
        self.update_objects(_poses, _pdbstrings)
        if self._displayed:
            self.viewer.update()

    def show_viewer(self):
        self.viewer.show()
        self._displayed = True


@attr.s(kw_only=True, slots=False)
class NGLviewViewer(ViewerBase):
    def setup(self):
        self.nglview = self._maybe_import_backend()
        self.viewer = self.nglview.widget.NGLWidget()

    def add_objects(self, _pose=None, _pdbstring=None):
        for _model in range(len(_pose)):
            p = _pose[_model]
            if p is not None:
                structure = self.nglview.adaptor.RosettaStructure(p)
            else:
                s = _pdbstring[_model]
                structure = self.nglview.adaptor.TextStructure(s, ext="pdb")
            self.viewer.add_structure(structure)

    def remove_objects(self):
        self.viewer.clear()
        for component_id in self.viewer._ngl_component_ids:
            self.viewer.remove_component(component_id)

    def update_viewer(self, _pose=None, _pdbstring=None):
        self.update_objects(_pose=_pose, _pdbstring=_pdbstring)

    def show_viewer(self):
        self.viewer.display(gui=True, style="ngl")
        self.viewer._ipython_display_()


@attr.s(kw_only=True, slots=False)
class PyMOLViewer(ViewerBase):
    def __attrs_pre_init__(self):
        self.pymol = self._maybe_import_backend()

        raise NotImplementedError(
            f"{self.__class__.__name__} is not currently supported."
        )

    def update_viewer(self, _pose, _pdbstring):
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
class SetupViewer(SetupBase):
    packed_and_poses_and_pdbs = attr.ib(
        type=Union[PackedPose, Pose, Iterable[Union[PackedPose, Pose]], None],
        default=None,
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
            widgets=self.widgets,
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
