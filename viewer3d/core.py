import attr
import collections
import logging
import os
import pyrosetta.distributed.io as io

from ipywidgets.widgets import Widget
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Iterable, List, Optional, Tuple, Union

from viewer3d.base import ViewerBase, Base3D
from viewer3d.config import BACKENDS
from viewer3d.converters import _to_poses_pdbstrings


_logger = logging.getLogger("viewer3d.core")


@attr.s(kw_only=True, slots=False)
class Py3DmolViewer(ViewerBase):
    def setup(self):
        self.py3Dmol = self._maybe_import_backend()
        self.viewer = self.py3Dmol.view(
            width=self.window_size[0],
            height=self.window_size[1],
        )

    def add_object(self, _poses, _pdbstrings, _model):
        _pose = _poses[_model]
        if _pose is not None:
            _pdbstring = io.to_pdbstring(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
        self.viewer.addModel(_pdbstring, "pdb")
        self.apply_modules(_pose, _pdbstring, _model)
        if self._displayed:
            self.viewer.update()

    def add_objects(self, _poses, _pdbstrings, _model):
        if _model is None:
            for _m in range(len(_poses)):
                self.add_object(_poses, _pdbstrings, _m)
        elif isinstance(_model, int):
            self.add_object(_poses, _pdbstrings, _model)

    def remove_objects(self, _model):
        self.viewer.removeAllShapes()
        self.viewer.removeAllSurfaces()
        self.viewer.removeAllLabels()
        if _model is None:
            self.viewer.removeAllModels()
        elif isinstance(_model, int):
            self.viewer.removeModel(_model)
        if self._displayed:
            self.viewer.update()

    def show_viewer(self):
        self.viewer.show()


@attr.s(kw_only=True, slots=False)
class NGLviewViewer(ViewerBase):
    def setup(self):
        self.nglview = self._maybe_import_backend()
        self.viewer = self.nglview.widget.NGLWidget()

    def set_window_size(self):
        """Resize the NGLWidget window."""
        self.viewer._remote_call(
            "setSize",
            targe="Widget",
            args=[f"{self.window_size[0]}px", f"{self.window_size[1]}px"],
        )

    def add_object(self, _poses, _pdbstrings, _model):
        _pose = _poses[_model]
        if _pose is not None:
            structure = self.nglview.adaptor.RosettaStructure(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
            structure = self.nglview.adaptor.TextStructure(_pdbstring, ext="pdb")
        self.viewer._load_data(structure)
        self.viewer._ngl_component_ids.append(structure.id)
        self.viewer._update_component_auto_completion()

    def add_objects(self, _poses, _pdbstrings, _model):
        _model_range = range(len(_poses))
        if _model is None:
            for _m in _model_range:
                self.add_object(_poses, _pdbstrings, _m)
            for _m in _model_range:
                self.apply_modules(_poses[_m], _pdbstrings[_m], _m)
        elif isinstance(_model, int):
            self.add_object(_poses, _pdbstrings, _model)
            self.apply_modules(_poses[_model], _pdbstrings[_model], _model)

    def remove_objects(self, _model):
        component_ids = self.viewer._ngl_component_ids
        if _model is None:
            for component_id in component_ids:
                component_index = component_ids.index(component_id)
                self.viewer.remove_component(component_id)
                self.viewer.clear(component=component_index)
        elif isinstance(_model, int):
            for component_id in component_ids:
                component_index = component_ids.index(component_id)
                if component_index == _model:
                    self.viewer.remove_component(component_id)
                    self.viewer.clear(component=component_index)
                    break

    def show_viewer(self):
        self.viewer.display(gui=self.gui, style="ngl")
        self.set_window_size()
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
class SetupViewer(Base3D):
    packed_and_poses_and_pdbs = attr.ib(
        type=Optional[Union[PackedPose, Pose, Iterable[Union[PackedPose, Pose]]]],
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
            n_decoys=self.n_decoys,
            widgets=self.widgets,
            auto_show=self.auto_show,
            backend=self.backend,
            gui=self.gui,
        )

    def initialize_viewer(self):
        if self.backend == BACKENDS[0]:
            if self.gui:
                _logger.debug(f"GUI is not supported for `{self.backend}` backend.")
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
    gui=None,
    auto_show=None,
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
        gui=gui,
        auto_show=auto_show,
    ).initialize_viewer()

    return viewer
