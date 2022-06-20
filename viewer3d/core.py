import attr
import logging
import pyrosetta.distributed.io as io

from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Generic, Iterable, List, Optional, TypeVar, Union

from viewer3d.base import ViewerBase
from viewer3d.initialization import InitBase
from viewer3d.config import BACKENDS
from viewer3d.converters import _to_poses_pdbstrings


_logger: logging.Logger = logging.getLogger("viewer3d.core")
ViewerType = TypeVar("ViewerType", bound=ViewerBase)


@attr.s(kw_only=True, slots=False)
class Py3DmolViewer(ViewerBase):
    def setup(self) -> None:
        self.py3Dmol = self._maybe_import_backend()
        self.viewer = self.py3Dmol.view(
            width=self.window_size[0],
            height=self.window_size[1],
        )

    def add_object(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        _pose = _poses[_model]
        if _pose is not None:
            _pdbstring = io.to_pdbstring(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
        self.viewer.addModel(_pdbstring, "pdb")
        self.apply_modules(_pose, _pdbstring, _model)

    def add_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        if _model is None:
            for _m in range(len(_poses)):
                self.add_object(_poses, _pdbstrings, _m)
        elif isinstance(_model, int):
            self.add_object(_poses, _pdbstrings, _model)

    def remove_objects(self, _model: Optional[int]) -> None:
        self.viewer.removeAllShapes()
        self.viewer.removeAllSurfaces()
        self.viewer.removeAllLabels()
        if _model is None:
            self.viewer.removeAllModels()
        elif isinstance(_model, int):
            self.viewer.removeModel(_model)

    def set_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        self.remove_objects(_model)
        self.add_objects(_poses, _pdbstrings, _model)

    def update(self) -> None:
        self.viewer.update()

    def show_viewer(self) -> None:
        self.viewer.show()


@attr.s(kw_only=True, slots=False)
class NGLViewViewer(ViewerBase):
    def setup(self) -> None:
        self.nglview = self._maybe_import_backend()
        self.viewer = self.nglview.widget.NGLWidget()

    def set_window_size(self) -> None:
        """Resize the NGLWidget window."""
        self.viewer._remote_call(
            "setSize",
            targe="Widget",
            args=[f"{self.window_size[0]}px", f"{self.window_size[1]}px"],
        )

    def add_object(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        _pose = _poses[_model]
        if _pose is not None:
            structure = self.nglview.adaptor.RosettaStructure(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
            structure = self.nglview.adaptor.TextStructure(_pdbstring, ext="pdb")
        self.viewer._load_data(structure)
        self.viewer._ngl_component_ids.append(structure.id)
        self.viewer._update_component_auto_completion()

    def add_model(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        if _model is None:
            for _m in range(len(_poses)):
                self.add_object(_poses, _pdbstrings, _m)
        elif isinstance(_model, int):
            self.add_object(_poses, _pdbstrings, _model)

    def apply_to_model(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        if _model is None:
            for _m in range(len(_poses)):
                self.apply_modules(_poses[_m], _pdbstrings[_m], _m)
        elif isinstance(_model, int):
            self.apply_modules(_poses[_model], _pdbstrings[_model], _model)

    def add_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        self.add_model(_poses, _pdbstrings, _model)
        self.apply_to_model(_poses, _pdbstrings, _model)

    def get_remove_component_ids(self, _model: Optional[int]) -> List[str]:
        component_ids = self.viewer._ngl_component_ids
        remove_component_ids = []
        if _model is None:
            for component_id in component_ids:
                remove_component_ids.append(component_id)
        elif isinstance(_model, int):
            for component_id in component_ids:
                component_index = component_ids.index(component_id)
                if component_index == _model:
                    remove_component_ids.append(component_id)
                    break

        return remove_component_ids

    def remove_objects(self, _model: Optional[int]) -> None:
        for component_id in self.get_remove_component_ids(_model):
            self.viewer.remove_component(component_id)

    def remove_component_ids(self, component_ids: List[str]) -> None:
        for component_id in component_ids:
            self.viewer.remove_component(component_id)

    def set_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        remove_component_ids = self.get_remove_component_ids(_model)
        self.add_model(_poses, _pdbstrings, _model)
        self.remove_component_ids(remove_component_ids)
        self.apply_to_model(_poses, _pdbstrings, _model)

    def update(self) -> None:
        pass

    def show_viewer(self) -> None:
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
class SetupViewer(InitBase):
    packed_and_poses_and_pdbs = attr.ib(
        type=Optional[Union[PackedPose, Pose, Iterable[Union[PackedPose, Pose]]]],
        default=None,
    )

    def __attrs_post_init__(self) -> None:
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
            widgets=self.widgets,
            auto_show=self.auto_show,
            backend=self.backend,
            gui=self.gui,
        )

    def initialize_viewer(self) -> Generic[ViewerType]:
        if self.backend == BACKENDS[0]:
            if self.gui:
                _logger.debug(f"GUI is not supported for `{self.backend}` backend.")
            viewer = Py3DmolViewer(**self.viewer_kwargs)
        elif self.backend == BACKENDS[1]:
            viewer = NGLViewViewer(**self.viewer_kwargs)
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
) -> Generic[ViewerType]:
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
