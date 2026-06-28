__author__ = "Jason C. Klima"

import attr
import logging
import pyrosetta.distributed.io as io
import sys
import time

from IPython.display import (
    Javascript,
    display,
)
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose

from viewer3d.base import Viewer
from viewer3d.config import BACKENDS
from viewer3d.converters import _to_poses_pdbstrings
from viewer3d.initialization import InitBase
from viewer3d.modules.base import ModuleBase
from viewer3d.type_defs import (
    DefaultDict,
    GenericViewer,
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.core")


@attr.define(kw_only=True, slots=False, frozen=False)
class Py3DmolViewer(Viewer):
    """Viewer object for the `py3Dmol` backend."""

    def setup(self) -> None:
        self.py3Dmol = self._maybe_import_backend()
        self.viewer = self.py3Dmol.view(
            width=self.window_size[0],
            height=self.window_size[1],
        )

    def add_object(self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]) -> None:
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

    def png(self, **kwargs) -> None:
        if self._in_notebook() and self._displayed:
            self.viewer.png()
        else:
            _logger.error("The `Viewer` object must be displayed in a Jupyter notebook.")

    def save_png(self, filename: str, **kwargs) -> None:
        if self._in_notebook() and self._displayed:
            div_id = f"3dmolviewer_{self.viewer.uniqueid}"
            display(Javascript(f"""
            const container = document.getElementById("{div_id}");
            if (!container) {{
                console.log("Container not found");
            }} else {{
                const canvas = container.querySelector("canvas");
                if (!canvas) {{
                    console.log("Canvas not found");
                }} else {{
                    const dataURL = canvas.toDataURL("image/png");

                    const a = document.createElement("a");
                    a.href = dataURL;
                    a.download = "{filename}";
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);
                }}
            }}
            """))
            _logger.info(f"Saved PNG file: '{filename}'")
        else:
            _logger.error(
                "Cannot save PNG file: the `Viewer` object must be displayed in a Jupyter notebook."
            )


@attr.define(kw_only=True, slots=False, frozen=False)
class NGLViewViewer(Viewer):
    """Viewer object for the `nglview` backend."""

    def setup(self) -> None:
        self.nglview = self._maybe_import_backend()
        self.viewer = self.nglview.widget.NGLWidget()
        self._component_ids: List[str] = []

    def set_window_size(self) -> None:
        """Resize the NGLWidget window."""
        self.viewer._remote_call(
            "setSize",
            targe="Widget",
            args=[f"{self.window_size[0]}px", f"{self.window_size[1]}px"],
        )

    def add_object(self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]) -> None:
        _pose = _poses[_model]
        if _pose is not None:
            structure = self.nglview.adaptor.RosettaStructure(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
            structure = self.nglview.adaptor.TextStructure(_pdbstring, ext="pdb")
        component = self.viewer.add_component(
            structure,
            default_representation=False,
        )
        if _model not in set(range(len(self._component_ids))):
            self._component_ids.append(component.id)
        else:
            self._component_ids.insert(_model, component.id)

    def add_model(self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]) -> None:
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

    def remove_objects(self, _model: Optional[int]) -> None:
        if _model is None:
            for component_id in self._component_ids:
                self.viewer.remove_component(component_id)
            self._component_ids.clear()
        else:
            component_id = self._component_ids.pop(_model)
            self.viewer.remove_component(component_id)

    def set_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        self.remove_objects(_model)
        self.add_objects(_poses, _pdbstrings, _model)

    def update(self) -> None:
        pass

    def show_viewer(self) -> None:
        self.viewer.display(gui=self.gui, style="ngl")
        self.set_window_size()
        self.viewer._ipython_display_()

    def png(
        self,
        frame: Optional[int] = None,
        factor: int = 4,
        antialias: bool = True,
        trim: bool = False,
        transparent: bool = False,
        **kwargs,
    ) -> None:
        if self._in_notebook() and self._displayed:
            image = self.viewer.render_image(
                frame=frame,
                factor=factor,
                antialias=antialias,
                trim=trim,
                transparent=transparent,
            )
            display(image)
        else:
            _logger.error("The `Viewer` object must be displayed in a Jupyter notebook.")

    def save_png(
        self,
        filename: str,
        factor: int = 4,
        antialias: bool = True,
        trim: bool = False,
        transparent: bool = False,
        **kwargs,
    ) -> None:
        if self._in_notebook() and self._displayed:
            self.viewer.download_image(
                filename=filename,
                factor=factor,
                antialias=antialias,
                trim=trim,
                transparent=transparent,
            )
            _logger.info(f"Saved PNG file: '{filename}'")
        else:
            _logger.error(
                "Cannot save PNG file: the `Viewer` object must be displayed in a Jupyter notebook."
            )


@attr.define(kw_only=True, slots=False, frozen=False)
class PyMOLViewer(Viewer):
    """Viewer object for the `pymol` backend."""

    def wait_for_gui(self) -> None:
        timeout = 10.0  # seconds
        start_time = time.time()
        elapsed_time = time.time() - start_time
        while elapsed_time < timeout:
            try:
                self.viewer.test()
                self.viewer.delete("all")
                break
            except:
                elapsed_time = time.time() - start_time
        else:
            raise RuntimeError("Launching PyMOL GUI exceeded timeout!")

    def launch_pymol(self) -> None:
        command_line = " ".join(
            [
                "pymol",
                "-R",
                "-k",
                "-q",
                f"-W {int(self.window_size[0])}",
                f"-H {int(self.window_size[1])}",
            ]
        )
        sys.modules["subprocess"].Popen(command_line, shell=True)
        self.viewer = sys.modules["xmlrpc.client"].ServerProxy("http://localhost:9123")
        self.wait_for_gui()

    def setup(self) -> None:
        self.pymol = self._maybe_import_backend()

    def add_object(self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]) -> None:
        _pose = _poses[_model]
        if _pose is not None:
            _pdbstring = io.to_pdbstring(_pose)
        else:
            _pdbstring = _pdbstrings[_model]
        self.viewer.load_raw(_pdbstring, "pdb", str(_model))
        self.apply_modules(_pose, _pdbstring, _model)

    def add_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        if hasattr(self, "viewer"):
            if _model is None:
                for _m in range(len(_poses)):
                    self.add_object(_poses, _pdbstrings, _m)
            elif isinstance(_model, int):
                self.add_object(_poses, _pdbstrings, _model)

    def remove_objects(self, _model: Optional[int]) -> None:
        if hasattr(self, "viewer"):
            if _model is None:
                self.viewer.delete("all")
            elif isinstance(_model, int):
                self.viewer.delete(str(_model))

    def set_objects(
        self, _poses: List[Pose], _pdbstrings: List[str], _model: Optional[int]
    ) -> None:
        self.remove_objects(_model)
        self.add_objects(_poses, _pdbstrings, _model)

    def update(self) -> None:
        pass

    def show_viewer(self):
        """Display PyMOLViewer."""
        self.launch_pymol()

    def png(
        self,
        width: int = 4000,
        height: int = 3000,
        ray: bool = True,
        ray_trace_mode: int = 0,
        ray_shadows: bool = False,
        opaque_background: bool = True,
        **kwargs,
    ) -> None:
        if self._in_notebook() and self._displayed:
            self.viewer.set("ray_trace_mode", ray_trace_mode)
            self.viewer.set("ray_shadows", int(ray_shadows))
            self.viewer.set("opaque_background", int(opaque_background))
            if ray:
                self.viewer.ray(width, height)
        else:
            _logger.error("The `Viewer` object must be displayed.")

    def save_png(
        self,
        filename: str,
        width: int = 2000,
        height: int = 1000,
        dpi: int = 300,
        ray: bool = True,
        ray_trace_mode: int = 0,
        ray_shadows: bool = False,
        opaque_background: bool = True,
        **kwargs,
    ) -> None:
        if self._in_notebook() and self._displayed:
            self.viewer.set("ray_trace_mode", ray_trace_mode)
            self.viewer.set("ray_shadows", int(ray_shadows))
            self.viewer.set("opaque_background", int(opaque_background))
            self.viewer.png(filename, width, height, dpi, int(ray))
            _logger.info(f"Saved PNG file: '{filename}'")
        else:
            _logger.error("Cannot save PNG file: the `Viewer` object must be displayed.")


@attr.define(kw_only=True, slots=False, frozen=False)
class SetupViewer(InitBase):
    """Initialize a `Viewer` object with the user-provided arguments."""

    packed_and_poses_and_pdbs: Optional[
        Union[PackedPose, Pose, str, Iterable[Union[PackedPose, Pose, str]]]
    ] = attr.field(default=None)

    def __attrs_post_init__(self) -> None:
        self.poses: DefaultDict[int, List[Optional[Pose]]]
        self.pdbstrings: DefaultDict[int, List[Optional[str]]]
        self.poses, self.pdbstrings = _to_poses_pdbstrings(self.packed_and_poses_and_pdbs)
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

    def initialize_viewer(self) -> GenericViewer:
        if self.backend == BACKENDS[0]:
            if self.gui:
                _logger.info(f"GUI is not supported for `{self.backend}` backend.")
            viewer = Py3DmolViewer(**self.viewer_kwargs)
        elif self.backend == BACKENDS[1]:
            viewer = NGLViewViewer(**self.viewer_kwargs)
        elif self.backend == BACKENDS[2]:
            if self.gui:
                _logger.info(f"GUI is not supported for `{self.backend}` backend.")
            viewer = PyMOLViewer(**self.viewer_kwargs)

        return viewer


def init(
    packed_and_poses_and_pdbs: Optional[
        Union[PackedPose, Pose, str, Iterable[Union[PackedPose, Pose, str]]]
    ] = None,
    window_size: Optional[Tuple[Union[int, float], Union[int, float]]] = (1200, 800),
    modules: Optional[List[ModuleBase]] = [],
    delay: Optional[Union[int, float]] = 0.0,
    continuous_update: Optional[bool] = False,
    backend: Optional[Union[int, str]] = 0,
    gui: Optional[bool] = False,
    auto_show: Optional[bool] = False,
) -> GenericViewer:
    """
    Initialize a `Py3DmolViewer`, `NGLViewViewer`, or a `PyMOLViewer` object.

    Args:
        packed_and_poses_and_pdbs: An optional `PackedPose`, `Pose`, a `str` of a PDB string or
            a filesystem path to a ".pdb" file, or an iterable of these objects.

        window_size: An optional `list` or `tuple` of `int` or `float` values for the
            `(width, height)` dimensions of the displayed window screen size.

        modules: An optional `list` of instantiated visualization modules to apply upon
            changing amongst `Pose` or PDB string `str` objects with the interactive slider
            widget, matching `viewer3d.set*` object names.

        delay: An optional `float` or `int` time delay in seconds before rendering the
            visualization in a Jupyter notebook, which is useful to prevent overburdening the
            Jupyter notebook client if `for` looping over quick modifications to a `Pose`, and
            must be `>=0`.

        continuous_update: A `bool` object. When using the interactive slider widget, `False`
            restricts rendering to mouse button release events.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization. The currently supported backends are `py3Dmol` (`0`), `nglview`
            (`1`), and `pymol` (`2`).

        gui: a `bool` object only supported by the `nglview` backend. If `True`, then show
            the NGLView graphical user interface.

        auto_show: a `bool` object. If `True`, then automatically run the `show` method on the
            returned `*Viewer` instance upon calling this `init` function.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
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
