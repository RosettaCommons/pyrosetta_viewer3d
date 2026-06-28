__author__ = "Jason C. Klima"

import attr
import logging
import time

from pyrosetta import Pose

from viewer3d.base3d import Base3D
from viewer3d.initialization import InitBase, BACKENDS
from viewer3d.pose import PoseBase
from viewer3d.type_defs import (
    List,
    NoReturn,
    Optional,
)
from viewer3d.widgets import WidgetsBase

_logger: logging.Logger = logging.getLogger("viewer3d.base")


@attr.define(kw_only=False, slots=False, frozen=False)
class Viewer(Base3D, InitBase, PoseBase, WidgetsBase):
    _displayed: bool = attr.field(default=False, init=False)

    def __attrs_post_init__(self) -> None:
        """Post-initialization setup."""
        self.setup()
        self._maybe_setup_colab()
        if self.auto_show:
            self.show()

    def apply_setZoomTo(self, _pose: Pose, _pdbstring: str, _model: int) -> None:
        """Apply `setZoomTo` module to model."""
        func = getattr(self._setZoomTo, f"apply_{self.backend}")
        self.viewer = func(
            self.viewer,
            _pose,
            _pdbstring,
            _model,
        )

    def apply_modules(self, _pose: Pose, _pdbstring: str, _model: int) -> None:
        """Apply visualization modules to model."""
        if not self._displayed:
            self.apply_setZoomTo(_pose, _pdbstring, _model)
        for _module in self.modules:
            func = getattr(_module, f"apply_{self.backend}")
            self.viewer = func(self.viewer, _pose, _pdbstring, _model)

    def update_objects(
        self,
        _poses: List[Pose],
        _pdbstrings: List[str],
        _model: Optional[int],
        _add_objects: bool,
        _remove_objects: bool,
    ):
        """Setup Viewer in Jupyter notebook."""
        assert len(_poses) == len(
            _pdbstrings
        ), "Number of `Pose` objects and PDB `str` objects must be equal."
        if self.delay:
            time.sleep(self.delay)
        if _remove_objects and _add_objects:
            self.set_objects(_poses, _pdbstrings, _model)
        elif _remove_objects and not _add_objects:
            self.remove_objects(_model)
        elif _add_objects and not _remove_objects:
            self.add_objects(_poses, _pdbstrings, _model)
        if self._displayed:
            self.update()

    def update_viewer(
        self,
        index: Optional[int] = None,
        model: Optional[int] = None,
        add_objects: bool = True,
        remove_objects: bool = True,
    ) -> Optional[NoReturn]:
        """
        Update Viewer in Jupyter notebook.

        Args:
            index: an optional `int` object representing the poses or pdbstrings index
                to update. If `None`, then update all poses or pdbstrings in the displayed
                index.
                Default: `None`
            model: an optional `int` object representing the model in the poses or pdbstrings
                index to update. If `None`, then update all models in the poses or pdbstrings
                index to update.
                Default: `None`
            add_objects: an optional `bool` object. If `True`, then add objects to the viewer.
                If `False`, then do not add objects to the viewer.
                Default: `True`
            remove_objects: an optional `bool` object. If `True`, then remove objects from the
                viewer. If `False`, then do not remove objects from the viewer.
                Default: `True`

        Raises:
            `IndexError` if index does not exist.

        Returns:
            `None`
        """
        if index is None:
            index = self.get_decoy_widget_index()
        if index in self.poses.keys():
            if hasattr(self, "decoy_widget"):
                if self._current_index != index:
                    self.decoy_widget.children[0].value = index
                    return
            self.update_objects(
                self.poses[index],
                self.pdbstrings[index],
                model,
                add_objects,
                remove_objects,
            )
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )

    def _initialize_viewer_once(self) -> None:
        if not getattr(self, "_initialized_once", False):
            self._initialized_once = True
            self._current_index = None  # Ensure first call runs
            self.update_decoy(index=0)

    def _invalidate_viewer_state(self):
        if hasattr(self, "_initialized_once"):
            self._initialized_once = False
        if hasattr(self, "_current_index"):
            self._current_index = None

    def show(self, force=False) -> None:
        """Display Viewer in Jupyter notebook."""
        if Base3D._in_notebook() or force:
            if self.backend == BACKENDS[2]:
                self.show_viewer()
                self._initialize_viewer_once()
                self._clear_output()
                self.display_widgets()
            else:
                self._initialize_viewer_once()
                self._clear_output()
                self._toggle_window(self.window_size)
                self.display_widgets()
                self.show_viewer()
            self._toggle_scrolling()
            self._displayed: bool = True
