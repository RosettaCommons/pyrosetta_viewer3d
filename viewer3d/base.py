import attr
import logging
import time

from pyrosetta import Pose
from typing import List, NoReturn, Optional

from viewer3d.base3d import Base3D, expand_notebook
from viewer3d.initialization import InitBase, BACKENDS
from viewer3d.pose import PoseBase
from viewer3d.widgets import WidgetsBase


_logger: logging.Logger = logging.getLogger("viewer3d.base")


@attr.s(kw_only=False, slots=False)
class ViewerBase(Base3D, InitBase, PoseBase, WidgetsBase):
    _displayed = attr.ib(type=bool, default=False, init=False)

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
        """Setup Viewer in Jupyter notebook. """
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
                self.decoy_widget.children[0].value = index
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

    def show(self) -> None:
        """Display Viewer in Jupyter notebook."""
        if self._in_notebook():
            if self.backend == BACKENDS[2]:
                self.show_viewer()
                self.display_widgets()
                # self.update_viewer()  # TODO: display decoy widget for indices
            else:
                self._clear_output()
                self._toggle_window(self.window_size)
                self.display_widgets()
                self.show_viewer()
                self._toggle_scrolling()
            self._displayed = True
