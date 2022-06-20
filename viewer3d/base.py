import attr
import logging
import time

from typing import Optional

from viewer3d.base3d import Base3D, expand_notebook
from viewer3d.initialization import InitBase
from viewer3d.pose import PoseBase
from viewer3d.widgets import WidgetsBase


_logger: logging.Logger = logging.getLogger("viewer3d.base")


@attr.s(kw_only=False, slots=False)
class ViewerBase(Base3D, InitBase, PoseBase, WidgetsBase):
    _displayed = attr.ib(type=bool, default=False, init=False)

    def __attrs_post_init__(self):
        self.setup()
        self._maybe_setup_colab()
        if self.auto_show:
            self.show()

    def apply_setZoomTo(self, _pose, _pdbstring, _model):
        func = getattr(self._setZoomTo, f"apply_{self.backend}")
        self.viewer = func(
            self.viewer,
            _pose,
            _pdbstring,
            _model,
        )

    def apply_modules(self, _pose, _pdbstring, _model):
        if not self._displayed:
            self.apply_setZoomTo(_pose, _pdbstring, _model)
        for _module in self.modules:
            func = getattr(_module, f"apply_{self.backend}")
            self.viewer = func(self.viewer, _pose, _pdbstring, _model)

    def update_objects(
        self,
        _poses,
        _pdbstrings,
        _model,
        _add_objects,
        _remove_objects,
    ):
        """
        Setup Viewer in Jupyter notebook.

        Args:
            _model: if `None`, then update all models. If `int`, then update a single model.
        """
        assert len(_poses) == len(
            _pdbstrings
        ), "Number of `Pose` objects and PDB `str` objects must be equal."
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
    ):
        if index is None:
            index = self.get_decoy_widget_index()
        if index in self.poses.keys():
            self.update_objects(
                self.poses[index],
                self.pdbstrings[index],
                model,
                add_objects,
                remove_objects,
            )
            if hasattr(self, "decoy_widget"):
                self.decoy_widget.children[0].value = index
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )

    def show(self):
        """Display Viewer in Jupyter notebook."""
        if self._in_notebook():
            self._clear_output()
            self._toggle_window(self.window_size)
            self.display_widgets()
            self.show_viewer()
            self._toggle_scrolling()
            self._displayed = True
