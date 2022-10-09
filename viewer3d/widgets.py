import attr
import logging
import sys

from ipywidgets import interactive, Image, IntSlider, Widget
from IPython.display import display
from typing import Any, Dict, List, NoReturn, Union

from viewer3d.converters import _to_widgets
from viewer3d.config import BACKENDS, COLORBAR_ATTR

_logger: logging.Logger = logging.getLogger("viewer3d.widgets")


@attr.s(kw_only=False, slots=False)
class WidgetsBase:
    def get_decoy_widget(self) -> Widget:
        # self.update_decoy(index=0) # TODO double check if needed
        return interactive(
            self.update_decoy,
            index=IntSlider(
                min=0,
                max=self.get_n_decoys() - 1,
                value=0,
                description="Decoys",
                continuous_update=self.continuous_update,
            ),
        )

    def get_widgets(self) -> List[Widget]:
        self.decoy_widget = self.get_decoy_widget()
        _widgets = self.widgets.copy()
        if self.backend != BACKENDS[2] and hasattr(self.viewer, COLORBAR_ATTR):
            _value = getattr(self.viewer, COLORBAR_ATTR)
            _image = Image(value=_value)
            _widgets.insert(0, _image)
        if self.get_n_decoys() > 1:
            _widgets.insert(0, self.decoy_widget)

        return _widgets

    def get_n_decoys(self) -> Union[int, NoReturn]:
        assert len(self.poses.keys()) == len(
            self.pdbstrings.keys()
        ), "Number of `Pose` objects and PDB `str` objects must be equal."
        return len(self.poses.keys())

    def update_decoy(self, index: int) -> None:
        self.update_viewer(index=index)

    def get_decoy_widget_index(self) -> int:
        if hasattr(self, "decoy_widget"):
            kwargs = self.decoy_widget.kwargs
        else:
            kwargs = {}
        index = kwargs["index"] if kwargs else 0
        return index

    def get_widgets_dict(self) -> Dict[str, Any]:
        """
        Return the widgets dictionary of descriptions and values.

        Returns:
            A `dict` object with each widget's "description" attribute as a key
                and "value" attribute as a value if they have been set.
        """
        return {
            widget.description: widget.value
            for widget in self.widgets
            if all(hasattr(widget, attr) for attr in ("description", "value"))
        }

    def set_widgets(self, widgets: Any) -> None:
        """
        Set the Viewer instance widgets.

        Args:
            widgets: a `Widget` object or an iterable of `Widget` objects.

        Returns:
            `None`
        """
        self.widgets = _to_widgets(widgets)

    def display_widgets(self) -> None:
        widgets = self.get_widgets()
        if widgets:
            display(*widgets)

    def _maybe_setup_colab(self) -> None:
        if "google.colab" in sys.modules:
            sys.modules["google.colab"].output.enable_custom_widget_manager()
