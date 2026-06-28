__author__ = "Jason C. Klima"

import attr
import math
import logging
import os
import sys

try:
    from IPython.core.display import display as core_display, HTML
except ImportError:
    from IPython.display import display as core_display, HTML
from IPython.display import clear_output

from viewer3d.config import _import_backend
from viewer3d.exceptions import ViewerImportError
from viewer3d.modules.base import ModuleBase
from viewer3d.type_defs import (
    Any,
    List,
    ModuleType,
    NoReturn,
    Self,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.base3d")


@attr.s(kw_only=False, slots=False, frozen=False)
class Base3D:
    def _maybe_import_backend(self) -> ModuleType:
        if self.backend not in sys.modules:
            try:
                _import_backend(self.backend)
            except ImportError:
                raise ViewerImportError(self.backend)

        return sys.modules[self.backend]

    def __add__(self, other: ModuleBase) -> Self:
        self.modules += [other]
        self._invalidate_viewer_state()
        return self

    def __radd__(self, other: ModuleBase) -> Self:
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __call__(self) -> None:
        self.show()

    def _clear_output(self) -> None:
        try:
            _logger.debug("IPython.display clearing Jupyter notebook cell output.")
            clear_output(wait=True)
        except NameError as e:
            _logger.debug(e)

    def _toggle_scrolling(self) -> None:
        try:
            _logger.debug("IPython.core.display toggling scrolling in Jupyter notebook cell.")
            core_display(HTML("<script>$('.output_scroll').removeClass('output_scroll')</script>"))
        except NameError as e:
            _logger.debug(e)

    def _toggle_window(self, _window_size: Tuple[int, int]) -> None:
        try:
            _logger.debug("IPython.core.display toggling cell window area in Jupyter notebook.")
            HTML("""<style>
                    .output_wrapper, .output {
                        height:auto !important;
                        max-height:%ipx;
                    }
                    .output_scroll {
                        box-shadow:none !important;
                        webkit-box-shadow:none !important;
                    }
                    </style>
            """ % math.ceil(_window_size[1]))
        except NameError as e:
            _logger.debug(e)

    @staticmethod
    def _in_notebook() -> bool:
        if os.environ.get("PAPERMILL_EXECUTION"):  # Non-interactive state
            return False
        else:
            try:
                get_ipython()
                return True
            except:
                return False

    def add(self, other: ModuleBase) -> Self:
        """Add a module to the Viewer instance."""
        return self.__add__(other)

    def reinit(self) -> None:
        """
        Subtract all modules from the Viewer instance. Alias of the `clear_modules` method.
        """
        _logger.warning("The 'reinit' method is deprecated. Please use the `clear` method instead.")
        self.clear_modules()

    def clear(self) -> None:
        """
        Subtract all modules from the Viewer instance. Alias of the `clear_modules` method.
        """
        self.clear_modules()

    def clear_modules(self) -> None:
        """Subtract all modules from the Viewer instance."""
        self.modules = []
        self._invalidate_viewer_state()

    def get_n_models(self, index: int) -> Union[NoReturn, int]:
        """Get the number of models at the poses index."""
        assert len(self.poses[index]) == len(
            self.pdbstrings[index]
        ), "Number of `Pose` objects and PDB `str` objects must be equal."
        return len(self.poses[index])

    def get_last_model(self, index: int) -> int:
        """Get last model index."""
        return self.get_n_models(index) - 1

    def get_modules(self) -> List[ModuleBase]:
        """
        Get a copy of the Viewer instance modules.

        Returns:
            A `list` object of `ModuleBase` objects.
        """
        return self.modules.copy()

    def set_modules(self, modules: Any) -> None:
        """
        Set the Viewer instance modules.

        Args:
            modules: a Viewer module or iterable of Viewer instance modules.

        Returns:
            `None`
        """
        self.modules = ModuleBase._to_modules(modules)
        self._invalidate_viewer_state()

    def reset(self) -> None:
        """Delete Viewer instance attributes."""
        self.poses = None
        self.pdbstrings = None
        self.window_size = None
        self.modules = []
        self.delay = None
        self.continuous_update = None
        self.widgets = []
        self._invalidate_viewer_state()


def expand_notebook() -> None:
    """Expand a Jupyter notebook cell to maximum width."""
    if Base3D._in_notebook():
        try:
            _logger.debug("IPython.core.display expanding Jupyter notebook cell width.")
            core_display(HTML("<style>.container { width:100% !important; }</style>"))
        except NameError:
            _logger.exception("IPython.core.display module not imported.")
