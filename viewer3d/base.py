import attr
import math
import logging
import sys

from ipywidgets import interact, IntSlider

try:
    from IPython.core.display import display, HTML
    from IPython.display import clear_output
except ImportError:
    _logger.error("IPython.core.display or IPython.display module cannot be imported.")
from typing import Generic, Tuple, TypeVar

from viewer3d.config import _import_backend
from viewer3d.exceptions import ViewerImportError


_logger = logging.getLogger("viewer3d.base")


@attr.s(kw_only=False, slots=False, frozen=False)
class ViewerBase:
    def __attrs_pre_init__(self):
        self._toggle_scrolling()

    def _maybe_import_backend(self):
        if self.backend not in sys.modules:
            try:
                _import_backend(self.backend)
            except ImportError:
                raise ViewerImportError(self.backend)

        return sys.modules[self.backend]

    def view(self, i=0):
        self.update(self.poses[i], self.pdbstrings[i])

    def setup(self):
        if self.n_decoys > 1:
            s_widget = IntSlider(
                min=0,
                max=self.n_decoys - 1,
                description="Decoys",
                continuous_update=self.continuous_update,
            )
            interact(self.view, i=s_widget)
        else:
            self.view()

    def __add__(self, other):
        self.modules += [other]

        return self

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __call__(self):
        return self.show()

    def _clear_output(self):
        try:
            _logger.debug("IPython.display clearing Jupyter notebook cell output.")
            clear_output(wait=True)
        except NameError as e:
            _logger.debug(e)

    def _toggle_scrolling(self):
        try:
            _logger.debug(
                "IPython.core.display toggling scrolling in Jupyter notebook cell."
            )
            display(
                HTML(
                    "<script>$('.output_scroll').removeClass('output_scroll')</script>"
                )
            )
        except NameError as e:
            _logger.debug(e)

    def _toggle_window(self, _window_size):
        try:
            _logger.debug(
                "IPython.core.display toggling cell window area in Jupyter notebook."
            )
            HTML(
                """<style>
                    .output_wrapper, .output {
                        height:auto !important;
                        max-height:%ipx;
                    }
                    .output_scroll {
                        box-shadow:none !important;
                        webkit-box-shadow:none !important;
                    }
                    </style>
            """
                % math.ceil(_window_size[1])
            )
        except NameError as e:
            _logger.debug(e)

    def add(self, other):
        """Add a module to the Viewer instance."""
        return self.__add__(other)

    def reinit(self):
        """Subtract all modules from the Viewer instance."""
        self.modules = []

    def reset(self):
        """Delete Viewer instance attributes."""
        self.poses = None
        self.pdbstrings = None
        self.window_size = None
        self.modules = None
        self.delay = None
        self.continuous_update = None


def expand_notebook():
    """Expand Jupyter notebook cell to maximum width."""
    try:
        _logger.debug("IPython.core.display expanding Jupyter notebook cell width.")
        display(HTML("<style>.container { width:100% !important; }</style>"))
    except NameError:
        _logger.exception("IPython.core.display module not imported.")
