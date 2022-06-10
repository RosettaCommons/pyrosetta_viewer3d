import attr
import collections
import math
import logging
import sys
import time

from ipywidgets import interactive, IntSlider
from ipywidgets.widgets import Widget
from IPython.display import display
from IPython.core.display import display as core_display, HTML
from IPython.display import clear_output
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import Generic, Iterable, List, Optional, Tuple, TypeVar, Union

from viewer3d.config import _import_backend, BACKENDS
from viewer3d.converters import _to_float, _to_widgets
from viewer3d.exceptions import ViewerImportError
from viewer3d.modules import ModuleBase
from viewer3d.tracer import silence_tracer
from viewer3d.validators import (
    _validate_add_pose,
    _validate_int_float,
    _validate_window_size,
)


_logger = logging.getLogger("viewer3d.base")


@attr.s(kw_only=False, slots=False)
class Base3D:
    poses = attr.ib(type=Iterable[Pose], default=None)
    pdbstrings = attr.ib(type=Iterable[str], default=None)
    n_decoys = attr.ib(type=int, default=None)
    window_size = attr.ib(
        type=Tuple[Union[int, float], Union[int, float]],
        default=None,
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of((int, float)),
                iterable_validator=attr.validators.instance_of(
                    collections.abc.Iterable
                ),
            ),
            _validate_window_size,
        ],
        converter=attr.converters.default_if_none(default=(1200, 800)),
    )
    modules = attr.ib(
        type=list,
        default=None,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(ModuleBase),
            iterable_validator=attr.validators.instance_of(list),
        ),
        converter=attr.converters.default_if_none(default=[]),
    )
    delay = attr.ib(
        type=float,
        default=None,
        validator=_validate_int_float,
        converter=attr.converters.pipe(
            attr.converters.default_if_none(0.25), _to_float
        ),
    )
    continuous_update = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    widgets = attr.ib(
        type=Optional[List[Widget]],
        default=None,
        validator=attr.validators.optional(
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of(Widget),
                iterable_validator=attr.validators.instance_of(list),
            )
        ),
        converter=_to_widgets,
    )
    backend = attr.ib(
        type=str,
        default=None,
        validator=[attr.validators.instance_of(str), attr.validators.in_(BACKENDS)],
        converter=attr.converters.default_if_none(default=BACKENDS[0]),
    )

    def _maybe_import_backend(self):
        if self.backend not in sys.modules:
            try:
                _import_backend(self.backend)
            except ImportError:
                raise ViewerImportError(self.backend)

        return sys.modules[self.backend]

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
            core_display(
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

    def _in_notebook(self):
        try:
            get_ipython()
            _in_notebook = True
        except:
            _in_notebook = False
        finally:
            return _in_notebook

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


@attr.s(kw_only=False, slots=False)
class WidgetsBase:
    decoy_widget = attr.ib(
        default=attr.Factory(
            lambda self: interactive(
                self.update_decoy,
                index=IntSlider(
                    min=0,
                    max=self.n_decoys - 1,
                    description="Decoys",
                    continuous_update=self.continuous_update,
                ),
            ),
            takes_self=True,
        )
    )

    def get_widgets(self):
        _widgets = self.widgets.copy()
        if self.n_decoys > 1:
            _widgets.insert(0, self.decoy_widget)
        return _widgets

    def set_widgets(self, obj):
        self.widgets = _to_widgets(obj)


@attr.s(kw_only=False, slots=False)
class ViewerBase(Base3D, WidgetsBase):
    def __attrs_post_init__(self):
        self.setup()

    @_validate_add_pose
    def add_pose(self, pose=None, index=None):
        if index is None:
            kwargs = self.decoy_widget.kwargs
            index = kwargs["index"] if kwargs else 0
        self.poses[index].append(pose)
        self.update_viewer(self.poses[index])

    @silence_tracer
    def apply_modules(self, _pose=None, _pdbstring=None):
        for _model in range(len(_pose)):
            for _module in self.modules:
                func = getattr(_module, f"apply_{self.backend}")
                self.viewer = func(
                    self.viewer,
                    _pose[_model],
                    _pdbstring[_model],
                    _model,
                )

    def update_objects(self, _pose=None, _pdbstring=None):
        """Setup Viewer in Jupyter notebook."""
        self.remove_objects()
        self.add_objects(_pose=_pose, _pdbstring=_pdbstring)
        self.apply_modules(_pose=_pose, _pdbstring=_pdbstring)

    def update_decoy(self, index=0):
        time.sleep(self.delay)
        self.update_viewer(self.poses[index], self.pdbstrings[index])

    def show(self):
        """Display Viewer in Jupyter notebook."""
        if self._in_notebook():
            self._toggle_scrolling()
            self._toggle_window(self.window_size)
            display(*self.get_widgets())
            self.show_viewer()


def expand_notebook():
    """Expand Jupyter notebook cell to maximum width."""
    try:
        _logger.debug("IPython.core.display expanding Jupyter notebook cell width.")
        core_display(HTML("<style>.container { width:100% !important; }</style>"))
    except NameError:
        _logger.exception("IPython.core.display module not imported.")
