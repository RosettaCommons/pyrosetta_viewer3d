import attr
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
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from typing import Generic, List, Optional, Tuple, TypeVar, Union

from viewer3d.config import BACKENDS
from viewer3d.config import _import_backend
from viewer3d.converters import _to_widgets
from viewer3d.exceptions import ViewerImportError


_logger = logging.getLogger("viewer3d.base")


@attr.s(kw_only=False, slots=False, frozen=False)
class WidgetsBase:
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

    def set_widgets(self, obj):
        self.widgets = _to_widgets(obj)

    def get_widgets(self):
        _widgets = self.widgets.copy()
        if self.n_decoys > 1:
            _widgets.insert(0, self.decoy_widget)
        return _widgets


@attr.s(kw_only=False, slots=False, frozen=False)
class ViewerBase(WidgetsBase):
    poses = attr.ib(type=Pose, default=None)
    pdbstrings = attr.ib(type=PackedPose, default=None)
    n_decoys = attr.ib(type=int, default=None)
    window_size = attr.ib(type=Tuple[Union[int, float]], default=None)
    modules = attr.ib(type=list, default=None)
    delay = attr.ib(type=float, default=None)
    continuous_update = attr.ib(type=bool, default=None)
    backend = attr.ib(type=str, default=None)

    def __attrs_post_init__(self):
        self.setup()
        self.decoy_widget = interactive(
            self.update_decoy,
            index=IntSlider(
                min=0,
                max=self.n_decoys - 1,
                description="Decoys",
                continuous_update=self.continuous_update,
            ),
        )

    def _maybe_import_backend(self):
        if self.backend not in sys.modules:
            try:
                _import_backend(self.backend)
            except ImportError:
                raise ViewerImportError(self.backend)

        return sys.modules[self.backend]

    def add_pose(self, pose=None, new_chain=False):
        assert isinstance(
            pose, Pose
        ), f"Object must be of type `Pose`. Received: {type(pose)}"
        kwargs = self.decoy_widget.kwargs
        index = kwargs["index"] if kwargs else 0
        self.poses[index] = self.poses[index].clone()
        append_pose_to_pose(self.poses[index], pose, new_chain=new_chain)
        self.update_viewer(self.poses[index])

    def apply_modules(self, _pose=None, _pdbstring=None):
        for _module in self.modules:
            func = getattr(_module, f"apply_{self.backend}")
            self.viewer = func(
                self.viewer,
                _pose,
                _pdbstring,
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
        self._toggle_scrolling()
        self._toggle_window(self.window_size)
        display(*self.get_widgets())
        self.show_viewer()

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
        core_display(HTML("<style>.container { width:100% !important; }</style>"))
    except NameError:
        _logger.exception("IPython.core.display module not imported.")
