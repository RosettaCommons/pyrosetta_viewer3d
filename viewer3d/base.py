import attr
import collections
import math
import logging
import pyrosetta.distributed.io as io
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
from typing import Generic, Iterable, List, Optional, Tuple, TypeVar, Union

from viewer3d.config import _import_backend, BACKENDS
from viewer3d.converters import _to_backend, _to_float, _to_widgets
from viewer3d.exceptions import ViewerImportError
from viewer3d.modules import ModuleBase, setZoomTo
from viewer3d.validators import _validate_int_float, _validate_window_size


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
            attr.converters.default_if_none(default=0.0), _to_float
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
        converter=[attr.converters.default_if_none(default=0), _to_backend],
    )
    auto_show = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    gui = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    _setZoomTo = attr.ib(
        type=ModuleBase,
        default=setZoomTo(),
        validator=attr.validators.instance_of(ModuleBase),
        init=False,
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
class PoseBase:
    def add_pose(self, pose, index=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        self.poses[index].append(pose)
        self.pdbstrings[index].append(None)
        if update_viewer:
            model = len(self.poses[index]) - 1
            self.add_objects(self.poses[index], self.pdbstrings[index], _model=model)

    def add_pdbstring(self, pdbstring, index=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        self.poses[index].append(None)
        self.pdbstrings[index].append(pdbstring)
        if update_viewer:
            model = len(self.poses[index]) - 1
            self.add_objects(self.poses[index], self.pdbstrings[index], _model=model)

    def remove_pose(self, index=None, model=None, update_viewer=True):
        self.remove_pdbstring(index=index, model=model, update_viewer=update_viewer)

    def remove_pdbstring(self, index=None, model=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = -1
        if len(self.pdbstrings[index]) > 0:
            self.poses[index].pop(model)
            self.pdbstrings[index].pop(model)
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes are empty at index `{index}`."
            )
        if update_viewer:
            self.remove_objects(model)

    def update_pose(self, pose, index=None, model=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.poses[index]))):
            model = 0
        if index in self.poses.keys():
            self.poses[index][model] = pose
            self.pdbstrings[index][model] = None
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_objects(self.poses[index], self.pdbstrings[index], _model=model)

    def update_pdbstring(self, pdbstring, index=None, model=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = 0
        if index in self.pdbstrings.keys():
            self.poses[index][model] = None
            self.pdbstrings[index][model] = pdbstring
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_objects(self.poses[index], self.pdbstrings[index], _model=model)

    def update_poses(self, poses, index=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        assert isinstance(poses, list)
        for pose in poses:
            assert isinstance(pose, Pose)
        self.poses[index] = poses
        self.pdbstrings[index] = [None] * len(poses)
        if update_viewer:
            self.update_viewer(index=index)

    def update_pdbstrings(self, pdbstrings, index=None, update_viewer=True):
        if index is None:
            index = self.get_decoy_widget_index()
        assert isinstance(pdbstrings, list)
        for pdbstring in pdbstrings:
            assert isinstance(pdbstring, str)
        self.poses[index] = [None] * len(pdbstrings)
        self.pdbstrings[index] = pdbstrings
        if update_viewer:
            self.update_viewer(index=index)


@attr.s(kw_only=False, slots=False)
class WidgetsBase:
    decoy_widget = attr.ib(
        default=attr.Factory(
            lambda self: interactive(
                self.update_viewer,
                index=IntSlider(
                    min=0,
                    max=self.n_decoys - 1,
                    description="Decoys",
                    continuous_update=self.continuous_update,
                ),
            ),
            takes_self=True,
        ),
    )

    def get_decoy_widget_index(self):
        kwargs = self.decoy_widget.kwargs
        index = kwargs["index"] if kwargs else 0
        return index

    def get_widgets(self):
        _widgets = self.widgets.copy()
        if self.n_decoys > 1:
            _widgets.insert(0, self.decoy_widget)
        else:
            self.update_viewer()
        return _widgets

    def get_widgets_dict(self):
        return {
            widget.description: widget.value
            for widget in self.widgets
            if all(hasattr(widget, attr) for attr in ("description", "value"))
        }

    def set_widgets(self, obj):
        self.widgets = _to_widgets(obj)

    def update_viewer(self, index: Optional[int] = None):
        time.sleep(self.delay)
        if index is None:
            index = self.get_decoy_widget_index()
        if index in self.poses.keys():
            self.update_objects(self.poses[index], self.pdbstrings[index])

    def _maybe_setup_colab(self):
        if "google.colab" in sys.modules:
            sys.modules["google.colab"].output.enable_custom_widget_manager()


@attr.s(kw_only=False, slots=False)
class ViewerBase(Base3D, PoseBase, WidgetsBase):
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

    def update_objects(self, _poses, _pdbstrings, _model: Optional[int] = None):
        """
        Setup Viewer in Jupyter notebook.

        Args:
            _models: if `None`, then update all models. If `int`, then update a single model.
        """
        assert len(_poses) == len(
            _pdbstrings
        ), "Number of `Pose` objects and PDB `str` objects must be equal."
        self.remove_objects(_model)
        self.add_objects(_poses, _pdbstrings, _model)

    def display_widgets(self):
        widgets = self.get_widgets()
        if widgets:
            display(*widgets)

    def set_modules(self, obj):
        self.modules = ModuleBase._to_modules(obj)

    def show(self):
        """Display Viewer in Jupyter notebook."""
        if self._in_notebook():
            self._clear_output()
            self._toggle_window(self.window_size)
            self.display_widgets()
            self.show_viewer()
            self._toggle_scrolling()
            self._displayed = True


def expand_notebook():
    """Expand Jupyter notebook cell to maximum width."""
    try:
        _logger.debug("IPython.core.display expanding Jupyter notebook cell width.")
        core_display(HTML("<style>.container { width:100% !important; }</style>"))
    except NameError:
        _logger.exception("IPython.core.display module not imported.")
