import attr
import collections
import logging

from ipywidgets import Widget
from pyrosetta import Pose
from typing import Iterable, List, Optional, Tuple, Union

from viewer3d.config import BACKENDS
from viewer3d.converters import _to_backend, _to_float, _to_widgets
from viewer3d.modules import ModuleBase, setZoomTo
from viewer3d.validators import (
    _validate_int_float,
    _validate_window_size,
)


_logger: logging.Logger = logging.getLogger("viewer3d.initialization")


@attr.s(kw_only=False, slots=False)
class InitBase:
    poses = attr.ib(type=Iterable[Pose], default=None)
    pdbstrings = attr.ib(type=Iterable[str], default=None)
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
