__author__ = "Jason C. Klima"

import attr
import collections
import logging

from ipywidgets import Widget
from pyrosetta import Pose

from viewer3d.config import BACKENDS
from viewer3d.converters import (
    _to_backend,
    _to_float,
    _to_widgets,
)
from viewer3d.modules import setZoomTo
from viewer3d.modules.base import ModuleBase
from viewer3d.type_defs import (
    DefaultDict,
    List,
    Optional,
    Tuple,
    Union,
)
from viewer3d.validators import (
    _validate_int_float,
    _validate_window_size,
)

_logger: logging.Logger = logging.getLogger("viewer3d.initialization")


@attr.define(kw_only=False, slots=False, frozen=False)
class InitBase:
    poses: DefaultDict[int, List[Optional[Pose]]] = attr.field(
        default=collections.defaultdict(list),
    )
    pdbstrings: DefaultDict[int, List[Optional[str]]] = attr.field(
        default=collections.defaultdict(list),
    )
    window_size: Tuple[Union[int, float], Union[int, float]] = attr.field(
        default=None,
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of((int, float)),
                iterable_validator=attr.validators.instance_of(collections.abc.Iterable),
            ),
            _validate_window_size,
        ],
        converter=attr.converters.default_if_none(default=(1200, 800)),
    )
    modules: List[ModuleBase] = attr.field(
        default=None,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(ModuleBase),
            iterable_validator=attr.validators.instance_of(list),
        ),
        converter=attr.converters.default_if_none(default=[]),
    )
    delay: float = attr.field(
        default=None,
        validator=_validate_int_float,
        converter=attr.converters.pipe(attr.converters.default_if_none(default=0.0), _to_float),
    )
    continuous_update: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    widgets: List[Widget] = attr.field(
        default=None,
        validator=attr.validators.optional(
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of(Widget),
                iterable_validator=attr.validators.instance_of(list),
            )
        ),
        converter=_to_widgets,
    )
    backend: str = attr.field(
        default=None,
        validator=[attr.validators.instance_of(str), attr.validators.in_(BACKENDS)],
        converter=[attr.converters.default_if_none(default=0), _to_backend],
    )
    auto_show: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    gui: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    _setZoomTo: ModuleBase = attr.field(
        default=setZoomTo(),
        validator=attr.validators.instance_of(ModuleBase),
        init=False,
    )
