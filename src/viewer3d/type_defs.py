__author__ = "Jason C. Klima"

import sys

if sys.version_info[:2] < (3, 9):
    from typing import Callable
else:
    from collections.abc import Callable

if sys.version_info[:2] < (3, 10):
    from typing_extensions import (
        Concatenate,
        ParamSpec,
    )
else:
    from typing import (
        Concatenate,
        ParamSpec,
    )

if sys.version_info[:2] < (3, 11):
    from typing_extensions import Self
else:
    from typing import Self

from types import ModuleType
from typing import (
    Any,
    DefaultDict,
    Dict,
    Generator,
    Generic,
    Iterable,
    List,
    NoReturn,
    Optional,
    OrderedDict,
    Tuple,
    TypeVar,
    Union,
    cast,
)

CallableType = TypeVar("T", bound=Callable[..., Any])
ViewerType = TypeVar("ViewerType")
GenericViewer = Generic[ViewerType]
