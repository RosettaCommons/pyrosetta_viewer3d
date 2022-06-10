import pyrosetta
import pyrosetta.distributed

from functools import wraps
from pyrosetta.rosetta.basic.options import get_integer_option, set_integer_option
from typing import (
    Any,
    Callable,
    TypeVar,
    cast,
)


T = TypeVar("T", bound=Callable[..., Any])


@pyrosetta.distributed.requires_init
def silence_tracer(func: T) -> T:
    """Silence PyRosetta tracer output."""

    @wraps(func)
    def wrapper(self, **kwargs):
        option = "out:level"
        user_out_level = get_integer_option(option)
        set_integer_option(option, 100)
        func(self, **kwargs)
        set_integer_option(option, user_out_level)

    return cast(T, wrapper)
