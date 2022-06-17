import pyrosetta.distributed

from functools import wraps
from typing import (
    Any,
    Callable,
    TypeVar,
    cast,
)


T = TypeVar("T", bound=Callable[..., Any])


def requires_init(func: T) -> T:
    @wraps(func)
    def wrapper(*args, **kwargs):
        init_kwargs = dict(
            options="",
            extra_options="-out:level 100",
            set_logging_handler="logging",
            notebook=None,
            silent=True,
        )
        pyrosetta.distributed.maybe_init(**init_kwargs)

        return func(*args, **kwargs)

    return cast(T, wrapper)
