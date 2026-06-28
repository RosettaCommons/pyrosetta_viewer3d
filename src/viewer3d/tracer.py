__author__ = "Jason C. Klima"

import pyrosetta.distributed

from functools import wraps

from viewer3d.type_defs import (
    Any,
    CallableType,
    cast,
)


def requires_init(func: CallableType) -> CallableType:
    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        init_kwargs = dict(
            options="",
            extra_options="-out:level 100",
            set_logging_handler="logging",
            notebook=None,
            silent=True,
        )
        pyrosetta.distributed.maybe_init(**init_kwargs)

        return func(*args, **kwargs)

    return cast(CallableType, wrapper)
