from functools import wraps
from pyrosetta import Pose
from typing import (
    Any,
    Iterable,
    Callable,
    NoReturn,
    Optional,
    TypeVar,
    Union,
    cast,
)

A = TypeVar("A", bound=Callable[..., Any])


def _validate_int(self, attribute: str, value: int) -> Optional[NoReturn]:
    """Validate that integers are greater than or equal to 1."""
    if value < 1:
        raise ValueError(
            f"`{attribute}` must be a positive integer greater than or equal to 1."
        )


def _validate_float(
    self, attribute: str, value: Union[float, int]
) -> Optional[NoReturn]:
    """Validate that floats are greater than or equal to 0.0"""
    msg = f"`{attribute}` must be a positive `float` or `int` value greater than or equal to 0."
    try:
        float(value)
    except:
        raise ValueError(msg)
    if value < 0.0:
        raise ValueError(msg)


def _validate_int_float(self, attribute: str, value: Any) -> Optional[NoReturn]:
    """Validate that floats are greater than or equal to 0.0 and integers are greater than or equal to 1."""
    if isinstance(value, int):
        _validate_int(self, attribute, value)
    elif isinstance(value, float):
        _validate_float(self, attribute, value)


def _validate_window_size(self, attribute: str, value: Iterable) -> Optional[NoReturn]:
    """Validate that the 'window_size' argument parameter is an iterable of length 2."""
    assert (
        len(value) == 2
    ), "Input argument 'window_size' must be an iterable of length 2."
    for v in value:
        _validate_int_float(self, attribute, v)


def _validate_add_pose(func: A) -> A:
    _func_name = func.__name__

    @wraps(func)
    def wrapper(self, pose=None, index=None, update_viewer=None):
        if not isinstance(pose, Pose):
            raise TypeError(
                f"The `{_func_name}` 'pose' keyword argument parameter must be of type `Pose`. Received: {type(pose)}"
            )
        if index is not None:
            if not isinstance(index, int):
                raise TypeError(
                    f"The `{_func_name}` 'index' keyword argument parameter must be of type `int`. Received: {type(index)}"
                )
            if index not in self.poses.keys():
                raise IndexError(
                    f"The `{_func_name}` 'index' keyword argument parameter `{index}` does not exist in `ViewerBase.poses` object."
                )
        return func(self, pose=pose, index=index)

    return cast(A, wrapper)
