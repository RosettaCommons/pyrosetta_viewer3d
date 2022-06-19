import logging

from decorator import decorator
from functools import wraps
from typing import Any, Callable, Iterable, NoReturn, Optional, TypeVar, Union, cast


_logger = logging.getLogger("viewer3d.validators")

V = TypeVar("V", bound=Callable[..., Any])


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


@decorator
def requires_show(func, self, *args, **kwargs):
    update_viewer = args[-1]
    if update_viewer and not self._displayed:
        _class_name = self.__class__.__name__
        _func_name = func.__name__
        _logger.warning(
            f"The `{_class_name}.{_func_name}` method requires calling `{_class_name}.show`."
        )
        self.show()
    return func(self, *args, **kwargs)
