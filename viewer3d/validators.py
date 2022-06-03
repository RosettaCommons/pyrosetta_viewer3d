from pyrosetta.distributed.cluster.validators import _validate_int, _validate_float

from typing import Iterable, NoReturn, Optional, Union


def _validate_int_float(
    self, attribute: str, value: Union[int, float]
) -> Optional[NoReturn]:
    if isinstance(value, int):
        _validate_int(self, attribute, value)
    elif isinstance(value, float):
        _validate_float(self, attribute, value)


def _validate_window_size(self, attribute: str, value: Iterable) -> Optional[NoReturn]:
    assert (
        len(value) == 2
    ), "Input argument 'window_size' must be an iterable of length 2."
    for v in value:
        _validate_int_float(self, attribute, v)
