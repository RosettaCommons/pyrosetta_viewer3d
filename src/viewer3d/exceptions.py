__author__ = "Jason C. Klima"

import attr
import logging
import os

from viewer3d.config import (
    BACKENDS,
    URLS,
)
from viewer3d.type_defs import Any

_logger: logging.Logger = logging.getLogger("viewer3d.exceptions")


@attr.define(kw_only=False, slots=True, frozen=True)
class ModuleNotImplementedError(NotImplementedError):
    """Exception raised upon implementing backends."""

    def __init__(self, class_name: str, backend: str) -> None:
        super().__init__(f"{class_name} is not supported for `{backend}` backend.")


@attr.define(kw_only=False, slots=True, frozen=True)
class ViewerImportError(ImportError):
    """Exception raised upon importing backends."""

    def __init__(self, backend: str) -> None:
        _backend_urls = dict(zip(BACKENDS, URLS))
        super().__init__(
            f"Using the '{backend}' backend requires the third-party package `{backend}`.{os.linesep}"
            + "Please install the package into your python environment. "
            + f"For installation instructions, visit:{os.linesep}"
            + f"{_backend_urls[backend]}{os.linesep}"
        )


@attr.define(kw_only=False, slots=True, frozen=True)
class ViewerInputError(Exception):
    """Exception raised for errors with the input argument `packed_and_poses_and_pdbs`."""

    def __init__(self, obj: Any) -> None:
        super().__init__(
            " ".join(
                "Input argument 'packed_and_poses_and_pdbs' must be an instance of \
                `Pose`, `PackedPose`, a valid filesystem path to a PDB file, a PDB string, \
                or an iterable of these objects. Input argument 'packed_and_poses_and_pdbs' \
                was invoked with: {0}".format(
                    obj
                ).split()
            )
        )
