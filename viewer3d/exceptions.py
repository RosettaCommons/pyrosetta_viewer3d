import os

from typing import Any

from viewer3d.config import BACKENDS, URLS


class ModuleNotImplementedError(NotImplementedError):
    """Exception raised upon implementing backends."""

    def __init__(self, class_name: str, backend: str):
        super().__init__(f"{class_name} is not supported for `{backend}` backend.")


class ViewerImportError(ImportError):
    """Exception raised upon importing backends."""

    def __init__(self, backend: str):
        _backend_urls = dict(zip(BACKENDS, URLS))
        super().__init__(
            f"Using the '{backend}' backend requires the third-party package `{backend}`.{os.linesep}"
            + "Please install the package into your python environment. "
            + f"For installation instructions, visit:{os.linesep}"
            + f"{_backend_urls[backend]}{os.linesep}"
        )


class ViewerInputError(Exception):
    """Exception raised for errors with the input argument `packed_and_poses_and_pdbs`."""

    def __init__(self, obj: Any):
        super().__init__(
            " ".join(
                "Input argument 'packed_and_poses_and_pdbs' should be an instance of \
                pyrosetta.rosetta.core.pose.Pose, pyrosetta.distributed.packed_pose.core.PackedPose, \
                or a valid path string to a .pdb file, or a list, set, or tuple of these objects. \
                Input argument 'packed_and_poses_and_pdbs' was invoked with: {0}".format(
                    obj
                ).split()
            )
        )
