import os


class ModuleInputError(Exception):
    """Exception raised for errors with the input argument `residue_selector`."""

    def __init__(self, obj):
        super().__init__(
            " ".join(
                "Input 'residue_selector' argument should be an instance of \
                pyrosetta.rosetta.core.select.residue_selector.ResidueSelector. \
                Input argument 'residue_selector' was invoked with: {0}".format(
                    obj
                ).split()
            )
        )


class ModuleNotImplementedError(NotImplementedError):
    """Exception raised upon implementing backends."""

    def __init__(self, class_name, backend):
        super().__init__(f"{class_name} is not supported for `{backend}` backend.")


class ViewerImportError(ImportError):
    """Exception raised upon importing backends."""

    def __init__(self, backend, url):
        super().__init__(
            f"Using the '{backend}' backend requires the third-party package `{backend}`.{os.linesep}"
            + "Please install the package into your python environment. "
            + f"For installation instructions, visit:{os.linesep}"
            + f"{url}{os.linesep}"
        )


class ViewerInputError(Exception):
    """Exception raised for errors with the input argument `packed_and_poses_and_pdbs`."""

    def __init__(self, obj):
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
