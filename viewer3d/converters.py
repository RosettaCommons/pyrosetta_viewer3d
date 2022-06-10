import collections
import logging
import pyrosetta
import pyrosetta.distributed.io as io
import os

from functools import singledispatch
from ipywidgets.widgets import Widget
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.pose.full_model_info import (
    get_res_num_from_pdb_info,
    get_chains_from_pdb_info,
)
from pyrosetta.rosetta.core.select import get_residues_from_subset
from typing import List

from viewer3d.exceptions import ViewerInputError


_logger = logging.getLogger("viewer3d.converters")


def _to_poses_pdbstrings(packed_and_poses_and_pdbs):
    @singledispatch
    def to_pose(obj):
        raise ViewerInputError(obj)

    to_pose.register(type(None), lambda obj: None)
    to_pose.register(PackedPose, lambda obj: io.to_pose(obj))
    to_pose.register(Pose, lambda obj: obj)
    to_pose.register(str, lambda obj: None)

    @singledispatch
    def to_pdbstring(obj):
        raise ViewerInputError(obj)

    @to_pdbstring.register(PackedPose)
    @to_pdbstring.register(Pose)
    def _(obj):
        return io.to_pdbstring(obj)

    @to_pdbstring.register(str)
    def _(obj):
        if not os.path.isfile(obj):
            raise ViewerInputError(obj)
        else:
            with open(obj, "r") as f:
                return f.read()

    if isinstance(
        packed_and_poses_and_pdbs, collections.abc.Iterable
    ) and not isinstance(packed_and_poses_and_pdbs, (Pose, PackedPose)):
        poses, pdbstrings = map(
            list,
            zip(*[(to_pose(p), to_pdbstring(p)) for p in packed_and_poses_and_pdbs]),
        )
    else:
        poses = [to_pose(packed_and_poses_and_pdbs)]
        pdbstrings = [to_pdbstring(packed_and_poses_and_pdbs)]

    return poses, pdbstrings


@singledispatch
def _to_float(obj):
    try:
        return float(obj)
    except ValueError:
        raise ValueError(
            "Input argument 'delay' should be an instance of float that is >= 0. Setting 'delay' to default."
        )


_to_float.register(int, lambda obj: float(obj))
_to_float.register(float, lambda obj: obj)


def _to_0_if_le_0(obj):
    return 1e-10 if isinstance(obj, (float, int)) and obj <= 0 else obj


def _to_1_if_gt_1(obj):
    return 1 if isinstance(obj, (float, int)) and obj > 1 else obj


def _to_widgets(objs) -> List[Widget]:
    @singledispatch
    def _to_widget(obj):
        raise ValueError(
            "The 'widgets' viewer attribute must be an instance of `Widget` "
            + f"or an iterable of `Widget` objects. Received: {type(obj)}"
        )

    _to_widget.register(Widget, lambda obj: obj)

    if objs is None:
        _widgets = []
    elif isinstance(objs, collections.abc.Iterable):
        _widgets = list(map(_to_widget, objs))
    else:
        _widgets = [_to_widget(objs)]

    return _widgets


def _pose_to_residue_chain_tuples(pose, residue_selector, logger=_logger):
    """
    Given a `Pose` object and `ResidueSelector` object, return a `tuple` of `list`s containing
    PDB residue numbers and chain IDs for the selection.
    """

    pdb_numbering = list(
        zip(get_res_num_from_pdb_info(pose), get_chains_from_pdb_info(pose))
    )
    residues_from_subset = list(get_residues_from_subset(residue_selector.apply(pose)))
    residue_chain_tuples = [pdb_numbering[i - 1] for i in residues_from_subset]

    if len(residue_chain_tuples) == 0:
        logger.info(
            "ResidueSelector {0} is empty and did not select any residues!".format(
                residue_selector
            )
        )
        return [], []
    else:
        return map(list, zip(*residue_chain_tuples))


def _pdbstring_to_pose(pdbstring, class_name, logger=_logger):
    """Convert pdbstring to a `Pose` with logging."""
    logger.info(
        " ".join(
            "{0} requires pyrosetta.rosetta.core.pose.Pose object but given input .pdb file. \
        Now instantiating pyrosetta.rosetta.core.pose.Pose object from input .pdb file. \
        For faster performance, either input pyrosetta.rosetta.core.pose.Pose \
        or pyrosetta.distributed.packed_pose.core.PackedPose objects to viewer3d.init, \
        or do not add {0} objects that require a pyrosetta.rosetta.core.pose.Pose object.  \
        ".format(
                class_name
            ).split()
        )
    )
    return io.to_pose(io.pose_from_pdbstring(pdbstring))
