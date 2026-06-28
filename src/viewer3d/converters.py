__author__ = "Jason C. Klima"

import collections
import logging
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy
import pyrosetta.distributed.io as io
import os
import requests
import time

from functools import singledispatch
from io import StringIO
from ipywidgets.widgets import Widget
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.pose.full_model_info import (
    get_res_num_from_pdb_info,
    get_chains_from_pdb_info,
)
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import ResidueSelector
from pyrosetta.toolbox.rcsb import download_from_web

from viewer3d.config import BACKENDS
from viewer3d.exceptions import ViewerInputError
from viewer3d.type_defs import (
    Any,
    DefaultDict,
    List,
    NoReturn,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.converters")


def _py3Dmol_to_nglview_style(style: str) -> str:
    if style == "stick":
        style = "licorice"
    elif style == "sphere":
        style = "spacefill"
    elif style == "cross":
        style = "point"
    elif style == "line":
        style = "line"
    return style


def _py3Dmol_to_pymol_style(style: str) -> str:
    if style == "stick":
        style = "sticks"
    elif style == "sphere":
        style = "spheres"
    elif style == "cross":
        style = "dots"
    elif style == "line":
        style = "lines"
    return style


def _is_pdbstring(obj: str) -> bool:
    if not obj:
        return False

    _pdb_record_prefixes = (
        "ATOM",
        "HETATM",
        "HEADER",
        "TITLE",
        "MODEL",
        "ENDMDL",
        "TER",
        "END",
    )
    for line in obj.splitlines():
        if any(line.startswith(prefix) for prefix in _pdb_record_prefixes):
            return True

    return False


def _to_poses_pdbstrings(
    packed_and_poses_and_pdbs: Any,
) -> Tuple[
    DefaultDict[int, List[Optional[Pose]]],
    DefaultDict[int, List[Optional[str]]],
]:
    @singledispatch
    def to_pose(obj: Any):
        raise ViewerInputError(obj)

    to_pose.register(type(None), lambda obj: None)
    to_pose.register(PackedPose, lambda obj: io.to_pose(obj))
    to_pose.register(Pose, lambda obj: obj)
    to_pose.register(str, lambda obj: None)

    @singledispatch
    def to_pdbstring(obj: Any) -> NoReturn:
        raise ViewerInputError(obj)

    @to_pdbstring.register(PackedPose)
    @to_pdbstring.register(Pose)
    def _(obj: Union[Pose, PackedPose]) -> None:
        return None

    @to_pdbstring.register(str)
    def _(obj: str) -> Union[NoReturn, str]:
        if os.path.isfile(obj):
            with open(obj, "r") as f:
                return f.read()
        elif _is_pdbstring(obj):
            return obj
        else:
            raise ViewerInputError(obj)

    to_pdbstring.register(type(None), lambda obj: None)

    def to_dict(objs: List[Any]) -> DefaultDict[int, List[Any]]:
        d = collections.defaultdict(list)
        for i, obj in enumerate(objs):
            d[i].append(obj)
        return d

    def remove_none(poses: List[Pose], pdbstrings: List[str]) -> Tuple[List[Pose], List[str]]:
        """Remove `NoneType` objects from models."""
        assert len(poses.keys()) == len(pdbstrings.keys())
        for index in range(len(poses.keys())):
            assert len(poses[index]) == len(pdbstrings[index])
            for model in range(len(poses[index])):
                if all(p[index][model] is None for p in (poses, pdbstrings)):
                    poses[index].pop(model)
                    pdbstrings[index].pop(model)

        return poses, pdbstrings

    if isinstance(packed_and_poses_and_pdbs, collections.abc.Iterable) and not isinstance(
        packed_and_poses_and_pdbs, (Pose, PackedPose, str)
    ):
        poses, pdbstrings = map(
            to_dict,
            map(
                list,
                zip(
                    *map(
                        lambda p: (to_pose(p), to_pdbstring(p)),
                        packed_and_poses_and_pdbs,
                    )
                ),
            ),
        )
    else:
        poses = to_dict([to_pose(packed_and_poses_and_pdbs)])
        pdbstrings = to_dict([to_pdbstring(packed_and_poses_and_pdbs)])

    poses, pdbstrings = remove_none(poses, pdbstrings)

    return poses, pdbstrings


@singledispatch
def _to_float(obj: Any) -> float:
    try:
        return float(obj)
    except ValueError:
        raise ValueError(
            "Input argument 'delay' should be an instance of float that is >= 0. Setting 'delay' to default."
        )


_to_float.register(int, lambda obj: float(obj))
_to_float.register(float, lambda obj: obj)


def _to_hex(obj: Any) -> Any:
    if isinstance(obj, str) and obj.startswith("#"):
        # return int(obj.replace("#", ""), 16) + int("0x200", 16)
        return int(obj.replace("#", ""), 16)
    else:
        return obj


def _to_0_if_le_0(obj: Any) -> Any:
    return 1e-10 if isinstance(obj, (float, int)) and obj <= 0 else obj


def _to_1_if_gt_1(obj: Any) -> Any:
    return 1 if isinstance(obj, (float, int)) and obj > 1 else obj


def _to_backend(obj: Any) -> str:
    if isinstance(obj, int):
        try:
            backend = BACKENDS[obj]
        except IndexError:
            raise IndexError(f"Backend index must be in: {tuple(range(len(BACKENDS)))}")
    elif isinstance(obj, str):
        if obj not in BACKENDS:
            raise ValueError(f"Backend must be in: {BACKENDS}")
        backend = obj
    else:
        raise TypeError(f"Backend must be an `int` or `str` object. Received: {type(obj)}")

    return backend


@singledispatch
def _int_to_str(obj):
    return obj


@_int_to_str.register(int)
def _to_str(obj):
    return f"#{obj:06x}".upper()


@singledispatch
def _int32_to_str(obj):
    return obj


@_int32_to_str.register(int)
def _to_str(obj):
    return f"#{obj:06X}"


@singledispatch
def _hex_to_rgb(obj):
    return obj


@_hex_to_rgb.register(str)
def _to_rgb(obj) -> List[float]:
    v = obj.lstrip("#")
    I = 255.0
    if len(v) == 6:  # RGB
        r, g, b = (int(v[i : i + 2], 16) for i in (0, 2, 4))
        return [r / I, g / I, b / I]

    if len(v) == 8:  # RGBA
        r, g, b, _a = (int(v[i : i + 2], 16) for i in (0, 2, 4, 6))
        return [r / I, g / I, b / I]

    raise ValueError(f"Invalid hex color: {obj}")


@_hex_to_rgb.register(int)
def _from_int(obj) -> List[float]:
    return _hex_to_rgb.dispatch(str)(_int32_to_str(obj))


def _to_int_color(color):
    if isinstance(color, int):
        return color
    if isinstance(color, str):
        color = color.strip()
        if color.startswith("#"):
            return int(color[1:], 16)
        if color.startswith("0x"):
            return int(color, 16)
    raise ValueError(f"Unsupported color format: {color}")


def _to_widgets(objs) -> List[Widget]:
    @singledispatch
    def _to_widget(obj: Any) -> NoReturn:
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


def _pose_to_residue_chain_tuples(
    pose, residue_selector: ResidueSelector, logger: logging.Logger = _logger
) -> Tuple[List[int], List[str]]:
    """
    Given a `Pose` object and `ResidueSelector` object, return a `tuple` of `list`s containing
    PDB residue numbers and chain IDs for the selection.
    """
    pdb_numbering = list(zip(get_res_num_from_pdb_info(pose), get_chains_from_pdb_info(pose)))
    residues_from_subset = list(get_residues_from_subset(residue_selector.apply(pose)))
    residue_chain_tuples = [pdb_numbering[i - 1] for i in residues_from_subset]

    if len(residue_chain_tuples) == 0:
        logger.info(
            "ResidueSelector {0} is empty and did not select any residues!".format(residue_selector)
        )
        return [], []
    else:
        return map(list, zip(*residue_chain_tuples))


def _get_nglview_selection(
    pose,
    residue_selector: ResidueSelector,
    show_hydrogens: bool = False,
    nbr_atom_only: bool = False,
    logger: logging.Logger = _logger,
) -> str:
    def _get_nbr_atom(rc: Tuple[str, str]) -> str:
        _residue = pose.residue(pose.pdb_info().pdb2pose(rc[1], int(rc[0])))
        _nbr_atom_index = _residue.type().nbr_atom()
        return _residue.atom_name(_nbr_atom_index).strip()

    def _format_residue(rc: Tuple[str, str]) -> str:
        if show_hydrogens:
            if nbr_atom_only:
                return f"({rc[0]}:{rc[1]}.{_get_nbr_atom(rc)})"
            else:
                return f"({rc[0]}:{rc[1]})"
        else:
            if nbr_atom_only:
                return f"({rc[0]}:{rc[1]}.{_get_nbr_atom(rc)} and not _H)"
            else:
                return f"({rc[0]}:{rc[1]} and not _H)"

    resi, chain = _pose_to_residue_chain_tuples(pose, residue_selector, logger=_logger)
    selection = " or ".join(map(_format_residue, zip(resi, chain)))

    return selection


def _get_pymol_selection(
    pose,
    residue_selector: ResidueSelector,
    show_hydrogens: bool = False,
    logger: logging.Logger = _logger,
) -> str:
    resi, chain = _pose_to_residue_chain_tuples(pose, residue_selector, logger=_logger)
    residue_chain_tuples = list(zip(map(str, resi), chain))
    selection = " or ".join(
        [f"(chain {_chain} and resi {_resi})" for _resi, _chain in residue_chain_tuples]
    )
    if not selection:
        return ""

    if not show_hydrogens:
        selection = f"({selection}) and not elem h"

    return selection


def _pdbstring_to_pose(pdbstring, class_name, logger=_logger):
    """Convert pdbstring to a `Pose` with logging."""
    logger.info(
        " ".join(
            "{0} requires `pyrosetta.rosetta.core.pose.Pose` object but given input '.pdb' file. \
        Now instantiating `pyrosetta.rosetta.core.pose.Pose` object from input '.pdb' file. \
        For faster performance, either input `pyrosetta.rosetta.core.pose.Pose` \
        or `pyrosetta.distributed.packed_pose.core.PackedPose` objects to `viewer3d.init`, \
        or do not add {0} objects that require a `pyrosetta.rosetta.core.pose.Pose` object.  \
        ".format(class_name).split()
        )
    )

    return io.to_pose(io.pose_from_pdbstring(pdbstring))


def _get_residue_chain_tuple(pose: Pose, res: int) -> Tuple[str, str]:
    residue, chain = map(lambda x: x.strip(), pose.pdb_info().pose2pdb(res).split())
    return residue, chain


def atom_array_to_pdbstring(atom_array: "AtomArray") -> str:
    """
    A helper function to convert a `biotite` `AtomArray` object into a PDB string.

    Args:
        atom_array: The input `AtomArray` object to convert.

    Raises:
        ImportError: If the `biotite` package is not installed.

    Returns:
        A `str` object representing the `AtomArray` object.
    """
    try:
        from biotite.structure.io.pdb import PDBFile
    except ImportError as ex:
        raise ImportError(
            f"{type(ex).__name__}: Please install the 'biotite' package "
            "into your virtual environment, then try again."
        ) from ex

    buffer = StringIO()
    pdb = PDBFile()
    pdb.set_structure(atom_array)
    pdb.write(buffer)

    return buffer.getvalue()


def atom_array_to_pose(atom_array: "AtomArray") -> Pose:
    """
    A helper function to convert a `biotite` `AtomArray` object into a PyRosetta `Pose` object.

    PyRosetta must first be initialized and loaded with any Rosetta topology files and/or
    Rosetta patch files to correctly represent the biomolecular structure in memory, otherwise
    PyRosetta will be initialized with default Rosetta command-line options.

    Args:
        atom_array: The input `AtomArray` object to convert.

    Returns:
        A `Pose` object representing the `AtomArray` object.
    """
    pdbstring = atom_array_to_pdbstring(atom_array)
    packed_pose = io.pose_from_pdbstring(pdbstring)

    return packed_pose.pose


def pdbstring_from_alphafold_id(
    alphafold_id: str, version: str = "v1", verbose: bool = True
) -> str:
    """
    A helper function to download a PDB file of an AlphaFold ID from the AlphaFold Protein
    Structure Database and return the PDB string.

    Args:
        alphafold_id: The input AlphaFold ID. For example, `"AF-0000000013593985"` or
            `"AF-A0A485P7I0-F1"`.

        version: The AlphaFold ID version. For example, `"v1"` or `"v6"`.

        verbose: Whether to print timing information for the PDB file download.

    Raises:
        ValueError: If the PDB file cannot be downloaded.

    Returns:
        A `str` object representing the AlphaFold ID.
    """
    url = f"https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_{version}.pdb"
    if verbose:
        print(f"Starting download: '{url}'")
    t0 = time.perf_counter()
    pdb = requests.get(url)
    dt = time.perf_counter() - t0
    if verbose:
        print(f"Finished download in {dt:.2f} seconds.")
    if pdb.ok:
        return pdb.text
    else:
        raise ValueError(f"Could not download PDB file: '{url}'")


def pose_from_alphafold_id(alphafold_id: str, version: str = "v1", verbose: bool = True) -> Pose:
    """
    A helper function to download a PDB file of an AlphaFold ID from the AlphaFold Protein
    Structure Database and return a `Pose` object from it.

    PyRosetta must first be initialized and loaded with any Rosetta topology files and/or
    Rosetta patch files to correctly represent the biomolecular structure in memory, otherwise
    PyRosetta will be initialized with default Rosetta command-line options.

    Args:
        alphafold_id: The input AlphaFold ID. For example, `"AF-0000000013593985"` or
            `"AF-A0A485P7I0-F1"`.

        version: The AlphaFold ID version. For example, `"v1"` or `"v6"`.

        verbose: Whether to print timing information for the PDB file download.

    Returns:
        A `Pose` object representing the AlphaFold ID.
    """
    pdbstring = pdbstring_from_alphafold_id(alphafold_id, version=version, verbose=verbose)

    return io.pose_from_pdbstring(pdbstring).pose


def pdbstring_from_pdb_id(pdb_id: str, verbose: bool = True) -> str:
    """
    A helper function to download a PDB file of a PDB accession number from RCSB and return
    the PDB string.

    Args:
        pdb_id: The input PDB accession number.

        verbose: Whether to print timing information for the PDB file download.

    Returns:
        A `str` object representing the PDB accession number.
    """
    url = f"http://files.rcsb.org/download/{pdb_id}.pdb"
    if verbose:
        print(f"Starting download: '{url}'")
    t0 = time.perf_counter()
    with download_from_web(url) as pdb_data:
        pdbstring = "\n".join(pdb_data)
    dt = time.perf_counter() - t0
    if verbose:
        print(f"Finished download in {dt:.2f} seconds.")

    return pdbstring


def pose_from_pdb_id(pdb_id: str, verbose: bool = True) -> Pose:
    """
    A helper function to download a PDB file of a PDB accession number from RCSB and return a
    `Pose` object from it.

    PyRosetta must first be initialized and loaded with any Rosetta topology files and/or
    Rosetta patch files to correctly represent the biomolecular structure in memory, otherwise
    PyRosetta will be initialized with default Rosetta command-line options.

    Args:
        pdb_id: The input PDB accession number.

        verbose: Whether to print timing information for the PDB file download.

    Returns:
        A `Pose` object representing the PDB accession number.
    """
    pdbstring = pdbstring_from_pdb_id(pdb_id, verbose=verbose)

    return io.pose_from_pdbstring(pdbstring).pose


def get_matplotlib_cmap(name: str, num: int = 256) -> List[str]:
    """
    A helper function to return a palette from a named `matplotlib` colormap.

    Args:
        name: The colormap name (e.g., `"viridis"`).

        num: The number of samples to generate. Must be a non-negative integer.

    Returns:
        A `Pose` object representing the AlphaFold ID.
    """
    if not isinstance(name, str):
        raise TypeError("The `name` argument value must be a `str` object.")
    if not isinstance(num, int):
        raise TypeError("The `num` keyword argument value must be an `int` object.")
    if num < 0:
        raise ValueError("The `num` keyword argument value must be non-negative.")

    cmap = plt.get_cmap(name)

    return [mcolors.to_hex(cmap(x)) for x in numpy.linspace(0, 1, num)]
