import attr
import itertools
import logging
import pyrosetta
import pyrosetta.distributed.io as io
import sys

from PIL import ImageColor
from pyrosetta.rosetta.core.conformation import is_disulfide_bond
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector,
)
from pyrosetta.rosetta.core.id import AtomID
from typing import Dict, Optional, Tuple, Union


from viewer3d.config import BACKENDS
from viewer3d.converters import (
    _get_nglview_selection,
    _get_residue_chain_tuple,
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
    _to_0_if_le_0,
    _to_1_if_gt_1,
)
from viewer3d.exceptions import ModuleNotImplementedError
from viewer3d.tracer import requires_init


_logger = logging.getLogger("viewer3d.modules")


@attr.s(kw_only=False, slots=True)
class ModuleBase:
    pass


@attr.s(kw_only=False, slots=True)
class setBackgroundColor(ModuleBase):
    """
    Set Viewer background color with either Hexcode or standard colors.

    Parameters
    ----------
    first : optional
        `color`

        Hexcode literal (e.g. 0xffffffff) or `str` indicating a standard color (e.g. "black").
        Default: 0xffffffff

    Returns
    -------
    A Viewer instance.
    """

    color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default=0xFFFFFFFF),
    )

    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        viewer.setBackgroundColor(self.color)

        return viewer

    def apply_nglview(self, viewer, pose, pdbstring, model):
        viewer.background = self.color

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setDisulfides(ModuleBase):
    """
    Display disulfide bonds according to `pyrosetta.rosetta.core.conformation.is_disulfide_bond()`
    for all combinations of cysteine residues in each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "black").
        Default: "gold"

    second : optional
        `radius`

        `float` or `int` indicating the radius of the stick connecting the atoms participating in each disulfide bond.
        Default: 0.5

    Returns
    -------
    A Viewer instance.
    """

    color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default="gold"),
    )
    radius = attr.ib(
        default=None,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=attr.converters.default_if_none(default=0.1),
    )
    sulfur_atom_name = attr.ib(
        default="SG",
        type=str,
        validator=attr.validators.instance_of(str),
        init=False,
    )

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        cys_res = [i for i, aa in enumerate(pose.sequence(), start=1) if aa == "C"]
        for (i, j) in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_xyz = pose.xyz(
                    AtomID(pose.residue(i).atom_index(self.sulfur_atom_name), i)
                )
                j_xyz = pose.xyz(
                    AtomID(pose.residue(j).atom_index(self.sulfur_atom_name), j)
                )
                viewer.addCylinder(
                    {
                        "radius": self.radius,
                        "color": self.color,
                        "fromCap": True,
                        "toCap": True,
                        "start": {"x": i_xyz[0], "y": i_xyz[1], "z": i_xyz[2]},
                        "end": {"x": j_xyz[0], "y": j_xyz[1], "z": j_xyz[2]},
                    }
                )

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        cys_res = [i for i, aa in enumerate(pose.sequence(), start=1) if aa == "C"]
        selection_disulfides = []
        for (i, j) in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_res, i_chain = _get_residue_chain_tuple(pose, i)
                j_res, j_chain = _get_residue_chain_tuple(pose, j)
                i_sele = f"{i_res}:{i_chain}.{self.sulfur_atom_name}"
                j_sele = f"{j_res}:{j_chain}.{self.sulfur_atom_name}"
                selection_disulfides.append([i_sele, j_sele])
        selection = " or ".join([f"({s[0]} or {s[1]})" for s in selection_disulfides])
        viewer.add_representation(
            repr_type="ball+stick",
            selection=selection,
            color=self.color,
            radius=self.radius,
            component=model,
        )
        viewer.add_distance(
            atom_pair=selection_disulfides,
            color=self.color,
            radius=self.radius,
            label_visible=False,
        )

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setHydrogenBonds(ModuleBase):
    """
    Display hydrogen bonds according to `pyrosetta.rosetta.core.pose.Pose.get_hbonds()`
    in each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "yellow").
        Default: "black"

    second : optional
        `dashed`

        `True` or `False` to show hydrogen bonds as dashed lines.
        If `True`, then option `radius` must be `None`.
        If `False`, then option `radius` must be specified.
        Default: True

    third : optional
        `radius`

        `float` or `int` indicating the radius of the solid (non-dashed) stick connecting the atoms participating
        in each hydrogen bond. If set, this automatically sets the option `dashed` to `False`.
        Default: None

    Returns
    -------
    A Viewer instance.
    """

    color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default="black"),
    )
    dashed = attr.ib(
        default=True,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    radius = attr.ib(
        default=None,
        type=Optional[Union[float, int]],
        validator=attr.validators.instance_of((float, int, type(None))),
        converter=_to_0_if_le_0,
    )

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        hbond_set = pose.get_hbonds()
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_xyz = pose.residue(r.don_res()).xyz(r.don_hatm())
                    acc_xyz = pose.residue(r.acc_res()).xyz(r.acc_atm())
                    if self.radius:
                        if self.dashed:
                            _logger.warning(
                                " ".join(
                                    "setHydrogenBonds argument 'radius' cannot be set with argument 'dashed' set to True. \
                                Setting argument 'dashed' to False.".split()
                                )
                            )
                        viewer.addCylinder(
                            {
                                "radius": self.radius,
                                "color": self.color,
                                "fromCap": True,
                                "toCap": True,
                                "start": {
                                    "x": don_xyz[0],
                                    "y": don_xyz[1],
                                    "z": don_xyz[2],
                                },
                                "end": {
                                    "x": acc_xyz[0],
                                    "y": acc_xyz[1],
                                    "z": acc_xyz[2],
                                },
                            }
                        )
                    else:
                        viewer.addLine(
                            {
                                "dashed": self.dashed,
                                "color": self.color,
                                "start": {
                                    "x": don_xyz[0],
                                    "y": don_xyz[1],
                                    "z": don_xyz[2],
                                },
                                "end": {
                                    "x": acc_xyz[0],
                                    "y": acc_xyz[1],
                                    "z": acc_xyz[2],
                                },
                            }
                        )

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        hbond_set = pose.get_hbonds()
        selection_hbonds = []
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_res = r.don_res()
                    don_hatm_name = (
                        pose.residue(don_res).atom_name(r.don_hatm()).strip()
                    )
                    don_residue, don_chain = _get_residue_chain_tuple(pose, don_res)
                    don_sele = f"{don_residue}:{don_chain}.{don_hatm_name}"
                    acc_res = r.acc_res()
                    acc_atm_name = pose.residue(acc_res).atom_name(r.acc_atm()).strip()
                    acc_residue, acc_chain = _get_residue_chain_tuple(pose, acc_res)
                    acc_sele = f"{acc_residue}:{acc_chain}.{acc_atm_name}"
                    selection_hbonds.append([don_sele, acc_sele])

        viewer.add_distance(
            atom_pair=selection_hbonds,
            color=self.color,
            label_visible=False,
        )

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setHydrogens(ModuleBase):
    """
    Show all or only polar hydrogen atoms in each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `color`

        `str` indicating a standard color (e.g. "grey").
        Default: "white"

    second : optional
        `radius`

        `float` or `int` indicating the radius of the hydrogen atom stick represnetations.
        Default: 0.05

    third : optional
        `polar_only`

        `True` or `False`. `True` to show only polar hydrogen atoms, and `False` to show all hydrogen atoms.
        Default: False

    Returns
    -------
    A Viewer instance.
    """

    color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default="white"),
    )
    color_rgb = attr.ib(
        default=attr.Factory(
            lambda self: ImageColor.getcolor(self.color, "RGB"), takes_self=True
        ),
        type=Tuple[int, int, int],
        validator=attr.validators.instance_of(tuple),
        init=False,
    )
    radius = attr.ib(
        default=0.05,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    polar_only = attr.ib(
        default=None,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    residue_selector = attr.ib(
        default=None,
        type=Optional[ResidueSelector],
        validator=attr.validators.optional(
            attr.validators.instance_of(ResidueSelector)
        ),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )

    def _addCylinder(self, _viewer, i_xyz, j_xyz):
        _viewer.addCylinder(
            {
                "radius": self.radius,
                "color": self.color,
                "fromCap": True,
                "toCap": True,
                "start": {"x": i_xyz[0], "y": i_xyz[1], "z": i_xyz[2]},
                "end": {"x": j_xyz[0], "y": j_xyz[1], "z": j_xyz[2]},
            }
        )
        return _viewer

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(
            pose, self.residue_selector, logger=_logger
        )
        residue_chain_tuples = list(zip(map(str, resi), chain))
        if pose.is_fullatom():
            for i in range(1, pose.total_residue() + 1):
                residue_chain_tuple = tuple(_get_residue_chain_tuple(pose, i))
                if residue_chain_tuple in residue_chain_tuples:
                    r = pose.residue(i)
                    h_begin = r.attached_H_begin()
                    h_end = r.attached_H_end()
                    for h in range(1, len(h_begin) + 1):
                        i_index = h_begin[h]
                        j_index = h_end[h]
                        if all(q != 0 for q in [i_index, j_index]):
                            i_xyz = r.atom(h).xyz()
                            for j in range(i_index, j_index + 1):
                                if self.polar_only:
                                    if r.atom_is_polar_hydrogen(j):
                                        j_xyz = r.atom(j).xyz()
                                        viewer = self._addCylinder(viewer, i_xyz, j_xyz)
                                else:
                                    j_xyz = r.atom(j).xyz()
                                    viewer = self._addCylinder(viewer, i_xyz, j_xyz)

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(
            pose, self.residue_selector, logger=_logger
        )
        residue_chain_tuples = list(zip(map(str, resi), chain))
        if pose.is_fullatom():
            selection = []
            for i in range(1, pose.total_residue() + 1):
                residue, chain = _get_residue_chain_tuple(pose, i)
                if (residue, chain) in residue_chain_tuples:
                    r = pose.residue(i)
                    h_begin = r.attached_H_begin()
                    h_end = r.attached_H_end()
                    for h in range(1, len(h_begin) + 1):
                        i_index = h_begin[h]
                        j_index = h_end[h]
                        if all(q != 0 for q in [i_index, j_index]):
                            i_name = r.atom_name(h).strip()
                            i_sele = f"{residue}:{chain}.{i_name}"
                            for j in range(i_index, j_index + 1):
                                if self.polar_only:
                                    if r.atom_is_polar_hydrogen(j):
                                        j_name = r.atom_name(j).strip()
                                        j_sele = f"{residue}:{chain}.{j_name}"
                                        sele = f"({i_sele} or {j_sele})"
                                        selection.append(sele)
                                else:
                                    j_name = r.atom_name(j).strip()
                                    j_sele = f"{residue}:{chain}.{j_name}"
                                    sele = f"{i_sele} or {j_sele}"
                                    selection.append(sele)

            selection_hydrogens = " or ".join(selection)
            viewer.add_representation(
                repr_type="line",
                selection=selection_hydrogens,
                color=self.color,
                radius=self.radius,
                component=model,
            )

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=True, slots=True)
class setStyle(ModuleBase):
    """
    Show and color cartoon, and/or show heavy atoms with provided style, color and radius for each initialized
    `.pdb` file, `Pose` or `PackedPose` object. If the `residue_selector` argument is provided, apply styles
    only to the selected residues. If the `command` argument is provided, override all other arguments and pass
    `py3Dmol.view.setStyle()` commands to the Viewer.

    Parameters
    ----------
    first : optional
        `residue_selector`

        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` on which to apply the style(s).
        Default: None

    second : optional
        `cartoon`

        `True` or `False` to show cartoon representation.
        Default: True

    third : optional
        `cartoon_color`

        Hexcode literal (e.g. 0xAF10AB) or `str` indicating a standard color (e.g. "grey") for the cartoon representation.
        If "spectrum", apply reversed color gradient based on residue numbers. The option `cartoon` must also be set to `True`.
        Default: "spectrum"
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorSpec

    fourth : optional
        `style`

        `str` indicating a representation style of heavy atoms, choosing from either "line", "cross", "stick", or "sphere".
        Default: "stick"

    fifth : optional
        `colorscheme`

        `str` indicating the color scheme for heavy atoms represented by the `style` option. Options include:
            A lower-case standard color followed by "Carbon" (e.g. "orangeCarbon")
            "ssPyMOL": PyMol secondary colorscheme
            "ssJmol": Jmol secondary colorscheme
            "Jmol": Jmol primary colorscheme
            "default": default colorscheme
            "amino": amino acid colorscheme
            "shapely": shapely protien colorscheme
            "nucleic": nucleic acid colorscheme
            "chain": standard chain colorscheme
            "chainHetatm": chain Hetatm colorscheme
        Default: "blackCarbon"
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec

    sixth : optional
        `radius`

        `float` or `int` indicating the radius of the heavy atoms represented by the `style` option.
        Default: 0.1

    seventh : optional
        `label`

        `True` or `False` to show labels next to residues selected by the `residue_selector` option.
        Default: True

    eighth : optional
        `label_fontsize`

        `int` or `float` indicating the font size of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: 12

    ninth : optional
        `label_background`

        `True` or `False` to show the background of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: False

    tenth : optional
        `label_fontcolor`

        `str` indicating a standard font color (e.g. "grey") for label text next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: "black"

    eleventh : optional
        `command`

        `dict` or `tuple` of `dict`s of `py3Dmol.view.setStyle()` commands. If specified, this option overrides all other options.
        Default: None
        Example:
            command = {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "greyCarbon", "radius": 0.15}}
            view = viewer.init(poses) + viewer.setStyle(command=command)
            view.show()

    Returns
    -------
    A Viewer instance.
    """

    residue_selector = attr.ib(
        default=None,
        type=Optional[ResidueSelector],
        validator=attr.validators.optional(
            attr.validators.instance_of(ResidueSelector)
        ),
    )
    cartoon = attr.ib(
        default=True,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    cartoon_color = attr.ib(
        default=None,
        type=Optional[Union[str, int]],
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
    )
    cartoon_radius = attr.ib(
        default=None,
        type=Optional[Union[float, int]],
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=attr.converters.default_if_none(default=0.1),
    )
    cartoon_opacity = attr.ib(
        default=1,
        type=Union[float, int],
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=[_to_0_if_le_0, _to_1_if_gt_1],
    )
    style = attr.ib(
        default="stick",
        type=str,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("line", "cross", "stick", "sphere", "ball+stick")),
        ],
    )
    colorscheme = attr.ib(
        default=None,
        type=str,
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(default="blackCarbon"),
    )
    radius = attr.ib(
        default=0.1,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    label = attr.ib(
        default=True,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    label_fontsize = attr.ib(
        default=12,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    label_background = attr.ib(
        default=None,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    label_fontcolor = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default="black"),
    )
    command = attr.ib(
        default=None,
        type=Optional[Union[Tuple[Dict], Dict]],
        validator=attr.validators.optional(attr.validators.instance_of((tuple, dict))),
    )
    show_hydrogens = attr.ib(
        default=None,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if self.show_hydrogens:
            _logger.warning(
                "The 'show_hydrogens' attribute is not supported for the `py3Dmol` backend. "
                "Please add `viewer3d.setHydrogens()` instead."
            )
        if self.cartoon_color is None:
            self.cartoon_color = "spectrum"  # Set default

        if self.command:
            if isinstance(self.command, tuple):
                viewer.setStyle(*self.command)
            elif isinstance(self.command, dict):
                viewer.setStyle(self.command)
        else:
            if self.residue_selector:
                if pose is None:
                    pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

                resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

                if (not resi) and (not chain):
                    pass
                else:
                    if self.cartoon:
                        viewer.setStyle(
                            {"model": model, "resi": resi, "chain": chain},
                            {
                                "cartoon": {"color": self.cartoon_color},
                                self.style: {
                                    "colorscheme": self.colorscheme,
                                    "radius": self.radius,
                                },
                            },
                        )
                    else:

                        viewer.setStyle(
                            {"model": model, "resi": resi, "chain": chain},
                            {
                                self.style: {
                                    "colorscheme": self.colorscheme,
                                    "radius": self.radius,
                                }
                            },
                        )
                    if self.label:
                        viewer.addResLabels(
                            {"model": model, "resi": resi, "chain": chain},
                            {
                                "fontSize": self.label_fontsize,
                                "showBackground": self.label_background,
                                "fontColor": self.label_fontcolor,
                            },
                        )
            else:
                if self.cartoon:
                    viewer.setStyle(
                        {"model": model},
                        {
                            "cartoon": {"color": self.cartoon_color},
                            self.style: {
                                "colorscheme": self.colorscheme,
                                "radius": self.radius,
                            },
                        },
                    )
                else:
                    viewer.remove_cartoon()
                    viewer.setStyle(
                        {"model": model},
                        {
                            self.style: {
                                "colorscheme": self.colorscheme,
                                "radius": self.radius,
                            }
                        },
                    )

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        # Set defaults
        if self.cartoon_color is None:
            self.cartoon_color = "atomindex"
        if self.colorscheme == "blackCarbon":
            self.colorscheme = "element"
        if self.style == "stick":
            self.style = "licorice"
        elif self.style == "sphere":
            self.style = "spacefill"
        elif self.style == "cross":
            self.style = "point"
        elif self.style == "line":
            self.style = "line"
        default_selection = "*" if self.show_hydrogens else "not hydrogen"

        if self.command:
            _logger.warning(
                "The 'command' attribute is not supported for the `nglview` backend."
            )
        if self.residue_selector is not None:
            if pose is None:
                pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

            selection = _get_nglview_selection(
                pose, self.residue_selector, logger=_logger
            )
            selection_hydrogens = (
                f"({selection}) and not hydrogen" if self.show_hydrogens else selection
            )

            if not selection:
                pass
            else:
                viewer.add_representation(
                    repr_type=self.style,
                    selection=selection_hydrogens,
                    color=self.colorscheme,
                    radius=self.radius,
                    component=model,
                )
                if self.cartoon:
                    viewer.remove_cartoon(component=model)
                    viewer.add_representation(
                        repr_type="cartoon",
                        selection=selection,
                        color=self.cartoon_color,
                        radius=self.cartoon_radius,
                        opacity=self.cartoon_opacity,
                        component=model,
                    )
                if self.label:
                    viewer.add_representation(
                        repr_type="label",
                        labelType="res",
                        selection=selection_hydrogens,
                        showBorder=self.label_background,
                        borderColor="gray",
                        showBackground=self.label_background,
                        backgroundColor=self.colorscheme,
                        component=model,
                    )
        else:
            viewer.add_representation(
                repr_type=self.style,
                selection=default_selection,
                color=self.colorscheme,
                radius=self.radius,
                component=model,
            )
            if self.cartoon:
                viewer.remove_cartoon(component=model)
                viewer.add_representation(
                    repr_type="cartoon",
                    selection=default_selection,
                    color=self.cartoon_color,
                    radius=self.cartoon_radius,
                    opacity=self.cartoon_opacity,
                    component=model,
                )

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setSurface(ModuleBase):
    """
    Show the specified surface for each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `residue_selector`

        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` to select residues
        on which to apply the surface.
        Default: None

    second : optional
        `surface_type`

        `str` indicating surface type to be displayed. py3Dmol supports the following options:
            "VDW": Van der Waals surface
            "MS": Molecular surface
            "SES": Solvent excluded surface
            "SAS": Solvent accessible surface
            "AV": High quality molecular surface (only supported for `nglview` backend)
        Default: "VDW"

    third : optional
        `opacity`

        `float` or `int` between 0 and 1 for opacity of the displayed surface. Not currently supported for `nglview` backend.
        Default: 0.5

    fourth : optional
        `color`

        `str` indicating a standard color (e.g. "grey") of the surface to be displayed.
        Either `color` or `colorscheme` may be specified, where `colorscheme` overrides `color`.
        Default: None

    fifth : optional
        `colorscheme`

        `str` indicating the color scheme of the surface to be displayed.
        Either `color` or `colorscheme` may be specified, where `colorscheme` overrides `color`.
        Options include:
            A lower-case standard color followed by "Carbon" (e.g. "yellowCarbon")
            "ssPyMOL": PyMol secondary colorscheme
            "ssJmol": Jmol secondary colorscheme
            "Jmol": Jmol primary colorscheme
            "default": default colorscheme
            "amino": amino acid colorscheme
            "shapely": shapely protien colorscheme
            "nucleic": nucleic acid colorscheme
            "chain": standard chain colorscheme
            "chainHetatm": chain Hetatm colorscheme
        Default: None
        Reference: https://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec
        Not currently supported for `nglview` backend.

    Returns
    -------
    A Viewer instance.
    """

    residue_selector = attr.ib(
        default=None,
        type=ResidueSelector,
        validator=attr.validators.instance_of(ResidueSelector),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )
    surface_type = attr.ib(
        default="VDW",
        type=str,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("VDW", "MS", "SAS", "SES", "AV")),
        ],
    )
    opacity = attr.ib(
        default=0.5,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=[_to_0_if_le_0, _to_1_if_gt_1],
    )
    color = attr.ib(
        default=None,
        type=Optional[Union[str, int]],
        validator=attr.validators.instance_of((str, int, type(None))),
    )
    colorscheme = attr.ib(
        default=None,
        type=Optional[str],
        validator=attr.validators.instance_of((str, type(None))),
    )

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if self.surface_type == "AV":
            raise NotImplementedError(
                "The surface type 'AV' is not supported by the `py3Dmol` backend."
            )

        py3Dmol = sys.modules["py3Dmol"]
        surface_types_dict: Dict[str, int] = {
            "VDW": py3Dmol.VDW,
            "MS": py3Dmol.MS,
            "SAS": py3Dmol.SAS,
            "SES": py3Dmol.SES,
        }
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            if self.colorscheme:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity, "colorscheme": self.colorscheme},
                    {"model": model, "resi": resi, "chain": chain},
                )
            elif self.color:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity, "color": self.color},
                    {"model": model, "resi": resi, "chain": chain},
                )
            else:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity},
                    {"model": model, "resi": resi, "chain": chain},
                )

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        surface_types_dict: Dict[str, str] = {
            "VDW": "vws",
            "MS": "ms",
            "SAS": "sas",
            "SES": "ses",
            "AV": "av",
        }
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_nglview_selection(pose, self.residue_selector, logger=_logger)
        if not selection:
            pass
        else:
            viewer.add_representation(
                repr_type="surface",
                selection=selection,
                surfaceType=surface_types_dict[self.surface_type],
                color=self.color,
                component=model,
            )

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setZoom(ModuleBase):
    """
    Set the zoom magnification factor of each initialized `.pdb` file, `Pose` or `PackedPose` object.
    For the `py3Dmol` backend, values >1 zoom in, and values <1 zoom out.
    For the `nglview` backend, values >0 zoom in, and values <0 zoom out.

    Parameters
    ----------
    first : optional
        `factor`

        `float` or `int` indicating the zoom magnification factor.
        Default: 2

    Returns
    -------
    A Viewer instance.
    """

    factor = attr.ib(
        default=2,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
    )

    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        viewer.zoom(self.factor)
        return viewer

    def apply_nglview(self, viewer, pose, pdbstring, model):
        viewer.control.zoom(self.factor)
        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setZoomTo(ModuleBase):
    """
    Zoom to a `ResidueSelector` in each initialized `.pdb` file, `Pose` or `PackedPose` object.

    Parameters
    ----------
    first : optional
        `residue_selector`

        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` into which to zoom.
        Default: None

    Returns
    -------
    A Viewer instance.
    """

    residue_selector = attr.ib(
        default=None,
        type=ResidueSelector,
        validator=attr.validators.instance_of(ResidueSelector),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )

    @requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            viewer.zoomTo({"model": model}, {"resi": resi, "chain": chain})

        return viewer

    @requires_init
    def apply_nglview(self, viewer, pose, pdbstring, model):
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_nglview_selection(pose, self.residue_selector, logger=_logger)
        if not selection:
            viewer.center(selection="*", component=model)
        else:
            viewer.center(selection=selection, component=model)

        return viewer

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])
