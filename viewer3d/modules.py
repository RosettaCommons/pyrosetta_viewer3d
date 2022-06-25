import attr
import bokeh.palettes
import collections
import copy
import itertools
import logging
import numpy
import pyrosetta
import pyrosetta.distributed.io as io
import sys
import uuid

from functools import singledispatch
from pyrosetta import Pose
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.conformation import Residue, is_disulfide_bond
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector,
)

from typing import (
    Dict,
    Generic,
    Iterable,
    List,
    NoReturn,
    Optional,
    OrderedDict,
    Tuple,
    TypeVar,
    Union,
)

from viewer3d.colors import default_element_colors
from viewer3d.config import BACKENDS
from viewer3d.converters import (
    _get_nglview_selection,
    _get_residue_chain_tuple,
    _int_to_str,
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
    _py3Dmol_to_nglview_style,
    _to_hex,
    _to_0_if_le_0,
    _to_1_if_gt_1,
)
from viewer3d.exceptions import ModuleNotImplementedError
from viewer3d.tracer import requires_init


_logger: logging.Logger = logging.getLogger("viewer3d.modules")
ViewerType = TypeVar("ViewerType")


@attr.s(kw_only=False, slots=True)
class ModuleBase:
    @staticmethod
    def _to_modules(objs) -> List["ModuleBase"]:
        @singledispatch
        def _to_module(obj):
            raise ValueError(
                "The 'modules' viewer attribute must be an instance of `ModuleBase` "
                + f"or an iterable of `ModuleBase` objects. Received: {type(obj)}"
            )

        _to_module.register(ModuleBase, lambda obj: obj)

        if objs is None:
            _modules = []
        elif isinstance(objs, collections.abc.Iterable):
            _modules = list(map(_to_module, objs))
        else:
            _modules = [_to_module(objs)]

        return _modules

    def add_selection_scheme(
        self, name: str, selection_scheme: List[List[Union[str, int]]]
    ) -> None:
        cm = sys.modules["nglview"].color.ColormakerRegistry
        cm.add_selection_scheme(name, selection_scheme)

    def add_element_selection_scheme(self, name: str) -> None:
        """
        Add element-based selection scheme to NGLView ColormakerRegistry
        based on '<color>Carbon' naming system or an arbitrary `str` object.
        """
        _default_element_colors = copy.deepcopy(default_element_colors)
        if name.endswith("Carbon"):
            _default_element_colors["C"] = name.replace("Carbon", "")
        else:
            _default_element_colors["C"] = name
        _selection_scheme = [
            [_int_to_str(color), f"_{element}"]
            for (element, color) in _default_element_colors.items()
        ]
        self.add_selection_scheme(name, _selection_scheme)

    # def add_per_residue_element_selection_scheme(
    #     self, name: str, color: str, residue_chain_tuples: List[Tuple[str, str]]
    # ) -> None:
    #     """
    #     Add per-residue, per-element selection scheme to NGLView ColormakerRegistry.
    #     """
    #     _default_element_colors = copy.deepcopy(default_element_colors)
    #     _default_element_colors["C"] = color
    #     _selection_scheme = [
    #         [_int_to_str(color), f"{resi}:{chain}_{element}"]
    #         for (element, color) in _default_element_colors.items()
    #         for (resi, chain) in residue_chain_tuples
    #     ]
    #     cm = sys.modules["nglview"].color.ColormakerRegistry
    #     cm.add_selection_scheme(name, _selection_scheme)


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
        converter=[attr.converters.default_if_none(default=0xFFFFFFFF), _to_hex],
    )

    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        viewer.setBackgroundColor(self.color)

        return viewer

    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        viewer.background = self.color

        return viewer

    def apply_pymol(self) -> NoReturn:
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
        converter=[attr.converters.default_if_none(default="gold"), _to_hex],
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
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        cys_res = [i for i, aa in enumerate(pose.sequence(), start=1) if aa == "C"]
        selection_disulfides = []
        for (i, j) in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_res, i_chain = _get_residue_chain_tuple(pose, i)
                j_res, j_chain = _get_residue_chain_tuple(pose, j)
                i_sele = f"{i_res}:{i_chain}.{self.sulfur_atom_name}"
                j_sele = f"{j_res}:{j_chain}.{self.sulfur_atom_name}"
                selection_disulfides.append([i_sele, j_sele])
        if selection_disulfides:
            viewer.add_distance(
                atom_pair=selection_disulfides,
                color=self.color,
                radius=self.radius,
                label_visible=False,
            )
            selection = " or ".join(
                [f"({s[0]} or {s[1]})" for s in selection_disulfides]
            )
            viewer.add_representation(
                repr_type="ball+stick",
                selection=selection,
                color=self.color,
                radius=self.radius,
                component=model,
            )

        return viewer

    def apply_pymol(self) -> NoReturn:
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
        converter=[attr.converters.default_if_none(default="black"), _to_hex],
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
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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

    def apply_pymol(self) -> NoReturn:
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

    fourth : optional
        `residue_selector`

        An instance of `pyrosetta.rosetta.core.select.residue_selector.ResidueSelector` on which to apply the style(s).
        Default: None

    Returns
    -------
    A Viewer instance.
    """

    color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="white"), _to_hex],
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
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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

    def apply_pymol(self) -> NoReturn:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=True, slots=True)
class setPerResidueRealMetric(ModuleBase):
    """
    Show and color residues by a per-residue real metric in the `pose.scores` dictionary.

    first: required
        `scoretype`

        `str` object matching a scoretype name. The scoretype name must contain a
        trailing underscore followed by either a residue number (rosetta numbering) or
        pdb number (pdb numbering), per the default output of `PerResidueRealMetrics`.
        Examples:
            Scoretype "res_energy_18A" would be selected by `scoretype='res_energy'`
            Scoretype "custom_type_atomic_clashes_4" would be selected by `scoretype='atomic_clashes'`

    second: optional
        `vmin`

        `float` or `int` object representing the minimum value for color mapping
        to scoretype values. If `None`, set 'vmin' to the minimum scoretype value.
        Default: None

    third: optional
        `vmax`

        `float` or `int` object representing the maximum value for color mapping
        to scoretype values. If `None`, set 'vmin' to the maximum scoretype value.
        Default: None

    fourth: optional
        `palette`

        An iterable of `str` (or `int`) objects representing a color map.
        Default: `bokeh.palettes.Greens256`

    fifth: optional
        `log`

        `None` for linear color mapping of 'vmin' to 'vmax'. If an `int`
        or `float` object is provided, map colors spaced evenly on a log
        scale with the base provided.

    sixth: optional:
        `style`

        `str` indicating a representation style of heavy atoms, choosing from
        either "line", "cross", "stick", or "sphere".
        Default: "stick"

    seventh : optional
        `radius`

        `float` or `int` indicating the radius of the heavy atoms represented by
        the `style` option.
        Default: 0.1

    eighth : optional
        `show_hydrogens`

        `bool` object for the `nglview` backend to show all hydrogens.
        Defualt: False

    ninth : optional
        `bonds`

        `str` object for the `nglview` backend to show double bonds.
        Available options are: "off", "symmetric", or "offset"
        Defualt: "symmetric"

    Returns
    -------
    A Viewer instance.
    """

    scoretype = attr.ib(
        type=str,
        validator=attr.validators.instance_of(str),
    )
    vmin = attr.ib(
        default=None,
        type=Union[float, int],
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
    )
    vmax = attr.ib(
        default=None,
        type=Union[float, int],
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
    )
    palette = attr.ib(
        default=None,
        type=Iterable[Union[str, int]],
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of((str, int)),
            iterable_validator=attr.validators.instance_of(collections.abc.Iterable),
        ),
        converter=attr.converters.default_if_none(default=bokeh.palettes.Greens256),
    )
    log = attr.ib(
        default=None,
        type=Optional[Union[float, int]],
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=_to_0_if_le_0,
    )
    style = attr.ib(
        default="stick",
        type=str,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("line", "cross", "stick", "sphere", "ball+stick")),
        ],
    )
    radius = attr.ib(
        default=0.1,
        type=Union[float, int],
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    show_hydrogens = attr.ib(
        default=None,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    bonds = attr.ib(
        default=None,
        type=bool,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("off", "symmetric", "offset")),
        ],
        converter=attr.converters.default_if_none(default="symmetric"),
    )

    def set_vmin_vmax(self, pose: Pose) -> None:
        values = [
            value
            for (scoretype, value) in pose.scores.items()
            if self.scoretype in scoretype
        ]
        if not values:
            raise ValueError(
                f"Scoretype matching '{self.scoretype}' not found in `pose.scores` keys."
            )
        if self.vmin is None:
            self.vmin = min(values)
        if self.vmax is None:
            self.vmax = max(values)

    def get_elements_from_residue(self, residue: Residue) -> List[str]:
        residue_type = residue.type()
        elements = [
            residue_type.element(atom).name.upper()
            for atom in range(1, residue.natoms() + 1)
        ]
        return elements

    def get_nearest_value_from_keys(
        self, _dict: Dict[float, str], _value: float
    ) -> float:
        index, value = min(enumerate(_dict.keys()), key=lambda x: abs(x[1] - _value))

        return value

    def get_palette_value_dict(self) -> OrderedDict[float, str]:
        if self.log is not None:
            _space = numpy.logspace(
                self.vmin, self.vmax, len(self.palette), base=self.log
            )
        else:
            _space = numpy.linspace(self.vmin, self.vmax, len(self.palette))
        _palette_value_dict = collections.OrderedDict(zip(_space, self.palette))

        return _palette_value_dict

    @requires_init
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        self.set_vmin_vmax(pose)
        _palette_value_dict = self.get_palette_value_dict()
        _selection_scheme = []
        for scoretype, value in pose.scores.items():
            _r = scoretype.split("_")[-1]
            if _r.isdigit():
                resi, chain = _get_residue_chain_tuple(pose, int(_r))
            else:
                resi, chain = _r[:-1], _r[-1]
            _nearest_value = self.get_nearest_value_from_keys(
                _palette_value_dict, value
            )
            _C_color = _palette_value_dict[_nearest_value]
            # TODO Set element-based color scheme
            viewer.setStyle(
                {"model": model, "resi": resi, "chain": chain},
                {
                    self.style: {
                        "color": _C_color,
                        "radius": self.radius,
                    }
                },
            )

        return viewer

    @requires_init
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        self.style = _py3Dmol_to_nglview_style(self.style)

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        self.set_vmin_vmax(pose)
        _palette_value_dict = self.get_palette_value_dict()
        _default_element_colors = copy.deepcopy(default_element_colors)
        _selection_scheme = []
        for scoretype, value in pose.scores.items():
            if self.scoretype in scoretype:
                _r = scoretype.split("_")[-1]
                if _r.isdigit():
                    _residue, _chain = _get_residue_chain_tuple(pose, int(_r))
                else:
                    _residue, _chain = _r[:-1], _r[-1]
                _resnum = pose.pdb_info().pdb2pose(_chain, int(_residue))
                _elements_from_residue = self.get_elements_from_residue(
                    pose.residue(_resnum)
                )
                _nearest_value = self.get_nearest_value_from_keys(
                    _palette_value_dict, value
                )
                _C_color = _palette_value_dict[_nearest_value]
                _default_element_colors["C"] = _C_color
                for (_element, _element_color) in _default_element_colors.items():
                    if _element in _elements_from_residue:
                        _selection_scheme.append(
                            [
                                _int_to_str(_element_color),
                                f"{_residue}:{_chain} and _{_element}",
                            ]
                        )
        _selection_name = f"{self.scoretype}_{uuid.uuid4().hex}"
        self.add_selection_scheme(_selection_name, _selection_scheme)
        _default_selection = "*" if self.show_hydrogens else "not hydrogen"
        viewer.add_representation(
            repr_type=self.style,
            selection=_default_selection,
            color=_selection_name,
            radius=self.radius,
            multipleBond=self.bonds,
            component=model,
        )

        return viewer

    def apply_pymol(self) -> NoReturn:
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
        `cartoon_radius`

        Set the cartoon radius for the `nglview` backend.

    fifth : optional
        `cartoon_opacity`

        Set the cartoon opacity for the `nglview` backend.

    sixth : optional
        `style`

        `str` indicating a representation style of heavy atoms, choosing from either "line", "cross", "stick", or "sphere".
        Default: "stick"

    seventh : optional
        `colorscheme`

        For the `py3Dmol` backend, a `str` indicating the color scheme for
        heavy atoms represented by the `style` option.
        Options include:
            A lowercased standard color optionally followed by "Carbon" (e.g. "orangeCarbon")
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

    For the `nglview` backend, a `str` indicating the color scheme for
    heavy atoms represented by the `style` option.
    Options include:
        A lowercased standard color.
        "atomindex"
        "bfactor"
        "chainid"
        "chainindex"
        "chainname"
        "densityfit"
        "electrostatic"
        "element"
        "entityindex"
        "entitytype"
        "geoquality"
        "hydrophobicity"
        "modelindex"
        "moleculetype"
        "occupancy"
        "random"
        "residueindex"
        "resname"
        "sstruc"
        "uniform"
        "value"
        "volume"
    Default: "element"
    Reference: https://nglviewer.org/ngl/api/manual/coloring.html

    eighth : optional
        `radius`

        `float` or `int` indicating the radius of the heavy atoms represented by the `style` option.
        Default: 0.1

    ninth : optional
        `label`

        `True` or `False` to show labels next to residues selected by the `residue_selector` option.
        Default: True

    tenth : optional
        `label_fontsize`

        `int` or `float` indicating the font size of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: 12

    eleventh : optional
        `label_background`

        `True` or `False` to show the background of labels next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: False

    twelfth : optional
        `label_fontcolor`

        `str` indicating a standard font color (e.g. "grey") for label text next to residues selected by the `residue_selector` option,
        only if `label` is `True`.
        Default: "black"

    thirteenth : optional
        `command`

        `dict` or `tuple` of `dict`s of `py3Dmol.view.setStyle()` commands. If specified, this option overrides all other options.
        Default: None
        Example:
            command = {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "greyCarbon", "radius": 0.15}}
            view = viewer.init(poses) + viewer.setStyle(command=command)
            view.show()

    fourteenth : optional
        `show_hydrogens`

        `bool` object for the `nglview` backend to show all hydrogens.
        Defualt: False

    fifteenth : optional
        `bonds`

        `str` object for the `nglview` backend to show double bonds.
        Options are: "off", "symmetric", "offset"
        Defualt: "symmetric"

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
        converter=_to_hex,
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
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="blackCarbon"), _to_hex],
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
        converter=[attr.converters.default_if_none(default="black"), _to_hex],
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
    bonds = attr.ib(
        default=None,
        type=bool,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("off", "symmetric", "offset")),
        ],
        converter=attr.converters.default_if_none(default="symmetric"),
    )

    @requires_init
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        _colorscheme = "colorscheme" if isinstance(self.colorscheme, str) else "color"
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
                                    _colorscheme: self.colorscheme,
                                    "radius": self.radius,
                                },
                            },
                        )
                    else:
                        viewer.setStyle(
                            {"model": model, "resi": resi, "chain": chain},
                            {
                                self.style: {
                                    _colorscheme: self.colorscheme,
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
                                _colorscheme: self.colorscheme,
                                "radius": self.radius,
                            },
                        },
                    )
                else:
                    viewer.setStyle(
                        {"model": model},
                        {
                            self.style: {
                                _colorscheme: self.colorscheme,
                                "radius": self.radius,
                            }
                        },
                    )

        return viewer

    @requires_init
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        # Set defaults
        if self.cartoon_color is None:
            self.cartoon_color = "atomindex"
        if isinstance(self.colorscheme, str) and self.colorscheme.endswith("Carbon"):
            self.add_element_selection_scheme(self.colorscheme)
        self.style = _py3Dmol_to_nglview_style(self.style)
        default_selection = "*" if self.show_hydrogens else "not hydrogen"

        if self.command:
            _logger.warning(
                "The 'command' attribute is not supported for the `nglview` backend."
            )
        if self.residue_selector is not None:
            if pose is None:
                pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

            selection = _get_nglview_selection(
                pose,
                self.residue_selector,
                show_hydrogens=self.show_hydrogens,
                logger=_logger,
            )

            if not selection:
                pass
            else:
                if self.radius > 1e-10:
                    viewer.add_representation(
                        repr_type=self.style,
                        selection=selection,
                        color=self.colorscheme,
                        radius=self.radius,
                        multipleBond=self.bonds,
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
                        selection=selection,  # selection_hydrogens,
                        showBorder=self.label_background,
                        borderColor="gray",
                        showBackground=self.label_background,
                        backgroundColor=self.colorscheme,
                        component=model,
                    )
        else:
            if self.radius > 1e-10:
                viewer.add_representation(
                    repr_type=self.style,
                    selection=default_selection,
                    color=self.colorscheme,
                    radius=self.radius,
                    multipleBond=self.bonds,
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

    def apply_pymol(self) -> NoReturn:
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
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
        converter=_to_hex,
    )
    colorscheme = attr.ib(
        default=None,
        type=Optional[Union[str, int]],
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
        converter=_to_hex,
    )

    @requires_init
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
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
            if self.colorscheme is not None:
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
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        surface_types_dict: Dict[str, str] = {
            "VDW": "vws",
            "MS": "ms",
            "SAS": "sas",
            "SES": "ses",
            "AV": "av",
        }
        for _colorscheme in (self.color, self.colorscheme):
            if isinstance(_colorscheme, str) and _colorscheme.endswith("Carbon"):
                self.add_element_selection_scheme(_colorscheme)

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_nglview_selection(
            pose, self.residue_selector, show_hydrogens=True, logger=_logger
        )
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

        return viewer

    def apply_pymol(self) -> NoReturn:
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

    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        viewer.zoom(self.factor)
        return viewer

    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        viewer.control.zoom(self.factor)
        return viewer

    def apply_pymol(self) -> NoReturn:
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
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            viewer.zoomTo({"model": model, "resi": resi, "chain": chain})

        return viewer

    @requires_init
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_nglview_selection(
            pose, self.residue_selector, show_hydrogens=True, logger=_logger
        )
        if not selection:
            viewer.center(selection="*", component=model)
        else:
            viewer.center(selection=selection, component=model)

        return viewer

    def apply_pymol(self) -> NoReturn:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=True, slots=True, frozen=True)
class setTemplate(ModuleBase):
    """Template class for developing new visualization modules."""

    @requires_init
    def apply_py3Dmol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[0])

    @requires_init
    def apply_nglview(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

    @requires_init
    def apply_pymol(
        self, viewer: Generic[ViewerType], pose: Pose, pdbstring: str, model: int
    ) -> Generic[ViewerType]:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])
