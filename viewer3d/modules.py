import attr
import logging
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys

from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector,
)

from viewer3d.config import BACKENDS
from viewer3d.converters import (
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
    _to_0_if_le_0,
    _to_1_if_gt_1,
)
from viewer3d.exceptions import ModuleNotImplementedError

from typing import Dict, Optional, Tuple, Union

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

    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        viewer.setBackgroundColor(self.color)

        return viewer

    def apply_nglview(self, viewer, pose, pdbstring, **kwargs):
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
        converter=attr.converters.default_if_none(default=0.5),
    )

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        if pose is not None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        cys_res = []
        for i, aa in enumerate(pose.sequence(), start=1):
            if aa == "C":
                cys_res.append(i)
        for i in cys_res:
            for j in cys_res:
                if pyrosetta.rosetta.core.conformation.is_disulfide_bond(
                    pose.conformation(), i, j
                ):
                    i_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(
                            pose.residue(i).atom_index("SG"), i
                        )
                    )
                    j_xyz = pose.xyz(
                        pyrosetta.rosetta.core.id.AtomID(
                            pose.residue(j).atom_index("SG"), j
                        )
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

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

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

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        if pose is not None:
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

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

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

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        if pose is not None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        if pose.is_fullatom():
            for i in range(1, pose.total_residue() + 1):
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

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

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
        validator=attr.validators.instance_of((ResidueSelector, type(None))),
    )
    cartoon = attr.ib(
        default=True,
        type=bool,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    cartoon_color = attr.ib(
        default=None,
        type=Union[str, int],
        validator=attr.validators.instance_of((str, int)),
        converter=attr.converters.default_if_none(default="spectrum"),
    )
    style = attr.ib(
        default="stick",
        type=str,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("line", "cross", "stick", "sphere")),
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
        validator=attr.validators.instance_of((tuple, dict, type(None))),
    )

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        if self.command:
            if isinstance(self.command, tuple):
                viewer.setStyle(*self.command)
            elif isinstance(self.command, dict):
                viewer.setStyle(self.command)
        else:
            if self.residue_selector:
                if pose is not None:
                    pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

                resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

                if (not resi) and (not chain):
                    pass
                else:
                    if self.cartoon:
                        viewer.setStyle(
                            {"resi": resi, "chain": chain},
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
                            {"resi": resi, "chain": chain},
                            {
                                self.style: {
                                    "colorscheme": self.colorscheme,
                                    "radius": self.radius,
                                }
                            },
                        )
                    if self.label:
                        viewer.addResLabels(
                            {"resi": resi, "chain": chain},
                            {
                                "fontSize": self.label_fontsize,
                                "showBackground": self.label_background,
                                "fontColor": self.label_fontcolor,
                            },
                        )
            else:
                if self.cartoon:
                    viewer.setStyle(
                        {
                            "cartoon": {"color": self.cartoon_color},
                            self.style: {
                                "colorscheme": self.colorscheme,
                                "radius": self.radius,
                            },
                        }
                    )
                else:
                    viewer.setStyle(
                        {
                            self.style: {
                                "colorscheme": self.colorscheme,
                                "radius": self.radius,
                            }
                        }
                    )

        return viewer

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

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
        Default: "VDW"

    third : optional
        `opacity`

        `float` or `int` between 0 and 1 for opacity of the displayed surface.
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
            attr.validators.in_(("VDW", "MS", "SES", "SAS")),
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

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, surface_types_dict=None, **kwargs):
        if pose is not None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            if self.colorscheme:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity, "colorscheme": self.colorscheme},
                    {"resi": resi, "chain": chain},
                )
            elif self.color:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity, "color": self.color},
                    {"resi": resi, "chain": chain},
                )
            else:
                viewer.addSurface(
                    surface_types_dict[self.surface_type],
                    {"opacity": self.opacity},
                    {"resi": resi, "chain": chain},
                )

        return viewer

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])


@attr.s(kw_only=False, slots=True)
class setZoom(ModuleBase):
    """
    Set the zoom magnification factor of each initialized `.pdb` file, `Pose` or `PackedPose` object.
    Values >1 zoom in, and values <1 zoom out.

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

    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        viewer.zoom(self.factor)
        return viewer

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

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

    @pyrosetta.distributed.requires_init
    def apply_py3Dmol(self, viewer, pose, pdbstring, **kwargs):
        if pose is not None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector)

        if (not resi) and (not chain):
            pass
        else:
            viewer.zoomTo({"resi": resi, "chain": chain})

        return viewer

    def apply_nglview(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

    def apply_pymol(self):
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])
