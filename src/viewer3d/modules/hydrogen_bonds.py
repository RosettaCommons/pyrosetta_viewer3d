__author__ = "Jason C. Klima"

import attr
import logging

from pyrosetta import Pose

from viewer3d.converters import (
    _get_residue_chain_tuple,
    _pdbstring_to_pose,
    _to_hex,
    _to_0_if_le_0,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    GenericViewer,
    Optional,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.hydrogen_bonds")


@attr.define(kw_only=False, slots=True, frozen=True)
class setHydrogenBonds(ModuleBase):
    """
    Display hydrogen bonds according to `pyrosetta.Pose.get_hbonds` in each decoy.

    Attributes:
        color: A `str` indicating a standard color (e.g., `"yellow"`).

            Default: `"black"`

        dashed: A `bool` to show hydrogen bonds as dashed lines. If `True`, then the option
            `radius` must be `None`. If `False`, then the option `radius` must be specified.

            Default: `True`

        radius: A `float` or `int` indicating the radius of the solid (non-dashed) stick
            connecting the atoms participating in each hydrogen bond. If set, this
            automatically sets the option `dashed` to `False`.

            Default: `None`
    """

    color: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="black"), _to_hex],
    )
    dashed: bool = attr.field(
        default=True,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    radius: Optional[Union[float, int]] = attr.field(
        default=None,
        validator=attr.validators.instance_of((float, int, type(None))),
        converter=_to_0_if_le_0,
    )

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _color = "#000000" if self.color == "black" else self.color
        hbond_set = pose.get_hbonds()
        selection_hbonds = []
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_res = r.don_res()
                    don_hatm_name = pose.residue(don_res).atom_name(r.don_hatm()).strip()
                    don_residue, don_chain = _get_residue_chain_tuple(pose, don_res)
                    don_sele = f"{don_residue}:{don_chain}.{don_hatm_name}"
                    acc_res = r.acc_res()
                    acc_atm_name = pose.residue(acc_res).atom_name(r.acc_atm()).strip()
                    acc_residue, acc_chain = _get_residue_chain_tuple(pose, acc_res)
                    acc_sele = f"{acc_residue}:{acc_chain}.{acc_atm_name}"
                    selection_hbonds.append([don_sele, acc_sele])

        viewer.add_distance(
            atom_pair=selection_hbonds,
            color=_color,
            label_visible=False,
        )

        return viewer

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        hbond_id = 0
        hbond_objs = []
        hbond_set = pose.get_hbonds()
        for i in range(1, pose.total_residue() + 1):
            res_hbonds = hbond_set.residue_hbonds(i, False)
            if res_hbonds:
                for j in range(1, len(res_hbonds) + 1):
                    r = res_hbonds[j]
                    don_res = r.don_res()
                    don_hatm_name = pose.residue(don_res).atom_name(r.don_hatm()).strip()
                    don_residue, don_chain = _get_residue_chain_tuple(pose, don_res)
                    don_sele = f"(obj {model} and chain {don_chain} and resi {don_residue} and name {don_hatm_name})"
                    acc_res = r.acc_res()
                    acc_atm_name = pose.residue(acc_res).atom_name(r.acc_atm()).strip()
                    acc_residue, acc_chain = _get_residue_chain_tuple(pose, acc_res)
                    acc_sele = f"(obj {model} and chain {acc_chain} and resi {acc_residue} and name {acc_atm_name})"
                    hbond_obj = f"{model}_hbond_{hbond_id}"
                    viewer.distance(hbond_obj, don_sele, acc_sele, 9999.9, 0)
                    hbond_objs.append(hbond_obj)
                    hbond_id += 1

        if hbond_objs:
            group_name = f"{model}_hbonds"
            viewer.group(group_name, " ".join(hbond_objs))
            with self.out:
                if self.radius:
                    viewer.do(f"set dash_radius, {self.radius}, {group_name}")
                viewer.set("dash_color", self.color, group_name)
            viewer.hide("labels", group_name)

        return viewer
