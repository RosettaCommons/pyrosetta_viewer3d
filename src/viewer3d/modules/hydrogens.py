__author__ = "Jason C. Klima"

import attr
import logging

from pyrosetta import Pose
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector,
)

from viewer3d.converters import (
    _get_residue_chain_tuple,
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
    _to_hex,
    _to_0_if_le_0,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    GenericViewer,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.hydrogens")


@attr.define(kw_only=False, slots=True, frozen=True)
class setHydrogens(ModuleBase):
    """
    Show all or only polar hydrogen atoms in each decoy.

    Attributes:
        color: A `str` indicating a standard color (e.g., `"grey"`).

            Default: `"white"`

        radius: A `float` or `int` indicating the radius of the hydrogen atom stick
            representations. Note that for the `pymol` backend, this parameter controls the
            `"set_h_scale"` setting, which typically has a value of `0.4`.

            Default: `0.05`

        polar_only: A `bool` object to show only polar hydrogen atoms. If `True`, then show only
            polar hydrogen atoms, and if `False`, then show all hydrogen atoms.

            Default: `False`

        residue_selector: An instance of `ResidueSelector` on which to apply the style(s).

            Default: `None`
    """

    color: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="white"), _to_hex],
    )
    radius: Union[float, int] = attr.field(
        default=0.05,
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    polar_only: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    residue_selector: Optional[ResidueSelector] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(ResidueSelector)),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )

    def _addCylinder(
        self,
        _viewer: GenericViewer,
        i_xyz: Tuple[float, float, float],
        j_xyz: Tuple[float, float, float],
    ):
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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector, logger=_logger)
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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _color = "#000000" if self.color == "black" else self.color
        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector, logger=_logger)
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
                                    sele = f"({i_sele} or {j_sele})"
                                    selection.append(sele)

            selection_hydrogens = " or ".join(selection)
            viewer.add_representation(
                repr_type="line",
                selection=selection_hydrogens,
                color=_color,
                radius=self.radius,
                component=model,
            )

        return viewer

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        resi, chain = _pose_to_residue_chain_tuples(pose, self.residue_selector, logger=_logger)
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
                            i_sele = f"(obj {model} and chain {chain} and resi {residue} and name {i_name})"
                            for j in range(i_index, j_index + 1):
                                if self.polar_only:
                                    if r.atom_is_polar_hydrogen(j):
                                        j_name = r.atom_name(j).strip()
                                        j_sele = f"(obj {model} and chain {chain} and resi {residue} and name {j_name})"
                                        sele = f"({i_sele} or {j_sele})"
                                        selection.append(sele)
                                else:
                                    j_name = r.atom_name(j).strip()
                                    j_sele = f"(obj {model} and chain {chain} and resi {residue} and name {j_name})"
                                    sele = f"({i_sele} or {j_sele})"
                                    selection.append(sele)

        selection_hydrogens = " or ".join(selection)
        viewer.show("sticks", selection_hydrogens)
        viewer.color(self.color, f"({selection_hydrogens}) and (elem h)")
        viewer.set("stick_h_scale", self.radius)

        return viewer
