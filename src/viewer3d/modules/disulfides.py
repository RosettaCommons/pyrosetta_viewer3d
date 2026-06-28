__author__ = "Jason C. Klima"

import attr
import itertools
import logging

from pyrosetta import Pose
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.conformation import is_disulfide_bond

from viewer3d.converters import (
    _get_residue_chain_tuple,
    _pdbstring_to_pose,
    _to_hex,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    GenericViewer,
    List,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.disulfides")


@attr.define(kw_only=False, slots=True, frozen=True)
class setDisulfides(ModuleBase):
    """
    Display disulfide bonds according to `pyrosetta.rosetta.core.conformation.is_disulfide_bond`
    for all combinations of cysteine residues in each decoy.

    Attributes:
        color: A `str` indicating a standard color (e.g., `"black"`).

            Default: `"gold"`

        radius: A `float` or `int` indicating the radius of the stick connecting the atoms
            participating in each disulfide bond.

            Default: `0.5`
    """

    color: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="gold"), _to_hex],
    )
    radius: Union[float, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((float, int)),
        converter=attr.converters.default_if_none(default=0.1),
    )
    sulfur_atom_name: str = attr.field(
        default="SG",
        validator=attr.validators.instance_of(str),
        init=False,
    )

    def get_cys_res(self, pose: Pose) -> List[int]:
        return [i for i, aa in enumerate(pose.sequence(), start=1) if aa == "C"]

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        cys_res = self.get_cys_res(pose)
        for i, j in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_xyz = pose.xyz(AtomID(pose.residue(i).atom_index(self.sulfur_atom_name), i))
                j_xyz = pose.xyz(AtomID(pose.residue(j).atom_index(self.sulfur_atom_name), j))
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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _color = "#000000" if self.color == "black" else self.color
        cys_res = self.get_cys_res(pose)
        selection_disulfides = []
        for i, j in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_res, i_chain = _get_residue_chain_tuple(pose, i)
                j_res, j_chain = _get_residue_chain_tuple(pose, j)
                i_sele = f"{i_res}:{i_chain}.{self.sulfur_atom_name}"
                j_sele = f"{j_res}:{j_chain}.{self.sulfur_atom_name}"
                selection_disulfides.append([i_sele, j_sele])
        if selection_disulfides:
            viewer.add_distance(
                atom_pair=selection_disulfides,
                color=_color,
                radius=self.radius,
                label_visible=False,
            )
            selection = " or ".join([f"({s[0]} or {s[1]})" for s in selection_disulfides])
            viewer.add_representation(
                repr_type="ball+stick",
                selection=selection,
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

        cys_res = self.get_cys_res(pose)
        for i, j in itertools.product(cys_res, repeat=2):
            if is_disulfide_bond(pose.conformation(), i, j):
                i_res, i_chain = _get_residue_chain_tuple(pose, i)
                j_res, j_chain = _get_residue_chain_tuple(pose, j)
                i_sele = f"(obj {model} and chain {i_chain} and resi {i_res} and name {self.sulfur_atom_name})"
                j_sele = f"(obj {model} and chain {j_chain} and resi {j_res} and name {self.sulfur_atom_name})"
                viewer.bond(i_sele, j_sele)
                viewer.set_bond("stick_color", self.color, i_sele, j_sele)
                viewer.set_bond("stick_radius", self.radius, i_sele, j_sele)

        return viewer
