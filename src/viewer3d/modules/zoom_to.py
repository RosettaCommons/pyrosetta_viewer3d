__author__ = "Jason C. Klima"

import attr
import logging

from pyrosetta import Pose
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueSelector,
    TrueResidueSelector,
)

from viewer3d.converters import (
    _get_nglview_selection,
    _get_pymol_selection,
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import GenericViewer

_logger: logging.Logger = logging.getLogger("viewer3d.modules.zoom_to")


@attr.define(kw_only=False, slots=True, frozen=True)
class setZoomTo(ModuleBase):
    """
    Zoom to a provided `ResidueSelector` in each decoy.

    Attributes:
        residue_selector: An instance of `ResidueSelector` into which to zoom.

            Default: `TrueResidueSelector()`
    """

    residue_selector: ResidueSelector = attr.field(
        default=None,
        validator=attr.validators.instance_of(ResidueSelector),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
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

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_pymol_selection(
            pose,
            self.residue_selector,
            logger=_logger,
        )
        if not selection:
            viewer.orient(f"obj {model}")
        else:
            viewer.orient(f"obj {model} and ({selection})")

        return viewer
