__author__ = "Jason C. Klima"

import attr
import logging
import pyrosetta.distributed.io as io

from pyrosetta import Pose

from viewer3d.modules.base import ModuleBase
from viewer3d.type_defs import (
    GenericViewer,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.zoom")


@attr.define(kw_only=False, slots=True, frozen=True)
class setZoom(ModuleBase):
    """
    Set the zoom magnification factor of each decoy. For the `py3Dmol` backend, values `>1`
    zoom in, and values `<1` zoom out. For the `nglview` backend, values `>0` zoom in, and
    values `<0` zoom out. For the `pymol` backend, values `<0` zoom in, and values `>0` zoom
    out.

    Attributes:
        factor: A `float` or `int` indicating the zoom magnification factor.

            Default: `2`
    """

    factor: Union[int, float] = attr.field(
        default=2,
        validator=attr.validators.instance_of((float, int)),
    )

    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        viewer.zoom(self.factor)
        return viewer

    def apply_nglview(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        viewer.control.zoom(self.factor)
        return viewer

    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        viewer.zoom(model, self.factor)
        return viewer
