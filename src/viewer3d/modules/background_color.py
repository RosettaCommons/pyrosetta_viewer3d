__author__ = "Jason C. Klima"

import attr
import logging

from pyrosetta import Pose

from viewer3d.converters import (
    _hex_to_rgb,
    _int32_to_str,
    _to_hex,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.type_defs import (
    GenericViewer,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.background_color")


@attr.define(kw_only=False, slots=True, frozen=True)
class setBackgroundColor(ModuleBase):
    """
    Set the background color with either hexcode or standard colors.

    Attributes:
        color: Hexcode literal (e.g., `0xffffffff`) or `str` indicating a standard color (e.g.,
            `"black"`).

            Default: `0xffffffff`
    """

    color: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[
            attr.converters.default_if_none(default=0xFFFFFFFF),
            _to_hex,
        ],
    )

    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        viewer.setBackgroundColor(self.color)

        return viewer

    def apply_nglview(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        viewer.background = f"#{self.color:06x}" if isinstance(self.color, int) else self.color

        return viewer

    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        name = _int32_to_str(self.color)
        if viewer.get_color_index(name) == -1:
            rgb = _hex_to_rgb(name)
            with self.out:
                viewer.do(f"set_color {name}, {rgb}")
        viewer.bg_color(name)

        return viewer
