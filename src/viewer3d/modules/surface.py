__author__ = "Jason C. Klima"

import attr
import logging
import sys

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
    _to_hex,
    _to_0_if_le_0,
    _to_1_if_gt_1,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Dict,
    GenericViewer,
    Optional,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.surface")


@attr.define(kw_only=False, slots=True, frozen=True)
class setSurface(ModuleBase):
    """
    Show the specified surface for each decoy.

    Attributes:
        residue_selector: An instance of `ResidueSelector` to select residues on which to apply
            the surface.

            Default: `None`

        surface_type: A `str` indicating surface type to be displayed. The `py3Dmol` and
            `nglview` backends support the following options:

                - "VDW": Van der Waals surface
                - "MS": Molecular surface
                - "SES": Solvent excluded surface
                - "SAS": Solvent accessible surface
                - "AV": High quality molecular surface
                    (only supported by `nglview` backend)

            The `pymol` backend currently supports the `"SAS"` option to show solvent-accessible
            surface, otherwise by default it shows the solvent-excluded surface. Additional
            surface properties may be set by adding the `viewer3d.setStyle(command=...)`
            visualization module or further configured within the PyMOL GUI.

            Default: `"VDW"`

        opacity: A `float` or `int` between `0` and `1` for opacity of the displayed surface.
            Not currently supported for `nglview` backend.

            Default: `0.5`

        color: A `str` indicating a standard color (e.g., `"grey"`) of the surface to be
            displayed, or a Hexcode literal (e.g., `0xffffffff`). Either `color` or `colorscheme`
            arguments may be provided, where `colorscheme` overrides `color`.

            Default: `None`

        colorscheme: A `str` indicating the color scheme of the surface to be displayed. Either
            `color` or `colorscheme` arguments may be provided, where `colorscheme` overrides
            `color`. Not currently supported for `nglview` or `pymol` backends.
            For the `py3Dmol` backend, options include:

                - A lowercased standard color optionally
                    followed by "Carbon" (e.g.,
                    `"orangeCarbon"`).
                - "ssPyMOL": PyMol secondary colorscheme
                - "ssJmol": Jmol secondary colorscheme
                - "Jmol": Jmol primary colorscheme
                - "default": default colorscheme
                - "amino": amino acid colorscheme
                - "shapely": shapely protien colorscheme
                - "nucleic": nucleic acid colorscheme
                - "chain": standard chain colorscheme
                - "chainHetatm": chain Hetatm colorscheme

            Reference: [3dmol builtinColorSchemes](https://3dmol.org/doc/global.html#builtinColorSchemes)

            Default: `None`
    """

    residue_selector: ResidueSelector = attr.field(
        default=None,
        validator=attr.validators.instance_of(ResidueSelector),
        converter=attr.converters.default_if_none(default=TrueResidueSelector()),
    )
    surface_type: str = attr.field(
        default="VDW",
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("VDW", "MS", "SAS", "SES", "AV")),
        ],
    )
    opacity: Union[float, int] = attr.field(
        default=0.5,
        validator=attr.validators.instance_of((float, int)),
        converter=[_to_0_if_le_0, _to_1_if_gt_1],
    )
    color: Optional[Union[str, int]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
        converter=_to_hex,
    )
    colorscheme: Optional[Union[str, int]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
        converter=_to_hex,
    )

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
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
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        surface_types_dict: Dict[str, str] = {
            "VDW": "vws",
            "MS": "ms",
            "SAS": "sas",
            "SES": "ses",
            "AV": "av",
        }

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _color = "#000000" if self.color == "black" else self.color
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
                color=_color,
                opacity=self.opacity,
                component=model,
            )

        return viewer

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        selection = _get_pymol_selection(
            pose, self.residue_selector, show_hydrogens=True, logger=_logger
        )
        if not selection:
            pass
        else:
            _selection = f"obj {model} and {selection}"
            with self.out:
                viewer.do("set surface_mode, 3")
            viewer.show("surface", _selection)
            viewer.set("transparency", 1 - self.opacity, _selection)
            if isinstance(self.color, str):
                with self.out:
                    viewer.do(f"set surface_color, {self.color}, {_selection}")
            elif isinstance(self.color, int):
                name = _int32_to_str(self.color)
                if viewer.get_color_index(name) == -1:
                    rgb = _hex_to_rgb(name)
                    with self.out:
                        viewer.do(f"set_color {name}, {rgb}")
                with self.out:
                    viewer.do(f"set surface_color, {name}, {_selection}")
            if self.surface_type == "SAS":
                viewer.set("surface_solvent", True)
            else:
                viewer.set("surface_solvent", False)
                _logger.warning(
                    "Disabling the 'surface_solvent' setting in PyMOL since 'SAS' "
                    "was not passed to the 'surface_type' argument value."
                )

        return viewer
