__author__ = "Jason C. Klima"

import attr
import logging
import uuid

from pyrosetta import Pose
from pyrosetta.rosetta.core.select.residue_selector import ResidueSelector

from viewer3d.colors import default_element_colors
from viewer3d.converters import (
    _get_nglview_selection,
    _get_pymol_selection,
    _int_to_str,
    _pdbstring_to_pose,
    _pose_to_residue_chain_tuples,
    _py3Dmol_to_nglview_style,
    _py3Dmol_to_pymol_style,
    _to_hex,
    _to_int_color,
    _to_0_if_le_0,
    _to_1_if_gt_1,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Any,
    Dict,
    GenericViewer,
    Optional,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.style")


@attr.define(kw_only=True, slots=True, frozen=True)
class setStyle(ModuleBase):
    """
    Show and color cartoon, and/or show heavy atoms with provided style, color and radius for
    each decoy. If the `residue_selector` argument is provided, apply styles only to the
    selected residues. If the `command` argument is provided, override all other arguments and
    pass arbitrary `py3Dmol.view.setStyle` commands (for the `py3Dmol` backend) or PyMOL
    commands (for the `pymol` backend) to the visualization.

    Attributes:
        residue_selector: An instance of `ResidueSelector` on which to apply the style(s).

            Default: `None`

        cartoon: A `bool` to show cartoon representation.

            Default: `True`

        cartoon_color: A hexcode literal (e.g., `0xAF10AB`) or `str` indicating a standard
            color (e.g., `"grey"`) for the cartoon representation. If `"spectrum"`, apply
            reversed color gradient based on residue numbers. The option `cartoon` must also be
            set to `True`. Reference: [3dmol ColorSpec](https://3dmol.org/doc/global.html#ColorSpec)

            Default: `"spectrum"`

        cartoon_radius: Set the cartoon radius for the `nglview` and `pymol` backends.

            Default: `0.1`

        cartoon_opacity: Set the cartoon opacity for the `nglview` and `pymol` backends.

            Default: `1`

        style: A `str` indicating a representation style of heavy atoms, choosing from either
            `"line"`, `"cross"`, `"stick"`, or `"sphere"`.

            Default: `"stick"`

        colorscheme: A hexcode literal (e.g., `0xAF10AB`) or a `str` object. If a `str` object
            is provided for the `py3Dmol` backend, the `str` indicates the color scheme for
            heavy atoms represented by the `style` option. Options include:

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

            If a `str` object is provided for the `nglview` backend, the `str` indicates the
            color scheme for heavy atoms represented by the `style` option. Options include:

                - A lowercased standard color optionally
                    followed by "Carbon" (e.g.,
                    `"orangeCarbon"`). Standard colors
                    that are supported by the `py3Dmol`
                    backend may or may not be supported
                    by the `nglview` backend.
                - "atomindex"
                - "bfactor"
                - "chainid"
                - "chainindex"
                - "chainname"
                - "densityfit"
                - "electrostatic"
                - "element"
                - "entityindex"
                - "entitytype"
                - "geoquality"
                - "hydrophobicity"
                - "modelindex"
                - "moleculetype"
                - "occupancy"
                - "random"
                - "residueindex"
                - "resname"
                - "sstruc"
                - "uniform"
                - "value"
                - "volume"

            Reference: [NGLView Coloring](https://nglviewer.org/ngl/api/manual/coloring.html)

            If a `str` object is provided for the `pymol` backend, the `str` indicates the
            color scheme for heavy atoms represented by the `style` option. Options include:

                - A lowercased standard color optionally
                    followed by "Carbon" (e.g., "red" or
                    `"orangeCarbon"`). Standard colors
                    that are supported by the `py3Dmol`
                    backend may or may not be supported
                    by the `pymol` backend.
                - A standard PyMOL color (e.g., `"grey"`)

            Default: `"blackCarbon"`

        radius: A `float` or `int` indicating the radius of the heavy atoms represented by the
            `style` option.

            Default: `0.1`

        label: A `bool` to show labels next to residues selected by the `residue_selector`
            option.

            Default: `True`

        label_fontsize: A `int` or `float` indicating the font size of labels next to residues
            selected by the `residue_selector` option, only if `label` is `True`.

            Default: `12`

        label_background: A `bool` to show the background of labels next to residues selected
            by the `residue_selector` option, only if `label` is `True`.

            Default: `False`

        label_fontcolor: A `str` indicating a standard font color (e.g., `"grey"`) for label
            text next to residues selected by the `residue_selector` option, only if `label`
            is `True`.

            Default: `"black"`

        command: For the `py3Dmol` backend, a `dict` or `tuple` of `dict` objects of arbitrary
            `py3Dmol.view.setStyle` commands. For the `pymol` backend, a `tuple` of `str`
            objects of arbitrary PyMOL commands. If specified with the `py3Dmol` or `pymol`
            backend, then this option overrides all other options. This option is not supported
            and ignored by the `nglview` backend.

            Example for `py3Dmol` backend:

                >>> cmds = (
                ...   {"hetflag": True},
                ...   {
                ...     "stick": {
                ...       "singleBond": True,
                ...       "colorscheme": "greyCarbon",
                ...       "radius": 0.2,
                ...     }
                ...   }
                ... )
                >>> v = viewer3d.init(poses, backend=0)
                >>> v += viewer3d.setStyle(command=cmds)
                >>> v.show()

            Example for `pymol` backend:

                >>> cmds = (
                ...   "set cartoon_oval_width, 1.7",
                ...   "set cartoon_rect_width, 1.7",
                ...   "set cartoon_loop_radius, 1.3",
                ... )
                >>> v = viewer3d.init(poses, backend=2)
                >>> v += viewer3d.setStyle(command=cmds)
                >>> v.show()

            Default: `None`

        show_hydrogens: A `bool` object for the `nglview` and `pymol` backends to show all
            hydrogens.

            Defualt: `False`

        bonds: A `str` object for the `nglview` and `pymol` backends to show double bonds.
            Options are: `"off"`, `"symmetric"`, `"offset"`.

            Defualt: `"symmetric"`
    """

    residue_selector: Optional[ResidueSelector] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(ResidueSelector)),
    )
    cartoon: bool = attr.field(
        default=True,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    cartoon_color: Optional[Union[str, int]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((str, int))),
        converter=_to_hex,
    )
    cartoon_radius: Optional[Union[float, int]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=attr.converters.default_if_none(default=0.1),
    )
    cartoon_opacity: Union[float, int] = attr.field(
        default=1,
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=[_to_0_if_le_0, _to_1_if_gt_1],
    )
    style: str = attr.field(
        default="stick",
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("line", "cross", "stick", "sphere", "ball+stick")),
        ],
    )
    colorscheme: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="blackCarbon"), _to_hex],
    )
    radius: Union[float, int] = attr.field(
        default=0.1,
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    label: bool = attr.field(
        default=True,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    label_fontsize: Union[float, int] = attr.field(
        default=12,
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
    )
    label_background: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    label_fontcolor: Union[str, int] = attr.field(
        default=None,
        validator=attr.validators.instance_of((str, int)),
        converter=[attr.converters.default_if_none(default="black"), _to_hex],
    )
    command: Optional[Union[Tuple[Dict[Any, Any]], Dict[Any, Any], Tuple[str]]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((tuple, dict))),
    )
    show_hydrogens: bool = attr.field(
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    bonds: bool = attr.field(
        default=None,
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("off", "symmetric", "offset")),
        ],
        converter=attr.converters.default_if_none(default="symmetric"),
    )

    def get_colorscheme_name(self):
        if isinstance(self.colorscheme, str) and self.colorscheme.endswith("Carbon"):
            _element_colors = default_element_colors.copy()
            _colorscheme_name = uuid.uuid4().hex
            _C_color = self.colorscheme[: -len("Carbon")]
            if _C_color == "black":
                _C_color = "#000000"
            _element_colors["C"] = _C_color
            _selection_scheme = [
                [_int_to_str(_element_color), f"_{_element}"]
                for _element, _element_color in _element_colors.items()
            ]
            self.add_selection_scheme(_colorscheme_name, _selection_scheme)
        elif isinstance(self.colorscheme, int):
            _element_colors = default_element_colors.copy()
            _colorscheme_name = uuid.uuid4().hex
            _element_colors["C"] = self.colorscheme
            _func_str = self.get_element_colors_func_str(_element_colors)
            self.add_scheme_func(_colorscheme_name, _func_str)
        else:
            _colorscheme_name = self.colorscheme

        return _colorscheme_name

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        _style = self.style
        _colorscheme = "colorscheme" if isinstance(self.colorscheme, str) else "color"
        if self.show_hydrogens:
            _logger.warning(
                "The 'show_hydrogens' attribute is not supported for the `py3Dmol` backend. "
                "Please add `viewer3d.setHydrogens()` instead."
            )
        _cartoon_color = "spectrum" if self.cartoon_color is None else self.cartoon_color

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
                                "cartoon": {
                                    "color": _cartoon_color,
                                    "opacity": self.cartoon_opacity,
                                },
                                _style: {
                                    _colorscheme: self.colorscheme,
                                    "radius": self.radius,
                                },
                            },
                        )
                    else:
                        viewer.setStyle(
                            {"model": model, "resi": resi, "chain": chain},
                            {
                                _style: {
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
                            "cartoon": {"color": _cartoon_color, "opacity": self.cartoon_opacity},
                            _style: {
                                _colorscheme: self.colorscheme,
                                "radius": self.radius,
                            },
                        },
                    )
                else:
                    viewer.setStyle(
                        {"model": model},
                        {
                            _style: {
                                _colorscheme: self.colorscheme,
                                "radius": self.radius,
                            }
                        },
                    )

        return viewer

    @requires_init
    def apply_nglview(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        # Set defaults
        _cartoon_color = "atomindex" if self.cartoon_color is None else self.cartoon_color
        if _cartoon_color == "black":
            _cartoon_color = "#000000"
        _style = _py3Dmol_to_nglview_style(self.style)
        default_selection = "*" if self.show_hydrogens else "not hydrogen"
        if self.command:
            _logger.warning("The 'command' attribute is not supported for the `nglview` backend.")
        else:
            if self.residue_selector:
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
                            repr_type=_style,
                            selection=selection,
                            color=self.get_colorscheme_name(),
                            radius=self.radius,
                            multipleBond=self.bonds,
                            component=model,
                        )
                    if self.cartoon:
                        viewer.add_representation(
                            repr_type="cartoon",
                            selection=selection,
                            color=_cartoon_color,
                            radius=self.cartoon_radius,
                            opacity=self.cartoon_opacity,
                            component=model,
                        )
                    if self.label:
                        nbr_res_selection = _get_nglview_selection(
                            pose,
                            self.residue_selector,
                            show_hydrogens=self.show_hydrogens,
                            nbr_atom_only=True,
                            logger=_logger,
                        )
                        viewer.add_representation(
                            repr_type="label",
                            labelType="res",
                            selection=nbr_res_selection,
                            showBorder=self.label_background,
                            borderColor="gray",
                            showBackground=self.label_background,
                            backgroundColor=self.colorscheme,  # Default to white if unrecognized
                            colorValue=self.label_fontcolor,
                            component=model,
                        )
            else:
                if self.radius > 1e-10:
                    viewer.add_representation(
                        repr_type=_style,
                        selection=default_selection,
                        color=self.get_colorscheme_name(),
                        radius=self.radius,
                        multipleBond=self.bonds,
                        component=model,
                    )
                if self.cartoon:
                    viewer.add_representation(
                        repr_type="cartoon",
                        selection=default_selection,
                        color=_cartoon_color,
                        radius=self.cartoon_radius,
                        opacity=self.cartoon_opacity,
                        component=model,
                    )

        return viewer

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        # Set defaults
        if isinstance(self.colorscheme, str) and self.colorscheme.endswith("Carbon"):
            _colorscheme = self.colorscheme[: -len("Carbon")]
        else:
            _colorscheme = self.colorscheme
        _style = _py3Dmol_to_pymol_style(self.style)
        default_selection = f"obj {model}" if self.show_hydrogens else f"obj {model} and not elem h"

        if self.command:
            _msg = (
                "The 'command' attribute must be a `tuple` of `str` objects to be used "
                "by the `pymol` backend. Skipping object of type: {0}"
            )
            if isinstance(self.command, tuple):
                for obj in self.command:
                    if isinstance(obj, str):
                        with self.out:
                            viewer.do(obj)
                    else:
                        _logger.warning(_msg.format(type(obj)))
            else:
                _logger.warning(_msg.format(type(obj)))
        else:
            if self.residue_selector:
                if pose is None:
                    pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

                selection = _get_pymol_selection(
                    pose,
                    self.residue_selector,
                    show_hydrogens=self.show_hydrogens,
                    logger=_logger,
                )
                if not selection:
                    pass
                else:
                    selection = f"obj {model} and ({selection})"
                    if self.cartoon:
                        if self.cartoon_color is None:
                            viewer.spectrum("index", "rainbow", selection)
                        else:
                            viewer = self.apply_pymol_cartoon_color(
                                viewer, self.cartoon_color, selection
                            )
                        with self.out:
                            viewer.color("atomic", "not elem C")
                        viewer.set("cartoon_transparency", 1 - self.cartoon_opacity)
                        with self.out:
                            viewer.set("cartoon_loop_radius", self.cartoon_radius)
                            viewer.set("cartoon_rect_width", self.cartoon_radius)
                            viewer.set("cartoon_oval_width", self.cartoon_radius)
                    if self.radius > 1e-10:
                        viewer.show(_style, selection)
                        if _style == "sticks":
                            viewer.set("stick_radius", self.radius, selection)
                        elif _style == "spheres":
                            viewer.set("sphere_scale", self.radius, selection)
                        elif _style == "dots":
                            viewer.set("dot_width", self.radius, selection)
                        elif _style == "lines":
                            viewer.set("line_width", self.radius, selection)
                        viewer = self.apply_pymol_color(viewer, _colorscheme, selection)
                        with self.out:
                            viewer.color("atomic", "not elem C")
                    if self.label:
                        resi, chain = _pose_to_residue_chain_tuples(
                            pose, self.residue_selector, logger=_logger
                        )
                        for _resi, _chain in zip(resi, chain):
                            _resnum = pose.pdb_info().pdb2pose(_chain, int(_resi))
                            _residue = pose.residue(_resnum)
                            _nbr_atom_index = _residue.type().nbr_atom()
                            _nbr_atom_name = _residue.atom_name(_nbr_atom_index).strip()
                            _label_selection = f"(obj {model} and chain {_chain} and resi {_resi} and name {_nbr_atom_name})"
                            _label = f"{_resi}-{_chain}"
                            with self.out:
                                viewer.do(f"label {_label_selection}, '{_label}'")
                            viewer.set("label_color", self.label_fontcolor, _label_selection)
                            viewer.set("label_size", self.label_fontsize)
                            if self.label_background:
                                viewer.set("label_bg_color", "white")
                                viewer.set("label_bg_outline", 1)
            else:
                if self.cartoon:
                    if self.cartoon_color is None:
                        viewer.spectrum("index", "rainbow", default_selection)
                    else:
                        viewer = self.apply_pymol_cartoon_color(
                            viewer, self.cartoon_color, default_selection
                        )
                    viewer.color("atomic", "not elem C")
                    viewer.set("cartoon_transparency", 1 - self.cartoon_opacity)
                    with self.out:
                        viewer.set("cartoon_loop_radius", self.cartoon_radius)
                        viewer.set("cartoon_rect_width", self.cartoon_radius)
                        viewer.set("cartoon_oval_width", self.cartoon_radius)
                if self.radius > 1e-10:
                    viewer.show("sticks", default_selection)
                    viewer.set("stick_radius", self.radius)
                    viewer = self.apply_pymol_color(viewer, _colorscheme, default_selection)
                    viewer.color("atomic", "not elem C")

            if self.bonds == "off":
                viewer.set("valence", 0)
            else:
                viewer.set("valence", 1)
                if self.bonds == "symmetric":
                    viewer.set("valence_mode", 1)
                elif self.bonds == "offset":
                    viewer.set("valence_mode", 0)

        return viewer
