__author__ = "Jason C. Klima"

import attr
import collections
import copy
import logging
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
import uuid

from bokeh.palettes import Greens256
from io import BytesIO
from pyrosetta import Pose
from pyrosetta.rosetta.core.conformation import Residue

from viewer3d.colors import default_element_colors
from viewer3d.config import COLORBAR_ATTR
from viewer3d.converters import (
    _get_residue_chain_tuple,
    _int_to_str,
    _pdbstring_to_pose,
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
    Generator,
    GenericViewer,
    Iterable,
    List,
    NoReturn,
    Optional,
    OrderedDict,
    Tuple,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.per_residue_real_metric")


@attr.define(kw_only=True, slots=True, frozen=True)
class setPerResidueRealMetric(ModuleBase):
    """
    Show and color residues by a per-residue real metric in the `Pose.cache` (or `Pose.scores`)
    dictionary.

    Attributes:
        scoretype: A required `str` object matching a scoretype name. The scoretype name must
            contain a trailing underscore followed by either a residue number (for Pose
            numbering) or a PDB number (for PDB numbering), per the default output of a
            `PerResidueRealMetric` scoring result. For example, scoretype `"res_energy_18A"`
            would be selected by `scoretype="res_energy"`, and scoretype
            `"custom_type_atomic_clashes_4"` would be selected by `scoretype="atomic_clashes"`.

        vmin: A `float` or `int` object representing the minimum scoretype value for the
            palette. If `None`, this value is automatically set to the minimum scoretype value.

            Default: `None`

        vmax: A `float` or `int` object representing the maximum scoretype value for the
            palette. If `None`, this value is automatically set to the maximum scoretype value.

            Default: `None`

        palette: An iterable of `str` or `int` objects representing a color map.

            Default: `bokeh.palettes.Greens256`

        log: `None` to map colors spaced evenly on a linear scale between `vmin` to `vmax`. If
            an `int` or `float` object is provided, map colors spaced evenly on a log scale
            with the base provided.

            Default: `None`

        style: A `str` indicating a representation style of heavy atoms, choosing from either
            `"line"`, `"cross"`, `"stick"`, or `"sphere"`.

            Default: `"stick"`

        radius: A `float` or `int` indicating the radius of the heavy atoms represented by the
            `style` option.

            Default: `0.1`

        show_hydrogens: A `bool` object for the `nglview` and `pymol` backends to show all
            hydrogens.

            Defualt: `False`

        bonds: A `str` object for the `nglview` and `pymol` backends to show double bonds.
            Available options are: `"off"`, `"symmetric"`, or `"offset"`.

            Defualt: `"symmetric"`

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

        colorbar: A `bool` to show the colorbar axis.

            Default: `True`

        colorbar_extremes: An interable of two `bool` objects representing whether to plot
            colorbar extremes for `vmin` and `vmax` options.

            Default: `(False, False)`

        colorbar_label: A `str` object to label the colorbar axis.

            Default: `None` uses the `scoretype` attribute.

        colorbar_fontsize: A positive `int` representing the colorbar label and tickmark label
            fontsize.

            Default: `20`

        colorbar_nticks: An `int` representing the number of tickmarks to show on the colorbar
            axis, automatically interpolated between the `vmin` and `vmax` attributes.

            Default: `11`

        colorbar_discrete_ticks: A `bool` to discretize the tick marks on the colorbar.

            Default: `False`
    """

    scoretype: str = attr.field(
        validator=attr.validators.instance_of(str),
    )
    vmin: Union[float, int] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
    )
    vmax: Union[float, int] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
    )
    palette: Iterable[Union[str, int]] = attr.field(
        default=None,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of((str, int)),
            iterable_validator=attr.validators.instance_of(collections.abc.Iterable),
        ),
        converter=attr.converters.default_if_none(default=Greens256),
    )
    log: Optional[Union[float, int]] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of((float, int))),
        converter=_to_0_if_le_0,
    )
    style: str = attr.field(
        default="stick",
        validator=[
            attr.validators.instance_of(str),
            attr.validators.in_(("line", "cross", "stick", "sphere", "ball+stick")),
        ],
    )
    radius: Union[float, int] = attr.field(
        default=0.1,
        validator=attr.validators.instance_of((float, int)),
        converter=_to_0_if_le_0,
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
    colorbar: bool = attr.field(
        default=True,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    colorbar_extremes: Tuple[bool, bool] = attr.field(
        default=None,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(bool),
            iterable_validator=attr.validators.instance_of(collections.abc.Iterable),
        ),
        converter=attr.converters.default_if_none(default=(False, False)),
    )
    colorbar_label: Optional[str] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )
    colorbar_fontsize: int = attr.field(
        default=None,
        validator=attr.validators.instance_of(int),
        converter=attr.converters.default_if_none(default=20),
    )
    colorbar_nticks: int = attr.field(
        default=None,
        validator=attr.validators.instance_of(int),
        converter=attr.converters.default_if_none(default=11),
    )
    colorbar_discrete_ticks: bool = attr.field(
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )

    @colorbar_extremes.validator
    def _len_2(self, attribute, value):
        if len(value) != 2:
            raise ValueError(f"The '{attribute}' attribute must be an interable of length 2.")

    @log.validator
    def _not_bool(self, attribute, value):
        if isinstance(value, bool):
            raise ValueError(
                f"The '{attribute}' attribute must be an `int` or `float` object, or `None`."
            )

    def _iter_matching_scores(self, pose: Pose) -> Generator[Tuple[str, Any], None, None]:
        if hasattr(pose, "cache"):
            scores = dict(pose.cache.metrics.per_residue_real.all)
        else:
            scores = dict(pose.scores)
        resnums = list(map(str, range(1, pose.size() + 1)))
        pdb_info = pose.pdb_info()
        for key, value in scores.items():
            if self.scoretype in key:
                _r = key.split("_")[-1]
                if _r.isdigit():
                    if _r in resnums:
                        yield key, value
                elif pdb_info.pdb2pose(_r[-1], int(_r[:-1])):
                    yield key, value

    def _parse_residue_chain(self, pose: Pose, scoretype: str) -> Tuple[str, str]:
        _r = scoretype.split("_")[-1]
        if _r.isdigit():
            return _get_residue_chain_tuple(pose, int(_r))
        else:
            return (_r[:-1], _r[-1])

    def get_vmin_vmax(
        self, pose: Pose, _scores: Dict[str, float]
    ) -> Tuple[Union[float, int], Union[float, int]]:
        values = list(_scores.values())
        if not values:
            _scores_accessor = "Pose.cache" if hasattr(pose, "cache") else "Pose.scores"
            raise ValueError(
                f"Scoretype matching '{self.scoretype}' not found in `{_scores_accessor}` keys."
            )
        vmin = min(values) if self.vmin is None else self.vmin
        vmax = max(values) if self.vmax is None else self.vmax

        return vmin, vmax

    def get_elements_from_residue(self, residue: Residue) -> List[str]:
        residue_type = residue.type()
        elements = [
            residue_type.element(atom).name.upper() for atom in range(1, residue.natoms() + 1)
        ]
        return elements

    def get_nearest_value_from_keys(self, _palette_keys: numpy.array, _value: float) -> float:
        idx = numpy.abs(_palette_keys - _value).argmin()
        return _palette_keys[idx]

    def get_active_element_colors(self, pose: Pose) -> Dict[str, int]:
        _unique_elements = set()
        for i in range(1, pose.size() + 1):
            res = pose.residue(i)
            for j in range(1, res.natoms() + 1):
                elem = res.atom_type(j).element().strip().upper()
                if elem != "C":
                    _unique_elements.add(elem)
        _element_colors = {
            elem: color
            for elem, color in default_element_colors.items()
            if elem in _unique_elements
        }

        return _element_colors

    def get_palette_value_dict(
        self, vmin: Union[float, int], vmax: Union[float, int]
    ) -> Union[OrderedDict[float, str], NoReturn]:
        if self.log is not None:
            if any(v <= 0 for v in (vmin, vmax)):
                raise ValueError(
                    "The 'vmin' and 'vmax' attributes must be >0 for logarithmic color mapping."
                )
            _space = numpy.logspace(
                math.log(vmin, self.log),
                math.log(vmax, self.log),
                num=len(self.palette),
                base=self.log,
            )
        else:
            _space = numpy.linspace(
                vmin,
                vmax,
                num=len(self.palette),
            )
        _palette_value_dict = collections.OrderedDict(zip(_space, self.palette))

        return _palette_value_dict

    def get_colorbar(self, space: numpy.ndarray) -> bytes:
        fig, ax = plt.subplots(figsize=(20, 1))
        fig.subplots_adjust(bottom=0.5)
        cmap = (matplotlib.colors.ListedColormap(self.palette)).with_extremes(
            over=self.palette[-1], under=self.palette[0]
        )
        indexes = numpy.round(numpy.linspace(0, len(space) - 1, self.colorbar_nticks)).astype(int)
        ticks = list(space[indexes])
        if self.colorbar_discrete_ticks:
            ticks = [int(round(t)) for t in ticks]
        norm = matplotlib.colors.BoundaryNorm(ticks, cmap.N)
        if self.colorbar_label is None:
            label = self.scoretype
        else:
            label = self.colorbar_label
        if self.colorbar_extremes[0] and self.colorbar_extremes[-1]:
            extend = "both"
        elif self.colorbar_extremes[0] and not self.colorbar_extremes[-1]:
            extend = "min"
        elif not self.colorbar_extremes[0] and self.colorbar_extremes[-1]:
            extend = "max"
        else:
            extend = "neither"
        cbar = fig.colorbar(
            matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
            cax=ax,
            ticks=ticks,
            extend=extend,
            spacing="uniform",
            orientation="horizontal",
        )
        cbar.ax.tick_params(labelsize=self.colorbar_fontsize)
        cbar.set_label(label=label, size=self.colorbar_fontsize)
        buffer = BytesIO()
        plt.savefig(buffer, format="png", bbox_inches="tight", dpi=300)
        buffer.seek(0)
        colorbar = buffer.read()
        plt.close()

        return colorbar

    def viewer_setattr_colorbar(self, viewer: GenericViewer, space: numpy.array) -> GenericViewer:
        colorbar = self.get_colorbar(space)
        setattr(viewer, COLORBAR_ATTR, colorbar)

        return viewer

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _scores = dict(self._iter_matching_scores(pose))
        _cartoon_color = "spectrum" if self.cartoon_color is None else self.cartoon_color
        _style = self.style
        _vmin, _vmax = self.get_vmin_vmax(pose, _scores)
        _palette_value_dict = self.get_palette_value_dict(_vmin, _vmax)
        _palette_keys = numpy.array(list(_palette_value_dict.keys()))

        for scoretype, value in _scores.items():
            resi, chain = self._parse_residue_chain(pose, scoretype)
            _nearest_value = self.get_nearest_value_from_keys(_palette_keys, value)
            _C_color = _palette_value_dict[_nearest_value]
            if self.cartoon:
                viewer.setStyle(
                    {"model": model, "resi": resi, "chain": chain},
                    {
                        "cartoon": {"color": _cartoon_color, "opacity": self.cartoon_opacity},
                        _style: {
                            "color": _C_color,
                            "radius": self.radius,
                        },
                    },
                )
            else:
                viewer.setStyle(
                    {"model": model, "resi": resi, "chain": chain},
                    {
                        _style: {
                            "color": _C_color,
                            "radius": self.radius,
                        },
                    },
                )

        _element_colors = self.get_active_element_colors(pose)
        for _elem, _elem_color in _element_colors.items():
            if self.cartoon:
                viewer.setStyle(
                    {
                        "model": model,
                        "elem": _elem,
                    },
                    {
                        "cartoon": {"color": _cartoon_color, "opacity": self.cartoon_opacity},
                        _style: {
                            _style: {"color": _elem_color},
                            "radius": self.radius,
                        },
                    },
                )
            else:
                viewer.setStyle(
                    {
                        "model": model,
                        "elem": _elem,
                    },
                    {
                        _style: {
                            _style: {"color": _elem_color},
                            "radius": self.radius,
                        }
                    },
                )

        if self.colorbar:
            viewer = self.viewer_setattr_colorbar(viewer, _palette_keys)

        return viewer

    @requires_init
    def apply_nglview(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _scores = dict(self._iter_matching_scores(pose))
        _cartoon_color = "atomindex" if self.cartoon_color is None else self.cartoon_color
        if _cartoon_color == "black":
            _cartoon_color = "#000000"
        _style = _py3Dmol_to_nglview_style(self.style)
        if self.cartoon:
            viewer.add_representation(
                repr_type="cartoon",
                selection="*",
                color=_cartoon_color,
                radius=self.cartoon_radius,
                opacity=self.cartoon_opacity,
                component=model,
            )
        _vmin, _vmax = self.get_vmin_vmax(pose, _scores)
        _palette_value_dict = self.get_palette_value_dict(_vmin, _vmax)
        _palette_keys = numpy.array(list(_palette_value_dict.keys()))
        _default_element_colors = copy.deepcopy(default_element_colors)

        _rules = []
        for scoretype, value in _scores.items():
            _residue, _chain = self._parse_residue_chain(pose, scoretype)
            _resnum = pose.pdb_info().pdb2pose(_chain, int(_residue))
            _elements_from_residue = self.get_elements_from_residue(pose.residue(_resnum))
            _nearest_value = self.get_nearest_value_from_keys(_palette_keys, value)
            _C_color = _palette_value_dict[_nearest_value]
            _element_colors = _default_element_colors.copy()
            _element_colors["C"] = _C_color
            for _element, _element_color in _element_colors.items():
                if _element in _elements_from_residue:
                    _rules.append((int(_residue), _chain, _element, _element_color))
        _selection_name = f"{self.scoretype}_{uuid.uuid4().hex}"
        self.add_scheme_func_from_rules(_selection_name, _rules)
        _default_selection = "*" if self.show_hydrogens else "not hydrogen"
        viewer.add_representation(
            repr_type=_style,
            selection=_default_selection,
            color=_selection_name,
            radius=self.radius,
            multipleBond=self.bonds,
            component=model,
        )

        if self.colorbar:
            viewer = self.viewer_setattr_colorbar(viewer, _palette_keys)

        return viewer

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:

        if pose is None:
            pose = _pdbstring_to_pose(pdbstring, self.__class__.__name__)

        _scores = dict(self._iter_matching_scores(pose))
        _style = _py3Dmol_to_pymol_style(self.style)
        default_selection = f"obj {model}" if self.show_hydrogens else f"obj {model} and not elem h"

        if self.cartoon:
            viewer.show("cartoon", default_selection)
            if self.cartoon_color is None:
                viewer.spectrum("index", "rainbow", default_selection)
            else:
                viewer = self.apply_pymol_cartoon_color(
                    viewer, self.cartoon_color, default_selection
                )
            viewer.set("cartoon_transparency", 1 - self.cartoon_opacity)
            with self.out:
                viewer.do(f"set cartoon_loop_radius, {self.cartoon_radius}")
                viewer.do(f"set cartoon_rect_width, {self.cartoon_radius}")
                viewer.do(f"set cartoon_oval_width, {self.cartoon_radius}")
        if self.radius > 1e-10:
            viewer.show(_style, default_selection)
            if _style == "sticks":
                viewer.set("stick_radius", self.radius, default_selection)
            elif _style == "spheres":
                viewer.set("sphere_scale", self.radius, default_selection)
            elif _style == "dots":
                viewer.set("dot_width", self.radius, default_selection)
            elif _style == "lines":
                viewer.set("line_width", self.radius, default_selection)

        _vmin, _vmax = self.get_vmin_vmax(pose, _scores)
        _palette_value_dict = self.get_palette_value_dict(_vmin, _vmax)
        _palette_keys = numpy.array(list(_palette_value_dict.keys()))
        for scoretype, value in _scores.items():
            resi, chain = self._parse_residue_chain(pose, scoretype)
            selection = f"(obj {model} and chain {chain} and resi {resi})"
            if not self.show_hydrogens:
                selection = f"{selection} and (not elem h)"
            _nearest_value = self.get_nearest_value_from_keys(_palette_keys, value)
            _color = _palette_value_dict[_nearest_value]
            viewer = self.apply_pymol_color(viewer, _color, selection)
        with self.out:
            viewer.do("color atomic, not elem C")

        if self.bonds == "off":
            viewer.set("valence", 0)
        else:
            viewer.set("valence", 1)
            if self.bonds == "symmetric":
                viewer.set("valence_mode", 1)
            elif self.bonds == "offset":
                viewer.set("valence_mode", 0)

        if self.colorbar:
            viewer = self.viewer_setattr_colorbar(viewer, _palette_keys)

        return viewer
