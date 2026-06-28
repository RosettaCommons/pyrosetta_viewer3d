__author__ = "Jason C. Klima"

import attr
import collections
import copy
import json
import logging
import sys

from functools import singledispatch
from ipywidgets.widgets import Output
from pyrosetta import Pose

from viewer3d.colors import default_element_colors
from viewer3d.config import BACKENDS
from viewer3d.converters import (
    _hex_to_rgb,
    _int_to_str,
    _int32_to_str,
    _to_int_color,
)
from viewer3d.exceptions import ModuleNotImplementedError
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Any,
    Dict,
    GenericViewer,
    List,
    Union,
)

_logger: logging.Logger = logging.getLogger("viewer3d.modules.base")


@attr.define(kw_only=False, slots=True, frozen=True)
class ModuleBase:

    out = attr.field(default=Output(), init=False)

    @staticmethod
    def _to_modules(objs: Any) -> List["ModuleBase"]:
        @singledispatch
        def _to_module(obj):
            raise ValueError(
                "The 'modules' viewer attribute must be an instance of `ModuleBase` "
                + f"or an iterable of `ModuleBase` objects. Received: {type(obj)}"
            )

        _to_module.register(ModuleBase, lambda obj: obj)

        if objs is None:
            _modules = []
        elif isinstance(objs, collections.abc.Iterable):
            _modules = list(map(_to_module, objs))
        else:
            _modules = [_to_module(objs)]

        return _modules

    def get_element_colors_func_str(self, element_colors: Dict[str, int]) -> str:
        func_lines = ["this.atomColor = function (atom) {"]
        for e, c in element_colors.items():
            func_lines.append(f"if (atom.element === '{e}') return 0x{c:06X};")
        func_lines.append("return 0xFFFFFF;")
        func_lines.append("};")

        return "\n".join(func_lines)

    def add_scheme_func(self, name: str, func_str: str) -> None:
        sys.modules["nglview"].color.ColormakerRegistry.add_scheme_func(name, func_str)

    def add_scheme_func_from_rules(self, selection_name: str, rules: List[Any]) -> None:
        color_map = {
            f"{chain}|{resi}|{element}": _to_int_color(color)
            for resi, chain, element, color in rules
        }
        js_map = json.dumps(color_map)
        func_str = f"""
        var colorMap = {js_map};
        this.atomColor = function (atom) {{
            var key = atom.chainname + "|" + atom.resno + "|" + atom.element;
            var c = colorMap[key];
            if (c !== undefined) return c;
            return 0xFFFFFF;
        }};
        """
        sys.modules["nglview"].color.ColormakerRegistry.add_scheme_func(selection_name, func_str)

    def add_selection_scheme(self, name: str, selection_scheme: List[List[str]]) -> None:

        cm = sys.modules["nglview"].color.ColormakerRegistry
        cm.add_selection_scheme(name, selection_scheme)

    def add_element_selection_scheme(self, name: str) -> None:
        """
        Add element-based selection scheme to NGLView apply_pymol_color
        based on '<color>Carbon' naming system or an arbitrary `str` object.
        """

        _default_element_colors = copy.deepcopy(default_element_colors)

        if name.endswith("Carbon"):
            _default_element_colors["C"] = name[: -len("Carbon")]
        else:
            _default_element_colors["C"] = name
        _selection_scheme = [
            [_int_to_str(color), f"_{element}"]
            for (element, color) in _default_element_colors.items()
        ]
        self.add_selection_scheme(name, _selection_scheme)

    def _maybe_set_pymol_color(self, viewer: GenericViewer, color: Union[int, str]) -> str:
        name = _int32_to_str(color)
        if viewer.get_color_index(name) == -1:
            rgb = _hex_to_rgb(name)
            with self.out:
                viewer.do(f"set_color {name}, {rgb}")

        return name

    def apply_pymol_color(
        self, viewer: GenericViewer, color: Union[int, str], selection: str
    ) -> GenericViewer:
        name = self._maybe_set_pymol_color(viewer, color)
        with self.out:
            viewer.do(f"color {name}, {selection}")

        return viewer

    def apply_pymol_cartoon_color(
        self, viewer: GenericViewer, color: Union[int, str], selection: str
    ) -> GenericViewer:
        name = self._maybe_set_pymol_color(viewer, color)
        viewer.set("cartoon_color", name, selection)

        return viewer


@attr.define(kw_only=True, slots=True, frozen=True)
class setTemplate(ModuleBase):
    """Template class for developing new visualization modules."""

    @requires_init
    def apply_py3Dmol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[0])

    @requires_init
    def apply_nglview(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[1])

    @requires_init
    def apply_pymol(
        self, viewer: GenericViewer, pose: Pose, pdbstring: str, model: int
    ) -> GenericViewer:
        raise ModuleNotImplementedError(self.__class__.name__, BACKENDS[2])
