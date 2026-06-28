__author__ = "Jason C. Klima"

import io
import matplotlib.pyplot as plt
import PIL
import logging

from ipywidgets.widgets import (
    FloatSlider,
    HBox,
    HTML,
    Image,
    VBox,
    Output,
)
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.select.residue_selector import LayerSelector
from viewer3d.converters import _to_backend
from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.modules import (
    setBackgroundColor,
    setDisulfides,
    setStyle,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    Iterable,
    List,
    Optional,
    Tuple,
    Union,
)

out = Output()
_logger: logging.Logger = logging.getLogger("viewer3d.presets.core_boundary_surface")


@requires_init
def coreBoundarySurface(
    packed_and_poses_and_pdbs: Optional[
        Union[PackedPose, Pose, str, Iterable[Union[PackedPose, Pose, str]]]
    ] = None,
    colorschemes: Tuple[Union[int, str], Union[int, str], Union[int, str]] = (
        0xF57900,
        0x00CC00,
        0x729FCF,
    ),
    cartoon_color: Union[int, str] = "white",
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = (1200, 800),
    continuous_update: bool = True,
    backend: Union[int, str] = 0,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Interactively visualize core, boundary, and surface layer residue selectors with cartoon
    representation using the sidechain neighbor-specific options. Reference:
    https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector

    Args:
        packed_and_poses_and_pdbs: A `PackedPose`, `Pose`, or `str` of PDB string or a valid
            filesystem path to a ".pdb" file, or an iterable of these objects.

        colorschemes: A `3`-`tuple` of `int` or `str` objects representing the colorschemes for
            the core, boundary, and surface layers, respectively. See the `colorscheme` keyword
            argument of the `viewer3d.setStyle` docstring for help on available colorschemes.

        cartoon_color: A `str` or `int` representing the cartoon color.

        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.

        continuous_update: a `bool` object. When using the interactive slider widget,
            `False` restricts rendering to mouse button release events.

        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if not isinstance(colorschemes, tuple) or not len(colorschemes) == 3:
        raise ValueError("The 'colorschemes' argument value must be a 3-tuple.")

    with out:
        core_selector = LayerSelector()
        core_selector.set_layers(True, False, False)
        core_selector.set_use_sc_neighbors(True)
        boundary_selector = LayerSelector()
        boundary_selector.set_layers(False, True, False)
        boundary_selector.set_use_sc_neighbors(True)
        surface_selector = LayerSelector()
        surface_selector.set_layers(False, False, True)
        surface_selector.set_use_sc_neighbors(True)

    backend = _to_backend(backend)
    modules = [
        setStyle(
            residue_selector=core_selector,
            cartoon=backend == "py3Dmol",
            cartoon_color=cartoon_color,
            colorscheme=colorschemes[0],
            style="stick",
            radius=0.25,
            label=False,
        ),
        setStyle(
            residue_selector=boundary_selector,
            cartoon=backend == "py3Dmol",
            cartoon_color=cartoon_color,
            colorscheme=colorschemes[1],
            style="stick",
            radius=0.25,
            label=False,
        ),
        setStyle(
            residue_selector=surface_selector,
            cartoon=backend == "py3Dmol",
            cartoon_color=cartoon_color,
            colorscheme=colorschemes[2],
            style="stick",
            radius=0.25,
            label=False,
        ),
        setDisulfides(radius=0.25),
        setBackgroundColor(color="white"),
    ]
    if backend != "py3Dmol":
        modules.append(
            setStyle(
                cartoon=True,
                cartoon_color=cartoon_color,
                cartoon_radius=0.25,
                cartoon_opacity=0.9,
                radius=0.0,
                label=False,
            )
        )

    # Reference: https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector
    # sc_neighbor_dist_exponent=(1.0 &Real): Alters the value of n, the distance exponent.
    # sc_neighbor_dist_midpoint=(9.0 &Real): Alters the value of m, the distance falloff midpoint.
    # sc_neighbor_angle_shift_factor=(0.5 &Real): Alters the value of a, the angular shift factor.
    # sc_neighbor_angle_exponent=(2.0 &Real): Alters the value of b, the angular sharpness value.
    # sc_neighbor_denominator=(1.0 &Real): Alters the value by which the overall expression is divided.

    angle_exponent = FloatSlider(
        min=-4,
        max=4,
        step=0.1,
        value=2,
        description="angle_exponent (b)",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    angle_shift_factor = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=0.5,
        description="angle_shift_factor (a)",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    dist_exponent = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=1,
        description="dist_exponent (n)",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    denominator = FloatSlider(
        min=0.1,
        max=10,
        step=0.1,
        value=1,
        description="denominator",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    dist_midpoint = FloatSlider(
        min=0,
        max=20,
        step=1,
        value=9,
        description="dist_midpoint (m)",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    core_cutoff = FloatSlider(
        min=0,
        max=10,
        step=0.1,
        value=5.2,
        description="core_cutoff",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    surface_cutoff = FloatSlider(
        min=0,
        max=10,
        step=0.1,
        value=2,
        description="surface_cutoff",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )

    view = init(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs,
        window_size=window_size,
        modules=modules,
        backend=backend,
    )

    def set_angle_exponent(angle_exponent):
        with out:
            core_selector.set_angle_exponent(angle_exponent.new)
            boundary_selector.set_angle_exponent(angle_exponent.new)
            surface_selector.set_angle_exponent(angle_exponent.new)
        view.update_viewer()

    def set_angle_shift_factor(angle_shift_factor):
        with out:
            core_selector.set_angle_shift_factor(angle_shift_factor.new)
            boundary_selector.set_angle_shift_factor(angle_shift_factor.new)
            surface_selector.set_angle_shift_factor(angle_shift_factor.new)
        view.update_viewer()

    def set_dist_exponent(dist_exponent):
        with out:
            core_selector.set_dist_exponent(dist_exponent.new)
            boundary_selector.set_dist_exponent(dist_exponent.new)
            surface_selector.set_dist_exponent(dist_exponent.new)
        view.update_viewer()

    def set_sc_neighbor_denominator(denominator):
        with out:
            core_selector.set_sc_neighbor_denominator(denominator.new)
            boundary_selector.set_sc_neighbor_denominator(denominator.new)
            surface_selector.set_sc_neighbor_denominator(denominator.new)
        view.update_viewer()

    def set_sc_neighbor_dist_midpoint(dist_midpoint):
        with out:
            core_selector.set_sc_neighbor_dist_midpoint(dist_midpoint.new)
            boundary_selector.set_sc_neighbor_dist_midpoint(dist_midpoint.new)
            surface_selector.set_sc_neighbor_dist_midpoint(dist_midpoint.new)
        view.update_viewer()

    def set_core_cutoff(core_cutoff):
        with out:
            core_selector.set_cutoffs(core=core_cutoff.new, surf=surface_cutoff.value)
            boundary_selector.set_cutoffs(core=core_cutoff.new, surf=surface_cutoff.value)
            surface_selector.set_cutoffs(core=core_cutoff.new, surf=surface_cutoff.value)
        view.update_viewer()

    def set_surface_cutoff(surface_cutoff):
        with out:
            core_selector.set_cutoffs(core=core_cutoff.value, surf=surface_cutoff.new)
            boundary_selector.set_cutoffs(core=core_cutoff.value, surf=surface_cutoff.new)
            surface_selector.set_cutoffs(core=core_cutoff.value, surf=surface_cutoff.new)
        view.update_viewer()

    angle_exponent.observe(set_angle_exponent, names="value")
    angle_shift_factor.observe(set_angle_shift_factor, names="value")
    dist_exponent.observe(set_dist_exponent, names="value")
    denominator.observe(set_sc_neighbor_denominator, names="value")
    dist_midpoint.observe(set_sc_neighbor_dist_midpoint, names="value")
    core_cutoff.observe(set_core_cutoff, names="value")
    surface_cutoff.observe(set_surface_cutoff, names="value")

    def latex_to_image_bytes(latex, fontsize=11, dpi=300):
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.patch.set_alpha(0)
        plt.text(0, 0, f"${latex}$", fontsize=fontsize)
        plt.axis("off")
        buf = io.BytesIO()
        plt.savefig(buf, format="png", dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        buf.seek(0)
        img = PIL.Image.open(buf)
        bbox = img.getbbox()
        if bbox:
            img = img.crop(bbox)
        out = io.BytesIO()
        img.save(out, format="PNG")
        return out.getvalue()

    advanced_labels = VBox(
        [
            HTML("<b>Advanced parameters:</b>"),
            HBox(
                [
                    HTML("distance factor ="),
                    Image(
                        value=latex_to_image_bytes(
                            r"\frac{1}{1 + \exp(n(d - m))}",
                            fontsize=11,
                        ),
                        format="png",
                    ),
                    HTML(
                        "where `d` is the distance of the neighbor from the residue CA, "
                        "`m` is the midpoint of the distance falloff, and `n` is a falloff "
                        "exponent factor that determines the sharpness of the distance falloff "
                        "(with higher values giving sharper falloff near the midpoint distance)."
                    ),
                ]
            ),
            HBox(
                [
                    HTML("angle factor ="),
                    Image(
                        value=latex_to_image_bytes(
                            r"\left(\frac{\cos(\theta) + a}{1 + a}\right)^b",
                            fontsize=10,
                        ),
                        format="png",
                    ),
                    HTML(
                        "where `θ` is the angle between the CA-CB vector and the CA-neighbor vector, "
                        "`a` is an offset factor that widens the cone somewhat, and `b` is an exponent "
                        "that determines the sharpness of the angular falloff "
                        "(with lower values resulting in a broader cone with a sharper edge falloff)."
                    ),
                ]
            ),
        ]
    )

    view.set_widgets(
        [
            core_cutoff,
            surface_cutoff,
            advanced_labels,
            dist_exponent,
            dist_midpoint,
            angle_exponent,
            angle_shift_factor,
            denominator,
        ]
    )
    view.update_viewer()

    return view
