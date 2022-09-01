"""
Display simple preset custom viewers for routine visualizations.
"""
import bokeh.palettes
import collections
import logging
import pyrosetta
import pyrosetta.distributed.io as io
import viewer3d

from IPython.display import display
from ipywidgets.widgets import (
    Button,
    Checkbox,
    Dropdown,
    FloatSlider,
    HBox,
    Image,
    IntSlider,
    Label,
    Output,
    Text,
    ToggleButtons,
    VBox,
    interact,
    interactive,
)
from pyrosetta import Pose
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.scoring.sasa import SasaMethodHPMode
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector,
    LayerSelector,
    NotResidueSelector,
    OrResidueSelector,
    ResiduePropertySelector,
    TrueResidueSelector,
)
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import (
    PerResidueClashMetric,
    PerResidueEnergyMetric,
    PerResidueSasaMetric,
)
from pyrosetta.rosetta.protocols.hbnet import UnsatSelector

from viewer3d.config import COLORBAR_ATTR
from viewer3d.converters import _to_backend, _to_poses_pdbstrings
from viewer3d.pose import apply_metric_to_poses
from viewer3d.tracer import requires_init


_logger: logging.Logger = logging.getLogger("viewer3d.presets")
out = Output()


@requires_init
def coreBoundarySurface(
    packed_and_poses_and_pdbs=None,
    window_size=(1200, 800),
    continuous_update=True,
    backend=0,
):
    """
    Interactively visualize core, boundary, and surface layer residue selectors
    with cartoon representation.

    Args:
        packed_and_poses_and_pdbs: An optional `PackedPose`, `Pose`, or `str` of a valid path
            to a .pdb file, or an iterable of these objects.
            Default: `None`
        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.
            Default: `(1200, 800)`
        continuous_update: a `bool` object. When using the interactive slider widget,
            `False` restricts rendering to mouse button release events.
            Default: `True`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization. The currently supported backends are 'py3Dmol' and 'nglview'.
            Default: `0` or `py3Dmol`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

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
        viewer3d.setStyle(
            residue_selector=core_selector,
            cartoon=True,
            cartoon_color="white",
            colorscheme=0xF57900,
            style="stick",
            radius=0.25,
            label=False,
        ),
        viewer3d.setStyle(
            residue_selector=boundary_selector,
            cartoon=True,
            cartoon_color="white",
            colorscheme=0x00CC00,
            style="stick",
            radius=0.25,
            label=False,
        ),
        viewer3d.setStyle(
            residue_selector=surface_selector,
            cartoon=True,
            cartoon_color="white",
            colorscheme=0x729FCF,
            style="stick",
            radius=0.25,
            label=False,
        ),
        viewer3d.setDisulfides(radius=0.25),
    ]

    angle_exponent = FloatSlider(
        min=-4,
        max=4,
        step=0.1,
        value=2,
        description="angle_exponent",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    angle_shift_factor = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=0.5,
        description="angle_shift_factor",
        continuous_update=continuous_update,
        style={"description_width": "initial"},
    )
    dist_exponent = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=1,
        description="dist_exponent",
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
        description="dist_midpoint",
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

    view = viewer3d.init(
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
            boundary_selector.set_cutoffs(
                core=core_cutoff.new, surf=surface_cutoff.value
            )
            surface_selector.set_cutoffs(
                core=core_cutoff.new, surf=surface_cutoff.value
            )
        view.update_viewer()

    def set_surface_cutoff(surface_cutoff):
        with out:
            core_selector.set_cutoffs(core=core_cutoff.value, surf=surface_cutoff.new)
            boundary_selector.set_cutoffs(
                core=core_cutoff.value, surf=surface_cutoff.new
            )
            surface_selector.set_cutoffs(
                core=core_cutoff.value, surf=surface_cutoff.new
            )
        view.update_viewer()

    angle_exponent.observe(set_angle_exponent, names="value")
    angle_shift_factor.observe(set_angle_shift_factor, names="value")
    dist_exponent.observe(set_dist_exponent, names="value")
    denominator.observe(set_sc_neighbor_denominator, names="value")
    dist_midpoint.observe(set_sc_neighbor_dist_midpoint, names="value")
    core_cutoff.observe(set_core_cutoff, names="value")
    surface_cutoff.observe(set_surface_cutoff, names="value")

    advanced_labels = map(
        Label,
        (
            "Advanced parameters:",
            "distance factor = 1 / (1 + exp( n*(d - m) ) ), where d is the distance of the neighbor from the residue CA, m is the midpoint of the distance falloff, and n is a falloff exponent factor that determines the sharpness of the distance falloff (with higher values giving sharper falloff near the midpoint distance).",
            "angle factor = ( (cos(theta)+a)/(1+a) )^b, where theta is the angle between the CA-CB vector and the CA-neighbor vector, a is an offset factor that widens the cone somewhat, and b is an exponent that determines the sharpness of the angular falloff (with lower values resulting in a broader cone with a sharper edge falloff).",
        ),
    )

    view.set_widgets(
        [
            core_cutoff,
            surface_cutoff,
            *advanced_labels,
            dist_exponent,
            dist_midpoint,
            angle_exponent,
            angle_shift_factor,
            denominator,
        ]
    )
    view.update_viewer()

    return view


@requires_init
def ligandsAndMetals(*args, **kwargs):
    """
    Display residues with `ResidueProperty.LIGAND` as 'brownCarbon' sticks with opaque surface,
    and `ResidueProperty.METAL` as 'chainHetatm' spheres, with 'spectrum' cartoon representation,
    disulfide bonds, polar hydrogens, and dashed hydrogen bonds.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    with out:
        metals_selector = ResiduePropertySelector(ResidueProperty.METAL)
        ligands_selector = ResiduePropertySelector(ResidueProperty.LIGAND)

    view = (
        viewer3d.init(*args, **kwargs)
        + viewer3d.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15)
        + viewer3d.setStyle(
            residue_selector=ligands_selector,
            style="stick",
            colorscheme="brownCarbon",
            radius=0.5,
            label=True,
        )
        + viewer3d.setStyle(
            residue_selector=metals_selector,
            style="sphere",
            colorscheme="chainHetatm",
            radius=1.5,
            label=True,
        )
        + viewer3d.setHydrogenBonds()
        + viewer3d.setDisulfides(radius=0.15)
        + viewer3d.setHydrogens(color="white", radius=0.033, polar_only=True)
        + viewer3d.setSurface(
            residue_selector=ligands_selector,
            surface_type="VDW",
            opacity=0.5,
            colorscheme="brownCarbon",
        )
        + viewer3d.setZoomTo(residue_selector=ligands_selector)
    )

    return view


@requires_init
def perResidueClashMetric(
    poses,
    vmin=0,
    vmax=10,
    log=None,
    palette=None,
    backend=1,
):
    """
    Score the input pose(s) with `PerResidueClashMetric` and color sidechains by
    per-residue clash score, with cartoon backbone, polar hydrogens, hydrogen bonds,
    and disulfide bonds also shown.

    Args:
        poses: a required `Pose` object or iterable of `Pose` objects to score and display.
        vmin: a `float` or `int` object representing the minimum clash score value for color map.
            If `None`, set 'vmin' to the minimum scoretype value.
            Default: `0`
        vmax: a `float` or `int` object representing the maximum clash score value for color map.
            If `None`, set 'vmin' to the maximum scoretype value.
            Default: `10`
        log: `None` to map colors spaced evenly on a linear scale between 'vmin' to 'vmax'.
            If an `int` or `float` object is provided, map colors spaced evenly on a log
            scale with the base provided.
            Default: `None`
        palette: an iterable of `str` (or `int`) objects representing a color map.
            Default: `list(reversed(bokeh.palettes.Reds256))`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.
            Default: `1` or `nglview`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    c = PerResidueClashMetric()
    c.set_output_as_pdb_nums(output_as_pdb_nums=True)
    c.set_residue_selector(TrueResidueSelector())
    c.set_secondary_residue_selector(TrueResidueSelector())
    c.set_soft_dampening(dampening=0.33)
    c.set_use_hydrogens(use_hydrogens=True)
    c.set_use_soft_clash(soft_clash_check=True)
    apply_metric_to_poses(c, poses)
    v = viewer3d.init(poses, backend=backend)
    if palette is None:
        palette = list(reversed(bokeh.palettes.Reds256))
    v += viewer3d.setPerResidueRealMetric(
        scoretype="atomic_clashes",
        colorbar_label="Per-Residue Clashes",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(False, True),
    )
    v += viewer3d.setHydrogens(polar_only=True, color="lightgray")
    v += viewer3d.setHydrogenBonds()
    v += viewer3d.setDisulfides()

    return v


@requires_init
def perResidueEnergyMetric(
    poses,
    scorefxn=None,
    vmin=(-5),
    vmax=5,
    log=None,
    palette=None,
    backend=1,
):
    """
    Score the input pose(s) with `PerResidueEnergyMetric` and color sidechains by
    per-residue total score, with cartoon backbone, polar hydrogens, hydrogen bonds,
    and disulfide bonds also shown.

    Args:
        poses: a required `Pose` object or iterable of `Pose` objects to score and display.
        scorefxn: an optional scorefunction to use.
            Default: 'ref2015'
        vmin: a `float` or `int` object representing the minimum energy value for color map.
            If `None`, set 'vmin' to the minimum scoretype value.
            Default: `-5`
        vmax: a `float` or `int` object representing the maximum energy value for color map.
            If `None`, set 'vmin' to the maximum scoretype value.
            Default: `5`
        log: `None` to map colors spaced evenly on a linear scale between 'vmin' to 'vmax'.
            If an `int` or `float` object is provided, map colors spaced evenly on a log
            scale with the base provided.
            Default: `None`
        palette: an iterable of `str` (or `int`) objects representing a color map.
            Default: `list(bokeh.palettes.Greens256) + list(reversed(bokeh.palettes.Reds256))`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.
            Default: `1` or `nglview`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")
    e = PerResidueEnergyMetric()
    e.set_scorefunction(scorefxn)
    e.set_output_as_pdb_nums(output_as_pdb_nums=True)
    apply_metric_to_poses(e, poses)
    v = viewer3d.init(poses, backend=backend)
    if palette is None:
        palette = list(bokeh.palettes.Greens256) + list(
            reversed(bokeh.palettes.Reds256)
        )
    scorefxn_name = scorefxn.get_name()
    weights_ext = ".wts"
    if scorefxn_name.endswith(weights_ext):
        scorefxn_name = scorefxn_name[: -len(weights_ext)]
    v += viewer3d.setPerResidueRealMetric(
        scoretype="res_energy",
        colorbar_label=f"Per-Residue Energy ({scorefxn_name})",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(True, True),
    )
    v += viewer3d.setHydrogens(polar_only=True, color="lightgray")
    v += viewer3d.setHydrogenBonds()
    v += viewer3d.setDisulfides()

    return v


@requires_init
def perResidueSasaMetric(
    poses,
    mode=0,
    vmin=None,
    vmax=None,
    log=None,
    palette=None,
    backend=1,
):
    """
    Score the input pose(s) with `PerResidueSasaMetric` and color sidechains by
    per-residue SASA score, with cartoon backbone, polar hydrogens, hydrogen bonds,
    and disulfide bonds also shown.

    Args:
        poses: a required `Pose` object or iterable of `Pose` objects to score and display.
        mode: an optional `int` object to set the SASA mode:
            `0`: all SASA.
            `1`: hydrophobic only.
            `2`: polar only.
            Default: `0`
        vmin: a `float` or `int` object representing the minimum SASA value for color map.
            If `None`, set 'vmin' to the minimum scoretype value.
            Default: `None`
        vmax: a `float` or `int` object representing the maximum SASA value for color map.
            If `None`, set 'vmin' to the maximum scoretype value.
            Default: `None`
        log: `None` to map colors spaced evenly on a linear scale between 'vmin' to 'vmax'.
            If an `int` or `float` object is provided, map colors spaced evenly on a log
            scale with the base provided.
            Default: `None`
        palette: an iterable of `str` (or `int`) objects representing a color map.
            Default: `bokeh.palettes.Viridis256`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.
            Default: `1` or `nglview`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if mode == 0:
        sasa_mode = SasaMethodHPMode.ALL_SASA
    elif mode == 1:
        sasa_mode = SasaMethodHPMode.HYDROPHOBIC_SASA
    elif mode == 2:
        sasa_mode = SasaMethodHPMode.POLAR_SASA
    else:
        raise ValueError(
            f"The 'mode' argument must be an `int` object in {list(range(3))}."
        )
    s = PerResidueSasaMetric()
    s.set_mode(sasa_mode)
    s.set_output_as_pdb_nums(True)
    s.set_residue_selector(TrueResidueSelector())
    apply_metric_to_poses(s, poses)
    v = viewer3d.init(poses, backend=backend)
    if palette is None:
        palette = list(bokeh.palettes.Viridis256)
    v += viewer3d.setPerResidueRealMetric(
        scoretype="res_sasa",
        colorbar_label="Per-Residue SASA (Å$^{2}$)",
        vmin=vmin,
        vmax=vmax,
        radius=0.2,
        log=log,
        palette=palette,
        colorbar_extremes=(False, True),
    )
    v += viewer3d.setHydrogens(polar_only=True, color="lightgray")
    v += viewer3d.setHydrogenBonds()
    v += viewer3d.setDisulfides()

    return v


@requires_init
def unsatSelector(
    pose,
    scorefxn=None,
    hbond_energy_cutoff=(-0.5),
    backend=0,
):
    """
    Visualize residues with unsatisfied backbone amine and backbone carbonyl hydrogen bonds.

    Residues with an unsatisfied backbone amine are colored blue, residues with
    an unsatisfied backbone carbonyl are colored red, and residues with both an
    unsatisfied backbone amine and unsatisfied backbone carbonyl are colored yellow.
    Residues with satisfied backbone hydrogen bonds are colored black. Cartoon
    representation is shown and hydrogen bonds are shown in black dashed lines.

    Args:
        pose: a required `Pose` object to display.
        scorefxn: an optional scorefunction to use.
            Default: 'ref2015'
        hbond_energy_cutoff: an optional energy cutoff for selecting unsatisfied
            hydrogen bonds.
            Default: `-0.5`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.
            Default: `0` or `py3Dmol`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()

    unsat_amine_selector = UnsatSelector()
    unsat_amine_selector.set_scorefxn(scorefxn)
    unsat_amine_selector.set_consider_mainchain_only(False)
    unsat_amine_selector.set_hbond_energy_cutoff(hbond_energy_cutoff)
    unsat_amine_selector.set_legacy(False)
    unsat_amine_selector.set_mode(False)
    with out:
        unsat_amine_selector.apply(pose)

    unsat_carbonyl_selector = UnsatSelector()
    unsat_carbonyl_selector.set_scorefxn(scorefxn)
    unsat_carbonyl_selector.set_consider_mainchain_only(False)
    unsat_carbonyl_selector.set_hbond_energy_cutoff(hbond_energy_cutoff)
    unsat_carbonyl_selector.set_legacy(False)
    unsat_carbonyl_selector.set_mode(True)
    with out:
        unsat_carbonyl_selector.apply(pose)

    unsat_carbonyl_and_amine_selector = AndResidueSelector(
        unsat_amine_selector, unsat_carbonyl_selector
    )
    unsat_selector = OrResidueSelector(unsat_amine_selector, unsat_carbonyl_selector)
    not_unsat_selector = NotResidueSelector(unsat_selector)

    view = (
        viewer3d.init(pose, backend=backend)
        + viewer3d.setBackgroundColor("white")
        + viewer3d.setStyle(style="stick", colorscheme="blackCarbon", radius=0.15)
        + viewer3d.setHydrogenBonds(color="black")
        + viewer3d.setDisulfides(radius=0.15)
        + viewer3d.setHydrogens(
            residue_selector=not_unsat_selector,
            color="gray",
            radius=0.033,
            polar_only=True,
        )
        + viewer3d.setStyle(
            residue_selector=unsat_amine_selector,
            style="stick",
            colorscheme="blueCarbon",
            radius=0.15,
            label=True,
        )
        + viewer3d.setHydrogens(
            residue_selector=unsat_amine_selector,
            color="white",
            radius=0.033,
            polar_only=True,
        )
        + viewer3d.setStyle(
            residue_selector=unsat_carbonyl_selector,
            style="stick",
            colorscheme="redCarbon",
            radius=0.15,
            label=True,
        )
        + viewer3d.setStyle(
            residue_selector=unsat_carbonyl_and_amine_selector,
            style="stick",
            colorscheme="yellowCarbon",
            radius=0.15,
            label=True,
        )
        + viewer3d.setHydrogens(
            residue_selector=unsat_carbonyl_selector,
            color="white",
            radius=0.033,
            polar_only=True,
        )
        + viewer3d.setZoomTo(residue_selector=unsat_selector)
    )

    return view


@requires_init
def makeBundle(
    modules=[],
    aa="VAL",
    num_helices=4,
    backend=0,
    window_size=None,
    continuous_update=True,
):
    """
    Interactively build and visualize a helical bundle parametrically with
    core, boundary, and surface layer residue selector representations.

    Args:
        modules: an optional `list` object containing instantiated visualization modules.
            Default: `[]`
        aa: an optional `str` object representing the 3-letter amino acid for
            the poly-XXX helical bundle.
            Default: 'VAL'
        num_helices: an `int` object representing the number of helices to generate.
            Default: `4`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization.
            Default: `0` or `py3Dmol`
        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.
            Default: `(1200, 800)`
        continuous_update: a `bool` object. When using the interactive widgets,
            `False` restricts rendering to mouse button release events.
            Default: `False`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Ajasja Ljubetic, Jason C. Klima"

    from pyrosetta.rosetta.protocols.helical_bundle import (
        BPC_delta_omega0,
        BPC_r0,
        BPC_invert_helix,
        MakeBundle,
    )
    from pyrosetta.rosetta.protocols.simple_moves import AddPDBInfoMover
    from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import (
        construct_poly_XXX_pose,
    )
    from pyrosetta.rosetta.utility import vector1_unsigned_long

    backend = _to_backend(backend)
    if not modules:
        with out:
            core_selector = LayerSelector()
            core_selector.set_layers(True, False, False)
            boundary_selector = LayerSelector()
            boundary_selector.set_layers(False, True, False)
            surface_selector = LayerSelector()
            surface_selector.set_layers(False, False, True)
        modules = [
            viewer3d.setStyle(
                residue_selector=core_selector,
                cartoon=False if backend == "nglview" else True,
                colorscheme=0xF57900,
                style="stick",
                radius=0.25,
                label=False,
            ),
            viewer3d.setStyle(
                residue_selector=boundary_selector,
                cartoon=False if backend == "nglview" else True,
                colorscheme=0x00CC00,
                style="stick",
                radius=0.25,
                label=False,
            ),
            viewer3d.setStyle(
                residue_selector=surface_selector,
                cartoon=False if backend == "nglview" else True,
                colorscheme=0x729FCF,
                style="stick",
                radius=0.25,
                label=False,
            ),
        ]
    if backend == "nglview":
        modules.append(
            viewer3d.setStyle(
                cartoon=True,
                radius=0,
                label=False,
            )
        )

    pose = pyrosetta.Pose()
    view = viewer3d.init(
        packed_and_poses_and_pdbs=pose,
        window_size=window_size,
        modules=modules,
        backend=backend,
    )
    with out:
        mb = MakeBundle()
    mb.set_reset_pose(True)
    mb.set_use_degrees(True)
    add_pdb_info_mover = AddPDBInfoMover()

    def make_poly_X(pose):
        positions = vector1_unsigned_long()
        for i in range(1, pose.size() + 1):
            positions.append(i)
        construct_poly_XXX_pose(
            aa=aa,
            pose=pose,
            positions=positions,
            keep_pro=False,
            keep_gly=False,
            keep_disulfide_cys=True,
        )

    def update_bundle():
        with out:
            mb.apply(pose)
        add_pdb_info_mover.apply(pose)
        make_poly_X(pose)
        view.update_viewer()

    def initialize_bundle():
        for i in range(1, num_helices + 1):
            mb.add_helix()
            mb.helix(i).set_helix_length(length.value)
            mb.helix(i).calculator_op().real_parameter(BPC_delta_omega0).set_value(
                360 / num_helices * (i - 1)
            )
            mb.helix(i).calculator_op().real_parameter(BPC_r0).set_value(
                r0.value
            )  # in angstrem
        update_bundle()

    def on_length_change(change):
        for i in range(1, num_helices + 1):
            if chosen_helix.value == "all" or chosen_helix.value == i:
                mb.helix(i).set_helix_length(length.value)
        update_bundle()

    def on_param_change(change):
        """Takes the name of the parameter from the change.owner.description and se"""
        for i in range(1, num_helices + 1):
            if chosen_helix.value == "all" or chosen_helix.value == i:
                param_enum = getattr(
                    pyrosetta.rosetta.protocols.helical_bundle,
                    f"BPC_{change.owner.description}",
                )
                mb.helix(i).calculator_op().real_parameter(param_enum).set_value(
                    change.new
                )
        update_bundle()

    def on_invert_change(change):
        for i in range(1, num_helices + 1):
            if chosen_helix.value == "all" or chosen_helix.value == i:
                mb.helix(i).calculator_op().boolean_parameter(
                    BPC_invert_helix
                ).set_value(bool(change.new))
        update_bundle()

    chosen_helix = ToggleButtons(
        options=["all"] + [i + 1 for i in range(num_helices)],
        description="chosen_helix",
    )
    r0 = FloatSlider(
        min=1,
        max=10,
        step=0.1,
        value=5,
        description="r0",
        continuous_update=continuous_update,
    )
    length = IntSlider(
        min=14,
        max=50,
        value=28,
        description="length",
        continuous_update=continuous_update,
    )
    omega0 = FloatSlider(
        min=-5,
        max=5,
        step=0.1,
        value=0,
        description="omega0",
        continuous_update=continuous_update,
    )
    delta_omega1 = FloatSlider(
        min=0,
        max=360,
        description="delta_omega1",
        style={"description_width": "initial"},
        continuous_update=continuous_update,
    )
    invert = Checkbox(value=False, description="invert")

    length.observe(on_length_change, names="value")
    r0.observe(on_param_change, names="value")
    omega0.observe(on_param_change, names="value")
    delta_omega1.observe(on_param_change, names="value")
    invert.observe(on_invert_change, names="value")

    save_button = Button(description="save PDB")
    save_edit = Text(value="bundle.pdb", description="filename")

    def save_pdb(sender):
        pose.dump_pdb(save_edit.value)

    save_button.on_click(save_pdb)
    save_box = HBox([save_button, save_edit], description="save_box")

    view.set_widgets(
        [
            chosen_helix,
            length,
            r0,
            omega0,
            delta_omega1,
            invert,
            save_box,
        ]
    )
    initialize_bundle()

    return view


@requires_init
def rosettaViewer(
    packed_and_poses_and_pdbs=None,
    window_size=(1200, 800),
    continuous_update=True,
    backend=1,
):
    """
    Interactively visualize the following `viewer3d` presets:
        0: perResidueEnergyMetric
        1: perResidueClashMetric
        2: perResidueSasaMetric
        3: unsatSelector
        4: ligandsAndMetals

    TODO:
        - If using `pyrosetta.pose_from_sequence`, `py3Dmol` may not show `perResidueClashMetric` correctly.

    Args:
        packed_and_poses_and_pdbs: An optional `PackedPose`, `Pose`, or `str` of a valid path
            to a .pdb file, or an iterable of these objects.
            Default: `None`
        window_size: an optional `list` or `tuple` of `int` or `float` values for the
            (width, height) dimensions of the displayed window screen size.
            Default: `(1200, 800)`
        continuous_update: a `bool` object. When using the interactive slider widget,
            `False` restricts rendering to mouse button release events.
            Default: `True`
        backend: an optional `str` or `int` object representing the backend to use for
            the visualization. The currently supported backends are 'py3Dmol' and 'nglview'.
            Default: `1` or `nglview`

    Returns:
        A `Py3DmolViewer` instance, a `NGLViewViewer` instance, or a `PyMOLViewer` instance.
    """
    __author__ = "Jason C. Klima"

    presets = (
        perResidueEnergyMetric,
        perResidueClashMetric,
        perResidueSasaMetric,
        unsatSelector,
        ligandsAndMetals,
    )  # Presets to display using an IntSlider widget
    _preset_scoretypes = ("res_energy", "atomic_clashes", "res_sasa")
    backend = _to_backend(backend)
    view = viewer3d.init(
        packed_and_poses_and_pdbs,
        window_size=window_size,
        continuous_update=continuous_update,
        backend=backend,
    )

    def on_preset(i):
        preset = presets[i]
        _poses, _pdbstrings = _to_poses_pdbstrings(packed_and_poses_and_pdbs)
        # Set up poses
        poses = []
        for (i, p) in _poses.items():
            poses.extend(p)
        # Clear scores
        for pose in poses:
            for scoretype in list(pose.scores.keys()):
                for preset_scoretype in _preset_scoretypes:
                    if preset_scoretype in scoretype:
                        try:
                            pose.scores.pop(scoretype)
                        except KeyError:
                            pass
                        finally:
                            break
        # Score
        v = preset(
            poses,
            backend=backend,
        )
        view.poses = v.poses
        view.pdbstrings = v.pdbstrings
        view.set_modules(v.get_modules())
        view.update_viewer()
        if hasattr(view.viewer, COLORBAR_ATTR):
            _value = getattr(view.viewer, COLORBAR_ATTR)
            display(Image(value=_value))
        else:
            display(Label(value=""))

    dropdown_options = [
        (preset.__name__, i) for (i, preset) in enumerate(presets, start=0)
    ]
    dropdown = Dropdown(
        options=dropdown_options,
        value=0,
        description="Preset",
    )
    preset_widget = interactive(
        on_preset,
        i=dropdown,
    )
    view.set_widgets(preset_widget)

    return view


@requires_init
def templatePreset(*args, **kwargs):
    """
    Add a description of the preset Viewer here
    """
    __author__ = ""
    view = viewer3d.init(*args, **kwargs)

    # Add custom Viewer commands here

    return view
