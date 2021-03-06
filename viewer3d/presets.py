"""
Display preset custom viewers for routine visualizations.
"""
import pyrosetta


import logging
import viewer3d

from ipywidgets.widgets import (
    Button,
    Checkbox,
    FloatSlider,
    IntSlider,
    HBox,
    Output,
    Text,
    ToggleButtons,
    VBox,
    Label,
)
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    LayerSelector,
    ResiduePropertySelector,
)
from viewer3d.converters import _to_backend
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

    # TODO description is currently not shown -- move it to dropdown in the future?
    advanced_label = Label(
        "Advanced parameters:",
        description="""
distance factor = 1 / (1 + exp( n*(d - m) ) ), where d is the distance of the neighbor from the residue CA, m is the midpoint of the distance falloff, and n is a falloff exponent factor that determines the sharpness of the distance falloff (with higher values giving sharper falloff near the midpoint distance).

angle factor = ( (cos(theta)+a)/(1+a) )^b, where theta is the angle between the CA-CB vector and the CA-neighbor vector, a is an offset factor that widens the cone somewhat, and b is an exponent that determines the sharpness of the angular falloff (with lower values resulting in a broader cone with a sharper edge falloff).
    """,
    )

    view.set_widgets(
        [
            core_cutoff,
            surface_cutoff,
            advanced_label,
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
def templatePreset(*args, **kwargs):
    """
    Add a description of the preset Viewer here
    """
    __author__ = ""
    view = viewer3d.init(*args, **kwargs)

    # Add custom Viewer commands here

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
    # TODO set crick parameters file
    #mb.set
    # TODO this is a bug in make bundle, because it does not expose the setters for residue name 
    # mb.residue_name(aa)
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
        step=0.05,
        value=0,
        description="omega0",
        continuous_update=continuous_update,
    )
    delta_omega1 = FloatSlider(
        min=-180,
        max=180,
        value=0,
        step=1,
        description="delta_omega1",
        style={"description_width": "initial"},
        continuous_update=continuous_update,
    )

    z0_offset = FloatSlider(
      min=-3,
      max=3,
      value=0,
      step=0.1,
      description="z0_offset",
      style={"description_width": "initial"},
      continuous_update=continuous_update,
    )

    z1_offset = FloatSlider(
      min=-3,
      max=3,
      value=0,
      step=0.1,
      description="z1_offset",
      style={"description_width": "initial"},
      continuous_update=continuous_update,
    )

    invert = Checkbox(value=False, description="invert")

    length.observe(on_length_change, names="value")
    r0.observe(on_param_change, names="value")
    omega0.observe(on_param_change, names="value")
    delta_omega1.observe(on_param_change, names="value")
    z0_offset.observe(on_param_change, names="value")
    z1_offset.observe(on_param_change, names="value")
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
            z0_offset,
            z1_offset,
            invert,
            save_box,
        ]
    )
    initialize_bundle()

    return view
