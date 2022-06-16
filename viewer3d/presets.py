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
)
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    LayerSelector,
    ResiduePropertySelector,
)
from viewer3d.converters import _to_backend
from viewer3d.tracer import requires_init


_logger = logging.getLogger("viewer3d.presets")
out = Output()


@requires_init
def coreBoundarySurface(
    packed_and_poses_and_pdbs=None,
    window_size=(1200, 800),
    continuous_update=True,
    backend="py3Dmol",
):
    """
    Display core residues as 'blackCarbon' sticks, boundary residues as 'greyCarbon' sticks, and surface residues
    as 'whiteCarbon' sticks, with 'spectrum' cartoon representation, using the default arguments in
    `pyrosetta.rosetta.core.select.residue_selector.LayerSelector()` to select layers.
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
            cartoon=True,
            radius=0,
            label=False,
        ),
        viewer3d.setStyle(
            residue_selector=core_selector,
            cartoon=False,
            colorscheme=0xF57900,
            style="stick",
            radius=0.25,
            label=False,
        ),
        viewer3d.setStyle(
            residue_selector=boundary_selector,
            cartoon=False,
            colorscheme=0x00CC00,
            style="stick",
            radius=0.25,
            label=False,
        ),
        viewer3d.setStyle(
            residue_selector=surface_selector,
            cartoon=False,
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
    )
    angle_shift_factor = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=0.5,
        description="angle_shift_factor",
        continuous_update=continuous_update,
    )
    dist_exponent = FloatSlider(
        min=-2,
        max=2,
        step=0.1,
        value=1,
        description="dist_exponent",
        continuous_update=continuous_update,
    )
    denominator = FloatSlider(
        min=0.1,
        max=10,
        step=0.1,
        value=1,
        description="denominator",
        continuous_update=continuous_update,
    )
    dist_midpoint = FloatSlider(
        min=0,
        max=20,
        step=1,
        value=9,
        description="dist_midpoint",
        continuous_update=continuous_update,
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

    angle_exponent.observe(set_angle_exponent, names="value")
    angle_shift_factor.observe(set_angle_shift_factor, names="value")
    dist_exponent.observe(set_dist_exponent, names="value")
    denominator.observe(set_sc_neighbor_denominator, names="value")
    dist_midpoint.observe(set_sc_neighbor_dist_midpoint, names="value")

    view.set_widgets(
        [angle_exponent, angle_shift_factor, dist_exponent, denominator, dist_midpoint]
    )
    view.update_viewer()

    return view


@requires_init
def ligandsAndMetals(*args, **kwargs):
    """
    Display residues with `ResidueProperty.LIGAND` as 'brownCarbon' sticks with opaque surface,
    and `ResidueProperty.METAL` as 'chainHetatm' spheres, with 'spectrum' cartoon representation,
    disulfide bonds, polar hydrogens, and dashed hydrogen bonds.
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
    backend="py3Dmol",
    window_size=None,
    continuous_update=True,
):
    """
    Add a description of the preset Viewer here
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
        view.update_pose(pose)

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
    save_edit = Text(value="bundle.pdb")

    def save_pdb(sender):
        pose.dump_pdb(save_edit.value)

    save_button.on_click(save_pdb)
    save_box = HBox([save_button, save_edit])

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
