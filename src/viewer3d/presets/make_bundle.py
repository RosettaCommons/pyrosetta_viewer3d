__author__ = "Ajasja Ljubetic, Jason C. Klima"

import logging
import pyrosetta

from ipywidgets.widgets import (
    Button,
    Checkbox,
    FloatSlider,
    HBox,
    IntSlider,
    Output,
    Text,
    ToggleButtons,
)

from pyrosetta.rosetta.core.select.residue_selector import LayerSelector

from viewer3d.converters import _to_backend
from viewer3d.core import (
    Py3DmolViewer,
    NGLViewViewer,
    PyMOLViewer,
    init,
)
from viewer3d.modules.base import ModuleBase
from viewer3d.modules import (
    setBackgroundColor,
    setStyle,
)
from viewer3d.tracer import requires_init
from viewer3d.type_defs import (
    List,
    Optional,
    Tuple,
    Union,
)

out = Output()
_logger: logging.Logger = logging.getLogger("viewer3d.presets.make_bundle")


@requires_init
def makeBundle(
    modules: List[ModuleBase] = [],
    aa: str = "VAL",
    num_helices: int = 4,
    backend: Union[int, str] = 0,
    window_size: Optional[
        Union[Tuple[Union[int, float], Union[int, float]], List[Union[int, float]]]
    ] = None,
    continuous_update: bool = True,
) -> Union[Py3DmolViewer, NGLViewViewer, PyMOLViewer]:
    """
    Interactively build and visualize a helical bundle parametrically with core, boundary,
    and surface layer residue selector representations.

    Args:
        modules: An optional `list` object containing instantiated visualization modules.

        aa: An optional `str` object representing the 3-letter amino acid for the poly-XXX
            helical bundle.

        num_helices: An `int` object representing the number of helices to generate.

        backend: An optional `str` or `int` object representing the backend to use for the
            visualization.

        window_size: An optional `list` or `tuple` of `int` or `float` values for the
            `(width, height)` dimensions of the displayed window screen size.

        continuous_update: A `bool` object. When using the interactive widgets, `False`
            restricts rendering to mouse button release events.

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
            setBackgroundColor(color="white"),
            setStyle(
                residue_selector=core_selector,
                cartoon=backend == "py3Dmol",
                cartoon_color=0x000000,
                colorscheme=0xF57900,
                style="stick",
                radius=0.25,
                label=False,
            ),
            setStyle(
                residue_selector=boundary_selector,
                cartoon=backend == "py3Dmol",
                cartoon_color=0x000000,
                colorscheme=0x00CC00,
                style="stick",
                radius=0.25,
                label=False,
            ),
            setStyle(
                residue_selector=surface_selector,
                cartoon=backend == "py3Dmol",
                cartoon_color=0x000000,
                colorscheme=0x729FCF,
                style="stick",
                radius=0.25,
                label=False,
            ),
        ]
        if backend != "py3Dmol":
            modules.append(
                setStyle(
                    cartoon=True,
                    cartoon_color=0x000000,
                    cartoon_radius=0.25,
                    cartoon_opacity=0.9,
                    radius=0.0,
                    label=False,
                )
            )

    pose = pyrosetta.Pose()
    view = init(
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
    # mb.set
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
            mb.helix(i).calculator_op().real_parameter(BPC_r0).set_value(r0.value)  # in angstrem
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
                mb.helix(i).calculator_op().real_parameter(param_enum).set_value(change.new)
        update_bundle()

    def on_invert_change(change):
        for i in range(1, num_helices + 1):
            if chosen_helix.value == "all" or chosen_helix.value == i:
                mb.helix(i).calculator_op().boolean_parameter(BPC_invert_helix).set_value(
                    bool(change.new)
                )
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
