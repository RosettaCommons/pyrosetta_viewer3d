"""
Display preset custom viewers for routine visualizations.
"""

import logging

from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    LayerSelector,
    ResiduePropertySelector,
)

_logger = logging.getLogger("viewer3d.presets")


def coreBoundarySurface(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Display core residues as 'blackCarbon' sticks, boundary residues as 'greyCarbon' sticks, and surface residues
    as 'whiteCarbon' sticks, with 'spectrum' cartoon representation, using the default arguments in
    `pyrosetta.rosetta.core.select.residue_selector.LayerSelector()` to select layers.
    """
    __author__ = "Jason C. Klima"

    core_selector = LayerSelector()
    core_selector.set_layers(True, False, False)
    boundary_selector = LayerSelector()
    boundary_selector.set_layers(False, True, False)
    surface_selector = LayerSelector()
    surface_selector.set_layers(False, False, True)

    view = viewer.init(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs
    )
    view.add(viewer.setStyle())
    view.add(
        viewer.setStyle(
            residue_selector=core_selector,
            style="stick",
            colorscheme="blackCarbon",
            radius=0.25,
            label=False,
        )
    )
    view.add(
        viewer.setStyle(
            residue_selector=boundary_selector,
            style="stick",
            colorscheme="greyCarbon",
            radius=0.25,
            label=False,
        )
    )
    view.add(
        viewer.setStyle(
            residue_selector=surface_selector,
            style="stick",
            colorscheme="whiteCarbon",
            radius=0.25,
            label=False,
        )
    )
    view.add(viewer.setDisulfides(radius=0.25))

    return view.show()


def ligandsAndMetals(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Display residues with `ResidueProperty.LIGAND` as 'brownCarbon' sticks with opaque surface,
    and `ResidueProperty.METAL` as 'chainHetatm' spheres, with 'spectrum' cartoon representation,
    disulfide bonds, polar hydrogens, and dashed hydrogen bonds.
    """
    __author__ = "Jason C. Klima"

    metals_selector = ResiduePropertySelector(ResidueProperty.METAL)
    ligands_selector = ResiduePropertySelector(ResidueProperty.LIGAND)

    view = (
        viewer.init(
            packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs
        )
        + viewer.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15)
        + viewer.setStyle(
            residue_selector=ligands_selector,
            style="stick",
            colorscheme="brownCarbon",
            radius=0.5,
            label=True,
        )
        + viewer.setStyle(
            residue_selector=metals_selector,
            style="sphere",
            colorscheme="chainHetatm",
            radius=1.5,
            label=True,
        )
        + viewer.setHydrogenBonds()
        + viewer.setDisulfides(radius=0.15)
        + viewer.setHydrogens(color="white", radius=0.033, polar_only=True)
        + viewer.setSurface(
            residue_selector=ligands_selector,
            surface_type="VDW",
            opacity=0.5,
            colorscheme="brownCarbon",
        )
        + viewer.setZoomTo(residue_selector=ligands_selector)
    )

    return view()


def templatePreset(packed_and_poses_and_pdbs=None, *args, **kwargs):
    """
    Add a description of the preset Viewer here
    """
    view = viewer.init(
        packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs
    )

    # Add custom Viewer commands here

    return view()
