import glob
import os
import pyrosetta
import pyrosetta.distributed.io as io
import tempfile
import unittest
import viewer3d

from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector,
    ResiduePropertySelector,
    SecondaryStructureSelector,
)


class TestViewer(unittest.TestCase):

    with tempfile.TemporaryDirectory() as workdir:

        def setUp(self, local_dir=workdir):

            if not os.path.isdir(local_dir):
                os.mkdir(local_dir)

            poses = [io.pose_from_sequence("TEST" * i) for i in range(1, 4)]
            for i, pose in enumerate(poses, start=1):
                with open(os.path.join(local_dir, "tmp_{0}.pdb".format(i)), "w") as f:
                    f.write(io.to_pdbstring(pose))

        def tearDown(self, local_dir=workdir):

            if os.path.isdir(local_dir):
                pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
                for pdbfile in pdbfiles:
                    os.remove(pdbfile)
                os.rmdir(local_dir)

        def test_viewer_with_pdbfiles(self, local_dir=workdir):

            pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
            viewer3d.presets.coreBoundarySurface(pdbfiles, continuous_update=True)
            viewer3d.presets.ligandsAndMetals(pdbfiles, window_size=(200.01, 200.01))
            view = viewer3d.init(pdbfiles, (1600, 400), delay=0.1234567890) \
                + viewer3d.setBackgroundColor("black") \
                + viewer3d.setStyle(style="line", colorscheme="blueCarbon")
            view.show()
            self.assertEqual(view.poses, [None] * len(pdbfiles))
            self.assertEqual(len(view.modules), 2)
            view.reinit()
            self.assertEqual(view.poses, [None] * len(pdbfiles))
            self.assertEqual(len(view.modules), 0)
            view.reset()
            self.assertIsNone(view.poses)
            self.assertIsNone(view.pdbstrings)

        def test_viewer_with_poses(self, local_dir=workdir):

            pdbfiles = glob.glob(os.path.join(local_dir, "*.pdb"))
            packed_poses = [io.pose_from_file(pdbfile) for pdbfile in pdbfiles]
            poses = [io.to_pose(p) for p in packed_poses]
            pose = poses[0]
            viewer3d.presets.coreBoundarySurface(packed_poses, delay=0)
            viewer3d.presets.ligandsAndMetals(packed_poses, continuous_update=True, window_size=(100., 100.))
            modules = [
                viewer3d.setBackgroundColor("grey"),
                viewer3d.setStyle(style="sphere", colorscheme="greenCarbon", radius=1.)
            ]
            view = sum([viewer3d.init(poses)] + modules)
            view()
            self.assertEqual(view.poses, poses)
            self.assertEqual(len(view.modules), 2)
            view.reinit()
            self.assertEqual(view.poses, poses)
            self.assertEqual(len(view.modules), 0)
            view.modules = modules
            self.assertEqual(len(view.modules), len(modules))
            view.reset()
            self.assertIsNone(view.poses)
            self.assertIsNone(view.pdbstrings)
            view = viewer3d.init(poses, modules=modules)
            self.assertEqual(len(view.modules), len(modules))
            view.reset()
            self.assertIsNone(view.modules)

            metals_selector = ResiduePropertySelector(ResidueProperty(31))
            ligands_selector = ResiduePropertySelector(ResidueProperty(2))
            view = viewer3d.init(poses, window_size=(800, 600)) \
                + viewer3d.setStyle() \
                + viewer3d.setStyle(residue_selector=ligands_selector, style="stick", colorscheme="magentaCarbon", radius=0.5) \
                + viewer3d.setStyle(residue_selector=metals_selector, style="sphere", colorscheme="chainHetatm", radius=1.5)
            view.reset()

            polar_residue_selector = ResiduePropertySelector(ResidueProperty(52))
            view = viewer3d.init(packed_poses)
            view.add(viewer3d.setStyle(radius=0.1))
            view.add(viewer3d.setStyle(residue_selector=polar_residue_selector, colorscheme="whiteCarbon", radius=0.25, label=False))
            view.add(viewer3d.setHydrogens(color="white", polar_only=True, radius=0.1))
            view.add(viewer3d.setHydrogenBonds(color="black"))
            view.add(viewer3d.setDisulfides(radius=0.1))
            view()
            view.reset()

            view = sum(
                [
                    viewer3d.init(poses),
                    viewer3d.setStyle(cartoon=False, style="sphere", radius=1.5, colorscheme="darkgreyCarbon"),
                    viewer3d.setZoom(factor=0.95)
                ]
            )
            view()
            view.reset()

            command_tuple = {"hetflag": True}, {"stick": {"singleBond": False, "colorscheme": "whiteCarbon", "radius": 0.25}}
            command_dict = {"hetflag": True}
            chA = ChainSelector("A")
            chB = ChainSelector("B")
            view = sum(
                [
                    viewer3d.init(poses),
                    viewer3d.setStyle(cartoon_color="lightgrey", radius=0.25),
                    viewer3d.setSurface(residue_selector=chA, colorscheme="greenCarbon", opacity=0.65, surface_type="VDW"),
                    viewer3d.setSurface(residue_selector=chB, color="blue", opacity=0.75, surface_type="SAS"),
                    viewer3d.setDisulfides(radius=0.25),
                    viewer3d.setZoom(factor=1.5),
                    viewer3d.setStyle(command=command_tuple),
                    viewer3d.setStyle(command=command_dict)
                ]
            )
            view()
            view.reset()

            helix_selector = SecondaryStructureSelector("H")
            sheet_selector = SecondaryStructureSelector("E")
            loop_selector = SecondaryStructureSelector("L")
            modules = [
                viewer3d.setBackgroundColor(color="black"),
                viewer3d.setStyle(residue_selector=helix_selector, cartoon_color="blue", label=False, radius=0),
                viewer3d.setStyle(residue_selector=sheet_selector, cartoon_color="red", label=False, radius=0),
                viewer3d.setStyle(residue_selector=loop_selector, cartoon_color="white", label=False, radius=0),
                viewer3d.setZoomTo(residue_selector=sheet_selector)
            ]
            viewer3d.init(poses, window_size=(1200, 600), modules=modules).show()

            view = viewer3d.init(pose, delay=0.15) \
                + viewer3d.setStyle(radius=0.1) \
                + viewer3d.setDisulfides(radius=0.1)
            backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
            minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            for _ in range(3):
                backrub.apply(pose)
                minimize.apply(pose)
                view.show()
            view.reset()

            def myCustomPreset(packed_and_poses_and_pdbs=None, *args, **kwargs):
                """
                Add a description of the preset Viewer here
                """
                # Add custrom ResidueSelectors
                metals_selector = ResiduePropertySelector(ResidueProperty(31))
                ligands_selector = ResiduePropertySelector(ResidueProperty(2))
                # Add custom Viewer commands
                view = viewer3d.init(packed_and_poses_and_pdbs=packed_and_poses_and_pdbs, *args, **kwargs) \
                    + viewer3d.setBackgroundColor("white") \
                    + viewer3d.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15) \
                    + viewer3d.setStyle(residue_selector=ligands_selector, style="stick", colorscheme="brownCarbon", radius=0.5, label=True) \
                    + viewer3d.setStyle(residue_selector=metals_selector, style="sphere", colorscheme="chainHetatm", radius=1.5, label=True) \
                    + viewer3d.setHydrogenBonds() \
                    + viewer3d.setDisulfides(radius=0.15) \
                    + viewer3d.setHydrogens(color="white", radius=0.033, polar_only=True) \
                    + viewer3d.setSurface(residue_selector=ligands_selector, surface_type="VDW", opacity=0.5, color="magenta") \
                    + viewer3d.setSurface(residue_selector=metals_selector, surface_type="VDW", opacity=0.5, color="magenta") \
                    + viewer3d.setZoomTo(residue_selector=ligands_selector)
                return view()

            myCustomPreset(pose)
