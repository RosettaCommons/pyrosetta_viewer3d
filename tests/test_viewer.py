__author__ = "Jason C. Klima"

import bokeh.palettes
import itertools
import glob
import logging
import os
import pyrosetta
import pyrosetta.distributed.io as io
import random
import tempfile
import unittest
import viewer3d as v3d

from io import StringIO
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector,
    ResiduePropertySelector,
    SecondaryStructureSelector,
    TrueResidueSelector,
)
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import (
    PerResidueClashMetric,
    PerResidueEnergyMetric,
    PerResidueSasaMetric,
)

from tests.utils import (
    has_biotite,
    has_nglview,
    has_py3Dmol,
    has_pymol,
)

if has_biotite():
    from biotite.structure.io.pdb import PDBFile


_logger: logging.Logger = logging.getLogger("viewer3d.tests.test_viewer")


class TestViewer(unittest.TestCase):
    """Smoke test for `viewer3d` library."""

    workdir = tempfile.TemporaryDirectory().name

    def setUp(self) -> None:
        random.seed(111)
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

        poses = [io.pose_from_sequence("TEST" * i) for i in range(1, 4)]
        for i, pose in enumerate(poses, start=1):
            with open(os.path.join(self.workdir, "tmp_{0}.pdb".format(i)), "w") as f:
                f.write(io.to_pdbstring(pose))

    def tearDown(self) -> None:
        if os.path.isdir(self.workdir):
            pdbfiles = glob.glob(os.path.join(self.workdir, "*.pdb"))
            for pdbfile in pdbfiles:
                os.remove(pdbfile)
            os.rmdir(self.workdir)

    @staticmethod
    def pose_to_atom_array(pose):
        pdbstring = io.to_pdbstring(pose)
        pdb = PDBFile.read(StringIO(pdbstring))

        return pdb.get_structure(model=1)

    def viewer_with_atom_arrays(self, backend=0) -> None:
        pdbfiles = glob.glob(os.path.join(self.workdir, "*.pdb"))
        atom_arrays = [
            TestViewer.pose_to_atom_array(pyrosetta.pose_from_file(pdbfile)) for pdbfile in pdbfiles
        ]
        pdbstrings = list(map(v3d.atom_array_to_pdbstring, atom_arrays))
        poses = list(map(v3d.atom_array_to_pdbstring, atom_arrays))
        v3d.presets.coreBoundarySurface(pdbstrings, continuous_update=True, backend=backend)
        v3d.presets.ligandsAndMetals(poses, window_size=(200.01, 200.01), backend=backend)
        view = (
            v3d.init(poses, (1600, 400), delay=0.1234567890, backend=backend)
            + v3d.setBackgroundColor("black")
            + v3d.setStyle(style="line", colorscheme="blueCarbon")
        )
        view.show()
        self.assertEqual(list(itertools.chain(*view.poses.values())), [None] * len(atom_arrays))
        self.assertEqual(len(view.modules), 2)
        view.clear()
        self.assertEqual(list(itertools.chain(*view.poses.values())), [None] * len(atom_arrays))
        self.assertEqual(len(view.modules), 0)
        view.reset()
        self.assertIsNone(view.poses)
        self.assertIsNone(view.pdbstrings)

    def viewer_with_pdbfiles(self, backend=0) -> None:
        pdbfiles = glob.glob(os.path.join(self.workdir, "*.pdb"))
        v3d.presets.coreBoundarySurface(pdbfiles, continuous_update=True, backend=backend)
        v3d.presets.ligandsAndMetals(pdbfiles, window_size=(200.01, 200.01), backend=backend)
        view = (
            v3d.init(pdbfiles, (1600, 400), delay=0.1234567890, backend=backend)
            + v3d.setBackgroundColor("black")
            + v3d.setStyle(style="line", colorscheme="blueCarbon")
        )
        view.show()
        self.assertEqual(list(itertools.chain(*view.poses.values())), [None] * len(pdbfiles))
        self.assertEqual(len(view.modules), 2)
        view.clear()
        self.assertEqual(list(itertools.chain(*view.poses.values())), [None] * len(pdbfiles))
        self.assertEqual(len(view.modules), 0)
        view.reset()
        self.assertIsNone(view.poses)
        self.assertIsNone(view.pdbstrings)

    def viewer_with_poses(self, backend=0) -> None:
        pdbfiles = glob.glob(os.path.join(self.workdir, "*.pdb"))
        packed_poses = [io.pose_from_file(pdbfile) for pdbfile in pdbfiles]
        poses = [io.to_pose(p) for p in packed_poses]
        pose = poses[0]
        v3d.presets.coreBoundarySurface(packed_poses, backend=backend)
        v3d.presets.ligandsAndMetals(
            packed_poses, continuous_update=True, window_size=(100.0, 100.0), backend=backend
        )
        modules = [
            v3d.setBackgroundColor("grey"),
            v3d.setStyle(style="sphere", colorscheme="greenCarbon", radius=1.0),
        ]
        view = sum([v3d.init(poses, backend=backend)] + modules)
        view()
        self.assertListEqual(list(itertools.chain(*view.poses.values())), poses)
        self.assertEqual(len(view.modules), 2)
        view.clear()
        self.assertListEqual(list(itertools.chain(*view.poses.values())), poses)
        self.assertEqual(len(view.modules), 0)
        view.modules = modules
        self.assertEqual(len(view.modules), len(modules))
        view.reset()
        self.assertIsNone(view.poses)
        self.assertIsNone(view.pdbstrings)
        view = v3d.init(poses, modules=modules, backend=backend)
        self.assertEqual(len(view.modules), len(modules))
        view.reset()
        self.assertEqual(view.modules, [])
        self.assertEqual(view.widgets, [])

        metals_selector = ResiduePropertySelector(ResidueProperty(31))
        ligands_selector = ResiduePropertySelector(ResidueProperty(2))
        view = (
            v3d.init(poses, window_size=(800, 600), backend=backend)
            + v3d.setStyle()
            + v3d.setStyle(
                residue_selector=ligands_selector,
                style="stick",
                colorscheme="magentaCarbon",
                radius=0.5,
            )
            + v3d.setStyle(
                residue_selector=metals_selector,
                style="sphere",
                colorscheme="chainHetatm",
                radius=1.5,
            )
        )
        view.reset()

        polar_residue_selector = ResiduePropertySelector(ResidueProperty(52))
        view = v3d.init(packed_poses, backend=backend)
        view.add(v3d.setStyle(radius=0.1))
        view.add(
            v3d.setStyle(
                residue_selector=polar_residue_selector,
                colorscheme="whiteCarbon",
                radius=0.25,
                label=False,
            )
        )
        view.add(v3d.setHydrogens(color="white", polar_only=True, radius=0.1))
        view.add(v3d.setHydrogenBonds(color="black"))
        view.add(v3d.setDisulfides(radius=0.1))
        view()
        view.reset()

        view = sum(
            [
                v3d.init(poses, backend=backend),
                v3d.setStyle(
                    cartoon=False,
                    style="sphere",
                    radius=1.5,
                    colorscheme="darkgreyCarbon",
                ),
                v3d.setZoom(factor=0.95),
            ]
        )
        view()
        view.reset()

        command_tuple = {"hetflag": True}, {
            "stick": {
                "singleBond": False,
                "colorscheme": "whiteCarbon",
                "radius": 0.25,
            }
        }
        command_dict = {"hetflag": True}
        chA = ChainSelector("A")
        chB = ChainSelector("B")
        view = sum(
            [
                v3d.init(poses, backend=backend),
                v3d.setStyle(cartoon_color="lightgrey", radius=0.25),
                v3d.setSurface(
                    residue_selector=chA,
                    colorscheme="greenCarbon",
                    opacity=0.65,
                    surface_type="VDW",
                ),
                v3d.setSurface(
                    residue_selector=chB,
                    color="blue",
                    opacity=0.75,
                    surface_type="SAS",
                ),
                v3d.setDisulfides(radius=0.25),
                v3d.setZoom(factor=1.5),
                v3d.setStyle(command=command_tuple),
                v3d.setStyle(command=command_dict),
            ]
        )
        view()
        view.reset()

        helix_selector = SecondaryStructureSelector("H")
        sheet_selector = SecondaryStructureSelector("E")
        loop_selector = SecondaryStructureSelector("L")
        modules = [
            v3d.setBackgroundColor(color="black"),
            v3d.setStyle(
                residue_selector=helix_selector,
                cartoon_color="blue",
                label=False,
                radius=0,
            ),
            v3d.setStyle(
                residue_selector=sheet_selector,
                cartoon_color="red",
                label=False,
                radius=0,
            ),
            v3d.setStyle(
                residue_selector=loop_selector,
                cartoon_color="white",
                label=False,
                radius=0,
            ),
            v3d.setZoomTo(residue_selector=sheet_selector),
        ]
        v3d.init(poses, window_size=(1200, 600), modules=modules, backend=backend).show()

        view = (
            v3d.init(pose, delay=0.15, backend=backend)
            + v3d.setStyle(radius=0.1)
            + v3d.setDisulfides(radius=0.1)
        )
        backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
        minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
        for _ in range(3):
            backrub.apply(pose)
            minimize.apply(pose)
            view.show()
        view.reset()

        scorefxn = pyrosetta.create_score_function("ref2015")
        e = PerResidueEnergyMetric()
        e.set_scorefunction(scorefxn)
        v = v3d.init(backend=backend, gui=True)
        palette = list(bokeh.palettes.Greens256) + list(reversed(bokeh.palettes.Reds256))
        v += v3d.setStyle(radius=0)
        v += v3d.setPerResidueRealMetric(
            scoretype="energy", vmin=-10, vmax=10, radius=0.2, log=10, palette=palette
        )
        v += v3d.setHydrogens(polar_only=True, color="lightgray")
        v += v3d.setHydrogenBonds()
        v += v3d.setDisulfides()
        for h in range(10):
            _pose = pose.clone()
            if hasattr(_pose, "cache"):
                _pose.cache.clear()
            else:
                _pose.scores.clear()
            for i in range(20):
                minimize.apply(_pose)
            e.apply(_pose)
            v.add_pose(_pose, index=h, update_viewer=False)
        v.show()

        viewer = v3d.presets.perResidueEnergyMetric(pose, backend=backend)
        viewer.show()
        viewer = v3d.presets.perResidueClashMetric(poses, backend=backend)
        viewer.show()
        viewer = v3d.presets.perResidueSasaMetric(pose, backend=backend)
        viewer.show()
        viewer = v3d.presets.unsatSelector(pose, backend=backend)
        viewer.show()
        viewer = v3d.presets.rosettaViewer(poses, backend=backend)
        viewer.show()

        def myCustomPreset(*args, **kwargs):
            """
            Add a description of the preset Viewer here
            """
            # Add custrom ResidueSelectors
            metals_selector = ResiduePropertySelector(ResidueProperty(31))
            ligands_selector = ResiduePropertySelector(ResidueProperty(2))
            # Add custom Viewer commands
            view = (
                v3d.init(*args, **kwargs)
                + v3d.setBackgroundColor("white")
                + v3d.setStyle(style="stick", colorscheme="lightgreyCarbon", radius=0.15)
                + v3d.setStyle(
                    residue_selector=ligands_selector,
                    style="stick",
                    colorscheme="brownCarbon",
                    radius=0.5,
                    label=True,
                )
                + v3d.setStyle(
                    residue_selector=metals_selector,
                    style="sphere",
                    colorscheme="chainHetatm",
                    radius=1.5,
                    label=True,
                )
                + v3d.setHydrogenBonds()
                + v3d.setDisulfides(radius=0.15)
                + v3d.setHydrogens(color="white", radius=0.033, polar_only=True)
                + v3d.setSurface(
                    residue_selector=ligands_selector,
                    surface_type="VDW",
                    opacity=0.5,
                    color="magenta",
                )
                + v3d.setSurface(
                    residue_selector=metals_selector,
                    surface_type="VDW",
                    opacity=0.5,
                    color="magenta",
                )
                + v3d.setZoomTo(residue_selector=ligands_selector)
            )
            return view()

        myCustomPreset(packed_and_poses_and_pdbs=pose, backend=backend)

        v = v3d.init(pose, delay=0, backend=backend)
        backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
        minimize = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
        v.set_modules([v3d.setStyle(), v3d.setDisulfides()])
        for h in range(5):
            for i in range(5):
                for j in range(5):
                    backrub.apply(pose)
                    minimize.apply(pose)
                v.add_pose(pose.clone(), index=h, update_viewer=False)
        v.show()

        v = v3d.init(pose, delay=0, backend=backend)
        v.set_modules(
            [
                v3d.setStyle(),
                v3d.setDisulfides(),
                v3d.setPerResidueRealMetric(
                    scoretype="atomic_clashes",
                    vmin=0,
                    vmax=8,
                    log=None,
                    cartoon=True,
                    cartoon_opacity=1,
                    show_hydrogens=False,
                    style="sphere",
                    radius=1,
                    palette=bokeh.palettes.Blues256,
                    colorbar=True,
                    colorbar_label="Per-Residue Clashes",
                    colorbar_extremes=(False, True),
                    colorbar_fontsize=12,
                    colorbar_nticks=20,
                ),
            ]
        )
        c = PerResidueClashMetric()
        c.set_output_as_pdb_nums(output_as_pdb_nums=True)
        c.set_residue_selector(TrueResidueSelector())
        c.set_secondary_residue_selector(TrueResidueSelector())
        c.set_soft_dampening(dampening=0.33)
        c.set_use_hydrogens(use_hydrogens=True)
        c.set_use_soft_clash(soft_clash_check=True)

        _pose = pose.clone()
        for h in range(5):
            for i in range(5):
                for j in range(5):
                    backrub.apply(pose)
                    minimize.apply(pose)
                if hasattr(_pose, "cache"):
                    _pose.cache.clear()
                else:
                    _pose.scores.clear()
                c.apply(_pose)
                e.apply(_pose)
                v.add_pose(_pose.clone(), index=h, update_viewer=False)
        v.show()
        v += v3d.setPerResidueRealMetric(
            scoretype="res_energy",
            vmin=-10,
            vmax=10,
            log=None,
            cartoon=True,
            cartoon_opacity=1,
            show_hydrogens=False,
            style="stick",
            radius=0.25,
            palette=bokeh.palettes.Reds256,
            colorbar=True,
            colorbar_label="Per-Residue Energy",
            colorbar_extremes=(True, True),
            colorbar_fontsize=8,
            colorbar_nticks=10,
        )

        _pose = pyrosetta.pose_from_sequence("TESTING")
        c.apply(_pose)
        v = v3d.init(_pose, backend=backend)
        v += v3d.setPerResidueRealMetric(
            scoretype="atomic_clashes",
            vmin=0,
            vmax=8,
            log=None,
            cartoon=True,
            cartoon_opacity=1,
            show_hydrogens=False,
            style="line",
            radius=1,
            palette=bokeh.palettes.Blues256,
            colorbar=True,
            colorbar_label="Per-Residue Clashes",
            colorbar_extremes=(False, True),
            colorbar_fontsize=12,
            colorbar_nticks=20,
        )
        v()
        v += v3d.setBackgroundColor("black")
        v()
        v += v3d.setStyle()
        v.show()
        v += v3d.setPerResidueRealMetric(scoretype="mystery_metric")
        if v._in_notebook():
            with self.assertRaises(ValueError):
                v.show()

        _pose = pyrosetta.pose_from_sequence("TEST/PER/RES/SASA")
        v = v3d.init(_pose, backend=backend, auto_show=True)
        s = PerResidueSasaMetric()
        s.apply(_pose)
        v += v3d.setPerResidueRealMetric(scoretype="res_sasa", style="cross")
        v.show()

        for rescale in (False, True):
            _pose = pyrosetta.pose_from_sequence("TEST/PLDDT")
            for res in range(1, _pose.size() + 1):
                _pose.pdb_info().temperature(
                    res=res,
                    atom_index=_pose.residue(res).atom_index("CA"),
                    t=random.uniform(0.0, 1.0 if rescale else 100.0),
                )
            v = v3d.alphaFoldPLDDT(_pose, rescale=rescale, backend=backend).show()

    @unittest.skipIf(not has_biotite(), "The 'biotite' package is not installed.")
    @unittest.skipIf(not has_py3Dmol(), "The 'py3Dmol' package is not installed.")
    def test_py3Dmol_with_atom_arrays(self):
        self.viewer_with_atom_arrays(backend=0)

    @unittest.skipIf(not has_py3Dmol(), "The 'py3Dmol' package is not installed.")
    def test_py3Dmol_with_pdbfiles(self):
        self.viewer_with_pdbfiles(backend=0)

    @unittest.skipIf(not has_py3Dmol(), "The 'py3Dmol' package is not installed.")
    def test_py3Dmol_with_poses(self):
        self.viewer_with_poses(backend=0)

    @unittest.skipIf(not has_biotite(), "The 'biotite' package is not installed.")
    @unittest.skipIf(not has_nglview(), "The 'nglview' package is not installed.")
    def test_nglview_with_atom_arrays(self):
        self.viewer_with_atom_arrays(backend=1)

    @unittest.skipIf(not has_nglview(), "The 'nglview' package is not installed.")
    def test_nglview_with_pdbfiles(self):
        self.viewer_with_pdbfiles(backend=1)

    @unittest.skipIf(not has_nglview(), "The 'nglview' package is not installed.")
    def test_nglview_with_poses(self):
        self.viewer_with_poses(backend=1)

    @unittest.skipIf(not has_biotite(), "The 'biotite' package is not installed.")
    @unittest.skipIf(not has_pymol(), "The 'pymol' package is not installed.")
    def test_pymol_with_atom_arrays(self):
        self.viewer_with_atom_arrays(backend=2)

    @unittest.skipIf(not has_pymol(), "The 'pymol' package is not installed.")
    def test_pymol_with_pdbfiles(self):
        self.viewer_with_pdbfiles(backend=2)

    @unittest.skipIf(not has_pymol(), "The 'pymol' package is not installed.")
    def test_pymol_with_poses(self):
        self.viewer_with_poses(backend=2)
