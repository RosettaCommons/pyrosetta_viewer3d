# viewer3d
Display `PackedPose` objects, `Pose` objects, or `.pdb` files within a Jupyter notebook and Google Colab.

# Description
The `viewer3d` macromolecular viewer quickly renders `.pdb` files, dynamically instantiating `Pose` objects if required for certain visualization modules (matching the name `viewer3d.set*`). So when adding visualization modules to the viewer or using presets, passing `Pose` or `PackedPose` objects to the viewer is suggested for quicker rendering. If a `Pose` object or `list`, `tuple`, or `set` of `Pose` objects are provided to the viewer, the `Pose`(s) pointer location(s) in memory remain fixed, and so the viewer can dynamically update upon `Pose` conformational changes by calling the `show()` method. The viewer applies visualization modules in the same order they are added (from left to right), so layering different styles (and `ResidueSelector`s) on top of one another becomes possible. The user must have already initialized PyRosetta providing `.params` files for any ligands and non-canonical residues in the input molecule(s), otherwise `pyrosetta.distributed` automatically initializes PyRosetta with default command line options.

# Contributing
Please open pull requests from a custom branch.

# PyPI Releases
1. Update version number in `setup.py`
2. `pip install --user --upgrade setuptools wheel twine`
3. `rm -rf dist/*`
4. `python setup.py sdist bdist_wheel`
5. `python -m twine upload dist/*`
