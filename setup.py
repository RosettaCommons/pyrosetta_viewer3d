from setuptools import setup

setup(
    name="viewer3d",
    version="1.0.0",
    description="Display PackedPose objects, Pose objects, or PDB files within a Jupyter notebook and Google Colab.",
    url="https://github.com/RosettaCommons/pyrosetta_viewer3d",
    author="Jason C. Klima",
    license="MIT",
    packages=[
        "attrs>=18.2.0",
        "ipywidgets>=7.5.1",
        "jupyter>=1.0.0",
        "py3Dmol>=1.8.0",
        "pyrosetta>=2022.21",
    ],
    zip_safe=False,
)
