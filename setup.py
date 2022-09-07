from setuptools import setup

setup(
    name="viewer3d",
    version="1.0.0",
    description="Display PackedPose objects, Pose objects, or PDB files within a Jupyter notebook and Google Colab.",
    url="https://github.com/RosettaCommons/pyrosetta_viewer3d",
    author="Jason C. Klima, Ajasja Ljubetic",
    license="MIT",
    packages=["viewer3d"],  # these are the files that get packaged up
    install_requires=[
        "attrs>=18.2.0",
        "bokeh>=2.4.2",
        "ipywidgets>=7.5.1",
        "jupyter>=1.0.0",
        "matplotlib>=3.5",
        "py3Dmol>=1.8.0",
        # TODO: add conditional dependencies, for example pip install pyrosetta_viewer3d[nglview]
        # "pyrosetta>=2022.21", # Pyrosetta is not on pip, so adding it as an explicit dependency would only cause problems.
    ],
    zip_safe=False,
)
