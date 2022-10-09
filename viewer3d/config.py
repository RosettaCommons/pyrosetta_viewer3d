from typing import Tuple

BACKENDS: Tuple[str, str, str] = ("py3Dmol", "nglview", "pymol")
URLS: Tuple[str, str, str] = (
    "https://pypi.org/project/py3Dmol/",
    "https://pypi.org/project/nglview/",
    "https://anaconda.org/schrodinger/pymol/",
)
COLORBAR_ATTR = "__viewer3d_colorbar__"


def _import_backend(backend: str) -> None:
    if backend == BACKENDS[0]:
        import py3Dmol
        import subprocess
        import xmlrpc.client
    elif backend == BACKENDS[1]:
        import nglview
    elif backend == BACKENDS[2]:
        import pymol
