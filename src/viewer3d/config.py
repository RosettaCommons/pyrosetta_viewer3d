__author__ = "Jason C. Klima"

import importlib

from viewer3d.type_defs import Tuple

BACKENDS: Tuple[str, str, str] = ("py3Dmol", "nglview", "pymol")
URLS: Tuple[str, str, str] = (
    "https://pypi.org/project/py3Dmol/",
    "https://pypi.org/project/nglview/",
    "https://pypi.org/project/pymol-open-source/",
)
COLORBAR_ATTR = "__viewer3d_colorbar__"


def _import_backend(backend: str) -> None:
    if backend not in BACKENDS:
        raise ValueError(f"Backend must be in {BACKENDS}. Received: {backend}")
    importlib.import_module(backend)
    if backend == BACKENDS[2]:
        importlib.import_module("subprocess")
        importlib.import_module("xmlrpc.client")
