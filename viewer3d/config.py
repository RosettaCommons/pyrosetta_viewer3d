from typing import Tuple

BACKENDS: Tuple[str, str, str] = ("py3Dmol", "nglview", "pymol")


def _import_backend(backend: str) -> None:
    if backend == BACKENDS[0]:
        import py3Dmol
    elif backend == BACKENDS[1]:
        import nglview
    elif backend == BACKENDS[2]:
        import pymol
