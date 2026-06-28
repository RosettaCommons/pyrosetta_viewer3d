__author__ = "Jason C. Klima"

import importlib
import os
import sys

from contextlib import contextmanager


def module_is_installed(module: str) -> bool:
    """Test if a module is already imported or installed."""

    if sys.modules.get(module):
        return True
    try:
        importlib.import_module(module)
        return True
    except ImportError:
        return False


def has_py3Dmol():
    return module_is_installed("py3Dmol")


def has_nglview():
    return module_is_installed("nglview")


def has_pymol():
    return module_is_installed("pymol")


def has_biotite():
    return module_is_installed("biotite")


@contextmanager
def set_temp_env(key, value):
    """Temporarily set an environment variable and restore the original."""

    old_value = os.environ.get(key)
    os.environ[key] = value
    try:
        yield
    finally:
        if old_value is None:
            del os.environ[key]
        else:
            os.environ[key] = old_value
