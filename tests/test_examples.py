__author__ = "Jason C. Klima"

try:
    import papermill as pm
except ImportError:
    print(
        "Testing the Jupyter notebook examples requires the 'papermill' package. "
        "Please install this package into your virtual environment to continue. "
        "For more information, please visit:\n"
        "https://pypi.org/project/papermill/\n"
    )
    raise

import glob
import logging
import tempfile
import unittest

from tests.utils import (
    has_nglview,
    has_py3Dmol,
    has_pymol,
    set_temp_env,
)

_logger: logging.Logger = logging.getLogger("viewer3d.tests.test_examples")


class TestNotebooks(unittest.TestCase):
    """Smoke test for `viewer3d` Jupyter notebook examples."""

    def run_notebook(self, input_path, backend):
        output_path = tempfile.NamedTemporaryFile(suffix=".ipynb").name
        with set_temp_env("PAPERMILL_EXECUTION", "1"):
            pm.execute_notebook(
                input_path=input_path,
                output_path=output_path,
                parameters={"backend": backend},
                kernel_name="python3",
            )

    @unittest.skipIf(not has_py3Dmol(), "The 'py3Dmol' package is not installed.")
    def test_py3Dmol(self):
        for notebook in glob.glob("examples/*.ipynb"):
            self.run_notebook(notebook, backend=0)

    @unittest.skipIf(not has_nglview(), "The 'nglview' package is not installed.")
    def test_nglview(self):
        for notebook in glob.glob("examples/*.ipynb"):
            self.run_notebook(notebook, backend=1)

    @unittest.skipIf(not has_pymol(), "The 'pymol' package is not installed.")
    def test_pymol(self):
        for notebook in glob.glob("examples/*.ipynb"):
            self.run_notebook(notebook, backend=2)
