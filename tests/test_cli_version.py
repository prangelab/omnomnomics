import contextlib
import io
import sys
import unittest
from unittest.mock import patch

from omnomnomics import __version__
from omnomnomics.cli import main


class CliVersionTests(unittest.TestCase):
    def test_top_level_version(self):
        output = io.StringIO()
        with patch.object(sys, "argv", ["omnomnomics", "--version"]):
            with contextlib.redirect_stdout(output):
                main()

        self.assertEqual(output.getvalue().strip(), f"omnomnomics {__version__}")
