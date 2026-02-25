"""
PyInstaller runtime hook for lv_explorer

Force-import stdlib modules that PyInstaller normally excludes.
"""

import sys
import unittest  # noqa: F401
import unittest.mock  # noqa: F401
import pyparsing.testing  # noqa: F401
