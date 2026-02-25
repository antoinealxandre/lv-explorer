# -*- mode: python ; coding: utf-8 -*-
# =============================================================================
#  LV Explorer — PyInstaller spec
#  Usage:  pyinstaller lv_explorer.spec
# =============================================================================

import sys
from pathlib import Path

ROOT = Path(SPECPATH)

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
a = Analysis(
    [str(ROOT / "main.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=[
        # Bundle all assets (3D bust model)
        (str(ROOT / "assets"), "assets"),
    ],
    hiddenimports=[
        # Standard library (commonly missed by PyInstaller)
        "unittest",
        "unittest.mock",
        "heapq",
        "json",
        "csv",
        "pathlib",
        "glob",
        "tempfile",
        "shutil",
        "subprocess",
        "re",
        "collections",
        "operator",
        "functools",
        "itertools",
        "warnings",
        "locale",
        "codecs",
        "dataclasses",
        # PyVista / VTK (large import tree)
        "vtkmodules",
        "vtkmodules.all",
        "vtkmodules.util",
        "vtkmodules.util.numpy_support",
        "pyvista",
        "pyvista.plotting",
        "pyvista.utilities",
        "pyvistaqt",
        # Qt
        "qtpy",
        "PyQt5",
        "PyQt5.QtCore",
        "PyQt5.QtGui",
        "PyQt5.QtWidgets",
        "PyQt5.sip",
        # Scipy
        "scipy.spatial",
        "scipy.spatial._ckdtree",
        "scipy.stats",
        "scipy.special",
        "scipy._lib.messagestream",
        # Numpy / matplotlib & dependencies
        "numpy",
        "matplotlib",
        "matplotlib.backends.backend_qt5agg",
        "matplotlib.rcsetup",
        "matplotlib._fontconfig_pattern",
        "pyparsing",
        "pyparsing.testing",
        "packaging",
        "kiwisolver",
        "cycler",
        "dateutil",
        "dateutil.tz",
        # Optional dashboard dependency
        "lv_explorer.core.dashboard_manager",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "test",
        "unittest",
        "distutils",
    ],
    noarchive=False,
    optimize=0,
)

# ---------------------------------------------------------------------------
# PYZ
# ---------------------------------------------------------------------------
pyz = PYZ(a.pure)

# ---------------------------------------------------------------------------
# EXE  (one-dir build — faster startup than onefile, easier to debug)
# ---------------------------------------------------------------------------
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="LVExplorer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,       # No terminal window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    # icon="assets/icon.ico",  # Uncomment once you have an .ico file
)

# ---------------------------------------------------------------------------
# COLLECT  — assemble dist/LVExplorer/
# ---------------------------------------------------------------------------
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="LVExplorer",
)
