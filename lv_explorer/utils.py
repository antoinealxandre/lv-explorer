"""
lv_explorer.utils — Shared utility helpers
"""

import os
import sys


def resource_path(relative_path: str) -> str:
    """Return the absolute path to a bundled resource.

    Works both during development (running from source) and when the app is
    frozen by PyInstaller (``sys._MEIPASS`` is set to the temp extraction dir).

    Parameters
    ----------
    relative_path:
        Path relative to the repository/bundle root, e.g. ``"assets/human_bust.obj"``.

    Returns
    -------
    str
        Absolute path that can be passed to file-reading functions.
    """
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        # PyInstaller bundle — resources are extracted to sys._MEIPASS
        base = sys._MEIPASS  # type: ignore[attr-defined]
    else:
        # Running from source — go up from lv_explorer/ to repo root
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    return os.path.join(base, relative_path)
