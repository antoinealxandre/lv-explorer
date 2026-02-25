#!/usr/bin/env python3
"""
LV Explorer — Main entry point
================================
Usage:
    python main.py
    python main.py --data-path "data/extract/PRAD18"
"""

import sys
import os

# Ensure the repo root is on sys.path so `lv_explorer` is importable
# whether the script is run directly or bundled by PyInstaller.
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from lv_explorer.run import main

if __name__ == "__main__":
    main()
