#!/usr/bin/env python3
"""
LV Explorer - Lanceur principal
================================

Usage:
    python run.py [--data-path PATH_TO_PATIENT_FOLDER]

Examples:
    python run.py
    python run.py --data-path "topographie ventriculare/data/extract/PRAD19"
"""

import sys
import os
import argparse

# Add parent directory to path to allow importing lv_explorer package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lv_explorer import LVExplorerApp
from qtpy import QtWidgets


def main():
    
    parser = argparse.ArgumentParser(
        description='LV Explorer - Ventricular Analysis Application'
    )
    parser.add_argument(
        '--data-path',
        type=str,
        default=None,
        help='Path to patient data folder (optional)'
    )
    
    args = parser.parse_args()
    
    app = QtWidgets.QApplication.instance()
    if not app:
        app = QtWidgets.QApplication(sys.argv)
    
    window = LVExplorerApp(data_path=args.data_path)
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()