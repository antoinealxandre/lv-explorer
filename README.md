# LV Explorer

> **Alpha v0.1.0** — Ventricular Topography Analysis & Visualization

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![PyVista](https://img.shields.io/badge/PyVista-0.43%2B-green.svg)](https://pyvista.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

LV Explorer is a desktop application for interactive 3D visualization and quantitative analysis of left ventricular (LV) geometry derived from cardiac CT and MR imaging. It is designed for use in clinical research settings, particularly for ventricular tachycardia (VT) substrate mapping and ablation planning.

---

## Getting Started

### Prerequisites

- Python ≥ 3.8
- Windows 10/11 (also tested on Ubuntu 22.04)

### Option A — Run from source

```bash
# 1. Clone the repo
git clone https://github.com/antoinealxandre/lv-explorer.git
cd lv-explorer

# 2. Create a virtual environment
python -m venv .venv
.venv\Scripts\activate          # Windows
# source .venv/bin/activate     # Linux / macOS

# 3. Install dependencies
pip install -r requirements.txt

# 4. Launch
python main.py
# or with a specific patient folder:
python main.py --data-path "data/extract/PRAD18"
```

### Option B — Standalone `.exe` (Windows)

Download the latest release from the [Releases](https://github.com/YOUR_USERNAME/lv-explorer/releases) page.

```
LVExplorer_v0.1.0_alpha_win64.zip
└── LVExplorer.exe    ← double-click to launch
```

No Python installation required.

---
