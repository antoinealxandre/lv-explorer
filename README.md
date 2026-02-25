# LV Explorer

> **Alpha v0.1.0** — Ventricular Topography Analysis & Visualization

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![PyVista](https://img.shields.io/badge/PyVista-0.43%2B-green.svg)](https://pyvista.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

LV Explorer is a desktop application for interactive 3D visualization and quantitative analysis of left ventricular (LV) geometry derived from cardiac CT and MR imaging. It is designed for use in clinical research settings, particularly for ventricular tachycardia (VT) substrate mapping and ablation planning.

---

## Features

- **3D Mesh Visualization** — interactive rendering of LV endo/epicardium, scars, coronaries and surrounding anatomy (PyVista / VTK)
- **Wall Thickness Maps** — multi-layer thickness analysis (1–5 mm isocontours)
- **Scar Segmentation** — dense scar, border zone, transmurality maps from LGE-MRI
- **Conduction Velocity Simulation** — sigmoidal CV model derived from wall thickness
- **AHA 17-Segment Model** — per-segment quantitative metrics
- **Cohort Analysis Dashboard** — multi-patient statistics, correlations, outlier detection
- **Orientation Widget** — anatomical landmark alignment with a 3D human bust

---

## Screenshots

> *(Add screenshots here)*

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

## Building the `.exe` yourself

```bash
pip install pyinstaller
scripts\build_exe.bat
```

The output will be placed in `dist\LVExplorer\LVExplorer.exe`.

See [docs/installation.md](docs/installation.md) for full build instructions.

---

## Project Structure

```
lv-explorer/
├── lv_explorer/          # Main Python package
│   ├── core/             # Data loading, visualization, orientation
│   ├── metrics/          # AHA metrics catalog
│   └── ui/               # Qt windows and dialogs
├── assets/               # 3D reference model (human bust .obj/.mtl)
├── data/                 # Patient data — NOT committed (see .gitignore)
├── cohort_results/       # Analysis output — NOT committed
├── docs/                 # Documentation
├── scripts/              # Build utilities
├── .github/workflows/    # CI/CD — automated .exe build
├── main.py               # Entry point
├── requirements.txt      # Python dependencies
└── pyproject.toml        # Package metadata
```

---

## Data Format

Patient folders are expected to follow the Carto/ADAS-3D export convention:

```
PRAD18/
├── LV ENDO (CT).vtk
├── LV EPI DIST MAP (CT).vtk
├── SCAR (LE).vtk
├── LV WT 1mm (CT).vtk
└── ...
```

Patient data is **excluded from version control** for privacy reasons. See [docs/installation.md](docs/installation.md) for the expected folder layout.

---

## Roadmap

- [ ] Exportable PDF reports per patient
- [ ] DICOM import support
- [ ] Probabilistic scar estimation without LGE-MRI
- [ ] REST API for remote batch processing

---

## Citation

If you use LV Explorer in your research, please cite:

```
@software{lv_explorer_2026,
  author = {Antoine [SURNAME]},
  title  = {LV Explorer — Ventricular Topography Analysis},
  year   = {2026},
  url    = {https://github.com/antoinealxandre/lv-explorer}
}
```

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

> ⚠️ This software is intended for **research purposes only** and is **not a certified medical device**.
