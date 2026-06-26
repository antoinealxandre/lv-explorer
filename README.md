# LV Explorer

**LV Explorer** is a desktop application for interactive 3D visualization and analysis of left ventricular geometry. It is designed for clinical research on ventricular tachycardia (VT) substrate mapping and ablation planning.

---

> ⚠️ **Research use only.** Not intended for clinical decision-making.

---

## Demo

https://private-user-images.githubusercontent.com/188148970/613584460-7838bbb7-6a52-4d15-bd0d-805400ec794f.mp4?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3ODI0NjQ3OTcsIm5iZiI6MTc4MjQ2NDQ5NywicGF0aCI6Ii8xODgxNDg5NzAvNjEzNTg0NDYwLTc4MzhiYmI3LTZhNTItNGQxNS1iZDBkLTgwNTQwMGVjNzk0Zi5tcDQ_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjYwNjI2JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI2MDYyNlQwOTAxMzdaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1iZmI5MzM3YzFlYjFmZmM3YmMyM2JhMjc2ZDc5ZGE0OTUyOTg0MGRhNWEwODAxZjUwM2ZhNDE2OGZmMjg3YzBiJlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZyZXNwb25zZS1jb250ZW50LXR5cGU9dmlkZW8lMkZtcDQifQ.492s62HdEYAjxuSYoJbgSkbdW0Je-_F-hXj9syr2w88

> _3D interactive visualization of left ventricular geometry with different metrics._

---

## Features

- **3D LV Rendering** — Interactive 3D visualisation of the left ventricle from `.vtk` files (via PyVista).
- **Wall Thickness & Topographic Metrics** — Analysis of wall thickness to characterise high-risk areas:
  - _Tissue Classification_: Automatic identification of healthy tissue (> 5 mm), the border zone (1–5 mm) and dense scar tissue (< 1 mm).
  - _Ciaccio’s Ratio_: Detection of abrupt variations in thickness (‘anatomical cliffs’).
  - _Conduction Channels_: Localisation of corridors of viable tissue traversing the scar.
  - _Advanced Metrics_: Calculation of the Laplacian of thickness and local entropy to measure tissue disorganisation.
- **Integrated Score** — Customisable combination of different metrics to predict the location of critical isthmuses.
- **Global Geometry** — Calculation of the sphericity index to assess the overall remodelling of the ventricular cavity.
- **Propagation Simulation** — Rapid simulation (in a matter of seconds) of electrical activation with isochrones displayed every 10 ms.
- **Statistical Dashboard** — Control panel for adjusting clinical thresholds in real time and exporting key statistics.

---

## Data

Patient data is provided by [inHEART](https://www.inheart.fr/) and exported in `.vtk` (CARTO format).

```
data/
└── <PATIENT_ID>/
    ├── LV ENDO DIST MAP THRES (CT).vtk     # LV epicardium to endocardium distance map thresholded
    ├── LV EPI DIST MAP (CT).vtk            # LV epicardium to endocardium distance map
    └── ...
```

> ⚠️ Patient data is **not included** in this repository (de-identified, subject to data sharing agreements with inHEART).

---

## Installation

**Requirements:** Python ≥ 3.10 · Windows 10/11 · macOS 12+ · Ubuntu 22.04+

### From source

```bash
# 1. Clone the repository
git clone https://github.com/antoinealxandre/lv-explorer.git
cd lv-explorer

# 2. Create a virtual environment
python -m venv .venv
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run
python main.py
```

To open a specific patient folder directly:

```bash
python main.py --data-path "data/PATIENT_ID"
```

### Standalone executable (Windows)

Download the latest `.zip` from the [Releases](https://github.com/antoinealxandre/lv-explorer/releases) page and run `LVExplorer.exe` — no Python required.

---

## Project Structure

```
lv-explorer/
├── lv_explorer/          # Core application package
│   ├── app.py            # Application entry point
│   ├── viewer.py         # 3D rendering and interaction
│   ├── loader.py         # VTK file loading and parsing
│   ├── analysis.py       # Voltage, scar, and geometry analysis
│   └── ui/               # UI components and panels
├── assets/               # Icons and static assets
├── config/               # Default configuration files
├── tests/                # Unit and integration tests
├── docs/                 # Additional documentation
├── main.py               # CLI entry point
├── pyproject.toml        # Project metadata and build config
├── requirements.txt      # Python dependencies
└── lv_explorer.spec      # PyInstaller build spec
```

---

## Dependencies

| Package                                                | Purpose                          |
| ------------------------------------------------------ | -------------------------------- |
| [PyVista](https://pyvista.org/)                        | 3D mesh rendering                |
| [VTK](https://vtk.org/)                                | File I/O and geometry processing |
| [NumPy](https://numpy.org/)                            | Numerical computations           |
| [PyQt5](https://riverbankcomputing.com/software/pyqt/) | Desktop UI framework             |

---

## Acknowledgements

This project was developed during an engineering internship in the **Rhythmology Department** at **Hôpital Louis Pradel, Hospices Civils de Lyon (HCL)**, under the supervision of **Dr. Geoffroy Ditac**.

Patient data and clinical expertise were provided by [inHEART](https://www.inheart.fr/).

If you use this software in your research, please cite it as:

> Antoine Alexandre. _LV Explorer: Interactive 3D Visualization of Left Ventricular Geometry for VT Substrate Mapping._ Engineering internship project, Service de Rythmologie, Hôpital Louis Pradel — HCL, 2026. https://github.com/antoinealxandre/lv-explorer
