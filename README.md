# LV Explorer

**LV Explorer** is a desktop application for interactive 3D visualization and analysis of left ventricular geometry. It is designed for clinical research on ventricular tachycardia (VT) substrate mapping and ablation planning.

Patient data is sourced from [inHEART](https://www.inheart.fr/) and exported in `.vtk` format from CARTO electroanatomical mapping systems.

---

> ⚠️ **Research use only.** Not intended for clinical decision-making.

---

## Demo


![LV Explorer Demo]([https://github.com/user-attachments/assets/7838bbb7-6a52-4d15-bd0d-805400ec794f](https://private-user-images.githubusercontent.com/188148970/613584460-7838bbb7-6a52-4d15-bd0d-805400ec794f.mp4?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3ODI0NjQ3OTcsIm5iZiI6MTc4MjQ2NDQ5NywicGF0aCI6Ii8xODgxNDg5NzAvNjEzNTg0NDYwLTc4MzhiYmI3LTZhNTItNGQxNS1iZDBkLTgwNTQwMGVjNzk0Zi5tcDQ_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjYwNjI2JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI2MDYyNlQwOTAxMzdaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1iZmI5MzM3YzFlYjFmZmM3YmMyM2JhMjc2ZDc5ZGE0OTUyOTg0MGRhNWEwODAxZjUwM2ZhNDE2OGZmMjg3YzBiJlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZyZXNwb25zZS1jb250ZW50LXR5cGU9dmlkZW8lMkZtcDQifQ.492s62HdEYAjxuSYoJbgSkbdW0Je-_F-hXj9syr2w88)

> _3D interactive visualization of left ventricular geometry with scar overlay and bipolar voltage mapping._

---

## Features

- **3D LV Rendering** — Interactive visualization of left ventricular meshes from `.vtk` files
- **Voltage Mapping** — Bipolar and unipolar electrogram overlays on the endocardial surface
- **Scar Segmentation** — Automatic detection and display of dense scar and border zone regions
- **Wall Motion Analysis** — Regional wall thickness and motion quantification
- **Multi-Patient Navigation** — Load and switch between patient datasets seamlessly
- **Export** — Save screenshots and annotated views for clinical reports

---

## Data

Patient data is provided by [inHEART](https://www.inheart.fr/) and exported from CARTO electroanatomical mapping systems in `.vtk` format.

```
data/
└── <PATIENT_ID>/
    ├── lv_mesh.vtk          # Left ventricular surface mesh
    ├── voltage_map.vtk      # Bipolar voltage map
    └── scar_overlay.vtk     # Scar segmentation
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

## Contributing

Issues and pull requests are welcome. Please open an issue first to discuss significant changes.

---

## Acknowledgements

Patient data and clinical expertise provided by [inHEART](https://www.inheart.fr/).
