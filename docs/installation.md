# Installation Guide

## Table of Contents

1. [Run from source](#run-from-source)
2. [Build the Windows .exe](#build-the-windows-exe)
3. [Data folder layout](#data-folder-layout)
4. [Troubleshooting](#troubleshooting)

---

## Run from source

### Requirements

| Dependency | Version |
|------------|---------|
| Python | ≥ 3.8 |
| numpy | ≥ 1.21 |
| scipy | ≥ 1.7 |
| pyvista | ≥ 0.43 |
| pyvistaqt | ≥ 0.11 |
| qtpy | ≥ 2.0 |
| PyQt5 | ≥ 5.15 |
| matplotlib | ≥ 3.5 |

### Steps

```bash
# Clone
git clone https://github.com/antoinealxandre/lv-explorer.git
cd lv-explorer

# Virtual environment (recommended)
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux / macOS

# Install
pip install -r requirements.txt

# Run
python main.py
python main.py --data-path "data/extract/PRAD18"
```

---

## Build the Windows .exe

### Prerequisites

```bash
pip install pyinstaller pyinstaller-hooks-contrib
```

> All other runtime dependencies must also be installed in the same environment.

### Build

```bat
scripts\build_exe.bat
```

Or directly:

```bash
pyinstaller lv_explorer.spec --noconfirm
```

### Output

```
dist\
└── LVExplorer\
    ├── LVExplorer.exe    ← launch this
    ├── assets\           ← 3D reference model (bundled automatically)
    └── ...               ← VTK / Qt DLLs
```

The script also produces `dist\LVExplorer_v0.1.0_alpha_win64.zip` ready for distribution.

### Adding an application icon

1. Convert a PNG to ICO: [https://convertio.co/png-ico/](https://convertio.co/png-ico/)
2. Place the result at `assets/icon.ico`
3. Uncomment the `icon=` line in `lv_explorer.spec`

---

## Data folder layout

LV Explorer expects one folder per patient with VTK mesh files exported from Carto / ADAS-3D:

```
data/
└── extract/
    └── PRAD18/               ← patient ID
        ├── LV ENDO (CT).vtk
        ├── LV EPI DIST MAP (CT).vtk
        ├── LV WT 1mm (CT).vtk
        ├── LV WT 2mm (CT).vtk
        ├── LV WT 3mm (CT).vtk
        ├── LV WT 4mm (CT).vtk
        ├── LV WT 5mm (CT).vtk
        ├── SCAR (LE).vtk              ← optional (LGE-MRI)
        ├── SCAR TRANSMURALITY 50 (LE).vtk
        ├── DENSE SCAR (LE).vtk
        └── ...
```

> **Patient data is excluded from version control.** Never commit the `data/` folder.

---

## Troubleshooting

### `ImportError: DLL load failed` (Windows)

Install the latest [Microsoft Visual C++ Redistributable](https://aka.ms/vs/17/release/vc_redist.x64.exe).

### Application opens then immediately closes

Run from a terminal to see error messages:

```bat
dist\LVExplorer\LVExplorer.exe
```

### `ModuleNotFoundError: No module named 'vtkmodules'`

Ensure PyVista and VTK are installed in the build environment:

```bash
pip install pyvista pyvistaqt pyinstaller-hooks-contrib
pyinstaller lv_explorer.spec --noconfirm
```

### Human bust not displayed

The `assets/human_bust.obj` file should be automatically bundled. If missing, the app falls back to a geometric approximation silently.
