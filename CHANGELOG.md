# Changelog

All notable changes to LV Explorer will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned

- PDF export per patient
- DICOM import support
- Probabilistic scar estimation (CT-only workflow)

---

## [0.1.0-alpha] — 2026-02-25

### Added

- Initial alpha release of LV Explorer
- 3D interactive visualization of LV endo/epicardium meshes (PyVista / VTK)
- Wall thickness maps (1–5 mm isocontours)
- Scar segmentation overlay (dense scar, border zone, transmurality from LGE-MRI)
- Conduction velocity simulation via sigmoidal CV model
- AHA 17-segment model with per-segment quantitative metrics
- Multi-view adaptive grid layout (up to 9 simultaneous views)
- Cohort analysis dashboard (descriptive stats, correlations, outlier detection)
- Human bust orientation widget for anatomical alignment
- CLI entry point (`python main.py --data-path <folder>`)
- Standalone Windows `.exe` build via PyInstaller

### Known Limitations

- Windows 10/11 only for the `.exe` distribution
- Requires LGE-MRI data for scar metrics (estimated fallback available)
- No DICOM import — VTK mesh format only
