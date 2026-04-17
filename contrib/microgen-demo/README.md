
# microgen-demo

**Open source contribution — PyQt/pydantic GUI for microstructure generation**  
Based on [Mérope](https://github.com/MarcJos/Merope) and [Combs](https://www.sna-mc2013.org/) APIs.

---

## Overview of the GUI

![GUI screenshot](1.png)

*Main configuration window (PyQt5).*

---

## Command-line demonstration

Example output of `make demo`:

![make demo - part 1](2.png)
![make demo - part 2](3.png)
![make demo - part 3](4.png)

*Full output: Merope script, Combs script, and JSON session.*

---

## What this contribution brings

- **pydantic v2 data model** → validation, JSON serialization, round‑trip.
- **PyQt5 / PySide6 GUI** → visual microstructure configuration.
- **Script generation** for Merope and Combs (real API).
- **Unit tests** (pytest) and code coverage.

---

## Quick install

```bash
git clone https://github.com/CodingYayaToure/microgen-demo
cd microgen-demo
make install          # pydantic + pytest
make install-qt5      # or : make install-qt6
```

---

## Usage

```bash
make demo             # CLI — prints JSON + generated scripts
make gui              # PyQt5/PySide6 GUI
make test             # 30+ unit tests (pytest)
make coverage         # HTML coverage report
make example          # Generate merope_beton.py + combs_beton.py
make example-laguerre # Generate polycrystal Laguerre script
```

---

## Model architecture

```
MicrostructureConfig          ← root, single source of truth
├── name, description, backend, seed
├── SimulationBox             ← Lx, Ly, Lz, shape (GlobalShape.Tore)
│   └── as_merope_L()         → [Lx, Ly, Lz] for algo.setLength(L)
├── inclusions : list[AnyInclusion]   ← discriminated union by "shape"
│   ├── SphereFamily          ← RSA/WP/AlgoBool sphere packing
│   │   ├── volume_fraction   → algo.setRadiusGenerator([[R, phi]], [phase])
│   │   └── n_spheres         → algo.fillMaxRSA_3D(n)
│   ├── MultiSphereFamily     ← multi‑radius distribution
│   │   └── to_merope_args()  → (desiredRPhi, tabPhases)
│   ├── LaguerreTessFamily    ← polycrystal via voro++
│   │   └── LaguerreTess_3D(L, n_grains)
│   └── EllipsoidFamily       ← PolyInclusions_3D
├── VoxelGrid                 ← CartesianGrid_3D
│   ├── nx, ny, nz            → powers of 2 (FFT tmfft/AMITEX)
│   └── as_merope_list()      → [nx, ny, nz]
├── volume_fraction_estimate() → real‑time validation in GUI
├── voxel_size_um()           → physical diagnostic
├── to_merope_script()        → full Merope Python script
├── to_combs_script()         → full Combs Python script
├── to_json() / from_json()   → .microgen session (SALOME save)
```

---

## Using the real Merope API

Generated scripts call the **documented** Python API of [Mérope](https://github.com/MarcJos/Merope).  
Example snippet from `to_merope_script()`:

```python
import merope
import merope.vox

L = [0.001, 0.001, 0.001]            # SimulationBox.as_merope_L()
shape = "Tore"                        # GlobalShape.TORE

algo = merope.RSA_3D()
algo.setLength(L)
algo.setShape(shape)
algo.setSeed(42)                      # MicrostructureConfig.seed
algo.setRadiusGenerator([[5e-05, 0.15]], [1])   # SphereFamily (vf)
output = algo.proceed()

spheres = merope.SphereInclusions_3D(output)
multi = merope.MultiInclusions_3D()
multi.setInclusions(spheres)

structure = merope.Structure_3D(multi)
vox = merope.vox.Voxellation_3D(structure)
vox.proceed([128, 128, 128])          # VoxelGrid.as_merope_list()
vox.write("microstructure.vtk")
```

> **Note**: This project does **not** require Merope to be installed. The generated scripts are syntactically correct; their execution requires a working Merope installation (not included).

---

## pydantic v2 features used

- `Field(gt=0, discriminator=...)` — declarative validation
- `field_validator` — single‑field validation
- `model_validator(mode="after")` — cross‑field validation
- `Annotated[Union[...], Field(discriminator="shape")]` — discriminated union for polymorphic inclusions
- `model_dump_json()` / `model_validate_json()` — perfect JSON round‑trip

---

## PyQt5 / PySide6 abstraction

Backend compatibility in **2 lines**:

```python
try:
    from PyQt5.QtWidgets import QApplication, ...
    from PyQt5.QtCore import pyqtSignal as Signal
    QT_BACKEND = "PyQt5"
except ImportError:
    from PySide6.QtWidgets import QApplication, ...
    from PySide6.QtCore import Signal
    QT_BACKEND = "PySide6"
```

All GUI code (widgets, layouts) is identical for both backends.

---

## Tests

```
tests/test_models.py  — 30+ tests organised by class
  TestSimulationBox         → volume, as_merope_L, invalid rejections
  TestSphereFamily          → vf/n modes, mutual exclusion, volumes
  TestMultiSphereFamily     → to_merope_args, total vf < 1
  TestVoxelGrid             → power of 2, warning
  TestMicrostructureConfig  → vf, voxel size, JSON roundtrip,
                              Merope scripts (setLength, setRadiusGenerator,
                              MultiInclusions_3D, Voxellation_3D, setSeed...)
  TestPerformance           → 1000 serialisations < 1s
```

---

## Limitations

- Actual execution of generated scripts requires a working **Mérope** installation (not included).
- Full integration into **SALOME** is not provided.
- Combs API has not been tested (no access to the library), but script generation follows its documentation.

---

## References

- [Mérope (MarcJos/Merope)](https://github.com/MarcJos/Merope)  
  M. Josien, *Mérope: A Microstructure Generator for Simulation of Heterogeneous Materials*, 2024, hal04486524
- [Combs](https://www.sna-mc2013.org/) — C. Bourcier et al., SNA+MC 2013
- [SALOME platform](https://www.salome-platform.org)
- [pydantic v2 docs](https://docs.pydantic.dev/latest/)

```

