# Run log

Newest entries at the bottom. Times are local.

## Phase 0: environment and verification (Fri Oct 9, 18:05-19:10)

| Check | Result |
|---|---|
| Branch | `cfd-pipeline` created from `main` |
| Python | numpy, scipy, matplotlib, pandas, pyvista, meshio, pypdf installed (Python 3.14) |
| ffmpeg | installed via winget |
| Docker | Docker Desktop started; `opencfd/openfoam-default:2406` pulled (fallback only) |
| META allowed folders | repo root added to `%LOCALAPPDATA%\meta-mcp\config.json` |
| AcuSolve 2019 | channel example runs on 2 and 16 cores; license checked out from HyperWorks units |
| AcuSolve features | 2019 Command Reference confirms SST, `REFERENCE_FRAME` (rotating frame), `MESH_MOTION` + `INTERFACE_SURFACE` (sliding mesh), `PERIODIC_BOUNDARY_CONDITION`, `multi_field = levelset` (free surface) |
| ANSA 24.1.6 batch | works via `cfd/scripts/run_ansa.ps1` (ANSA_BATCH license OK) |
| META 24.1.6 batch | works via `cfd/scripts/run_meta.ps1` (META_POST_BATCH license OK); reads AcuSolve results after `acuTrans -to ensight` and writes PNGs |

### Findings that change the plan

- **ANSA launcher bug:** `ansa64.bat` drops any argument containing `|`, so `-execscript "file.py|main()"` silently starts an idle ANSA with no arguments. `run_ansa.ps1` sets ANSA's environment itself and calls `ansa_win64.exe` directly.
- **No AcuSolve exporter in ANSA 24.** ANSA writes Fluent, CGNS, STAR-CD, Universal, Nastran and OpenFOAM. AcuSolve's `acuImport` reads Fluent, STAR-CD and I-DEAS universal. Route for Level 2/3: ANSA → Fluent `.msh` → `acuImport -from fluent`.
- **Level 1 mesh:** the idealized annulus is a perfectly structured hexahedral mesh, so it is generated directly in Python and written in AcuSolve's native format. This gives an exact, controllable mesh for the validation and convergence study. ANSA is used where its meshing matters (Levels 2 and 3).
- META reads results via EnSight Gold (`acuTrans -out -to ensight`). AcuSolve can also write `.h3d` for HyperView.

## Lost time: machine slept (Fri 18:38 to Sat 16:29)

Windows idle-slept about 30 minutes after setup started, which froze everything for roughly 22 hours. `cfd/scripts/keep_awake.ps1` now holds a "system required" request (the mechanism video players use) for 72 hours while the runs go. It changes no power settings and stops when the process ends.

## Phase 1a: laminar Couette validation (done)

Idealized annulus (R1 = 10.5 mm, R2 = 35 mm), 5 degree wedge with 8 cells around, slip top and bottom, Re = 50 with real water properties. Three radial meshes:

| Radial cells | Velocity RMS error (of wall speed) | Torque error |
|---|---|---|
| 40 | 3.1e-5 | +0.0166% |
| 80 | 8.5e-6 | +0.0047% |
| 160 | 2.7e-6 | +0.0017% |

Observed order 1.96 (second-order scheme), fine-grid GCI 0.0013%. A first attempt with only 1 cell around the wedge plateaued at a 0.05% torque offset because the walls were flat facets instead of arcs.

## Phase 1b: Taylor vortices (in progress)

- Independent linear-stability solver (`taylor_stability.py`, Chebyshev collocation) reproduces the literature value for radius ratio 0.5 (Re_c = 68.19). For this mug's ratio of 0.3: **Re_c = 73.22, axial wavelength 1.95 gaps**.
- Domain height set to exactly two critical wavelengths.
- Re = 150: 4 vortex cells as expected; measured growth rate 0.2011 vs 0.2037 from linear theory (1.3% off).
- Concurrent small runs were about 3x slower per step than a solo run. Disabling OpenMP pinning did not fix it. Re = 300 and 600 were stopped early (partial) to free cores for Level 2.

## Phase 2 meshing (pipeline working)

Final pipeline: gmsh/OpenCASCADE CAD and surface mesh (named groups) → ANSA batch mesh (orient, 6 prism layers with collapse at problem spots, Tetra Rapid volume) → Nastran → `nas_to_acusolve.py` (rotor/stator split by centroid, boundary faces matched to owning elements) → AcuSolve native files. Problems hit on the way:

- gmsh/OCC boolean of torus + spokes produces inverted solids at most spoke angles; fixed by subtracting all frother parts from the fluid in one multi-tool cut.
- ANSA's STEP import leaves a doubly periodic torus face unmeshed and reports flat faces 6-8% too small; switched from CAD import to importing gmsh's surface mesh.
- ANSA layers grow from the grey side; volumes are detected and oriented before layering.
- ANSA puts all layer prisms in one property and would not move them to the rotor property; the converter assigns regions by element centroid instead.

| Mesh | Frother surface size | Cells | Prisms |
|---|---|---|---|
| m0 | 0.35 mm | 505k | 106k |
| m1 | 0.18 mm | 1.45M | 396k |

A 1-thread check run of m1 reads cleanly; about 3 minutes per steady iteration on 1 thread.
