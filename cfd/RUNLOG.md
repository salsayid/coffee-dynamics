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

## Sat evening: Taylor sweep results, compute problems, and a change of approach

### Taylor sweep (done)

| Re | Growth rate, AcuSolve | Growth rate, linear theory | Saturated strength u_r/U |
|---|---|---|---|
| 60 | -0.035 | -0.083 | decays (fit skewed by noise floor) |
| 70 | -0.0165 | -0.0175 | decays |
| 76 | +0.0118 | +0.0140 | 0.033 |
| 80 | +0.0303 | +0.0325 | 0.052 |
| 150 | +0.2011 | +0.2037 | 0.123 |
| 300 | +0.305 | +0.312 | stopped early |
| 600 | +0.367 | +0.374 | stopped early |

Extrapolating strength squared to zero from Re = 76 and 80 gives **Re_c = 73.25 vs 73.22 from theory**.

### Why runs were slow, and the fix

- AcuSolve runs pin their OpenMP threads to the first cores by default, so concurrent jobs piled onto the same 4-8 cores while total CPU stayed at 11-15%. Explicit, non-overlapping `KMP_AFFINITY` core lists fixed it (axisymmetric startup went from 13.5 s to 2.2 s per step).
- **Core limit:** at times more than 16 solver threads were allocated, but measured machine CPU stayed between 11% and 34% (well under 16 cores' worth, i.e. 32 logical processors). From 18:20 on, allocation is held at 16 threads (8 + 4 + 4).
- The level-set auto strategy ran 20 mostly serial redistancing passes per step (31 s/step); capping it at 2-4 passes gives 5.7 s/step.

### Change of approach for the transient runs

The full 3D startup on the coarse 505k-cell mesh cost 66 s per time step (about 37 hours for 10 s). Because the torus head is axisymmetric apart from three thin spokes, the startup (Level 2b) and free surface (Level 3) now use a 5 degree wedge of the real tapered mug with the torus as a rotating wall (17k and 33k cells). The 3D steady rotating-frame run on the 1.45M-cell mesh continues and serves as the check on the wedge model's end state. Spokes and non-axisymmetric turbulence structure are not represented in the wedge runs.

### Early 3D steady result (not yet converged)

After 65 iterations the frother torque has settled at about -3.0 mN·m (3.8 W at 12,000 RPM, plausible for a handheld frother), but the mug walls only take +0.9 mN·m back. The coffee is still spinning up inside the steady solver; the run continues.

### META

Batch sessions (`meta_post64 -b -s`) load AcuSolve results via EnSight. Solid faces render white in batch mode, but cut planes show contours properly, so all volume renders use cut planes. `cfd/scripts/derive_fields.py` adds radial velocity, swirl velocity and speed to the EnSight case for META.
