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
