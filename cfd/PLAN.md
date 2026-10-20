# Coffee Dynamics: CFD Pipeline Plan

Unattended weekend run. ANSA builds and meshes the models, AcuSolve 2019 solves them (OpenFOAM in Docker is the automatic fallback), META post-processes, and Python handles validation plots and the README.

Progress is logged live in `cfd/RUNLOG.md`. Every run can be resumed from that log.

---

## 0. Decisions locked in

| Topic | Decision |
|---|---|
| Licenses / machine | ANSA, META and Altair licenses are OK to use on this machine |
| Solver | AcuSolve 2019 (`C:\Program Files\Altair\2019\hwcfdsolvers\acusolve\win64`) |
| Fallback | If AcuSolve fails at a step (license, cores, missing feature, unreadable results), that step switches to OpenFOAM in Docker and the reason is logged |
| Pre / post | ANSA 24.1.6 builds and meshes. META 24.1.6 post-processes. HyperMesh/AcuConsole only if AcuSolve input truly needs them |
| Cores | Up to 16 |
| Turbulence | RANS only (SST k-omega). Laminar for the low-Re Level 1 cases |
| Fluid | Water at 85°C |
| Scope | Levels 1, 2, 3, plus Level 4 prep |
| Frother | Coiled-spring whisk, dimensions **estimated** (section 2) |
| Frother position | Centered, head mid-plane 10 mm above the mug bottom |
| Mug | Linear taper: 35 mm inner radius at the rim, about 30 mm at the bottom (**estimated**) |
| Git | New branch `cfd-pipeline`, commit as I go, no push. Only scripts, case setups, figures and videos are tracked. Raw results are gitignored |
| Permissions granted | Start Docker and pull an OpenFOAM image; pip installs; install ffmpeg; add the repo to META's allowed folders; run ANSA scripts without per-step confirmation |
| Deliverables | META 3D renders and animations, validation plots, new README sections, side-by-side video template |

---

## 1. Physical inputs

| Quantity | Value | Source |
|---|---|---|
| Frother speed | 12,000 RPM = 1256.6 rad/s | README |
| Tip speed (R1 = 10.5 mm) | 13.19 m/s | computed |
| Water density, 85°C | 968.6 kg/m³ | standard tables |
| Kinematic viscosity | 3.3e-7 m²/s | README |
| Dynamic viscosity | 3.20e-4 Pa·s | ρν |
| Surface tension, water/air at 85°C | 0.0618 N/m | standard tables |
| Air density / viscosity (hot, humid) | 1.0 kg/m³ / 2.1e-5 Pa·s | standard tables |
| Gravity | 9.81 m/s², −z | |
| Liquid depth | 40 mm | README |

---

## 2. Geometry (built in ANSA from Python scripts)

All coordinates are in mm. The z axis is the mug axis, with z = 0 at the inside of the mug bottom.

### 2.1 Mug inner wall (fluid boundary)
- Truncated cone: r(z) = 30 + z/21, so r = 30 at z = 0, 31.90 at z = 40 (liquid surface) and 35 at z = 105 (rim).
- Flat bottom with no fillet.
- Only the fluid region is meshed. The mug wall itself is not modeled (no heat transfer).

### 2.2 Frother (estimated; every value is a parameter in `cfd/geometry/params.py`)
- **Coil ring (torus-shaped helical spring):** ring centerline radius 9.0 mm, coil outer diameter 3.0 mm (outer extent 10.5 mm, matching R1), wire diameter 0.4 mm, pitch 0.8 mm, about 70 turns.
- **Inner support ring:** wire diameter 0.8 mm at radius 9.0 mm, running through the coil.
- **Spokes:** 3 wires, 0.8 mm diameter, from the shaft to the support ring.
- **Shaft:** 3.5 mm diameter, vertical, from the head up through the liquid surface (and to z = 105 for Level 3).
- **Head mid-plane:** z = 10 mm.

### 2.3 Two frother variants
- **Variant A, envelope torus:** a solid torus (major radius 9.0, minor radius 1.5) plus spokes and shaft. Cheap and robust. Used for transient and free-surface runs.
- **Variant B, resolved coil:** the actual helical wire. High fidelity, steady runs only. Used to check how far off Variant A is.

### 2.4 Rotating region
- Cylinder r ≤ 13 mm, 6 ≤ z ≤ 14 mm, enclosing the head. Its faces are the interface between the rotating and stationary mesh.
- Shaft above the region: a rotating wall boundary (wall velocity = Ω × r).

### 2.5 Level 1 idealized domain
- Straight annulus: R1 = 10.5, R2 = 35 (matches the README derivation), gap d = 24.5.
- 5° wedge with periodic side faces. Axial length L = 4d = 98 mm, periodic top and bottom (room for two Taylor-vortex pairs).
- Structured hexahedral mesh. Three levels for mesh convergence: 40×160, 80×320, 160×640 (radial × axial).

---

## 3. Meshing (ANSA, headless)

- **Surface:** triangles. 0.10 mm on Variant B wire, 0.3 mm on the torus and spokes, 0.5 mm on the shaft, 1.0 mm on the mug walls, 0.4 mm band around the liquid surface (Level 3).
- **Boundary layers:** prism layers on all walls. First-layer height is chosen for y+ of about 30-60 with wall functions (about 20 µm on the frother, 60 µm on the mug). 5-8 layers, growth ratio 1.2.
- **Volume:** tetrahedra with size grading, with a refined zone over the rotating region.
- **Estimated sizes:** Level 2 Variant A 2-3 M cells, Variant B 6-10 M, Level 3 Variant A 4-6 M.
- **Quality gates** (scripted, logged): skewness < 0.9, no negative volumes, prism layers unbroken on the frother. Failing meshes are fixed and re-meshed before solving.
- **Boundary faces** become separate ANSA parts/properties: `frother_head`, `shaft`, `mug_wall`, `mug_bottom`, `top` (rigid lid or atmosphere), `mrf_interface`, `periodic_*`.
- **Export path, in order:** (1) ANSA's own AcuSolve output if this version has one; (2) a Nastran export plus `cfd/scripts/nas_to_acusolve.py`, which writes AcuSolve coordinate, connectivity and boundary-surface files; (3) ANSA's OpenFOAM output for the fallback solver.

---

## 4. Phases

### Phase 0: Environment and verification (about 2-3 h). Nothing proceeds until this passes.
1. Create branch `cfd-pipeline`. Add the `cfd/` folders (`geometry/ mesh/ cases/ scripts/ meta/ post/`) and `.gitignore` entries for meshes, results, `ACUSIM.DIR`, `processor*` and logs.
2. `pip install numpy scipy matplotlib pandas pyvista meshio` and add `requirements.txt`.
3. Install ffmpeg (`winget install Gyan.FFmpeg`).
4. Start Docker Desktop and pull `opencfd/openfoam-default` (fallback only).
5. Add the repo to META's allowed folders (`setup_client.ps1 -AllowedRoot`).
6. **AcuSolve check:** run a tiny laminar channel on 2 cores, then on 16 cores, to confirm the license checkout and core count. Read the 2019 docs to confirm SST, rotating frame, sliding mesh, periodic boundaries and level-set free surface. Whatever isn't available is routed to OpenFOAM.
7. **ANSA headless check:** run a script with `ansa64.bat -nogui` that builds a cylinder, meshes it and exports it. Confirm the export path from section 3.
8. **META check:** load the channel result (directly, or translated with `acuTrans` to EnSight), make a contour, save a PNG without the GUI, and also through `meta.ps1`. If META can't read it, fall back to PyVista renders.
9. Log the outcome of every check in `RUNLOG.md`.

### Phase 1, Level 1: Idealized check against the analytical profile (about 4-6 h)
- **1a, laminar steady at Re = 50** (Ω chosen from Re at the real viscosity). Compare u_θ(r) with u = Ar + B/r on all three meshes. Pass: L2 error below 0.5% on the medium mesh, plus a grid-convergence index. Also compare inner-wall torque with the exact value T = 4πμL Ω1 R1² R2² / (R2² − R1²).
- **1b, Taylor-vortex sweep:** transient laminar at Re = 50, 60, 70, 80, 90, 100, 150, 300, 1000. Measure the radial-velocity amplitude, then extrapolate amplitude² versus Re to zero to estimate the critical Re. Compare with linear-stability literature for η = 0.3. Render the vortex cells in META.
- **1c, the real Re (981k) with RANS** in the same geometry. Compare the turbulent mean profile with the laminar one (flat core, thin wall layers) and compare torque with published high-Re Taylor-Couette torque scaling. This shows quantitatively how far reality is from the laminar ground truth.
- **Figures:** `validation_profile.png`, `mesh_convergence.png`, `taylor_onset.png`, `laminar_vs_turbulent_profile.png`. META: Taylor-cell velocity slice, high-Re slice.

### Phase 2, Level 2: Real mug and frother, rigid flat liquid surface (about 6 h hands-on plus 20-30 h compute)
- **2a, steady rotating-frame runs (frozen rotor), SST:** Variant A on three meshes for convergence (torque, mean fluid spin); Variant B on one mesh. Compare A and B: if they differ by more than 15%, report that A is a coarse stand-in.
- **2b, spin-up from rest, Variant A:** time-accurate with a rotating frame, time step about 5e-4 s, t = 0 → 10 s (about 20k steps). If the measured cost per step projects past the budget, run to 5 s and fit an exponential to show where it's heading.
- **2c, rotating-mesh spot check:** about 0.1 s of a true sliding-mesh run, restarted from the 2b end state, to measure the frozen-rotor error in torque and fluid spin.
- **Saved over time (2b):** frother torque T(t) and power (the resistance felt in your hand); volume-averaged fluid angular velocity; relative velocity at the head; Re_rel(t) = (Ω1 − ω_fluid) R1 d / ν; volume-averaged turbulent kinetic energy and turbulent-to-molecular viscosity ratio. Field snapshots every 0.1 s.
- **Questions answered:** Does spin-up actually lower the frother-to-fluid Reynolds number, by how much, and on what time scale? Does torque drop the way your hand feels it? This tests the README's argument directly.
- **Figures:** `torque_vs_time.png`, `re_rel_vs_time.png`, `fluid_spinup.png`, `variantA_vs_B.png`. META: velocity on vertical and horizontal slices, streamlines, vortex iso-surfaces (Q-criterion, computed in PyVista and written back for META if the solver doesn't output it), spin-up animation.

### Phase 3, Level 3: Free surface and air entrainment (about 4 h hands-on plus 12-20 h compute)
- Domain up to the rim (z = 105) with air above the liquid at z = 40 and an open atmosphere boundary on top. Variant A head, rotating frame, SST.
- Interface tracking: AcuSolve level set if Phase 0 confirmed it, otherwise OpenFOAM `interFoam` with a rotating frame.
- t = 0 → 3 s (extended if the surface hasn't settled).
- **Measured:** dip depth over time; the final surface shape compared with a Rankine-vortex estimate from the Level 2 swirl; whether the dip reaches the head (air reaching the frother is where entrainment starts); air volume below z = 40 (entrainment measure).
- **Not modeled, stated plainly:** foam and its properties. Only the bubble-entrainment stage is captured.
- **Figures:** `surface_dip_vs_time.png`, `vortex_shape_vs_rankine.png`. META: free-surface render with velocity colors, cut-away animation.

### Phase 4, Level 4 prep: Real-coffee comparison (about 2 h)
- `cfd/scripts/video_compare.py`: an ffmpeg template that puts the real top-down slow-motion video next to the simulation's top surface, on a synced timeline with matching color scale and labels.
- `cfd/scripts/extract_surface_velocity.py`: exports the simulated surface velocity (z = 40 plane) as a grid, ready to compare with PIVlab output.
- `cfd/LEVEL4_INSTRUCTIONS.md`: filming setup (tracer, lighting, 240 fps, phone mount, scale reference), PIVlab settings, and how to run the comparison scripts.

### Phase 5, Post-processing pass in META (about 4-6 h, overlapping with compute)
- Saved META session files in `cfd/meta/*.ses` so every image can be reproduced.
- Consistent style: inferno colormap (matching the existing GIF), fixed legend ranges per level, the same camera angles across cases.
- Videos as MP4 (ffmpeg) plus small GIFs for the README in `figures/`.

### Phase 6, Write-up (about 2-3 h)
- New README sections in your voice: "Taking it to CFD", "Level 1: does the solver agree with my math?", "Taylor vortices", "My actual mug", "Spin-up: does the Reynolds number really drop?", "The vortex and the air", "What the simulation can't do (yet)", "Next: filming my coffee".
- Fix the earlier README issues: the η row label and value, the unscaled viscosity in the existing animation (disclose it), small number rounding. Add the mug taper to the mug table.
- Commit at the end of each phase with clear messages. No push.

---

## 5. Schedule (16 cores, weekend)

| Block | Hands-on | Compute |
|---|---|---|
| Phase 0 | 2-3 h | – |
| Phase 1 | 3 h | 2-4 h |
| Phase 2 geometry and mesh | 4-6 h | – |
| Phase 2 runs | 2 h | 20-30 h (Level 3 geometry and meshing happens in parallel) |
| Phase 3 runs | 2 h | 12-20 h |
| Phases 4-6 | 8-10 h | – |

Compute is the critical path. If projected costs blow past the budget, cuts happen in this order: Variant B mesh levels, then 2b end time (10 s → 5 s), then 3b end time, then the 2c length. Every cut is logged.

---

## 6. Stop / fallback rules

- AcuSolve step fails → retry once with a fix → switch that step to OpenFOAM in Docker → log the reason.
- Solution diverges → halve the time step or lower relaxation, start from a laminar or coarser solution, check mesh quality. After 3 failed attempts, mark it blocked and move on to independent work.
- META can't read a result → translate it (EnSight) → render with PyVista instead.
- License server unreachable → wait and retry every 15 min for 2 h, keep working on local-only tasks, log it.
- Never delete your existing files, force-push, or change git config.

---

## 7. Before you leave (needed for unattended running)

1. Turn off sleep and hibernate while plugged in, and pause Windows Update for the weekend. A reboot would kill the runs.
2. Stay logged in (locking is fine) and stay on the company network or VPN, since all licenses check out from network servers.
3. Close any ANSA, META or HyperWorks windows you have open so the license seats are free.
4. Optional: real frother measurements (coil diameter, wire thickness, head height) and the real inner bottom radius of the mug. Each one is a single-line change in `params.py`.
5. Long simulations run as background processes that keep going even if this chat stops. If I stop before finishing, tell me "continue from RUNLOG" and I'll pick up where the log ends.
