"""Level 2/3: write an AcuSolve case for the real mug + frother (MRF rotor, SST RANS)."""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import params as P  # noqa: E402

ROTATING_WALLS = ("frother_head", "shaft_rotor", "shaft_stator")
FIXED_WALLS = ("mug_wall", "mug_bottom")


def probes(mesh_dir, case_dir):
    """Nearest nodes to sample points: radial lines at the head plane and mid-depth, plus a vertical line."""
    crd = np.loadtxt(os.path.join(mesh_dir, "mesh.crd"))
    ids, xyz = crd[:, 0].astype(int), crd[:, 1:]
    out = {}
    for name, z in (("radial_head", P.HEAD_Z), ("radial_mid", 0.020), ("radial_top", 0.035)):
        pts = [(r, 0.0, z) for r in np.linspace(0.0005, P.mug_radius(z) - 0.0002, 60)]
        out[name] = pts
    out["vertical_r20"] = [(0.020, 0.0, z) for z in np.linspace(0.0005, P.LIQUID_H - 0.0005, 60)]
    sel = {}
    for name, pts in out.items():
        nn = []
        for p in pts:
            d = np.sum((xyz - np.array(p)) ** 2, axis=1)
            nn.append(int(ids[np.argmin(d)]))
        nn = list(dict.fromkeys(nn))
        np.savetxt(os.path.join(case_dir, f"probe_{name}.nbc"), nn, fmt="%d")
        sel[name] = nn
    return sel


def init_from(src_case, case_dir, mesh_dir):
    """Write nodal initial-condition files from the last output of another case on the same mesh."""
    sys.path.insert(0, HERE)
    import acu
    files = {}
    for var in ("velocity", "pressure", "kinetic_energy", "eddy_frequency", "eddy_viscosity"):
        ids, _, vals = acu.nodal(src_case, var, mesh_dir=mesh_dir)
        vals = vals.reshape(len(ids), -1)
        fn = os.path.join(case_dir, f"init_{var}.nic")
        np.savetxt(fn, np.column_stack([ids, vals]), fmt="%d" + " %.8e" * vals.shape[1])
        files[var] = os.path.basename(fn)
    return files


def write(case_dir, mesh_dir, transient=False, dt=5e-4, final_time=10.0, nsteps=400, out_every=200,
          tol=1e-3, omega=P.OMEGA, title="L2", restart=None, stagger=3):
    os.makedirs(case_dir, exist_ok=True)
    meta = json.load(open(os.path.join(mesh_dir, "mesh_meta.json")))
    rel = os.path.relpath(mesh_dir, case_dir).replace("\\", "/")
    L = []
    w = L.append
    w(f'ANALYSIS {{\n    title = "{title}"\n    type = {"transient" if transient else "steady"}\n}}\n')
    w("EQUATION {\n    flow = navier_stokes\n    turbulence = shear_stress_transport\n}\n")
    if transient:
        w(f"AUTO_SOLUTION_STRATEGY {{\n    max_time_steps = {nsteps}\n    final_time = {final_time}\n"
          f"    initial_time_increment = {dt}\n    auto_time_increment = off\n"
          f"    max_stagger_iterations = {stagger}\n    convergence_tolerance = {tol}\n}}\n")
    else:
        w(f"AUTO_SOLUTION_STRATEGY {{\n    max_time_steps = {nsteps}\n    convergence_tolerance = {tol}\n}}\n")
    w(f'DENSITY_MODEL( "water85" ) {{\n    density = {P.RHO}\n}}\n')
    w(f'VISCOSITY_MODEL( "water85" ) {{\n    viscosity = {P.MU}\n}}\n')
    w('MATERIAL_MODEL( "water85" ) {\n    density_model = "water85"\n    viscosity_model = "water85"\n}\n')
    w(f'REFERENCE_FRAME( "rotor" ) {{\n    centrifugal = on\n    coriolis = on\n'
      f'    rotation_center = {{ 0, 0, 0 }}\n    angular_velocity = {{ 0, 0, {omega} }}\n}}\n')
    w(f'COORDINATE {{\n    coordinates = Read( "{rel}/mesh.crd" )\n}}\n')
    for name, s in meta["element_sets"].items():
        rf = '\n    reference_frame = "rotor"' if s["region"] == "rotor" else ""
        w(f'ELEMENT_SET( "{name}" ) {{\n    elements = Read( "{rel}/{s["file"]}" )\n    shape = {s["shape"]}\n'
          f'    material_model = "water85"{rf}\n}}\n')
    for s in meta["surfaces"]:
        bc = s["bc"]
        tag = f'{bc}__{s["element_set"]}__{"tri" if s["shape"].startswith("three") else "quad"}'
        if bc in ROTATING_WALLS:
            extra = '    type = wall\n    reference_frame = "rotor"\n'
        elif bc in FIXED_WALLS:
            extra = "    type = wall\n"
        elif bc == "top":
            extra = "    type = slip\n"
        else:
            raise ValueError(bc)
        body = (f'    surfaces = Read( "{rel}/{s["file"]}" )\n    shape = {s["shape"]}\n'
                f'    element_set = "{s["element_set"]}"\n')
        w(f'SIMPLE_BOUNDARY_CONDITION( "{tag}" ) {{\n{body}{extra}}}\n')
        if bc in ROTATING_WALLS or bc in FIXED_WALLS:
            w(f'SURFACE_OUTPUT( "{tag}" ) {{\n{body}    integrated_output_frequency = 1\n}}\n')
    nut0 = 10 * P.NU
    k0 = 1e-6
    if restart:
        for var, fn in init_from(restart, case_dir, mesh_dir).items():
            w(f'NODAL_INITIAL_CONDITION( {var} ) {{\n    type = nodal_values\n    nodal_values = Read( "{fn}" )\n}}\n')
    else:
        w("NODAL_INITIAL_CONDITION( velocity ) {\n    default_values = { 0, 0, 0 }\n}\n")
        w("NODAL_INITIAL_CONDITION( pressure ) {\n    default_value = 0\n}\n")
        w(f"NODAL_INITIAL_CONDITION( kinetic_energy ) {{\n    default_value = {k0}\n}}\n")
        w(f"NODAL_INITIAL_CONDITION( eddy_frequency ) {{\n    default_value = {k0 / nut0:.6e}\n}}\n")
        w(f"NODAL_INITIAL_CONDITION( eddy_viscosity ) {{\n    default_value = {nut0:.6e}\n}}\n")
    w(f"NODAL_OUTPUT {{\n    output_frequency = {out_every}\n}}\n")
    w(f"RESTART_OUTPUT {{\n    output_frequency = {out_every}\n}}\n")
    for name in probes(mesh_dir, case_dir):
        w(f'TIME_HISTORY_OUTPUT( "{name}" ) {{\n    nodes = Read( "probe_{name}.nbc" )\n'
          f'    output_frequency = {1 if not transient else max(1, out_every // 20)}\n}}\n')
    w("RUN\n")
    open(os.path.join(case_dir, "case.inp"), "w").write("\n".join(L))
    open(os.path.join(case_dir, "Acusim.cnf"), "w").write("problem = case\ninput_file = case.inp\n")
    json.dump(dict(transient=transient, dt=dt, final_time=final_time, nsteps=nsteps, omega=omega,
                   mesh_dir=os.path.abspath(mesh_dir)), open(os.path.join(case_dir, "case_info.json"), "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("mesh_dir")
    ap.add_argument("--transient", action="store_true")
    ap.add_argument("--dt", type=float, default=5e-4)
    ap.add_argument("--final_time", type=float, default=10.0)
    ap.add_argument("--nsteps", type=int, default=400)
    ap.add_argument("--out_every", type=int, default=200)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--init_from", default=None, help="case dir whose last output seeds the initial fields")
    ap.add_argument("--stagger", type=int, default=3)
    a = ap.parse_args()
    write(a.case_dir, a.mesh_dir, a.transient, a.dt, a.final_time, a.nsteps, a.out_every, a.tol,
          restart=a.init_from, stagger=a.stagger)
    print("wrote", a.case_dir)
