"""Level 1: idealized Taylor-Couette annulus, structured hex wedge, written in AcuSolve native format.

Inner cylinder R1 rotates, outer cylinder R2 is fixed. Wedge sides use AcuSolve's axisymmetric
periodic condition; top and bottom are slip walls so the exact Couette solution applies.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "geometry"))
import params as P  # noqa: E402


def radial_grid(nr, r1, r2, first=None):
    """Uniform, or two-sided tanh clustering with first cell height `first` at both walls."""
    if first is None:
        return np.linspace(r1, r2, nr + 1)
    xi = np.linspace(0, 1, nr + 1)
    target = first / (r2 - r1)

    def dist(beta):
        s = 0.5 * (1 + np.tanh(beta * (2 * xi - 1)) / np.tanh(beta))
        return s

    lo, hi = 1e-6, 30.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if dist(mid)[1] > target:
            lo = mid
        else:
            hi = mid
    return r1 + (r2 - r1) * dist(0.5 * (lo + hi))


def build_mesh(case_dir, nr, nz, length, wedge_deg=5.0, nth=1, first=None, r1=P.R1, r2=P.R2_IDEAL):
    mdir = os.path.join(case_dir, "MESH.DIR")
    os.makedirs(mdir, exist_ok=True)
    r = radial_grid(nr, r1, r2, first)
    th = np.radians(np.linspace(0, wedge_deg, nth + 1))
    z = np.linspace(0, length, nz + 1)

    ni, nj, nk = nr + 1, nth + 1, nz + 1

    def nid(i, j, k):
        return 1 + i + ni * (j + nj * k)

    I, J, K = np.meshgrid(np.arange(ni), np.arange(nj), np.arange(nk), indexing="ij")
    ids = nid(I, J, K).ravel()
    R = r[I].ravel()
    T = th[J].ravel()
    Z = z[K].ravel()
    crd = np.column_stack([ids, R * np.cos(T), R * np.sin(T), Z])
    order = np.argsort(crd[:, 0])
    crd = crd[order]
    np.savetxt(os.path.join(mdir, "mesh.crd"), crd, fmt="%d %.12e %.12e %.12e")

    Ie, Je, Ke = np.meshgrid(np.arange(nr), np.arange(nth), np.arange(nz), indexing="ij")
    Ie, Je, Ke = Ie.ravel(), Je.ravel(), Ke.ravel()
    eid = 1 + np.arange(Ie.size)
    n1 = nid(Ie, Je, Ke)
    n2 = nid(Ie + 1, Je, Ke)
    n3 = nid(Ie + 1, Je + 1, Ke)
    n4 = nid(Ie, Je + 1, Ke)
    n5, n6, n7, n8 = (nid(Ie, Je, Ke + 1), nid(Ie + 1, Je, Ke + 1),
                      nid(Ie + 1, Je + 1, Ke + 1), nid(Ie, Je + 1, Ke + 1))
    cnn = np.column_stack([eid, n1, n2, n3, n4, n5, n6, n7, n8])
    np.savetxt(os.path.join(mdir, "mesh.cnn"), cnn, fmt="%d")

    # element ids follow the ravel order of (Ie, Je, Ke) built with indexing="ij"
    def eidx(i, j, k):
        return 1 + (i * nth + j) * nz + k

    surfs = {}
    sid = 1

    def add(name, rows):
        nonlocal sid
        arr = np.array(rows, dtype=np.int64)
        arr = np.column_stack([arr[:, 0], np.arange(sid, sid + len(arr)), arr[:, 1:]])
        sid += len(arr)
        surfs[name] = arr

    # node order gives outward normals (out of the fluid)
    add("inner", [(eidx(0, j, k), nid(0, j, k), nid(0, j, k + 1), nid(0, j + 1, k + 1), nid(0, j + 1, k))
                  for j in range(nth) for k in range(nz)])
    add("outer", [(eidx(nr - 1, j, k), nid(nr, j, k), nid(nr, j + 1, k), nid(nr, j + 1, k + 1), nid(nr, j, k + 1))
                  for j in range(nth) for k in range(nz)])
    add("bottom", [(eidx(i, j, 0), nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0))
                   for i in range(nr) for j in range(nth)])
    add("top", [(eidx(i, j, nz - 1), nid(i, j, nz), nid(i + 1, j, nz), nid(i + 1, j + 1, nz), nid(i, j + 1, nz))
                for i in range(nr) for j in range(nth)])
    for name, arr in surfs.items():
        np.savetxt(os.path.join(mdir, f"{name}.srf"), arr, fmt="%d")

    pairs = [(nid(i, 0, k), nid(i, nth, k)) for i in range(ni) for k in range(nk)]
    pbc = np.column_stack([np.arange(1, len(pairs) + 1), np.array(pairs)])
    np.savetxt(os.path.join(mdir, "wedge.pbc"), pbc, fmt="%d")

    # probe line at mid-gap on the j=0 face, all z, for radial-velocity amplitude
    imid = int(np.argmin(np.abs(r - 0.5 * (r1 + r2))))
    probe = np.array([nid(imid, 0, k) for k in range(nk)])
    np.savetxt(os.path.join(mdir, "probe_midgap.nbc"), probe, fmt="%d")
    radial = np.array([nid(i, 0, nk // 2) for i in range(ni)])
    np.savetxt(os.path.join(mdir, "probe_radial.nbc"), radial, fmt="%d")

    meta = dict(nr=nr, nz=nz, nth=nth, wedge_deg=wedge_deg, length=length, r1=r1, r2=r2,
                n_nodes=int(ni * nj * nk), n_elems=int(eid.size), first=first)
    with open(os.path.join(mdir, "mesh_info.txt"), "w") as f:
        for k_, v in meta.items():
            f.write(f"{k_} = {v}\n")
    return meta


def perturbed_initial(case_dir, omega, amp, r1=P.R1, r2=P.R2_IDEAL, nwaves=2, seed=0):
    """Couette profile plus a small axisymmetric (r,z) perturbation to trigger Taylor cells."""
    crd = np.loadtxt(os.path.join(case_dir, "MESH.DIR", "mesh.crd"))
    ids, x, y, z = crd.T
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    L = z.max() - z.min()
    ut = P.couette(r, omega, r1, r2)
    s = (r - r1) / (r2 - r1)
    shape = np.sin(np.pi * s) ** 2
    rng = np.random.default_rng(seed)
    ur = amp * omega * r1 * shape * (np.cos(2 * np.pi * nwaves * z / L) + 0.3 * rng.standard_normal(z.size))
    uz = np.zeros_like(z)
    vx = ur * np.cos(th) - ut * np.sin(th)
    vy = ur * np.sin(th) + ut * np.cos(th)
    out = np.column_stack([ids, vx, vy, uz])
    np.savetxt(os.path.join(case_dir, "MESH.DIR", "init_velocity.nic"), out, fmt="%d %.10e %.10e %.10e")


INP = """ANALYSIS {{
    title                         = "{title}"
    type                          = {atype}
}}

EQUATION {{
    flow                          = navier_stokes
    turbulence                    = {turb}
}}

AUTO_SOLUTION_STRATEGY {{
    max_time_steps                = {nsteps}
    {final_time_line}
    initial_time_increment        = {dt}
    auto_time_increment           = off
    convergence_tolerance         = {tol}
    {stagger_line}
}}

DENSITY_MODEL( "water85" ) {{
    density                       = {rho}
}}

VISCOSITY_MODEL( "water85" ) {{
    viscosity                     = {mu}
}}

MATERIAL_MODEL( "water85" ) {{
    density_model                 = "water85"
    viscosity_model               = "water85"
}}

REFERENCE_FRAME( "rotor" ) {{
    rotation_center               = {{ 0, 0, 0 }}
    angular_velocity              = {{ 0, 0, {omega} }}
}}

COORDINATE {{
    coordinates                   = Read( "MESH.DIR/mesh.crd" )
}}

ELEMENT_SET( "fluid" ) {{
    elements                      = Read( "MESH.DIR/mesh.cnn" )
    shape                         = eight_node_brick
    material_model                = "water85"
}}

SIMPLE_BOUNDARY_CONDITION( "inner" ) {{
    surfaces                      = Read( "MESH.DIR/inner.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    type                          = wall
    reference_frame               = "rotor"
}}

SIMPLE_BOUNDARY_CONDITION( "outer" ) {{
    surfaces                      = Read( "MESH.DIR/outer.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    type                          = wall
}}

SIMPLE_BOUNDARY_CONDITION( "top" ) {{
    surfaces                      = Read( "MESH.DIR/top.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    type                          = slip
}}

SIMPLE_BOUNDARY_CONDITION( "bottom" ) {{
    surfaces                      = Read( "MESH.DIR/bottom.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    type                          = slip
}}

PERIODIC_BOUNDARY_CONDITION( "wedge" ) {{
    variable                      = all
    type                          = axisymmetric
    nodal_pairs                   = Read( "MESH.DIR/wedge.pbc" )
    rotation_axis                 = {{ 0, 0, 0 ; 0, 0, 1 }}
}}

{init_block}

NODAL_OUTPUT {{
    output_frequency              = {out_freq}
}}

RESTART_OUTPUT {{
    output_frequency              = {out_freq}
}}

TIME_HISTORY_OUTPUT( "midgap" ) {{
    nodes                         = Read( "MESH.DIR/probe_midgap.nbc" )
    output_frequency              = {th_freq}
}}

SURFACE_OUTPUT( "inner" ) {{
    surfaces                      = Read( "MESH.DIR/inner.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    integrated_output_frequency   = 1
}}

SURFACE_OUTPUT( "outer" ) {{
    surfaces                      = Read( "MESH.DIR/outer.srf" )
    shape                         = four_node_quad
    element_set                   = "fluid"
    integrated_output_frequency   = 1
}}

RUN
"""


def write_inp(case_dir, title, omega, rho=P.RHO, mu=P.MU, transient=False, turb="none",
              dt=1e10, nsteps=200, final_time=None, tol=1e-6, out_freq=50, th_freq=5,
              init="couette", r1=P.R1, r2=P.R2_IDEAL):
    if init == "file":
        init_block = ('NODAL_INITIAL_CONDITION( velocity ) {\n'
                      '    type                          = nodal_values\n'
                      '    nodal_values                  = Read( "MESH.DIR/init_velocity.nic" )\n'
                      '}\n')
    else:
        init_block = ('NODAL_INITIAL_CONDITION( velocity ) {\n'
                      '    default_values                = { 0, 0, 0 }\n'
                      '}\n')
    init_block += ('NODAL_INITIAL_CONDITION( pressure ) {\n'
                   '    default_value                 = 0\n'
                   '}\n')
    if turb != "none":
        nut0 = 10 * mu / rho
        init_block += (f'NODAL_INITIAL_CONDITION( kinetic_energy ) {{\n'
                       f'    default_value                 = {(0.05 * omega * r1) ** 2:.6e}\n'
                       f'}}\n'
                       f'NODAL_INITIAL_CONDITION( eddy_frequency ) {{\n'
                       f'    default_value                 = {(0.05 * omega * r1) ** 2 / nut0:.6e}\n'
                       f'}}\n'
                       f'NODAL_INITIAL_CONDITION( eddy_viscosity ) {{\n'
                       f'    default_value                 = {nut0:.6e}\n'
                       f'}}\n')
    text = INP.format(
        title=title, atype="transient" if transient else "steady", turb=turb, nsteps=nsteps,
        final_time_line=f"final_time                    = {final_time}" if final_time else "",
        stagger_line="max_stagger_iterations        = 4" if transient else "",
        dt=dt, tol=tol, rho=rho, mu=mu, omega=omega, init_block=init_block,
        out_freq=out_freq, th_freq=th_freq)
    with open(os.path.join(case_dir, "case.inp"), "w") as f:
        f.write(text)
    with open(os.path.join(case_dir, "Acusim.cnf"), "w") as f:
        f.write("problem = case\ninput_file = case.inp\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--nr", type=int, default=80)
    ap.add_argument("--nz", type=int, default=4)
    ap.add_argument("--length", type=float, default=0.01)
    ap.add_argument("--first", type=float, default=None)
    ap.add_argument("--nth", type=int, default=1)
    ap.add_argument("--wedge", type=float, default=5.0)
    ap.add_argument("--re", type=float, required=True)
    ap.add_argument("--nu", type=float, default=P.NU, help="kinematic viscosity used to set omega")
    ap.add_argument("--transient", action="store_true")
    ap.add_argument("--turb", default="none")
    ap.add_argument("--dt", type=float, default=1e10)
    ap.add_argument("--nsteps", type=int, default=200)
    ap.add_argument("--final_time", type=float, default=None)
    ap.add_argument("--out_freq", type=int, default=50)
    ap.add_argument("--th_freq", type=int, default=5)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--perturb", type=float, default=0.0)
    a = ap.parse_args()

    os.makedirs(a.case_dir, exist_ok=True)
    info = build_mesh(a.case_dir, a.nr, a.nz, a.length, wedge_deg=a.wedge, nth=a.nth, first=a.first)
    rho = P.RHO
    mu = rho * a.nu
    omega = a.re * a.nu / (P.R1 * (P.R2_IDEAL - P.R1))
    init = "couette"
    if a.perturb > 0:
        perturbed_initial(a.case_dir, omega, a.perturb)
        init = "file"
    write_inp(a.case_dir, f"L1 Re={a.re:g}", omega, rho=rho, mu=mu, transient=a.transient, turb=a.turb,
              dt=a.dt, nsteps=a.nsteps, final_time=a.final_time, tol=a.tol, out_freq=a.out_freq,
              th_freq=a.th_freq, init=init)
    with open(os.path.join(a.case_dir, "case_info.txt"), "w") as f:
        f.write(f"re = {a.re}\nnu = {a.nu}\nrho = {rho}\nmu = {mu}\nomega = {omega}\n")
        for k, v in info.items():
            f.write(f"{k} = {v}\n")
    print(f"{a.case_dir}: {info['n_elems']} elems, {info['n_nodes']} nodes, omega={omega:.6g} rad/s")


if __name__ == "__main__":
    main()
