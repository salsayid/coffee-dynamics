"""Level 4: export the simulated surface swirl field on a Cartesian grid for comparison with PIVlab.

For the axisymmetric runs the swirl u_theta(r) near the free surface is spun out to a top-down (x, y)
grid, matching what a phone camera sees from above. Output: CSV (x, y, u, v in m and m/s) per time,
in the same layout PIVlab exports (columns x, y, u, v).
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import flowstats as F  # noqa: E402
import params as P  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=os.path.join(ROOT, "cfd", "cases", "L2b_axi_spinup"))
    ap.add_argument("--depth", type=float, default=0.002, help="sample this far below the still surface [m]")
    ap.add_argument("--grid", type=int, default=64)
    ap.add_argument("--out", default=os.path.join(ROOT, "cfd", "post", "surface_velocity"))
    a = ap.parse_args()

    info = json.load(open(os.path.join(a.case, "case_info.json")))
    mesh = F.Mesh(os.path.join(a.case, "MESH.DIR"))
    z0 = P.LIQUID_H - a.depth
    band = (np.abs(mesh.z - z0) < 6e-4) & (np.abs(mesh.th) < 1e-6)
    snaps = F.nodal_series(a.case, "velocity", every=info["out_every"])
    os.makedirs(a.out, exist_ok=True)
    R = P.mug_radius(z0)
    xs = np.linspace(-R, R, a.grid)
    X, Y = np.meshgrid(xs, xs)
    Rg = np.hypot(X, Y)
    Tg = np.arctan2(Y, X)
    for step, vel in snaps.items():
        t = step * info["dt"]
        r = mesh.r[band]
        ut = -vel[band, 0] * np.sin(mesh.th[band]) + vel[band, 1] * np.cos(mesh.th[band])
        o = np.argsort(r)
        u_theta = np.interp(Rg, r[o], ut[o], left=0.0, right=0.0)
        u_theta[Rg > R] = np.nan
        u = -u_theta * np.sin(Tg)
        v = u_theta * np.cos(Tg)
        fn = os.path.join(a.out, f"sim_surface_t{t:07.3f}s.csv")
        np.savetxt(fn, np.column_stack([X.ravel(), Y.ravel(), u.ravel(), v.ravel()]), delimiter=",",
                   header="x [m],y [m],u [m/s],v [m/s]", comments="")
    print(f"wrote {len(snaps)} frames to {a.out}")


if __name__ == "__main__":
    main()
