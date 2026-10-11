"""Add cylindrical velocity components (and speed) to an acuTrans EnSight case for META."""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ensight as E  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_file")
    ap.add_argument("--omega", type=float, default=None, help="if set, also write swirl relative to the frother")
    a = ap.parse_args()
    d = os.path.dirname(os.path.abspath(a.case_file))
    text = open(a.case_file).read()
    geo = os.path.join(d, text.split("model:")[1].split()[0])
    parts = E.read_geo(geo)
    steps, pad = E.case_steps(a.case_file)
    vel_pat = text.split("per node: velocity")[1].split()[0]
    base = os.path.join(d, os.path.dirname(vel_pat))
    stem = os.path.basename(vel_pat).split(".")[0]
    for s in steps:
        vel = E.read_var(os.path.join(d, vel_pat.replace("*" * pad, str(s).zfill(pad))), parts, vector=True)
        ur, ut, sp = [], [], []
        for part, v in zip(parts, vel):
            x, y = part["xyz"][:, 0], part["xyz"][:, 1]
            th = np.arctan2(y, x)
            ur.append(v[:, 0] * np.cos(th) + v[:, 1] * np.sin(th))
            ut.append(-v[:, 0] * np.sin(th) + v[:, 1] * np.cos(th))
            sp.append(np.linalg.norm(v, axis=1))
        num = str(s).zfill(pad)
        for name, vals in (("u_radial", ur), ("u_swirl", ut), ("speed", sp)):
            E.write_var(os.path.join(base, f"{stem}.{name}.{num}"), parts, vals, name)
    rel = os.path.dirname(vel_pat)
    for name in ("u_radial", "u_swirl", "speed"):
        E.add_to_case(a.case_file, name, f"{rel}/{stem}.{name}.{'*' * pad}")
    print(f"added u_radial, u_swirl, speed for {len(steps)} step(s) to {a.case_file}; parts: "
          + ", ".join(f"{p['id']}:{p['name']}({len(p['node_ids'])})" for p in parts))


if __name__ == "__main__":
    main()
