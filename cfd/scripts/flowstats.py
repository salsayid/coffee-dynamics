"""Integral flow quantities for the Level 2/3 cases: torques, mean swirl, frother-relative Reynolds number."""
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import acu  # noqa: E402
import params as P  # noqa: E402

ROTATING = ("frother_head", "shaft_rotor", "shaft_stator")
FIXED = ("mug_wall", "mug_bottom")


def surface_names(case):
    text = open(os.path.join(case, "case.inp")).read()
    return re.findall(r'SURFACE_OUTPUT\( "([^"]+)" \)', text)


def torques(case):
    """Per-step z-moment summed over rotating (frother) and fixed (mug) walls."""
    names = surface_names(case)
    srf = acu.surface_integrals(case)
    out = {}
    for i, name in enumerate(names, start=1):
        if i not in srf:
            continue
        a = srf[i]
        out[name] = (a[:, 0], a[:, 1], a[:, acu.OSI_COLS["moment"]][:, 2])
    steps = next(iter(out.values()))[0]
    times = next(iter(out.values()))[1]
    n = min(len(v[0]) for v in out.values())
    frother = sum(v[2][:n] for k, v in out.items() if k.split("__")[0] in ROTATING)
    mug = sum(v[2][:n] for k, v in out.items() if k.split("__")[0] in FIXED)
    return dict(step=steps[:n], time=times[:n], frother=frother, mug=mug)


class Mesh:
    def __init__(self, mesh_dir):
        crd = np.loadtxt(os.path.join(mesh_dir, "mesh.crd"))
        self.xyz = crd[:, 1:]
        meta = json.load(open(os.path.join(mesh_dir, "mesh_meta.json")))
        n = len(self.xyz)
        vol = np.zeros(n)
        rotor = np.zeros(n, bool)
        for name, s in meta["element_sets"].items():
            c = np.loadtxt(os.path.join(mesh_dir, s["file"]), dtype=np.int64)[:, 1:] - 1
            p = self.xyz[c]
            if s["shape"] == "four_node_tet":
                tets = [(0, 1, 2, 3)]
            elif s["shape"] == "six_node_wedge":
                tets = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
            else:
                tets = [(0, 1, 2, 4), (0, 2, 3, 4)]
            ev = np.zeros(len(c))
            for t in tets:
                a, b, cc, d = (p[:, k] for k in t)
                ev += np.abs(np.einsum("ij,ij->i", np.cross(b - a, cc - a), d - a)) / 6
            np.add.at(vol, c.ravel(), np.repeat(ev / c.shape[1], c.shape[1]))
            if s["region"] == "rotor":
                rotor[np.unique(c)] = True
        self.vol = vol
        self.rotor = rotor
        self.r = np.hypot(self.xyz[:, 0], self.xyz[:, 1])
        self.th = np.arctan2(self.xyz[:, 1], self.xyz[:, 0])
        self.z = self.xyz[:, 2]


def swirl_stats(mesh, vel, omega=P.OMEGA, liquid=None):
    """liquid: optional nodal weight (1 = water) for free-surface cases."""
    ut = -vel[:, 0] * np.sin(mesh.th) + vel[:, 1] * np.cos(mesh.th)
    w = mesh.vol * (liquid if liquid is not None else 1.0)
    I = np.sum(w * mesh.r**2)
    L = np.sum(w * mesh.r * ut)
    omega_bulk = L / I  # angular-momentum-weighted mean rotation rate
    # fluid rotation just outside the head, in a ring R1 < r < R1 + 5 mm, |z - head| < 3 mm
    ring = (mesh.r > P.R1) & (mesh.r < P.R1 + 0.005) & (np.abs(mesh.z - P.HEAD_Z) < 0.003)
    omega_ring = np.sum(w[ring] * ut[ring] / mesh.r[ring]) / np.sum(w[ring])
    gap = P.mug_radius(P.HEAD_Z) - P.R1
    re_rel = (omega - omega_ring) * P.R1 * gap / P.NU
    ke = 0.5 * P.RHO * np.sum(w * np.sum(vel**2, axis=1))
    return dict(omega_bulk=float(omega_bulk), omega_ring=float(omega_ring), re_rel=float(re_rel),
                re_static=float(omega * P.R1 * gap / P.NU), kinetic_energy=float(ke),
                ang_momentum=float(P.RHO * L))
