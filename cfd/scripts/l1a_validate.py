"""Level 1a: compare AcuSolve laminar Couette results with the exact solution on three meshes."""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import acu  # noqa: E402
import params as P  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG = os.path.join(ROOT, "figures")


def read_info(case):
    d = {}
    for line in open(os.path.join(case, "case_info.txt")):
        k, v = line.split("=", 1)
        d[k.strip()] = v.strip()
    return d


def analyse(case):
    info = read_info(case)
    omega, mu = float(info["omega"]), float(info["mu"])
    L, wedge = float(info["length"]), float(info["wedge_deg"])
    ids, xyz, vel = acu.nodal(case, "velocity")
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    th = np.arctan2(xyz[:, 1], xyz[:, 0])
    ut = -vel[:, 0] * np.sin(th) + vel[:, 1] * np.cos(th)
    exact = P.couette(r, omega)
    scale = omega * P.R1
    err = (ut - exact) / scale
    l2 = float(np.sqrt(np.mean(err**2)))
    linf = float(np.max(np.abs(err)))

    srf = acu.surface_integrals(case)
    mz = srf[1][-1, acu.OSI_COLS["moment"]][2]
    t_exact = P.couette_torque_per_length(omega, mu) * L * math.radians(wedge) / (2 * math.pi)
    torque_err = float((abs(mz) - t_exact) / t_exact)

    sel = (np.abs(th) < 1e-9) & (np.abs(xyz[:, 2] - xyz[:, 2].max() / 2) < 1e-9 + L / 8)
    order = np.argsort(r[sel])
    return dict(case=os.path.basename(case), nr=int(info["nr"]), omega=omega, l2=l2, linf=linf,
                torque_sim=float(abs(mz)), torque_exact=t_exact, torque_err=torque_err,
                r=r[sel][order].tolist(), ut=ut[sel][order].tolist())


def gci(f1, f2, f3, ratio=2.0):
    """Observed order and fine-grid GCI (Roache) for three systematically refined meshes, f1 finest."""
    p = math.log(abs((f3 - f2) / (f2 - f1))) / math.log(ratio)
    e21 = abs((f1 - f2) / f1)
    return p, 1.25 * e21 / (ratio**p - 1)


def main(cases):
    res = [analyse(c) for c in cases]
    res.sort(key=lambda d: d["nr"])
    for d in res:
        print(f"{d['case']}: nr={d['nr']}  L2={d['l2']:.3e}  Linf={d['linf']:.3e}  "
              f"torque err={d['torque_err']*100:+.4f}%")
    out = {"cases": [{k: v for k, v in d.items() if k not in ("r", "ut")} for d in res]}
    if len(res) >= 3:
        f3, f2, f1 = (d["torque_sim"] for d in res[:3])
        p, g = gci(f1, f2, f3)
        out["torque_order"] = p
        out["torque_gci_fine"] = g
        print(f"torque: observed order p={p:.2f}, fine-grid GCI={g*100:.4f}%")
    json.dump(out, open(os.path.join(ROOT, "cfd", "post", "l1a_validation.json"), "w"), indent=2)

    os.makedirs(FIG, exist_ok=True)
    fine = res[-1]
    omega = fine["omega"]
    rr = np.linspace(P.R1, P.R2_IDEAL, 400)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6))
    a.plot(rr * 1e3, P.couette(rr) , color="#2c5f8a", lw=2.2, label="Exact: $u_\\theta = Ar + B/r$")
    s = P.OMEGA / omega
    a.plot(np.array(fine["r"])[::4] * 1e3, np.array(fine["ut"])[::4] * s, "o", ms=5, mfc="none",
           color="#d1495b", label=f"AcuSolve ({fine['nr']} radial cells)")
    a.set_xlabel("Radial position r [mm]")
    a.set_ylabel("$u_\\theta$ [m/s] (scaled to 12,000 RPM)")
    a.set_title("Laminar profile: simulation vs. my hand derivation")
    a.legend(frameon=False)
    a.grid(alpha=0.2)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)

    nrs = np.array([d["nr"] for d in res])
    b.loglog(nrs, [d["l2"] for d in res], "o-", color="#2c5f8a", label="Velocity error (RMS)")
    b.loglog(nrs, [abs(d["torque_err"]) for d in res], "s-", color="#d1495b", label="Torque error")
    ref = abs(res[0]["torque_err"]) * (nrs[0] / nrs) ** 2
    b.loglog(nrs, ref, "--", color="gray", lw=1, label="2nd-order slope")
    b.set_xlabel("Radial cells")
    b.set_ylabel("Relative error")
    b.set_title("Mesh convergence")
    b.legend(frameon=False)
    b.grid(alpha=0.2, which="both")
    for sp in ("top", "right"):
        b.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "l1a_validation.png"), dpi=180, bbox_inches="tight")
    print("saved figures/l1a_validation.png")


if __name__ == "__main__":
    main(sys.argv[1:])
