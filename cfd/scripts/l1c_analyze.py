"""Level 1c: turbulent (SST RANS) vs laminar profile and torque at the real frother Reynolds number."""
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
CASE = os.path.join(ROOT, "cfd", "cases", "L1c_rans_Re981k")


def wendt_G(re, eta):
    """Wendt (1933) empirical torque law, G = T / (rho nu^2 L); fitted for 1e4 < Re < 1e5."""
    return 0.23 * eta**1.5 / (1 - eta) ** 1.75 * re**1.7


def main():
    info = {k.strip(): v.strip() for k, v in (l.split("=", 1) for l in open(os.path.join(CASE, "case_info.txt")))}
    omega, L, wedge = float(info["omega"]), float(info["length"]), float(info["wedge_deg"])
    re = float(info["re"])
    ids, xyz, v = acu.nodal(CASE, "velocity")
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    th = np.arctan2(xyz[:, 1], xyz[:, 0])
    ut = -v[:, 0] * np.sin(th) + v[:, 1] * np.cos(th)
    sel = (np.abs(th) < 1e-9) & (np.abs(xyz[:, 2] - xyz[:, 2].max() / 2) < 1e-6 + L / 4)
    o = np.argsort(r[sel])
    rr, uu = r[sel][o], ut[sel][o]
    U = omega * P.R1

    mz = acu.surface_integrals(CASE)[1][-1, acu.OSI_COLS["moment"]][2]
    t_full_per_len = abs(mz) * 360 / wedge / L
    t_lam_per_len = P.couette_torque_per_length(omega)
    eta = P.R1 / P.R2_IDEAL
    G_sim = t_full_per_len / (P.RHO * P.NU**2)
    G_w = wendt_G(re, eta)
    out = dict(re=re, midgap_u_over_U_turb=float(np.interp(0.5 * (P.R1 + P.R2_IDEAL), rr, uu) / U),
               midgap_u_over_U_lam=float(P.couette(0.5 * (P.R1 + P.R2_IDEAL), omega) / U),
               torque_per_length_turb=t_full_per_len, torque_per_length_lam=t_lam_per_len,
               torque_ratio_turb_over_lam=t_full_per_len / t_lam_per_len,
               G_sim=G_sim, G_wendt_extrapolated=G_w, sim_over_wendt=G_sim / G_w)
    json.dump(out, open(os.path.join(ROOT, "cfd", "post", "l1c_turbulent.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    xr = np.linspace(P.R1, P.R2_IDEAL, 300)
    ax.plot(xr * 1e3, P.couette(xr, omega) / U, color="#2c5f8a", lw=2, label="laminar (my derivation)")
    ax.plot(rr * 1e3, uu / U, color="#d1495b", lw=2, label=f"turbulent, SST RANS at Re = {re/1e3:.0f}k")
    ax.set_xlabel("radius [mm]")
    ax.set_ylabel("$u_\\theta$ / (frother tip speed)")
    ax.set_title(f"Same geometry, real Reynolds number: torque is {out['torque_ratio_turb_over_lam']:.0f}x the laminar value")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "figures", "l1c_turbulent_profile.png"), dpi=170, bbox_inches="tight")
    print("saved figures/l1c_turbulent_profile.png")


if __name__ == "__main__":
    main()
