"""Level 2b: spin-up from rest in the axisymmetric mug model. Does the frother-relative Re really drop?"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import flowstats as F  # noqa: E402
import params as P  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def main(case=os.path.join(ROOT, "cfd", "cases", "L2b_axi_spinup")):
    info = json.load(open(os.path.join(case, "case_info.json")))
    scale = 360.0 / info["wedge_deg"]
    dt = info["dt"]
    mesh = F.Mesh(os.path.join(case, "MESH.DIR"))
    tq = F.torques(case)
    t_tq = tq["time"]
    frother = -tq["frother"] * scale  # torque the motor must supply (positive)
    mug = tq["mug"] * scale

    snaps = F.nodal_series(case, "velocity", every=info["out_every"])
    rows = []
    profiles = {}
    zsel = np.abs(mesh.z - 0.020) < 4e-4
    for step, vel in snaps.items():
        s = F.swirl_stats(mesh, vel)
        s["time"] = step * dt
        rows.append(s)
        ut = -vel[:, 0] * np.sin(mesh.th) + vel[:, 1] * np.cos(mesh.th)
        sel = zsel & (np.abs(mesh.th) < 1e-6)
        o = np.argsort(mesh.r[sel])
        profiles[s["time"]] = (mesh.r[sel][o], ut[sel][o])
    t = np.array([r["time"] for r in rows])
    wb = np.array([r["omega_bulk"] for r in rows])
    wr = np.array([r["omega_ring"] for r in rows])
    rerel = np.array([r["re_rel"] for r in rows])
    re0 = rows[0]["re_static"]

    final = wb[-1]
    tau = float(np.interp(0.632 * final, wb, t)) if final > 0 and wb.max() >= 0.632 * final else float("nan")
    summary = dict(t_end=float(t[-1]), omega_bulk_final=float(final), omega_ring_final=float(wr[-1]),
                   re_rel_start=float(rerel[0]), re_rel_final=float(rerel[-1]), re_static=float(re0),
                   spinup_time_63=tau, torque_start=float(frother[min(5, len(frother) - 1)]),
                   torque_final=float(frother[-1]), mug_torque_final=float(mug[-1]),
                   power_final_W=float(frother[-1] * P.OMEGA))
    json.dump({"summary": summary, "rows": rows}, open(os.path.join(ROOT, "cfd", "post", "l2b_spinup.json"), "w"),
              indent=1)
    print(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    ax[0].plot(t_tq, frother * 1e3, color="#d1495b", lw=1.6, label="torque on frother (your hand)")
    ax[0].plot(t_tq, mug * 1e3, color="#2c5f8a", lw=1.6, label="torque passed to the mug wall")
    ax[0].set_xlabel("time since switching on [s]")
    ax[0].set_ylabel("torque [mN·m]")
    ax[0].set_title("Resistance felt in the hand")
    ax[0].legend(frameon=False)

    ax[1].plot(t, wb, color="#2c5f8a", lw=2, label="whole cup (angular-momentum average)")
    ax[1].plot(t, wr, color="#edae49", lw=2, label="ring just outside the frother")
    ax[1].axhline(P.OMEGA, color="gray", ls="--", lw=0.8)
    ax[1].text(t[-1], P.OMEGA * 0.97, "frother 1257 rad/s", ha="right", va="top", color="gray", fontsize=9)
    ax[1].set_xlabel("time since switching on [s]")
    ax[1].set_ylabel("rotation rate [rad/s]")
    ax[1].set_title("The coffee spins up")
    ax[1].legend(frameon=False, loc="center right")
    ax1b = ax[1].twinx()
    ax1b.plot(t, rerel / 1e3, color="#d1495b", lw=1.4, ls=":")
    ax1b.set_ylabel("frother-relative Re [thousands]", color="#d1495b")
    ax1b.tick_params(axis="y", colors="#d1495b")

    times = sorted(profiles)
    pick = [times[i] for i in np.unique(np.linspace(0, len(times) - 1, 6).astype(int))]
    cm = plt.cm.inferno
    for i, tt in enumerate(pick):
        r, u = profiles[tt]
        ax[2].plot(r * 1e3, u, color=cm(0.15 + 0.7 * i / max(1, len(pick) - 1)), lw=1.8, label=f"t = {tt:.1f} s")
    rr = np.linspace(P.R1, P.mug_radius(0.02), 200)
    ax[2].plot(rr * 1e3, P.couette(rr, r2=P.mug_radius(0.02)), "k--", lw=1, label="laminar Couette (README)")
    ax[2].set_xlabel("radius [mm]")
    ax[2].set_ylabel("swirl velocity $u_\\theta$ at mid-depth [m/s]")
    ax[2].set_title("Swirl profile developing")
    ax[2].legend(frameon=False, fontsize=8)
    for a in ax:
        a.grid(alpha=0.2)
    for a in (ax[0], ax[2]):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "figures", "l2b_spinup.png"), dpi=170, bbox_inches="tight")
    print("saved figures/l2b_spinup.png")


if __name__ == "__main__":
    main(*sys.argv[1:])
