"""Level 2 (3D, real frother with spokes): torque balance and swirl, compared with the axisymmetric wedge model."""
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
import acu  # noqa: E402
import flowstats as F  # noqa: E402
import params as P  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = os.path.join(ROOT, "cfd", "cases", "L2c_mrf_A_m1_pseudo")
STEADY = os.path.join(ROOT, "cfd", "cases", "L2a_mrf_A_m1")
MESH = os.path.join(ROOT, "cfd", "mesh", "acu_A_L2_m1")


def main():
    t1 = F.torques(STEADY)
    t2 = F.torques(CASE)
    mesh = F.Mesh(MESH)
    _, _, vel = acu.nodal(CASE, "velocity", mesh_dir=MESH)
    st = F.swirl_stats(mesh, vel)
    fr, mg = -t2["frother"][-1], t2["mug"][-1]
    out = dict(time=float(t2["time"][-1]), frother_torque=float(fr), mug_torque=float(mg),
               balance=float(mg / fr), power_W=float(fr * P.OMEGA), **st)
    axi = os.path.join(ROOT, "cfd", "post", "l2b_spinup.json")
    if os.path.exists(axi):
        a = json.load(open(axi))["summary"]
        out["axisymmetric_end_state"] = {k: a[k] for k in ("t_end", "omega_bulk_final", "omega_ring_final",
                                                            "re_rel_final", "torque_final", "mug_torque_final")}
    json.dump(out, open(os.path.join(ROOT, "cfd", "post", "l2c_3d.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))

    fig, ax = plt.subplots(figsize=(8, 4.4))
    n1 = len(t1["frother"])
    ax.plot(np.arange(n1), -t1["frother"] * 1e3, color="#d1495b", lw=1.2, alpha=0.5)
    ax.plot(np.arange(n1), t1["mug"] * 1e3, color="#2c5f8a", lw=1.2, alpha=0.5)
    x2 = n1 + np.arange(len(t2["frother"]))
    ax.plot(x2, -t2["frother"] * 1e3, color="#d1495b", lw=2, label="frother")
    ax.plot(x2, t2["mug"] * 1e3, color="#2c5f8a", lw=2, label="mug walls")
    ax.axvline(n1, color="gray", ls="--", lw=0.8)
    ax.text(n1 * 0.5, ax.get_ylim()[1] * 0.95, "steady solver\n(stalls)", ha="center", va="top", fontsize=9)
    ax.text(n1 + len(x2) * 0.5, ax.get_ylim()[1] * 0.95, "time-marching, 50 ms steps\n(torques balance)",
            ha="center", va="top", fontsize=9)
    ax.set_xlabel("solver step")
    ax.set_ylabel("torque [mN·m]")
    ax.set_title("3D model: torque on the frother vs. torque on the mug")
    ax.legend(frameon=False, loc="center right")
    ax.grid(alpha=0.2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "figures", "l2c_3d_torque.png"), dpi=170, bbox_inches="tight")
    print("saved figures/l2c_3d_torque.png")


if __name__ == "__main__":
    main()
