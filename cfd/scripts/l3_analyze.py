"""Level 3: free-surface vortex in the axisymmetric mug. When does air reach the frother?"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import flowstats as F  # noqa: E402
import params as P  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def section(mesh_dir, mesh):
    """Triangulation of the theta = 0 face (r, z) using the first face of each wedge/hex."""
    meta = json.load(open(os.path.join(mesh_dir, "mesh_meta.json")))
    tris = []
    n2 = mesh.xyz.shape[0] // 2
    for s in meta["element_sets"].values():
        c = np.loadtxt(os.path.join(mesh_dir, s["file"]), dtype=np.int64)[:, 1:] - 1
        base = c[:, :3] if s["shape"] == "six_node_wedge" else c[:, :4]
        base = np.where(base >= n2, base - n2, base)
        if base.shape[1] == 3:
            tris.append(base)
        else:
            tris += [base[:, [0, 1, 2]], base[:, [0, 2, 3]]]
    tris = np.vstack(tris)
    return mtri.Triangulation(mesh.r[:n2], mesh.z[:n2], tris), n2


def main(case=os.path.join(ROOT, "cfd", "cases", "L3_axi_freesurface")):
    info = json.load(open(os.path.join(case, "case_info.json")))
    scale = 360.0 / info["wedge_deg"]
    dt = info["dt"]
    mdir = os.path.join(case, "MESH.DIR")
    mesh = F.Mesh(mdir)
    tri, n2 = section(mdir, mesh)
    ls = F.nodal_series(case, "levelset", every=info["out_every"])
    wf = {k: (v[:, 0] if v.ndim > 1 else v) for k, v in
          F.nodal_series(case, "field", every=info["out_every"]).items()}  # column 0 = water
    vel = F.nodal_series(case, "velocity", every=info["out_every"])
    tq = F.torques(case)

    # nodes in a 0.6 mm shell around the torus surface: is the frother in water or in air?
    dtor = np.abs(np.hypot(mesh.r - P.COIL_RING_R, mesh.z - P.HEAD_Z) - P.COIL_D / 2)
    shell = dtor < 6e-4
    rows, shapes = [], {}
    for step in sorted(ls):
        phi = ls[step][:n2]
        t = step * dt
        fig = plt.figure()
        cs = plt.tricontour(tri, phi, levels=[0.0])
        plt.close(fig)
        segs = [p for p in cs.allsegs[0] if len(p) > 2]
        if not segs:
            continue
        main_seg = max(segs, key=lambda p: np.ptp(p[:, 0]))
        shapes[t] = main_seg
        allpts = np.vstack(segs)
        near_axis = allpts[allpts[:, 0] < 0.004]
        dip_bottom = float(near_axis[:, 1].min()) if len(near_axis) else float(main_seg[:, 1].min())
        rim_h = float(allpts[allpts[:, 0] > P.mug_radius(0.04) - 0.002][:, 1].max()) if np.any(
            allpts[:, 0] > P.mug_radius(0.04) - 0.002) else float("nan")
        water = wf[step]
        wet = float(np.sum(mesh.vol[shell] * water[shell]) / np.sum(mesh.vol[shell]))
        vol_water = float(np.sum(mesh.vol * water) * scale)
        st = F.swirl_stats(mesh, vel[step], liquid=water)
        rows.append(dict(time=t, dip_bottom_z=dip_bottom, dip_depth=P.LIQUID_H - dip_bottom, wall_height=rim_h,
                         frother_wetted=wet, water_volume=vol_water, omega_bulk=st["omega_bulk"]))
    v0 = rows[0]["water_volume"]
    for r in rows:
        r["water_volume_change"] = r["water_volume"] / v0 - 1
    reach = next((r["time"] for r in rows if r["dip_bottom_z"] < P.HEAD_Z + P.COIL_D / 2 + 0.0005), None)
    summary = dict(t_end=rows[-1]["time"], dip_depth_final_mm=rows[-1]["dip_depth"] * 1e3,
                   wall_rise_final_mm=(rows[-1]["wall_height"] - P.LIQUID_H) * 1e3,
                   time_air_reaches_frother=reach, frother_wetted_final=rows[-1]["frother_wetted"],
                   water_volume_drift=rows[-1]["water_volume_change"])
    json.dump({"summary": summary, "rows": rows}, open(os.path.join(ROOT, "cfd", "post", "l3_freesurface.json"), "w"),
              indent=1)
    print(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [1.1, 1, 1]})
    a = ax[0]
    zz = np.linspace(0, P.MUG_H, 50)
    for sgn in (1, -1):
        a.plot(sgn * P.mug_radius(zz) * 1e3, zz * 1e3, color="k", lw=1.5)
        th = np.linspace(0, 2 * np.pi, 80)
        a.fill(sgn * (P.COIL_RING_R + P.COIL_D / 2 * np.cos(th)) * 1e3, (P.HEAD_Z + P.COIL_D / 2 * np.sin(th)) * 1e3,
               color="#888888")
        a.fill([sgn * 0, sgn * P.SHAFT_D / 2 * 1e3, sgn * P.SHAFT_D / 2 * 1e3, 0],
               [P.HEAD_Z * 1e3 - 0.8, P.HEAD_Z * 1e3 - 0.8, P.MUG_H * 1e3, P.MUG_H * 1e3], color="#888888")
    a.plot([-P.MUG_R_BOTTOM * 1e3, P.MUG_R_BOTTOM * 1e3], [0, 0], color="k", lw=1.5)
    times = sorted(shapes)
    pick = [times[i] for i in np.unique(np.linspace(0, len(times) - 1, 6).astype(int))]
    cm = plt.cm.viridis
    for i, tt in enumerate(pick):
        p = shapes[tt]
        col = cm(i / max(1, len(pick) - 1))
        a.plot(p[:, 0] * 1e3, p[:, 1] * 1e3, color=col, lw=1.8, label=f"t = {tt:.2f} s")
        a.plot(-p[:, 0] * 1e3, p[:, 1] * 1e3, color=col, lw=1.8)
    a.set_aspect("equal")
    a.set_xlabel("x [mm]")
    a.set_ylabel("height [mm]")
    a.set_title("Coffee surface (cross-section of my mug)")
    a.legend(frameon=False, fontsize=8, loc="upper right")

    t = np.array([r["time"] for r in rows])
    ax[1].plot(t, [r["dip_depth"] * 1e3 for r in rows], color="#2c5f8a", lw=2, label="depth of the central dip")
    ax[1].axhline((P.LIQUID_H - P.HEAD_Z - P.COIL_D / 2) * 1e3, color="#d1495b", ls="--", lw=1,
                  label="top of the frother head")
    ax[1].set_xlabel("time since switching on [s]")
    ax[1].set_ylabel("dip depth below the still surface [mm]")
    ax[1].set_title("The vortex digs down to the frother")
    ax[1].legend(frameon=False)

    ax[2].plot(tq["time"], -tq["frother"] * scale * 1e3, color="#d1495b", lw=1.4, label="frother torque")
    ax[2].set_xlabel("time since switching on [s]")
    ax[2].set_ylabel("torque [mN·m]", color="#d1495b")
    ax2 = ax[2].twinx()
    ax2.plot(t, [r["frother_wetted"] * 100 for r in rows], color="#2c5f8a", lw=1.8)
    ax2.set_ylabel("frother surroundings that are water [%]", color="#2c5f8a")
    ax2.set_ylim(0, 105)
    ax[2].set_title("Air at the frother: resistance drops")
    for a in ax:
        a.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "figures", "l3_freesurface.png"), dpi=170, bbox_inches="tight")
    print("saved figures/l3_freesurface.png")


if __name__ == "__main__":
    main(*sys.argv[1:])
