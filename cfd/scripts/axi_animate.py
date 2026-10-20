"""Animate an axisymmetric case as a mirrored mug cross-section (speed in inferno, free surface in white)."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import flowstats as F  # noqa: E402
import params as P  # noqa: E402
from l3_analyze import section  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def mug_outline(ax, top):
    zz = np.linspace(0, top, 50)
    for s in (1, -1):
        ax.plot(s * P.mug_radius(zz) * 1e3, zz * 1e3, color="#cccccc", lw=2)
        th = np.linspace(0, 2 * np.pi, 60)
        ax.fill(s * (P.COIL_RING_R + P.COIL_D / 2 * np.cos(th)) * 1e3, (P.HEAD_Z + P.COIL_D / 2 * np.sin(th)) * 1e3,
                color="#9a9a9a", zorder=5)
    ax.fill([-P.SHAFT_D / 2 * 1e3, P.SHAFT_D / 2 * 1e3, P.SHAFT_D / 2 * 1e3, -P.SHAFT_D / 2 * 1e3],
            [P.HEAD_Z * 1e3 - 0.8, P.HEAD_Z * 1e3 - 0.8, top * 1e3, top * 1e3], color="#9a9a9a", zorder=5)
    ax.plot([-P.MUG_R_BOTTOM * 1e3, P.MUG_R_BOTTOM * 1e3], [0, 0], color="#cccccc", lw=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--vmax", type=float, default=4.0)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--every", type=int, default=None)
    ap.add_argument("--title", default="")
    a = ap.parse_args()

    info = json.load(open(os.path.join(a.case, "case_info.json")))
    mdir = os.path.join(a.case, "MESH.DIR")
    mesh = F.Mesh(mdir)
    tri, n2 = section(mdir, mesh)
    level = info["level"]
    top = P.LIQUID_H if level == 2 else P.MUG_H
    every = a.every or info["out_every"]
    vel = F.nodal_series(a.case, "velocity", every=every)
    # water volume fraction (column 0 of "field"); avoids relying on the level-set sign convention
    ls = ({k: (v[:, 0] if v.ndim > 1 else v) for k, v in F.nodal_series(a.case, "field", every=every).items()}
          if level == 3 else {})

    import matplotlib.tri as mtri
    xs = np.concatenate([tri.x, -tri.x])
    zs = np.concatenate([tri.y, tri.y])
    tris = np.vstack([tri.triangles, tri.triangles + n2])
    tri2 = mtri.Triangulation(xs * 1e3, zs * 1e3, tris)
    tmp = tempfile.mkdtemp()
    for k, (step, v) in enumerate(vel.items()):
        sp = np.linalg.norm(v[:n2], axis=1)
        s2 = np.concatenate([sp, sp])
        fig, ax = plt.subplots(figsize=(6.4, 6.4 if level == 2 else 9.0), facecolor="black")
        ax.set_facecolor("black")
        if level == 3:
            phi = ls[step][:n2]
            p2 = np.concatenate([phi, phi])
            s2 = np.where(p2 > 0.5, s2, np.nan)  # show the coffee only; air stays black
        ax.tripcolor(tri2, s2, cmap="inferno", vmin=0, vmax=a.vmax, shading="gouraud")
        if level == 3:
            ax.tricontour(tri2, p2, levels=[0.5], colors="white", linewidths=1.2)
        mug_outline(ax, top)
        ax.set_aspect("equal")
        ax.set_xlim(-37, 37)
        ax.set_ylim(-2, top * 1e3 + 2)
        ax.axis("off")
        t = step * info["dt"]
        ax.text(0.02, 0.98, f"t = {t:5.2f} s", transform=ax.transAxes, color="white", family="monospace",
                va="top", fontsize=13)
        if a.title:
            ax.text(0.98, 0.98, a.title, transform=ax.transAxes, color="#bbbbbb", ha="right", va="top", fontsize=10)
        sm = plt.cm.ScalarMappable(cmap="inferno", norm=plt.Normalize(0, a.vmax))
        cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.01)
        cb.set_label("speed [m/s]", color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        plt.setp(cb.ax.get_yticklabels(), color="white")
        fig.savefig(os.path.join(tmp, f"f{k:04d}.png"), dpi=a.dpi, facecolor="black", bbox_inches="tight")
        plt.close(fig)
    ff = shutil.which("ffmpeg")
    pal = os.path.join(tmp, "pal.png")
    subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%04d.png"),
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2,palettegen", pal], check=True)
    subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%04d.png"),
                    "-i", pal, "-lavfi", "scale=trunc(iw/2)*2:trunc(ih/2)*2[x];[x][1:v]paletteuse", a.out], check=True)
    mp4 = os.path.splitext(a.out)[0] + ".mp4"
    subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%04d.png"),
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(f"wrote {a.out} and {mp4} ({len(vel)} frames)")


if __name__ == "__main__":
    main()
