"""Level 4: side-by-side video of the real coffee (phone, top-down) and the simulated surface swirl.

  python cfd/scripts/video_compare.py --real my_coffee.mp4 --start 1.25 --out figures/compare.mp4

--start is the moment (seconds into your video) when the frother switches on, so both clips share
t = 0. The simulation side is rendered from the CSVs written by extract_surface_velocity.py, with
the same inferno colormap as the rest of the project.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def render_sim_frames(csv_dir, frame_dir, fps, duration, vmax):
    files = sorted(glob.glob(os.path.join(csv_dir, "sim_surface_t*.csv")))
    times = np.array([float(re.findall(r"t([\d.]+)s", f)[0]) for f in files])
    n = int(duration * fps)
    for i in range(n):
        t = i / fps
        f = files[int(np.argmin(np.abs(times - t)))]
        d = np.loadtxt(f, delimiter=",", skiprows=1)
        g = int(round(np.sqrt(len(d))))
        x, y = d[:, 0].reshape(g, g), d[:, 1].reshape(g, g)
        sp = np.hypot(d[:, 2], d[:, 3]).reshape(g, g)
        fig, ax = plt.subplots(figsize=(6, 6), facecolor="black")
        ax.set_facecolor("black")
        ax.pcolormesh(x * 1e3, y * 1e3, sp, cmap="inferno", vmin=0, vmax=vmax, shading="auto")
        s = 4
        ax.quiver(x[::s, ::s] * 1e3, y[::s, ::s] * 1e3, d[:, 2].reshape(g, g)[::s, ::s],
                  d[:, 3].reshape(g, g)[::s, ::s], color="white", alpha=0.6, scale=vmax * 25)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0.03, 0.96, f"simulation  t = {t:5.2f} s", transform=ax.transAxes, color="white",
                family="monospace", va="top")
        fig.savefig(os.path.join(frame_dir, f"f{i:05d}.png"), dpi=100, facecolor="black")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True, help="top-down phone video of the coffee")
    ap.add_argument("--start", type=float, default=0.0, help="seconds into the video when the frother starts")
    ap.add_argument("--duration", type=float, default=8.0)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--csv", default=os.path.join(ROOT, "cfd", "post", "surface_velocity"))
    ap.add_argument("--vmax", type=float, default=3.0, help="color scale top [m/s]")
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "compare.mp4"))
    a = ap.parse_args()
    ff = shutil.which("ffmpeg")
    if not ff:
        raise SystemExit("ffmpeg not found on PATH")
    tmp = tempfile.mkdtemp()
    render_sim_frames(a.csv, tmp, a.fps, a.duration, a.vmax)
    sim = os.path.join(tmp, "sim.mp4")
    subprocess.run([ff, "-y", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%05d.png"),
                    "-pix_fmt", "yuv420p", sim], check=True)
    # real clip: trim to t = 0, square-crop the center, scale to 600 px, then stack left/right
    flt = (f"[0:v]trim=start={a.start}:duration={a.duration},setpts=PTS-STARTPTS,fps={a.fps},"
           f"crop='min(iw,ih)':'min(iw,ih)',scale=600:600,drawtext=fontfile='C\\:/Windows/Fonts/arial.ttf':"
           f"text='my coffee':x=10:y=10:fontcolor=white:fontsize=24[l];[1:v]scale=600:600[r];"
           f"[l][r]hstack=inputs=2[v]")
    subprocess.run([ff, "-y", "-i", a.real, "-i", sim, "-filter_complex", flt, "-map", "[v]",
                    "-pix_fmt", "yuv420p", a.out], check=True)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
