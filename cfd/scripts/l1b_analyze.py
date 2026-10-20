"""Level 1b: Taylor-vortex onset from the AcuSolve sweep vs. linear stability theory."""
import glob
import json
import os
import re
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
from taylor_stability import growth  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
D = P.R2_IDEAL - P.R1
K_C = 3.215
RE_C_LINEAR = json.load(open(os.path.join(ROOT, "cfd", "post", "taylor_linear_stability.json")))["eta=0.3"]["Re_c"]


def inp_value(text, key):
    m = re.search(rf"{key}\s*=\s*([0-9.eE+-]+)", text)
    return float(m.group(1)) if m else None


def load(case):
    info = dict(l.split("=", 1) for l in open(os.path.join(case, "case_info.txt")))
    info = {k.strip(): v.strip() for k, v in info.items()}
    re_ = float(info["re"])
    omega = float(info["omega"])
    U = omega * P.R1
    text = open(os.path.join(case, "case.inp")).read()
    dt = inp_value(text, "initial_time_increment")
    thf = int(re.search(r'TIME_HISTORY_OUTPUT\( "midgap" \) \{[^}]*output_frequency\s*=\s*(\d+)', text).group(1))
    files = glob.glob(os.path.join(case, "case_set1_node*.oth"))
    if not files and glob.glob(os.path.join(case, "case.*.Log")):
        try:
            acu.trans(case, "-oth", "-othv", "velocity", "-to", "table")
        except RuntimeError as e:
            print(case, "acuTrans failed:", str(e)[-200:])
        files = glob.glob(os.path.join(case, "case_set1_node*.oth"))
    if not files:
        return None
    series = []
    for f in files:
        a = np.loadtxt(f, ndmin=2)
        series.append(a[:, 0])
    n = min(len(s) for s in series)
    ur = np.array([s[:n] for s in series])  # (nodes, time)
    t = np.arange(1, n + 1) * thf * dt
    amp = np.sqrt(2 * np.mean(ur**2, axis=0)) / U
    return dict(re=re_, t=t, amp=amp, U=U, dt=dt, tnd=D / U)


def growth_fit(d):
    a, t = d["amp"], d["t"]
    la = np.log(np.maximum(a, 1e-30))
    afin = a[-1]
    if afin > 5 * a[0]:  # grows: fit between 3x initial and 30% of saturation
        sel = (a > 3 * a[0]) & (a < 0.3 * afin)
    else:  # decays: skip the first tenth
        sel = t > 0.1 * t[-1]
    if sel.sum() < 5:
        return np.nan
    s = np.polyfit(t[sel], la[sel], 1)[0]
    return s * d["tnd"]  # nondimensional growth rate (time unit d/U)


def main():
    cases = sorted(glob.glob(os.path.join(ROOT, "cfd", "cases", "L1b_Re*")),
                   key=lambda c: float(re.findall(r"Re(\d+)", c)[0]))
    data = [d for d in (load(c) for c in cases) if d]
    rows = []
    for d in data:
        sig = growth_fit(d)
        last = d["amp"][-max(3, len(d["amp"]) // 10):]
        saturated = (np.ptp(last) / max(last.mean(), 1e-30)) < 0.03
        rows.append(dict(re=d["re"], amp_final=float(d["amp"][-1]), saturated=bool(saturated),
                         sigma_sim=float(sig), sigma_linear=float(growth(d["re"], K_C, 0.3)),
                         t_end=float(d["t"][-1])))
    for r in rows:
        print(f"Re={r['re']:6.0f}  A_final={r['amp_final']:.4f}  sat={r['saturated']}  "
              f"sigma sim={r['sigma_sim']:+.4f}  linear={r['sigma_linear']:+.4f}")

    sup = [r for r in rows if r["re"] > RE_C_LINEAR and r["re"] <= 100 and r["saturated"]]
    fit = None
    if len(sup) >= 2:
        x = np.array([r["re"] for r in sup])
        y = np.array([r["amp_final"] ** 2 for r in sup])
        m, c = np.polyfit(x, y, 1)
        fit = dict(slope=float(m), intercept=float(c), re_c_sim=float(-c / m), points=[r["re"] for r in sup])
        print(f"A^2 linear fit -> Re_c(sim) = {fit['re_c_sim']:.2f}   linear theory {RE_C_LINEAR:.2f}")
    json.dump({"rows": rows, "amp2_fit": fit, "re_c_linear": RE_C_LINEAR},
              open(os.path.join(ROOT, "cfd", "post", "l1b_taylor.json"), "w"), indent=2)

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cmap = plt.cm.viridis
    for i, d in enumerate(data):
        ax[0].semilogy(d["t"] / d["tnd"], d["amp"], color=cmap(i / max(1, len(data) - 1)), lw=1.6,
                       label=f"Re {d['re']:.0f}")
    ax[0].set_xlabel("time  [gap transit times $d/U$]")
    ax[0].set_ylabel("vortex strength  $u_r^{rms}/U$")
    ax[0].set_title("Vortex growth (or decay) from a tiny kick")
    ax[0].legend(frameon=False, fontsize=8, ncol=2)

    res = np.linspace(55, 160, 200)
    ax[1].plot(res, [growth(r, K_C, 0.3) for r in res], color="#2c5f8a", lw=2, label="Linear stability theory")
    ok = [r for r in rows if np.isfinite(r["sigma_sim"]) and r["re"] <= 160]
    ax[1].plot([r["re"] for r in ok], [r["sigma_sim"] for r in ok], "o", color="#d1495b", ms=7,
               label="AcuSolve")
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].axvline(RE_C_LINEAR, color="gray", ls="--", lw=0.8)
    ax[1].set_xlabel("Re")
    ax[1].set_ylabel("growth rate  [per $d/U$]")
    ax[1].set_title("Growth rate: grows above Re$_c$, decays below")
    ax[1].legend(frameon=False)

    good = [r for r in rows if r["saturated"] and r["amp_final"] > 1e-3]
    ax[2].plot([r["re"] for r in good], [r["amp_final"] ** 2 for r in good], "o", color="#d1495b", ms=7,
               label="AcuSolve (saturated)")
    if fit:
        xx = np.linspace(fit["re_c_sim"], max(fit["points"]) + 5, 50)
        ax[2].plot(xx, fit["slope"] * xx + fit["intercept"], "-", color="#d1495b", lw=1,
                   label=f"fit $\\to$ Re$_c$ = {fit['re_c_sim']:.1f}")
    ax[2].axvline(RE_C_LINEAR, color="#2c5f8a", ls="--", lw=1.2, label=f"theory Re$_c$ = {RE_C_LINEAR:.1f}")
    ax[2].set_xlabel("Re")
    ax[2].set_ylabel("(vortex strength)$^2$")
    ax[2].set_title("Onset: $A^2 \\propto$ (Re - Re$_c$)")
    ax[2].legend(frameon=False, fontsize=9)
    for a in ax:
        a.grid(alpha=0.2)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "figures", "l1b_taylor_onset.png"), dpi=170, bbox_inches="tight")
    print("saved figures/l1b_taylor_onset.png")


if __name__ == "__main__":
    main()
