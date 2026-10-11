"""Level 1b: set up and run the Taylor-vortex Reynolds sweep in the idealized annulus."""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import acu  # noqa: E402
import params as P  # noqa: E402
from taylor_stability import growth  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
K_C = 3.215
WAVELENGTH = 2 * 3.141592653589793 / K_C
D = P.R2_IDEAL - P.R1
LENGTH = 2 * WAVELENGTH * D


def plan(re):
    U = re * P.NU / D
    s = growth(re, K_C, 0.3)
    efold = abs(1 / s) * D / U
    tdiff = D * D / P.NU
    final = max(14 * efold, 0.5 * tdiff)
    # ~120 steps per e-fold keeps the implicit solve well conditioned (larger steps cost 20x GMRES iterations)
    dt = min(efold / 120, 2 * D / U)
    nsteps = int(final / dt) + 1
    return dt, final, nsteps


def setup(re, nr=40, nz=160):
    case = os.path.join(ROOT, "cfd", "cases", f"L1b_Re{int(re)}")
    dt, final, nsteps = plan(re)
    cmd = [sys.executable, os.path.join(HERE, "l1_case.py"), case, "--nr", str(nr), "--nz", str(nz),
           "--length", str(LENGTH), "--re", str(re), "--transient", "--dt", f"{dt:.6g}",
           "--nsteps", str(nsteps), "--final_time", f"{final:.6g}", "--out_freq", str(max(1, nsteps // 40)),
           "--th_freq", str(max(1, nsteps // 400)), "--perturb", "1e-3", "--tol", "1e-4"]
    subprocess.run(cmd, check=True)
    return case, nsteps


def run_case(re, np_):
    case, nsteps = setup(re)
    t0 = time.time()
    r = acu.run(case, np_)
    ok = "Total CPU" in open(os.path.join(case, "case.1.Log")).read()
    print(f"Re={re}: {'ok' if ok else 'FAILED'} in {time.time() - t0:.0f} s ({nsteps} steps)", flush=True)
    return re, ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("re", type=float, nargs="+")
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--parallel", type=int, default=4)
    a = ap.parse_args()
    with ThreadPoolExecutor(a.parallel) as ex:
        list(ex.map(lambda re: run_case(re, a.np), a.re))
