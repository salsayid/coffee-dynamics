"""Helpers for running AcuSolve and reading its output."""
import glob
import os
import subprocess

import numpy as np

ACU_BIN = r"C:\Program Files\Altair\2019\hwcfdsolvers\acusolve\win64\bin"

# column layout of acuTrans -osi table output
OSI_COLS = {"step": 0, "time": 1, "area": 2, "mass_flux": 3,
            "traction": slice(7, 10), "moment": slice(10, 13)}


def env():
    e = os.environ.copy()
    e["PATH"] = ACU_BIN + os.pathsep + e["PATH"]
    # concurrent runs otherwise pin their OpenMP threads to the same cores
    e["KMP_AFFINITY"] = "disabled"
    e["OMP_PROC_BIND"] = "false"
    return e


def run(case_dir, np_=16, extra=()):
    cmd = ["cmd", "/c", "acuRun.bat", "-np", str(np_), *extra]
    return subprocess.run(cmd, cwd=case_dir, env=env(), capture_output=True, text=True)


def trans(case_dir, *args):
    cmd = [os.path.join(ACU_BIN, "acuTrans.exe"), *args]
    r = subprocess.run(cmd, cwd=case_dir, env=env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout + r.stderr)
    return r.stdout


def nodal(case_dir, var="velocity", step="L", mesh_dir=None):
    """Return (ids, xyz, values) for one nodal variable at one time step."""
    for f in glob.glob(os.path.join(case_dir, "case_step*.out")):
        os.remove(f)
    trans(case_dir, "-out", "-to", "table", "-outv", var, "-ts", str(step))
    out = sorted(glob.glob(os.path.join(case_dir, "case_step*.out")))[-1]
    vals = np.loadtxt(out)
    crd = np.loadtxt(os.path.join(mesh_dir or os.path.join(case_dir, "MESH.DIR"), "mesh.crd"))
    return crd[:, 0].astype(int), crd[:, 1:4], vals


def surface_integrals(case_dir):
    """Dict surface-index -> array of rows from acuTrans -osi."""
    trans(case_dir, "-osi", "-osiv", "_all", "-to", "table")
    res = {}
    for f in sorted(glob.glob(os.path.join(case_dir, "case_srf*.osi"))):
        k = int(os.path.basename(f).split("srf")[1].split(".")[0])
        res[k] = np.atleast_2d(np.loadtxt(f))
    return res


def time_history(case_dir, var="velocity"):
    """Return (times, nodes, values[t, node, comp]) from TIME_HISTORY_OUTPUT."""
    for f in glob.glob(os.path.join(case_dir, "case_*.oth")):
        os.remove(f)
    trans(case_dir, "-oth", "-othv", var, "-to", "table")
    files = sorted(glob.glob(os.path.join(case_dir, "case_*.oth")))
    return files
