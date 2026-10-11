"""Generate and run META batch sessions that render cut-plane contours of an EnSight case to PNG.

Batch META here draws contours only on cut-plane sections (solid faces stay white), so every view
is a section. Each call renders one variable on one or more planes.
"""
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SES_DIR = os.path.join(ROOT, "cfd", "meta")
PS = os.path.join(os.environ.get("SystemRoot", r"C:\WINDOWS"), "System32", "WindowsPowerShell", "v1.0",
                  "powershell.exe")


def q(p):
    return '"' + p.replace("\\", "/") + '"'


def render(name, case_file, var, planes, view, png, srange=None, palette=None, size=(1400, 1000),
           hide_model=True, extra=()):
    """planes: list of (point xyz tuple, normal tuple)."""
    L = [f"read geom Ensight {q(case_file)}",
         f"read onlyfun Ensight {q(case_file)} all {var}",
         "options state 1",
         "grstyle scalarfringe onnode",
         "grstyle scalarfringe enable"]
    if palette:
        L.append(f"color fringebar scalarset {palette}")
    if srange:
        L.append(f"srange window act set {srange[0]},{srange[1]}")
    L.append(f"window size {size[0]},{size[1]}")
    for i, (p, n) in enumerate(planes):
        L.append(f"plane new cut{i} xyz {p[0]},{p[1]},{p[2]} {n[0]},{n[1]},{n[2]}")
    if planes:
        L.append("plane options onlysection enable all")
    if hide_model:
        L.append("erase all")
    L += list(extra)
    L.append(f"view default {view}")
    L.append("view center")
    L.append(f"write png {q(png)}")
    L.append("options quit")
    ses = os.path.join(SES_DIR, name + ".ses")
    os.makedirs(SES_DIR, exist_ok=True)
    open(ses, "w").write("\n".join(L) + "\n")
    r = subprocess.run([PS, "-ExecutionPolicy", "Bypass", "-File", os.path.join(ROOT, "cfd", "scripts", "run_meta.ps1"),
                        "-Session", ses, "-TimeoutMin", "10"], capture_output=True, text=True)
    ok = os.path.exists(png)
    print(f"{name}: {'ok' if ok else 'FAILED'} -> {png}")
    if not ok:
        print(r.stdout[-1500:])
    return ok


if __name__ == "__main__":
    case = sys.argv[1]
    render("test", case, "speed", [((0, 0, 0.010), (0, 0, 1))], "top", os.path.join(ROOT, "cfd", "post", "meta_test.png"))
