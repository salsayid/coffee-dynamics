"""Level 2/3 mesh pipeline: gmsh CAD + surface mesh -> ANSA layers + volume -> Nastran -> AcuSolve files."""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(__file__)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "cfd", "geometry"))
import build_cad  # noqa: E402
import mpar  # noqa: E402
import params as P  # noqa: E402

MESH = os.path.join(ROOT, "cfd", "mesh")

# named mesh levels: (frother surface size, max size, MRF zone size, first layer, layers) in mm
LEVELS = {
    "m0": (0.35, 2.5, 0.9, 0.04, 4),
    "m1": (0.18, 1.5, 0.5, 0.03, 6),
    "m2": (0.11, 1.0, 0.33, 0.02, 8),
    "mB": (0.10, 1.5, 0.5, 0.015, 5),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="A")
    ap.add_argument("--level", type=int, default=2)
    ap.add_argument("--mesh", default="m1", choices=sorted(LEVELS))
    ap.add_argument("--growth", type=float, default=1.25)
    ap.add_argument("--no-layers", action="store_true")
    a = ap.parse_args()

    h_head, h_max, h_mrf, first, nlay = LEVELS[a.mesh]
    tag = f"_{a.mesh}"
    name = f"mug_{a.variant}_L{a.level}{tag}"
    t0 = time.time()
    info = build_cad.build(a.variant, a.level, MESH, name_tag=tag, h_head=h_head, h_max=h_max, h_mrf=h_mrf)
    print(f"[{time.time()-t0:5.0f}s] gmsh surface mesh: {info['n_tri']} tris, min SICN {info['min_sicn']:.2f}",
          flush=True)

    layers = mpar.write("layers", os.path.join(MESH, name + "_layers.ansa_mpar"),
                        layers_first_height=f"{first}", layers_growth_factor=f"{a.growth}",
                        Requested_layers=f"{nlay}", Last_aspect="0.6",
                        Problematic_area_treatment="Collapse", retain_valid_layers="true")
    volume = mpar.write("volume", os.path.join(MESH, name + "_volume.ansa_mpar"),
                        create_volume_mesh="true", tetras_algorithm="Tetra Rapid",
                        tetras_max_growth_rate="1.2", tetras_max_elem_length=f"{h_max}")
    job = {
        "surface_nas": info["surface_mesh"],
        "pid": info["pid"],
        "layer_pids": ["frother_head", "shaft_rotor", "mug_wall", "mug_bottom"],
        "zone_type": {"mrf_side": "interior", "mrf_top": "interior", "mrf_bottom": "interior",
                      "liquid_level": "interior",
                      "top": "symmetry" if a.level == 2 else "pressure-outlet"},
        "layers": not a.no_layers,
        "layers_mpar": layers,
        "volume_mpar": volume,
        "rotor_box": {"r": P.MRF_R * 1000, "z0": P.MRF_Z0 * 1000, "z1": P.MRF_Z1 * 1000},
        "ansa_db": os.path.join(MESH, name + ".ansa"),
        "fluent_msh": os.path.join(MESH, name + ".msh"),
        "nastran": os.path.join(MESH, name + "_vol.nas"),
        "log": os.path.join(MESH, name + "_ansa.log"),
        "time_limit_min": 240,
    }
    json.dump(job, open(os.path.join(MESH, "ansa_job.json"), "w"), indent=1)
    ps = os.path.join(os.environ["SystemRoot"], "System32", "WindowsPowerShell", "v1.0", "powershell.exe")
    subprocess.run([ps, "-ExecutionPolicy", "Bypass", "-File", os.path.join(HERE, "run_ansa.ps1"), "-Script",
                    os.path.join(ROOT, "cfd", "geometry", "ansa_volume.py"), "-TimeoutMin", "300"],
                   capture_output=True, text=True)
    log = open(job["log"]).read()
    print(f"[{time.time()-t0:5.0f}s] ANSA finished\n{log}", flush=True)
    if "DONE" not in log:
        sys.exit("ANSA step failed")

    out = os.path.join(MESH, f"acu_{a.variant}_L{a.level}{tag}")
    subprocess.run([sys.executable, os.path.join(HERE, "nas_to_acusolve.py"), job["nastran"],
                    os.path.join(MESH, name + ".json"), out], check=True)
    print(f"[{time.time()-t0:5.0f}s] AcuSolve mesh in {out}", flush=True)


if __name__ == "__main__":
    main()
