"""ANSA batch script: gmsh surface mesh -> prism layers + tetra volume -> Fluent .msh.

Driven by cfd/mesh/ansa_job.json (written by cfd/scripts/make_mesh.py). Run with run_ansa.ps1.
"""
import json
import math
import traceback

import ansa
from ansa import base, batchmesh, constants, mesh, session

JOB = r"C:\Users\aalsayid\downloads\per\coffee-dynamics\cfd\mesh\ansa_job.json"


def log(f, msg):
    f.write(msg + "\n")
    f.flush()


def centroid(d, elem):
    c = base.Cog(elem)
    return c if c else (0.0, 0.0, 0.0)


def main():
    job = json.load(open(JOB))
    f = open(job["log"], "w")
    try:
        run(job, f)
    except Exception:
        log(f, "FAILED\n" + traceback.format_exc())
    f.close()
    session.Quit()


def run(job, f):
    session.New("discard")
    base.SetCurrentDeck(constants.FLUENT)
    d = base.CurrentDeck()
    base.InputNastran(job["surface_nas"], properties_id="keep-new")
    shells = base.CollectEntities(d, None, "SHELL")
    log(f, "imported shells: %d" % len(shells))

    props = {p._id: p for p in base.CollectEntities(d, None, "SHELL_PROPERTY")}
    by_name = {}
    for name, pid in job["pid"].items():
        p = props.get(pid)
        if p is None:
            log(f, "missing PID %d (%s)" % (pid, name))
            continue
        # undefined properties are silently dropped by OutputFluent
        base.SetEntityCardValues(d, p, {"Name": name, "ZONE_TYPE": job["zone_type"].get(name, "wall"),
                                        "DEFINED": "YES"})
        by_name[name] = p
    log(f, "properties: %s" % sorted(by_name))

    # layers grow from the grey side, so orient every wall towards the fluid first
    vols = mesh.VolumesDetect(1, return_volumes=True, include_facets=False, whole_db=True)
    log(f, "pre-layer volumes detected: %d" % (len(vols) if vols else 0))
    for v in vols or []:
        log(f, "  orient volume %d -> %s" % (v._id, mesh.OrientVolume(v)))
    if vols:
        mesh.VolumesDelete(vols)

    # boundary layers on no-slip walls only (not on the rigid lid / atmosphere or the MRF interface)
    wall_props = [by_name[n] for n in job["layer_pids"] if n in by_name]
    if wall_props and job.get("layers", True):
        ls = batchmesh.GetNewLayersScenario("layers", "PIDS")
        batchmesh.AddPartToMeshingScenario(wall_props, ls)
        for s in base.CollectEntities(d, None, "BATCH_MESH_LAYERS_SESSION"):
            log(f, "layers params read: %s" % batchmesh.ReadSessionMeshParams(s, job["layers_mpar"]))
        r = batchmesh.RunMeshingScenario(ls, job.get("time_limit_min", 120))
        log(f, "layers run: %s" % r)
        solids = base.CollectEntities(d, None, "SOLID")
        log(f, "solids after layers: %d" % len(solids))

    vs = batchmesh.GetNewVolumeScenario("volume", True, "one_part_for_volume", False)
    for s in base.CollectEntities(d, None, "BATCH_MESH_VOLUME_SESSION"):
        log(f, "volume params read: %s" % batchmesh.ReadSessionMeshParams(s, job["volume_mpar"]))
    r = batchmesh.RunMeshingScenario(vs, job.get("time_limit_min", 120))
    log(f, "volume run: %s" % r)

    solids = base.CollectEntities(d, None, "SOLID")
    log(f, "total solids: %d" % len(solids))
    vols = base.CollectEntities(d, None, "VOLUME")
    cards = [base.GetEntityCardValues(d, v, ["PID", "Volume", "Tetras", "Pentas", "Pyramids", "Grids"])
             for v in vols]
    for c in cards:
        log(f, "volume %s" % c)
    tet_vols = sorted([c for c in cards if c["Tetras"] > 0], key=lambda c: c["Volume"])
    rotor_pid, stator_pid = tet_vols[0]["PID"], tet_vols[-1]["PID"]
    sprops = {p._id: p for p in base.CollectEntities(d, None, "SOLID_PROPERTY")}
    base.SetEntityCardValues(d, sprops[rotor_pid], {"Name": "rotor", "ZONE_TYPE": "fluid", "DEFINED": "YES"})
    base.SetEntityCardValues(d, sprops[stator_pid], {"Name": "stator", "ZONE_TYPE": "fluid", "DEFINED": "YES"})

    # ANSA puts all layer prisms in one property; split them between the rotating and fixed regions
    rc = job["rotor_box"]
    moved = {"rotor": 0, "stator": 0}
    layer_pids = [c["PID"] for c in cards if c["PID"] not in (rotor_pid, stator_pid)]
    by_vol = {}
    for v in vols:
        pid = base.GetEntityCardValues(d, v, ["PID"])["PID"]
        if pid in layer_pids:
            by_vol[pid] = base.CollectEntities(d, v, "SOLID")
    if layer_pids and not any(by_vol.values()):
        by_vol = {pid: [] for pid in layer_pids}
        for e in base.CollectEntities(d, None, "SOLID"):
            pid = base.GetEntityCardValues(d, e, ["PID"])["PID"]
            if pid in by_vol:
                by_vol[pid].append(e)
    for c in cards:
        if c["PID"] not in layer_pids:
            continue
        els = by_vol.get(c["PID"], [])
        log(f, "layer property %d: %d solids" % (c["PID"], len(els)))
        to_rotor, to_stator = [], []
        for e in els:
            x, y, z = centroid(d, e)
            if math.hypot(x, y) < rc["r"] and rc["z0"] < z < rc["z1"]:
                to_rotor.append(e)
            else:
                to_stator.append(e)
        for e in to_rotor:
            base.SetEntityCardValues(d, e, {"PID": rotor_pid})
        for e in to_stator:
            base.SetEntityCardValues(d, e, {"PID": stator_pid})
        moved["rotor"] += len(to_rotor)
        moved["stator"] += len(to_stator)
    log(f, "layer prisms reassigned: %s" % moved)
    for pid, sp in sprops.items():
        n = base.GetEntityCardValues(d, sp, ["Name", "Num.Elem"])
        log(f, "solid property %d: %s" % (pid, n))

    base.SaveAs(job["ansa_db"])
    ok = base.OutputFluent(job["fluent_msh"], mode="all", format="ascii")
    log(f, "fluent output: %s -> %s" % (ok, job["fluent_msh"]))
    # Nastran keeps element->node connectivity and PIDs; our own converter writes AcuSolve files from it
    ok = base.OutputNastran(job["nastran"], mode="all", format="long", write_comments="off")
    log(f, "nastran output: %s -> %s" % (ok, job["nastran"]))
    log(f, "DONE")
