"""Convert the ANSA Nastran volume mesh (long format, mm) to AcuSolve native mesh files (m).

Elements are assigned to the rotating or fixed region by centroid, re-oriented to positive volume,
and every boundary shell is matched to its owning element face.
"""
import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "geometry"))
import params as P  # noqa: E402

NNODES = {"CTETRA": 4, "CPENTA": 6, "CPYRAM": 5, "CTRIA3": 3, "CQUAD4": 4}
SHAPE = {"CTETRA": "four_node_tet", "CPENTA": "six_node_wedge", "CPYRAM": "five_node_pyramid"}
SHORT = {"CTETRA": "tet", "CPENTA": "wedge", "CPYRAM": "pyr"}
# local faces per element type (node indices), any orientation
FACES = {
    "CTETRA": [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)],
    "CPENTA": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (0, 2, 5, 3)],
    "CPYRAM": [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (0, 3, 4)],
}


def read_nastran(path):
    def fields(line, long):
        if long:
            return [line[8 + 16 * i: 24 + 16 * i].strip() for i in range(4)]
        return [line[8 + 8 * i: 16 + 8 * i].strip() for i in range(8)]

    grids, elems, pnames = {}, defaultdict(list), {}
    cur, vals = None, []

    def flush():
        if cur is None:
            return
        name = cur
        if name == "GRID":
            grids[int(vals[0])] = (float(vals[2].replace("D", "E")), float(vals[3]), float(vals[4]))
        elif name in NNODES:
            n = NNODES[name]
            elems[name].append([int(vals[0]), int(vals[1])] + [int(v) for v in vals[2:2 + n]])

    with open(path) as fh:
        for line in fh:
            if line.startswith("$") or not line.strip():
                if line.startswith("$ANSA_NAME") or line.startswith("$HMNAME"):
                    pass
                continue
            if line[0] == "*" or line[0] == "+":
                vals += fields(line, True)
                continue
            flush()
            head = line[:8].strip()
            long = head.endswith("*")
            cur = head.rstrip("*")
            vals = fields(line, long) if cur in ("GRID",) or cur in NNODES else []
    flush()
    return grids, elems


def fix_float(s):
    s = s.strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:
        # Nastran shorthand like 1.5-3
        for i in range(1, len(s)):
            if s[i] in "+-" and s[i - 1] not in "eE":
                return float(s[:i] + "e" + s[i:])
        raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("nas")
    ap.add_argument("json_info", help="build_cad json with the pid map")
    ap.add_argument("out_dir")
    a = ap.parse_args()

    info = json.load(open(a.json_info))
    pid_name = {v: k for k, v in info["pid"].items()}

    print("reading", a.nas, flush=True)
    grids = {}
    elems = defaultdict(list)
    cur, vals = None, []

    def flush():
        if cur == "GRID":
            grids[int(vals[0])] = (fix_float(vals[2]), fix_float(vals[3]), fix_float(vals[4]))
        elif cur in NNODES:
            n = NNODES[cur]
            elems[cur].append([int(vals[0]), int(vals[1])] + [int(v) for v in vals[2:2 + n]])

    with open(a.nas) as fh:
        for line in fh:
            c0 = line[0]
            if c0 == "*":
                vals += [line[8 + 16 * i: 24 + 16 * i] for i in range(4)]
                continue
            if c0 == "$":
                continue
            if cur is not None:
                flush()
            head = line[:8].strip()
            cur = head.rstrip("*")
            if cur == "GRID" or cur in NNODES:
                vals = [line[8 + 16 * i: 24 + 16 * i] for i in range(4)] if head.endswith("*") else \
                       [line[8 + 8 * i: 16 + 8 * i] for i in range(8)]
            else:
                cur = None
    if cur is not None:
        flush()
    print(f"nodes {len(grids)}, " + ", ".join(f"{k} {len(v)}" for k, v in elems.items()), flush=True)

    ids = np.array(sorted(grids))
    xyz = np.array([grids[i] for i in ids]) * 1e-3  # mm -> m
    idx = np.full(ids.max() + 1, -1, dtype=np.int64)
    idx[ids] = np.arange(len(ids))

    os.makedirs(a.out_dir, exist_ok=True)
    crd = np.column_stack([np.arange(1, len(ids) + 1), xyz])
    np.savetxt(os.path.join(a.out_dir, "mesh.crd"), crd, fmt="%d %.12e %.12e %.12e")

    rc = (P.MRF_R, P.MRF_Z0, P.MRF_Z1)
    sets = {}
    face_lookup = {3: {}, 4: {}}
    eid_next = 1
    stats = {}
    for etype in ("CTETRA", "CPENTA", "CPYRAM"):
        if not elems.get(etype):
            continue
        arr = np.array(elems[etype], dtype=np.int64)
        conn = idx[arr[:, 2:]]
        p = xyz[conn]
        cen = p.mean(axis=1)
        if etype == "CTETRA":
            vol = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
            flip = vol < 0
            conn[flip] = conn[flip][:, [0, 2, 1, 3]]
        elif etype == "CPENTA":
            vol = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
            flip = vol < 0
            conn[flip] = conn[flip][:, [0, 2, 1, 3, 5, 4]]
        else:
            vol = np.einsum("ij,ij->i", np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]), p[:, 4] - p[:, 0])
            flip = vol < 0
            conn[flip] = conn[flip][:, [0, 3, 2, 1, 4]]
        r = np.hypot(cen[:, 0], cen[:, 1])
        rotor = (r < rc[0]) & (cen[:, 2] > rc[1]) & (cen[:, 2] < rc[2])
        stats[etype] = {"flipped": int(flip.sum()), "rotor": int(rotor.sum()), "stator": int((~rotor).sum())}
        for region, mask in (("rotor", rotor), ("stator", ~rotor)):
            sub = conn[mask]
            if len(sub) == 0:
                continue
            name = f"{region}_{SHORT[etype]}"
            eids = np.arange(eid_next, eid_next + len(sub))
            eid_next += len(sub)
            np.savetxt(os.path.join(a.out_dir, f"{name}.cnn"), np.column_stack([eids, sub + 1]), fmt="%d")
            sets[name] = {"shape": SHAPE[etype], "region": region, "count": int(len(sub)), "file": f"{name}.cnn"}
            for fl in FACES[etype]:
                fn = sub[:, fl]
                key = np.sort(fn, axis=1)
                k = key[:, 0] * (1 << 40) + key[:, 1] * (1 << 20) + key[:, 2]
                d = face_lookup[len(fl)]
                for kk, e, nodes, cc in zip(k.tolist(), eids.tolist(), fn.tolist(), (xyz[sub].mean(axis=1)).tolist()):
                    d.setdefault(kk, []).append((e, name, nodes, cc))
    print("elements:", json.dumps(stats), flush=True)

    surf = defaultdict(list)
    missing = defaultdict(int)
    for etype in ("CTRIA3", "CQUAD4"):
        for row in elems.get(etype, []):
            pid = row[1]
            bc = pid_name.get(pid, f"pid{pid}")
            if bc.startswith("mrf_") or bc == "liquid_level":
                continue
            nodes = idx[np.array(row[2:])].tolist()
            key = sorted(nodes)
            kk = key[0] * (1 << 40) + key[1] * (1 << 20) + key[2]
            owners = face_lookup[len(nodes)].get(kk, [])
            owners = [o for o in owners if sorted(o[2]) == key]
            if len(owners) != 1:
                missing[(bc, len(owners))] += 1
                continue
            e, setname, enodes, cen = owners[0]
            pts = xyz[nodes]
            nrm = np.cross(pts[1] - pts[0], pts[2] - pts[0])
            if np.dot(nrm, np.array(cen) - pts[0]) > 0:  # make the normal point out of the fluid element
                nodes = nodes[::-1]
            shape = "three_node_triangle" if len(nodes) == 3 else "four_node_quad"
            surf[(bc, setname, shape)].append([e] + [n + 1 for n in nodes])

    surfaces = []
    sid = 1
    for (bc, setname, shape), rows in sorted(surf.items()):
        arr = np.array(rows, dtype=np.int64)
        arr = np.column_stack([arr[:, 0], np.arange(sid, sid + len(arr)), arr[:, 1:]])
        sid += len(arr)
        fname = f"{bc}__{setname}__{'tri' if shape.startswith('three') else 'quad'}.srf"
        np.savetxt(os.path.join(a.out_dir, fname), arr, fmt="%d")
        surfaces.append({"bc": bc, "element_set": setname, "shape": shape, "count": int(len(arr)), "file": fname})
    meta = {"element_sets": sets, "surfaces": surfaces, "n_nodes": int(len(ids)),
            "missing_owner": {f"{k[0]}:{k[1]}": v for k, v in missing.items()}, "stats": stats}
    json.dump(meta, open(os.path.join(a.out_dir, "mesh_meta.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != "surfaces"}, indent=1))
    for s in surfaces:
        print(f"  {s['bc']:14s} {s['element_set']:14s} {s['shape']:20s} {s['count']}")


if __name__ == "__main__":
    main()
