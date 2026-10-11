"""Minimal EnSight Gold C-binary reader/writer for per-node variables (as written by acuTrans).

Used to add derived fields (radial/azimuthal velocity, vortex criteria, ...) to a case so META can draw them.
"""
import os
import re

import numpy as np

ELEM_NODES = {"point": 1, "bar2": 2, "tria3": 3, "quad4": 4, "tetra4": 4, "pyramid5": 5, "penta6": 6,
              "hexa8": 8, "bar3": 3, "tria6": 6, "quad8": 8, "tetra10": 10, "pyramid13": 13,
              "penta15": 15, "hexa20": 20}


class Reader:
    def __init__(self, path):
        self.b = open(path, "rb").read()
        self.p = 0

    def s(self):
        v = self.b[self.p:self.p + 80].decode("ascii", "replace").strip("\x00 ").strip()
        self.p += 80
        return v

    def i(self, n=1):
        v = np.frombuffer(self.b, dtype="<i4", count=n, offset=self.p)
        self.p += 4 * n
        return v if n > 1 else int(v[0])

    def f(self, n):
        v = np.frombuffer(self.b, dtype="<f4", count=n, offset=self.p)
        self.p += 4 * n
        return v

    def eof(self):
        return self.p >= len(self.b)


def read_geo(path):
    r = Reader(path)
    assert r.s().lower().startswith("c binary")
    r.s(); r.s()
    node_id = r.s().lower()
    elem_id = r.s().lower()
    nid_given = "given" in node_id or "ignore" in node_id
    eid_given = "given" in elem_id or "ignore" in elem_id
    if r.b[r.p:r.p + 6].decode().lower().startswith("extent"):
        r.s(); r.f(6)
    parts = []
    while not r.eof():
        tag = r.s()
        if not tag.lower().startswith("part"):
            break
        pid = r.i()
        name = r.s()
        r.s()  # coordinates
        n = r.i()
        ids = r.i(n) if nid_given else np.arange(1, n + 1)
        x, y, z = r.f(n), r.f(n), r.f(n)
        elems = {}
        while not r.eof():
            save = r.p
            et = r.s().lower()
            if et.startswith("part"):
                r.p = save
                break
            ne = r.i()
            if eid_given:
                r.i(ne)
            k = ELEM_NODES[et]
            elems[et] = r.i(ne * k).reshape(ne, k) if ne * k > 1 else np.array([[r.i()]])
        parts.append({"id": pid, "name": name, "node_ids": np.atleast_1d(ids),
                      "xyz": np.column_stack([x, y, z]), "elems": elems})
    return parts


def write_var(path, parts, values_by_part, desc, vector=False):
    """values_by_part: list of arrays (n,) or (n,3) aligned with parts."""
    with open(path, "wb") as fh:
        fh.write(desc.ljust(80)[:80].encode())
        for part, v in zip(parts, values_by_part):
            fh.write(b"part".ljust(80))
            fh.write(np.int32(part["id"]).tobytes())
            fh.write(b"coordinates".ljust(80))
            v = np.asarray(v, dtype="<f4")
            if vector:
                for c in range(3):
                    fh.write(np.ascontiguousarray(v[:, c]).tobytes())
            else:
                fh.write(v.tobytes())


def read_var(path, parts, vector=False):
    r = Reader(path)
    r.s()
    out = []
    for part in parts:
        r.s(); r.i(); r.s()
        n = len(part["node_ids"])
        if vector:
            out.append(np.column_stack([r.f(n), r.f(n), r.f(n)]))
        else:
            out.append(r.f(n))
    return out


def add_to_case(case_path, name, rel_pattern, vector=False):
    text = open(case_path).read()
    kind = "vector" if vector else "scalar"
    line = f"{kind} per node: {name} {rel_pattern}\n"
    if re.search(rf"per node:\s+{re.escape(name)}\s", text):
        text = re.sub(rf"^.*per node:\s+{re.escape(name)}\s.*\n", line, text, flags=re.M)
    else:
        text = text.replace("VARIABLE\n", "VARIABLE\n" + line)
    open(case_path, "w").write(text)


def case_steps(case_path):
    text = open(case_path).read()
    nums = re.search(r"filename numbers:\s*([\d\s]+?)\n\s*time values", text, re.S)
    nums = [int(x) for x in nums.group(1).split()] if nums else []
    pad = re.search(r"per node:\s+velocity\s+(\S+)", text).group(1).count("*")
    return nums, pad
