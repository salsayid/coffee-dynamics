"""Build the fluid-domain CAD for Levels 2 and 3 with gmsh/OpenCASCADE and export STEP (mm).

  variant A: frother head as a solid torus envelope + spokes + shaft
  variant B: resolved helical spring coil + support ring + spokes + shaft
  level 2: liquid only, flat rigid lid at the liquid surface
  level 3: liquid + air up to the mug rim
The rotating (MRF) cylinder is fragmented into its own volume so the interface is conformal.
"""
import argparse
import json
import math
import os
import sys

import gmsh

sys.path.insert(0, os.path.dirname(__file__))
import params as P  # noqa: E402

MM = 1000.0


def split_torus(occ, zc, ring, minor):
    """Torus as 8 patches (4 quadrants x upper/lower) so no face is periodic in both directions;
    ANSA's CAD import cannot mesh a doubly periodic toroidal face."""
    pieces = []
    big = ring + minor + 1.0
    for qx, qy in ((0, 0), (-1, 0), (-1, -1), (0, -1)):
        for qz in (0, -1):
            t = occ.addTorus(0, 0, zc, ring, minor)
            box = occ.addBox(qx * big, qy * big, zc + qz * (minor + 1.0), big, big, minor + 1.0)
            out, _ = occ.intersect([(3, t)], [(3, box)])
            pieces += out
    return pieces


def frother(occ, variant):
    zc = P.HEAD_Z * MM
    ring = P.COIL_RING_R * MM
    tools = []
    if variant == "A":
        tools += split_torus(occ, zc, ring, P.COIL_D / 2 * MM)
    else:
        tools += split_torus(occ, zc, ring, P.SUPPORT_WIRE_D / 2 * MM)
        a = (P.COIL_D - P.COIL_WIRE_D) / 2 * MM
        rw = P.COIL_WIRE_D / 2 * MM
        turns = int(round(2 * math.pi * ring / (P.COIL_PITCH * MM)))
        npts = turns * 16
        pts = []
        for i in range(npts):
            phi = 2 * math.pi * i / npts
            rr = ring + a * math.cos(turns * phi)
            pts.append(occ.addPoint(rr * math.cos(phi), rr * math.sin(phi), zc + a * math.sin(turns * phi)))
        pts.append(pts[0])
        spl = occ.addBSpline(pts)
        wire = occ.addWire([spl])
        # profile circle at phi=0, normal to the path tangent there
        tx, ty, tz = 0.0, ring + a, a * turns  # d/dphi of the path at phi=0 (unnormalised)
        n = math.sqrt(tx * tx + ty * ty + tz * tz)
        disk = occ.addDisk(ring + a, 0, zc, rw, rw, zAxis=[tx / n, ty / n, tz / n])
        pipe = occ.addPipe([(2, disk)], wire, "DiscreteTrihedron")
        tools += [e for e in pipe if e[0] == 3]
    r_shaft = P.SHAFT_D / 2 * MM
    spoke_r = P.SPOKE_D / 2 * MM
    for k in range(P.SPOKE_N):
        # offset keeps spokes off the torus parametric seam at phi = 0, which breaks the OCC boolean
        th = 2 * math.pi * k / P.SPOKE_N + math.radians(30)
        x0, y0 = 0.5 * r_shaft * math.cos(th), 0.5 * r_shaft * math.sin(th)
        L = ring - 0.5 * r_shaft
        tools.append((3, occ.addCylinder(x0, y0, zc, L * math.cos(th), L * math.sin(th), 0, spoke_r)))
    return tools, r_shaft, zc


def build(variant, level, out_dir, surface_mesh=True, name_tag="", h_head=None, h_max=1.5, h_mrf=0.5):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"mug_{variant}_L{level}")
    occ = gmsh.model.occ

    top = (P.LIQUID_H if level == 2 else P.MUG_H) * MM
    r_top = P.mug_radius(top / MM) * MM
    mug = occ.addCone(0, 0, 0, 0, 0, top, P.MUG_R_BOTTOM * MM, r_top)

    tools, r_shaft, zc = frother(occ, variant)
    # shaft starts slightly below the head plane so the spoke roots sit inside it (a small hub)
    hub = 0.8
    shaft = occ.addCylinder(0, 0, zc - hub, 0, 0, top - zc + hub + 1.0, r_shaft)
    tools.append((3, shaft))
    # fusing the head first breaks OCC (torus + spoke booleans fail at most angles); one multi-tool cut is robust
    v_mug = occ.getMass(3, mug)
    fluid, _ = occ.cut([(3, mug)], tools)
    occ.synchronize()
    v_fluid = sum(occ.getMass(*e) for e in fluid)

    mrf = occ.addCylinder(0, 0, P.MRF_Z0 * MM, 0, 0, (P.MRF_Z1 - P.MRF_Z0) * MM, P.MRF_R * MM)
    frag, fmap = occ.fragment(fluid, [(3, mrf)])
    keep = {e for e in fmap[0]}
    extra = [e for e in frag if e not in keep]
    if extra:
        occ.remove(extra, recursive=True)
    occ.synchronize()

    vols = gmsh.model.getEntities(3)
    info = {"variant": variant, "level": level, "units": "mm", "volumes": [],
            "frother_volume_removed_mm3": v_mug - v_fluid}
    for dim, tag in vols:
        cx, cy, cz = occ.getCenterOfMass(dim, tag)
        info["volumes"].append({"tag": tag, "volume_mm3": occ.getMass(dim, tag), "com": [cx, cy, cz]})
    name = f"mug_{variant}_L{level}{name_tag}"
    step = os.path.join(out_dir, name + ".step")
    gmsh.write(step)
    info["step"] = step
    info["faces"] = len(gmsh.model.getEntities(2))

    groups = classify(occ, level, top)
    info["groups"] = {k: len(v) for k, v in groups.items()}
    if surface_mesh:
        info.update(mesh_surfaces(groups, vols, out_dir, name, variant, h_head, h_max, h_mrf))
    json.dump(info, open(os.path.join(out_dir, name + ".json"), "w"), indent=2)
    gmsh.finalize()
    return info


def classify(occ, level, top):
    """Name every CAD face from its geometric type and position (mm)."""
    zc = P.HEAD_Z * MM
    mrf_r, z0, z1 = P.MRF_R * MM, P.MRF_Z0 * MM, P.MRF_Z1 * MM
    rs = P.SHAFT_D / 2 * MM
    g = {k: [] for k in ("frother_head", "shaft_rotor", "shaft_stator", "mug_wall", "mug_bottom",
                         "top", "mrf_side", "mrf_top", "mrf_bottom")}
    for dim, tag in gmsh.model.getEntities(2):
        typ = gmsh.model.getType(dim, tag)
        x, y, z = occ.getCenterOfMass(dim, tag)
        bb = occ.getBoundingBox(dim, tag)
        rmax = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]))
        if typ == "Cone":
            g["mug_wall"].append(tag)
        elif typ == "Plane" and abs(z) < 1e-6:
            g["mug_bottom"].append(tag)
        elif typ == "Plane" and abs(z - top) < 1e-6:
            g["top"].append(tag)
        elif typ == "Plane" and abs(z - z1) < 1e-6:
            g["mrf_top"].append(tag)
        elif typ == "Plane" and abs(z - z0) < 1e-6:
            g["mrf_bottom"].append(tag)
        elif typ == "Cylinder" and abs(rmax - mrf_r) < 1e-3:
            g["mrf_side"].append(tag)
        elif typ == "Cylinder" and abs(rmax - rs) < 1e-3 and bb[2] >= z1 - 1e-6:
            g["shaft_stator"].append(tag)
        elif typ == "Cylinder" and abs(rmax - rs) < 1e-3:
            g["shaft_rotor"].append(tag)
        else:
            g["frother_head"].append(tag)  # torus patches, spokes, hub cap, coil
    return {k: v for k, v in g.items() if v}


def mesh_surfaces(groups, vols, out_dir, name, variant, h_head=None, h_max=1.5, h_mrf=0.5):
    """Surface-only triangulation with physical groups, written as Nastran (PID = group id)."""
    head = groups["frother_head"] + groups.get("shaft_rotor", [])
    if h_head is None:
        h_head = 0.18 if variant == "A" else 0.12
    h_min = 0.06 if variant == "B" else min(0.12, h_head)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    f = gmsh.model.mesh.field
    dist = f.add("Distance")
    f.setNumbers(dist, "SurfacesList", head)
    f.setNumber(dist, "Sampling", 60)
    th = f.add("Threshold")
    f.setNumber(th, "InField", dist)
    f.setNumber(th, "SizeMin", h_head)
    f.setNumber(th, "SizeMax", h_max)
    f.setNumber(th, "DistMin", 0.3)
    f.setNumber(th, "DistMax", 6.0)
    mrf = f.add("Cylinder")
    f.setNumber(mrf, "Radius", P.MRF_R * MM + 1.0)
    f.setNumber(mrf, "VIn", h_mrf)
    f.setNumber(mrf, "VOut", h_max)
    f.setNumber(mrf, "XAxis", 0)
    f.setNumber(mrf, "YAxis", 0)
    f.setNumber(mrf, "ZAxis", (P.MRF_Z1 - P.MRF_Z0) * MM + 2.0)
    f.setNumber(mrf, "XCenter", 0)
    f.setNumber(mrf, "YCenter", 0)
    f.setNumber(mrf, "ZCenter", P.HEAD_Z * MM)
    mn = f.add("Min")
    f.setNumbers(mn, "FieldsList", [th, mrf])
    f.setAsBackgroundMesh(mn)

    ids = {}
    for i, (k, tags) in enumerate(sorted(groups.items()), start=1):
        gmsh.model.addPhysicalGroup(2, tags, i, name=k)
        ids[k] = i
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.SaveElementTagType", 2)
    gmsh.option.setNumber("Mesh.BdfFieldFormat", 2)
    bdf = os.path.join(out_dir, name + "_surf.nas")
    gmsh.write(bdf)
    q = gmsh.model.mesh.getElementQualities(gmsh.model.mesh.getElementsByType(2)[0], "minSICN")
    ntri = len(gmsh.model.mesh.getElementsByType(2)[0])
    return {"surface_mesh": bdf, "pid": ids, "n_tri": int(ntri), "min_sicn": float(min(q))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="A", choices=["A", "B"])
    ap.add_argument("--level", type=int, default=2, choices=[2, 3])
    a = ap.parse_args()
    out = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "mesh"))
    os.makedirs(out, exist_ok=True)
    print(json.dumps(build(a.variant, a.level, out), indent=1))
