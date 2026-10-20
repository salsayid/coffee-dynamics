"""Axisymmetric (5 degree wedge) model of the real mug + torus frother head for transient runs.

  level 2: spin-up from rest, flat rigid liquid surface (slip lid at 40 mm)
  level 3: water + air up to the rim with a level-set free surface and gravity

The torus head is axisymmetric, so it is a rotating wall in a fixed mesh (no MRF needed). The three
thin spokes are not represented. Geometry is built in the x-z plane (mm) and revolved by WEDGE_DEG.
"""
import argparse
import json
import math
import os
import sys

import gmsh
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
import params as P  # noqa: E402

MM = 1000.0
WEDGE_DEG = 5.0
R_AXIS = 0.05  # mm; tiny slip cylinder instead of a degenerate axis


def build_mesh(level, h_wall=0.12, h_far=0.8, first=0.02, nbl=8, out_dir=".", h_band=0.35):
    top = (P.LIQUID_H if level == 2 else P.MUG_H) * MM
    rb, rt = P.MUG_R_BOTTOM * MM, P.mug_radius(top / MM) * MM
    rs, zh = P.SHAFT_D / 2 * MM, P.HEAD_Z * MM - 0.8
    tc, tr, trad = P.COIL_RING_R * MM, P.HEAD_Z * MM, P.COIL_D / 2 * MM

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"axi_L{level}")
    g = gmsh.model.geo
    pts = [g.addPoint(*p, 0) for p in [(R_AXIS, 0), (rb, 0), (rt, top), (rs, top), (rs, zh), (R_AXIS, zh)]]
    pts = [g.addPoint(x, 0, z) for (x, z) in [(R_AXIS, 0), (rb, 0), (rt, top), (rs, top), (rs, zh), (R_AXIS, zh)]]
    names = ["bottom", "mug_wall", "top", "shaft", "hub", "axis"]
    lines = [g.addLine(pts[i], pts[(i + 1) % 6]) for i in range(6)]
    c = g.addPoint(tc, 0, tr)
    q = [g.addPoint(tc + trad * math.cos(a), 0, tr + trad * math.sin(a)) for a in (0, math.pi / 2, math.pi, 1.5 * math.pi)]
    arcs = [g.addCircleArc(q[i], c, q[(i + 1) % 4]) for i in range(4)]
    outer = g.addCurveLoop(lines)
    hole = g.addCurveLoop(arcs)
    surf = g.addPlaneSurface([outer, hole])
    g.synchronize()

    f = gmsh.model.mesh.field
    walls = [lines[0], lines[1], lines[3], lines[4]] + arcs
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", [lines[3], lines[4]] + arcs)
    f.setNumber(d, "Sampling", 200)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", h_wall)
    f.setNumber(t, "SizeMax", h_far)
    f.setNumber(t, "DistMin", 0.5)
    f.setNumber(t, "DistMax", 8.0)
    fl = [t]
    if level == 3:
        # free-surface band and the central column where the vortex dip forms
        band = f.add("Box")
        for k, v in dict(VIn=h_band, VOut=h_far, XMin=0, XMax=rt + 1, YMin=-1, YMax=1,
                         ZMin=P.LIQUID_H * MM - 12, ZMax=P.LIQUID_H * MM + 8, Thickness=4).items():
            f.setNumber(band, k, v)
        core = f.add("Box")
        for k, v in dict(VIn=h_band, VOut=h_far, XMin=0, XMax=10, YMin=-1, YMax=1,
                         ZMin=P.MRF_Z1 * MM, ZMax=P.LIQUID_H * MM + 8, Thickness=4).items():
            f.setNumber(core, k, v)
        fl += [band, core]
    mn = f.add("Min")
    f.setNumbers(mn, "FieldsList", fl)
    f.setAsBackgroundMesh(mn)
    bl = f.add("BoundaryLayer")
    f.setNumbers(bl, "CurvesList", walls)
    f.setNumbers(bl, "PointsList", [])
    f.setNumber(bl, "Size", first)
    f.setNumber(bl, "Ratio", 1.25)
    f.setNumber(bl, "Thickness", first * (1.25**nbl - 1) / 0.25)
    f.setNumber(bl, "Quads", 1)
    f.setAsBoundaryLayer(bl)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)

    # revolve the 2D mesh by hand: one layer of wedges/hexes
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    xz = ncoord.reshape(-1, 3)[:, [0, 2]]
    tag2i = {int(tg): i for i, tg in enumerate(ntags)}
    etypes, _, enodes = gmsh.model.mesh.getElements(2, surf)
    gmsh.finalize()

    n2 = len(xz)
    th = math.radians(WEDGE_DEG)
    xyz0 = np.column_stack([xz[:, 0], np.zeros(n2), xz[:, 1]])
    xyz1 = np.column_stack([xz[:, 0] * math.cos(th), xz[:, 0] * math.sin(th), xz[:, 1]])
    xyz = np.vstack([xyz0, xyz1]) / MM  # m
    elems = {}
    for et, en in zip(etypes, enodes):
        k = 3 if et == 2 else 4
        con = np.array([tag2i[int(v)] for v in en]).reshape(-1, k)
        full = np.hstack([con, con + n2])
        elems["six_node_wedge" if k == 3 else "eight_node_brick"] = full
    return xyz, elems, n2, dict(top=top / MM, rb=rb / MM, rt=rt / MM, rs=rs / MM, zh=zh / MM,
                                tc=tc / MM, tr=tr / MM, trad=trad / MM)


def orient(xyz, shape, con):
    p = xyz[con]
    if shape == "six_node_wedge":
        v = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
        con[v < 0] = con[v < 0][:, [0, 2, 1, 3, 5, 4]]
    else:
        v = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]), p[:, 4] - p[:, 0])
        con[v < 0] = con[v < 0][:, [0, 3, 2, 1, 4, 7, 6, 5]]
    return con


FACES = {"six_node_wedge": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (0, 2, 5, 3)],
         "eight_node_brick": [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (0, 3, 7, 4)]}


def classify_node_sets(xyz, geo):
    r = np.hypot(xyz[:, 0], xyz[:, 1])
    z = xyz[:, 2]
    tol = 2e-7
    return {
        "axis": np.abs(r - R_AXIS / MM) < tol,
        "bottom": np.abs(z) < tol,
        "top": np.abs(z - geo["top"]) < tol,
        "mug_wall": np.abs(r - (geo["rb"] + (geo["rt"] - geo["rb"]) * z / geo["top"])) < 1e-6,
        "shaft": (np.abs(r - geo["rs"]) < 1e-6) & (z >= geo["zh"] - tol),
        "hub": (np.abs(z - geo["zh"]) < tol) & (r <= geo["rs"] + tol),
        "torus": np.abs(np.hypot(r - geo["tc"], z - geo["tr"]) - geo["trad"]) < 2e-6,
    }


def write_mesh(case_dir, xyz, elems, n2, geo):
    mdir = os.path.join(case_dir, "MESH.DIR")
    os.makedirs(mdir, exist_ok=True)
    np.savetxt(os.path.join(mdir, "mesh.crd"), np.column_stack([np.arange(1, len(xyz) + 1), xyz]),
               fmt="%d %.12e %.12e %.12e")
    sets, eid = {}, 1
    allfaces = {}
    for shape, con in elems.items():
        con = orient(xyz, shape, con.copy())
        ids = np.arange(eid, eid + len(con))
        eid += len(con)
        name = "fluid_" + ("wedge" if shape.startswith("six") else "hex")
        np.savetxt(os.path.join(mdir, name + ".cnn"), np.column_stack([ids, con + 1]), fmt="%d")
        sets[name] = {"shape": shape, "file": name + ".cnn", "count": int(len(con)), "region": "stator"}
        for fl in FACES[shape]:
            fn = con[:, fl]
            for e, nodes in zip(ids.tolist(), fn.tolist()):
                key = tuple(sorted(nodes))
                allfaces.setdefault(key, []).append((e, name, nodes))
    nsets = classify_node_sets(xyz, geo)
    # periodic: node i on theta=0 pairs with node i + n2 on theta=WEDGE_DEG
    pbc = np.column_stack([np.arange(1, n2 + 1), np.arange(1, n2 + 1), np.arange(n2 + 1, 2 * n2 + 1)])
    np.savetxt(os.path.join(mdir, "wedge.pbc"), pbc, fmt="%d")
    surf = {}
    unassigned = 0
    for key, owners in allfaces.items():
        if len(owners) != 1:
            continue
        e, setname, nodes = owners[0]
        arr = np.array(nodes)
        if np.all(arr < n2) or np.all(arr >= n2):
            continue  # periodic sides
        bc = next((nm for nm, m in nsets.items() if np.all(m[arr])), None)
        if bc is None:
            unassigned += 1
            continue
        shape = "three_node_triangle" if len(nodes) == 3 else "four_node_quad"
        surf.setdefault((bc, setname, shape), []).append([e] + [n + 1 for n in nodes])
    surfaces, sid = [], 1
    for (bc, setname, shape), rows in sorted(surf.items()):
        a = np.array(rows)
        a = np.column_stack([a[:, 0], np.arange(sid, sid + len(a)), a[:, 1:]])
        sid += len(a)
        fn = f"{bc}__{setname}__{'tri' if shape.startswith('three') else 'quad'}.srf"
        np.savetxt(os.path.join(mdir, fn), a, fmt="%d")
        surfaces.append({"bc": bc, "element_set": setname, "shape": shape, "count": int(len(a)), "file": fn})
    meta = {"element_sets": sets, "surfaces": surfaces, "n_nodes": int(len(xyz)), "unassigned_faces": unassigned,
            "geo": geo, "wedge_deg": WEDGE_DEG}
    json.dump(meta, open(os.path.join(mdir, "mesh_meta.json"), "w"), indent=1)
    return meta


ROTATING = ("torus", "shaft", "hub")


def write_inp(case_dir, meta, level, dt, final_time, out_every, omega=P.OMEGA, stagger=2):
    L = []
    w = L.append
    w(f'ANALYSIS {{\n    title = "axisymmetric L{level}"\n    type = transient\n}}\n')
    eq = "EQUATION {\n    flow = navier_stokes\n    turbulence = shear_stress_transport\n"
    if level == 3:
        eq += f'    multi_field = levelset\n    fields = {{ "water", "air" }}\n'
        eq += f"    gravitational_acceleration = {{ 0, 0, {-P.G} }}\n"
    w(eq + "}\n")
    w(f"AUTO_SOLUTION_STRATEGY {{\n    max_time_steps = {int(round(final_time / dt)) + 1}\n"
      f"    final_time = {final_time}\n    initial_time_increment = {dt}\n    auto_time_increment = off\n"
      f"    max_stagger_iterations = {stagger}\n    convergence_tolerance = 1e-3\n}}\n")
    if level == 3:
        # the auto strategy runs 20 mostly-serial redistancing passes per step, which dominates run time
        w('STAGGER( "redistancing" ) {\n    equation = levelset_redistancing\n    min_stagger_iterations = 2\n'
          '    max_stagger_iterations = 4\n    convergence_tolerance = 0.01\n    lhs_update_frequency = 1\n}\n')
    w(f'DENSITY_MODEL( "water" ) {{\n    density = {P.RHO}\n}}\n')
    w(f'VISCOSITY_MODEL( "water" ) {{\n    viscosity = {P.MU}\n}}\n')
    w('MATERIAL_MODEL( "water" ) {\n    density_model = "water"\n    viscosity_model = "water"\n}\n')
    mat = "water"
    if level == 3:
        w(f'DENSITY_MODEL( "air" ) {{\n    density = {P.RHO_AIR}\n}}\n')
        w(f'VISCOSITY_MODEL( "air" ) {{\n    viscosity = {P.MU_AIR}\n}}\n')
        w('MATERIAL_MODEL( "air" ) {\n    density_model = "air"\n    viscosity_model = "air"\n}\n')
        w('FIELD( "water" ) {\n    material_model = "water"\n}\n')
        w('FIELD( "air" ) {\n    material_model = "air"\n}\n')
        w('FIELD_INTERACTION_MODEL( "free surface" ) {\n    type = levelset\n    fields_1 = { "air" }\n'
          '    fields_2 = { "water" }\n    thickness_type = auto\n}\n')
        w('MULTI_FIELD_MODEL( "water-air" ) {\n    field_interaction_models = { "free surface" }\n'
          '    viscosity_type = linear\n    immiscible_interface_type = smooth\n}\n')
        mat = "water-air"
    w(f'REFERENCE_FRAME( "rotor" ) {{\n    rotation_center = {{ 0, 0, 0 }}\n    angular_velocity = {{ 0, 0, {omega} }}\n}}\n')
    w('COORDINATE {\n    coordinates = Read( "MESH.DIR/mesh.crd" )\n}\n')
    matline = (f'    medium = multi_field\n    multi_field_model = "{mat}"\n' if level == 3
               else f'    material_model = "{mat}"\n')
    for name, s in meta["element_sets"].items():
        w(f'ELEMENT_SET( "{name}" ) {{\n    elements = Read( "MESH.DIR/{s["file"]}" )\n    shape = {s["shape"]}\n'
          f'{matline}}}\n')
    for s in meta["surfaces"]:
        bc = s["bc"]
        tag = f'{bc}__{s["element_set"]}__{"tri" if s["shape"].startswith("three") else "quad"}'
        body = (f'    surfaces = Read( "MESH.DIR/{s["file"]}" )\n    shape = {s["shape"]}\n'
                f'    element_set = "{s["element_set"]}"\n')
        if bc in ROTATING:
            extra = '    type = wall\n    reference_frame = "rotor"\n'
        elif bc in ("mug_wall", "bottom"):
            extra = "    type = wall\n"
        else:  # top lid, axis
            extra = "    type = slip\n"
        w(f'SIMPLE_BOUNDARY_CONDITION( "{tag}" ) {{\n{body}{extra}}}\n')
        if bc in ROTATING or bc in ("mug_wall", "bottom"):
            w(f'SURFACE_OUTPUT( "{tag}" ) {{\n{body}    integrated_output_frequency = 1\n}}\n')
    w('PERIODIC_BOUNDARY_CONDITION( "wedge" ) {\n    variable = all\n    type = axisymmetric\n'
      '    nodal_pairs = Read( "MESH.DIR/wedge.pbc" )\n    rotation_axis = { 0, 0, 0 ; 0, 0, 1 }\n}\n')
    nut0 = 10 * P.NU
    k0 = 1e-6
    w("NODAL_INITIAL_CONDITION( velocity ) {\n    default_values = { 0, 0, 0 }\n}\n")
    if level == 3:
        w('NODAL_INITIAL_CONDITION( pressure ) {\n    type = nodal_values\n    nodal_values = Read( "MESH.DIR/init_pressure.nic" )\n}\n')
        w('NODAL_INITIAL_CONDITION( field ) {\n    field = "water"\n    type = nodal_values\n'
          '    nodal_values = Read( "MESH.DIR/init_water.nic" )\n}\n')
    else:
        w("NODAL_INITIAL_CONDITION( pressure ) {\n    default_value = 0\n}\n")
    w(f"NODAL_INITIAL_CONDITION( kinetic_energy ) {{\n    default_value = {k0}\n}}\n")
    w(f"NODAL_INITIAL_CONDITION( eddy_frequency ) {{\n    default_value = {k0 / nut0:.6e}\n}}\n")
    w(f"NODAL_INITIAL_CONDITION( eddy_viscosity ) {{\n    default_value = {nut0:.6e}\n}}\n")
    w(f"NODAL_OUTPUT {{\n    output_frequency = {out_every}\n}}\n")
    w(f"RESTART_OUTPUT {{\n    output_frequency = {out_every * 5}\n}}\n")
    w("RUN\n")
    open(os.path.join(case_dir, "case.inp"), "w").write("\n".join(L))
    open(os.path.join(case_dir, "Acusim.cnf"), "w").write("problem = case\ninput_file = case.inp\n")


def initial_fields(case_dir, xyz):
    z = xyz[:, 2]
    zs = P.LIQUID_H
    water = (z < zs).astype(float)
    p = np.where(z < zs, P.RHO * P.G * (zs - z), P.RHO_AIR * P.G * (zs - z))
    ids = np.arange(1, len(z) + 1)
    np.savetxt(os.path.join(case_dir, "MESH.DIR", "init_water.nic"), np.column_stack([ids, water]), fmt="%d %.6f")
    np.savetxt(os.path.join(case_dir, "MESH.DIR", "init_pressure.nic"), np.column_stack([ids, p]), fmt="%d %.8e")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--level", type=int, default=2)
    ap.add_argument("--dt", type=float, default=2e-3)
    ap.add_argument("--final_time", type=float, default=10.0)
    ap.add_argument("--out_every", type=int, default=50)
    ap.add_argument("--h_wall", type=float, default=0.12)
    ap.add_argument("--h_far", type=float, default=0.8)
    ap.add_argument("--first", type=float, default=0.02)
    ap.add_argument("--stagger", type=int, default=2)
    ap.add_argument("--h_band", type=float, default=0.35)
    a = ap.parse_args()
    os.makedirs(a.case_dir, exist_ok=True)
    xyz, elems, n2, geo = build_mesh(a.level, a.h_wall, a.h_far, a.first, h_band=a.h_band)
    meta = write_mesh(a.case_dir, xyz, elems, n2, geo)
    if a.level == 3:
        initial_fields(a.case_dir, xyz)
    write_inp(a.case_dir, meta, a.level, a.dt, a.final_time, a.out_every, stagger=a.stagger)
    json.dump(dict(level=a.level, dt=a.dt, final_time=a.final_time, out_every=a.out_every, wedge_deg=WEDGE_DEG,
                   omega=P.OMEGA), open(os.path.join(a.case_dir, "case_info.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != "surfaces"}, indent=1))
    for s in meta["surfaces"]:
        print(f"  {s['bc']:10s} {s['element_set']:12s} {s['shape']:20s} {s['count']}")


if __name__ == "__main__":
    main()
