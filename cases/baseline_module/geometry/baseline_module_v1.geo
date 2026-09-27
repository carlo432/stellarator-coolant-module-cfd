// Baseline module v1
// Solver-neutral Gmsh geometry draft.
// Units: meters.
//
// Expanded heated channel:
// narrow inlet -> wider chamber -> narrow outlet.
//
// This file is a starting point. It has not been meshed here because Gmsh
// is not currently installed on this machine.

SetFactory("OpenCASCADE");

inlet_length = 0.10;
chamber_length = 0.20;
outlet_length = 0.10;

inlet_width = 0.02;
chamber_width = 0.08;
outlet_width = 0.02;
depth = 0.04;

lc = 0.01;

// 2D footprint points in x-y plane. Centerline is y = 0.
p1 = newp; Point(p1) = {0.00, -inlet_width/2, 0, lc};
p2 = newp; Point(p2) = {inlet_length, -inlet_width/2, 0, lc};
p3 = newp; Point(p3) = {inlet_length, -chamber_width/2, 0, lc};
p4 = newp; Point(p4) = {inlet_length + chamber_length, -chamber_width/2, 0, lc};
p5 = newp; Point(p5) = {inlet_length + chamber_length, -outlet_width/2, 0, lc};
p6 = newp; Point(p6) = {inlet_length + chamber_length + outlet_length, -outlet_width/2, 0, lc};
p7 = newp; Point(p7) = {inlet_length + chamber_length + outlet_length, outlet_width/2, 0, lc};
p8 = newp; Point(p8) = {inlet_length + chamber_length, outlet_width/2, 0, lc};
p9 = newp; Point(p9) = {inlet_length + chamber_length, chamber_width/2, 0, lc};
p10 = newp; Point(p10) = {inlet_length, chamber_width/2, 0, lc};
p11 = newp; Point(p11) = {inlet_length, inlet_width/2, 0, lc};
p12 = newp; Point(p12) = {0.00, inlet_width/2, 0, lc};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p5};
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p7};
l7 = newl; Line(l7) = {p7, p8};
l8 = newl; Line(l8) = {p8, p9};
l9 = newl; Line(l9) = {p9, p10};
l10 = newl; Line(l10) = {p10, p11};
l11 = newl; Line(l11) = {p11, p12};
l12 = newl; Line(l12) = {p12, p1};

loop1 = newll; Curve Loop(loop1) = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12};
s1 = news; Plane Surface(s1) = {loop1};

out[] = Extrude {0, 0, depth} {
  Surface{s1};
  Layers{8};
  Recombine;
};

// Extrusion return surfaces follow the original curve order:
// out[2] from l1, out[3] from l2, ..., out[13] from l12.
// Inlet cap is generated from l12; outlet cap is generated from l6.
Physical Surface("inlet") = {out[13]};
Physical Surface("outlet") = {out[7]};

// Heated wall choice: the wide lower chamber wall.
// Surface IDs may need verification in Gmsh after first generation.
Physical Surface("heated_wall") = {out[4]};

// Include the original source surface s1 and the extruded cap out[0], plus
// unheated side walls. Without s1, OpenFOAM receives undefined default faces.
Physical Surface("adiabatic_walls") = {s1, out[0], out[2], out[3], out[5], out[6], out[8], out[9], out[10], out[11], out[12]};
Physical Volume("fluid") = {out[1]};
