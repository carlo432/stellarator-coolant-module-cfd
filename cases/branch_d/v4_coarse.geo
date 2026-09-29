// Baseline module v4 offset-recirculation geometry
// Units: meters.
//
// Branch D geometry. This keeps the v3
// chamber, inlet, heated wall, depth, and patch names, but shifts the outlet to
// the upper side of the chamber so the lower heated wall is more likely to
// contain a low-flow/recirculation region. A backward-facing-step variant (v6)
// is the alternative geometry.

SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.0025;
Mesh.CharacteristicLengthMax = 0.0090;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

inlet_length = 0.10;
chamber_length = 0.20;
outlet_length = 0.10;

inlet_width = 0.02;
chamber_width = 0.08;
outlet_width = 0.02;
depth = 0.04;

lc_base = 0.0090;
lc_hot = 0.0025;

outlet_center_y = chamber_width/2 - outlet_width/2;
outlet_lower_y = outlet_center_y - outlet_width/2;
outlet_upper_y = outlet_center_y + outlet_width/2;

// 2D footprint points in x-y plane. The inlet stays centered at y = 0.
// The outlet is shifted upward but keeps the same area as the baseline outlet.
p1 = newp; Point(p1) = {0.00, -inlet_width/2, 0, lc_base};
p2 = newp; Point(p2) = {inlet_length, -inlet_width/2, 0, lc_base};
p3 = newp; Point(p3) = {inlet_length, -chamber_width/2, 0, lc_hot};
p4 = newp; Point(p4) = {inlet_length + chamber_length, -chamber_width/2, 0, lc_hot};
p5 = newp; Point(p5) = {inlet_length + chamber_length, outlet_lower_y, 0, lc_base};
p6 = newp; Point(p6) = {inlet_length + chamber_length + outlet_length, outlet_lower_y, 0, lc_base};
p7 = newp; Point(p7) = {inlet_length + chamber_length + outlet_length, outlet_upper_y, 0, lc_base};
p8 = newp; Point(p8) = {inlet_length + chamber_length, outlet_upper_y, 0, lc_base};
p9 = newp; Point(p9) = {inlet_length, chamber_width/2, 0, lc_base};
p10 = newp; Point(p10) = {inlet_length, inlet_width/2, 0, lc_base};
p11 = newp; Point(p11) = {0.00, inlet_width/2, 0, lc_base};

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
l11 = newl; Line(l11) = {p11, p1};

// Refine around the heated chamber wall before extrusion.
Field[1] = Distance;
Field[1].EdgesList = {l3};
Field[1].NNodesByEdge = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_hot;
Field[2].SizeMax = lc_base;
Field[2].DistMin = 0.004;
Field[2].DistMax = 0.025;
Background Field = 2;

loop1 = newll; Curve Loop(loop1) = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11};
s1 = news; Plane Surface(s1) = {loop1};

out[] = Extrude {0, 0, depth} {
  Surface{s1};
  Layers{16};
  Recombine;
};

// Extrusion return surfaces follow the curve order:
// out[2] from l1, out[3] from l2, ..., out[12] from l11.
// Inlet cap is generated from l11; outlet cap is generated from l6.
Physical Surface("inlet") = {out[12]};
Physical Surface("outlet") = {out[7]};

// Heated wall choice: same wide lower chamber wall as v3.
Physical Surface("heated_wall") = {out[4]};

// Include the original source surface s1 and the extruded cap out[0], plus
// unheated side walls. Without s1, OpenFOAM receives undefined default faces.
Physical Surface("adiabatic_walls") = {s1, out[0], out[2], out[3], out[5], out[6], out[8], out[9], out[10], out[11]};
Physical Volume("fluid") = {out[1]};
