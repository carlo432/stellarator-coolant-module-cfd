// Baseline module v5 — REAL 90-degree bend (Branch D, geometry-robustness test)
// Units: meters.
//
// Constant-width (0.02) channel making a sharp 90-deg turn. Flow enters +x along the
// horizontal leg, turns up the vertical leg, exits at top. The HEATED WALL is the INNER
// wall of the bend (x=0.10, y 0.02->0.22), where the flow separates just past the corner
// -> a recirculation/low-flow zone over the heated wall by a DIFFERENT mechanism than the
// v4 outlet-offset chamber. Tests whether recirculation->hotspot is geometry-robust.
//
// Invariants held vs v13/v14 + v22/v23: inlet 0.02 x 0.04 -> Re 10,000, mdot 1.80 kg/s,
// heated-wall area 0.20 x 0.04 = 0.0080 m^2, depth 0.04. Patch names: inlet/outlet/heated_wall/adiabatic_walls.

SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.0015;
Mesh.CharacteristicLengthMax = 0.0060;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

w       = 0.02;    // channel width
hlen    = 0.10;    // horizontal leg length (inlet face to inner corner)
vlen    = 0.20;    // vertical leg length (= heated wall length)
depth   = 0.04;    // z extrusion

lc_base = 0.0060;
lc_hot  = 0.0015;

// L-shaped footprint (CCW). Horizontal leg y:[0,0.02] x:[0,0.12]; vertical leg x:[0.10,0.12] y:[0.02,0.22].
p1 = newp; Point(p1) = {0.00,      0.00,      0, lc_base};  // inlet bottom-left
p2 = newp; Point(p2) = {hlen+w,    0.00,      0, lc_base};  // outer bottom corner
p3 = newp; Point(p3) = {hlen+w,    w+vlen,    0, lc_base};  // outer top (outlet right)
p4 = newp; Point(p4) = {hlen,      w+vlen,    0, lc_base};  // outlet left
p5 = newp; Point(p5) = {hlen,      w,         0, lc_hot };  // INNER corner (separation point)
p6 = newp; Point(p6) = {0.00,      w,         0, lc_base};  // inlet top-left

l1 = newl; Line(l1) = {p1, p2};   // bottom wall (adiabatic)
l2 = newl; Line(l2) = {p2, p3};   // outer wall (adiabatic)
l3 = newl; Line(l3) = {p3, p4};   // OUTLET (top)
l4 = newl; Line(l4) = {p4, p5};   // HEATED WALL (inner wall of bend, length 0.20)
l5 = newl; Line(l5) = {p5, p6};   // horizontal-leg top wall (adiabatic)
l6 = newl; Line(l6) = {p6, p1};   // INLET (x=0)

// Refine near the heated inner wall (l4) and the inner corner (l5) where separation occurs.
Field[1] = Distance;
Field[1].EdgesList = {l4, l5};
Field[1].NNodesByEdge = 250;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_hot;
Field[2].SizeMax = lc_base;
Field[2].DistMin = 0.004;
Field[2].DistMax = 0.025;
Background Field = 2;

loop1 = newll; Curve Loop(loop1) = {l1, l2, l3, l4, l5, l6};
s1 = news; Plane Surface(s1) = {loop1};

out[] = Extrude {0, 0, depth} { Surface{s1}; Layers{24}; Recombine; };
// l_k -> out[k+1]: l1->out[2] l2->out[3] l3->out[4] l4->out[5] l5->out[6] l6->out[7]
Physical Surface("inlet")  = {out[7]};       // l6
Physical Surface("outlet") = {out[4]};       // l3
Physical Surface("heated_wall") = {out[5]};  // l4
Physical Surface("adiabatic_walls") = {s1, out[0], out[2], out[3], out[6]};
Physical Volume("fluid") = {out[1]};
