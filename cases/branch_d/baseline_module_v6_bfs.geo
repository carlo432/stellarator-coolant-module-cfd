// Baseline module v6 — backward-facing step (BFS), predictive 4th point for the flushing-speed law.
// Units: meters. Inlet enters at the TOP of a height-0.04 channel; bottom drops by a step of 0.02 at
// x=0.10 -> strong, broad recirculation over the heated bottom wall -> VERY slow near-wall flow.
// The flushing-speed law predicts this should give a HIGH peak wall T. Test of predictiveness.
// Invariants: inlet 0.02 x 0.04 -> Re 10,000, mdot 1.80 kg/s, heated wall 0.20 x 0.04 = 0.0080 m^2, depth 0.04.

SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.0015;
Mesh.CharacteristicLengthMax = 0.0060;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

depth = 0.04; lc_base = 0.0060; lc_hot = 0.0015;

p1 = newp; Point(p1) = {0.00,  0.00, 0, lc_base};  // inlet duct bottom-left
p2 = newp; Point(p2) = {0.10,  0.00, 0, lc_hot };  // step top
p3 = newp; Point(p3) = {0.10, -0.02, 0, lc_hot };  // step bottom = channel floor start
p4 = newp; Point(p4) = {0.30, -0.02, 0, lc_hot };  // heated-wall end
p5 = newp; Point(p5) = {0.40, -0.02, 0, lc_base};  // channel floor end
p6 = newp; Point(p6) = {0.40,  0.02, 0, lc_base};  // outlet bottom-right
p7 = newp; Point(p7) = {0.00,  0.02, 0, lc_base};  // top-left

l1 = newl; Line(l1) = {p1, p2};  // inlet duct floor (adiabatic)
l2 = newl; Line(l2) = {p2, p3};  // backward-facing step (adiabatic)
l3 = newl; Line(l3) = {p3, p4};  // HEATED WALL (0.20)
l4 = newl; Line(l4) = {p4, p5};  // downstream floor (adiabatic)
l5 = newl; Line(l5) = {p5, p6};  // OUTLET
l6 = newl; Line(l6) = {p6, p7};  // ceiling (adiabatic)
l7 = newl; Line(l7) = {p7, p1};  // INLET

Field[1] = Distance; Field[1].EdgesList = {l2, l3}; Field[1].NNodesByEdge = 250;
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = lc_hot; Field[2].SizeMax = lc_base;
Field[2].DistMin = 0.004; Field[2].DistMax = 0.025; Background Field = 2;

loop1 = newll; Curve Loop(loop1) = {l1, l2, l3, l4, l5, l6, l7};
s1 = news; Plane Surface(s1) = {loop1};
out[] = Extrude {0, 0, depth} { Surface{s1}; Layers{24}; Recombine; };
// l_k -> out[k+1]
Physical Surface("inlet")  = {out[8]};       // l7
Physical Surface("outlet") = {out[6]};       // l5
Physical Surface("heated_wall") = {out[4]};  // l3
Physical Surface("adiabatic_walls") = {s1, out[0], out[2], out[3], out[5], out[7]};
Physical Volume("fluid") = {out[1]};
