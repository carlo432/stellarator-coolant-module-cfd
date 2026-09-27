// Stellarator liquid-blanket segment v1 toy geometry.
// Purpose: meshable planning geometry only. Do not treat as real stellarator CAD.
// Concept: curved/wavy first-wall segment with a local liquid blanket pocket, one inlet, one outlet.
// Units: meters.

SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.0015;
Mesh.CharacteristicLengthMax = 0.0060;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

depth = 0.040;
lc_hot = 0.0015;
lc_bulk = 0.0060;

// Coolant-facing first wall: deliberately non-axisymmetric/wavy local segment.
fw1 = newp; Point(fw1) = {0.000,  0.000, 0, lc_hot};
fw2 = newp; Point(fw2) = {0.040, -0.010, 0, lc_hot};
fw3 = newp; Point(fw3) = {0.100, -0.018, 0, lc_hot};
fw4 = newp; Point(fw4) = {0.165, -0.010, 0, lc_hot};
fw5 = newp; Point(fw5) = {0.225,  0.006, 0, lc_hot};
fw6 = newp; Point(fw6) = {0.285, -0.001, 0, lc_hot};
fw7 = newp; Point(fw7) = {0.325, -0.012, 0, lc_hot};

// Blanket back wall: offset from first wall with a mid-segment pocket/bulge.
bw1 = newp; Point(bw1) = {0.000,  0.040, 0, lc_bulk};
bw2 = newp; Point(bw2) = {0.040,  0.038, 0, lc_bulk};
bw3 = newp; Point(bw3) = {0.100,  0.045, 0, lc_bulk};
bw4 = newp; Point(bw4) = {0.165,  0.060, 0, lc_bulk};
bw5 = newp; Point(bw5) = {0.225,  0.062, 0, lc_bulk};
bw6 = newp; Point(bw6) = {0.285,  0.050, 0, lc_bulk};
bw7 = newp; Point(bw7) = {0.325,  0.038, 0, lc_bulk};

l_inlet = newl; Line(l_inlet) = {fw1, bw1};
l_back = newl; BSpline(l_back) = {bw1, bw2, bw3, bw4, bw5, bw6, bw7};
l_outlet = newl; Line(l_outlet) = {bw7, fw7};
l_firstwall = newl; BSpline(l_firstwall) = {fw1, fw2, fw3, fw4, fw5, fw6, fw7};

// Refine near the first wall and the pocket throat.
Field[1] = Distance;
Field[1].EdgesList = {l_firstwall};
Field[1].NNodesByEdge = 350;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_hot;
Field[2].SizeMax = lc_bulk;
Field[2].DistMin = 0.004;
Field[2].DistMax = 0.030;
Background Field = 2;

loop1 = newll; Curve Loop(loop1) = {l_inlet, l_back, l_outlet, -l_firstwall};
s1 = news; Plane Surface(s1) = {loop1};

out[] = Extrude {0, 0, depth} { Surface{s1}; Layers{24}; Recombine; };
// For this curve loop, lateral surfaces follow:
// out[2] inlet, out[3] blanket_back_wall, out[4] outlet, out[5] heated_first_wall.
Physical Surface("inlet") = {out[2]};
Physical Surface("outlet") = {out[4]};
Physical Surface("heated_first_wall") = {out[5]};
Physical Surface("blanket_back_wall") = {out[3]};
Physical Surface("side_walls") = {s1, out[0]};
Physical Volume("fluid") = {out[1]};
