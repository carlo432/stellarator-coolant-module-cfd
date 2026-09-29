#!/usr/bin/env python3
"""Generate a wall-resolved backward-facing-step (BFS) CHT blockMeshDict.
4 blocks: B1 inlet(upper) | B2 exp-upper | B3 exp-lower(recirc) | B4 solid slab.
Heated floor = y=0, x[0.1,0.4] (downstream of the step). Solid under it (y[-0.004,0]).
Wall-resolved: B3 graded to cluster cells at the heated interface (y=0) -> 1st-cell y+ ~0.6.
"""
import math

# geometry
x0,x1,x2 = 0.0, 0.10, 0.40
y0,y1,y2,y3 = -0.004, 0.0, 0.010, 0.020   # solid base / interface / step-top / channel-top
W = 0.04
# divisions
nx_in, nx_exp, nz = 40, 120, 16
ny_up = 25            # B1,B2  y[0.01,0.02]
ny_lo = 50            # B3     y[0,0.01], graded toward interface
ny_so = 8             # B4 solid y[-0.004,0]

# grading R for B3 so first cell (at the interface y=0) ~5e-5 m  (centre y+ ~0.6 at U~2.3)
def find_R(c1, L, N):
    lo,hi=1.0001,2.0
    for _ in range(200):
        k=(lo+hi)/2; c=L*(k-1)/(k**N-1)
        if c>c1: lo=k
        else: hi=k
    return k**(N-1)
R = find_R(5e-5, y2-y1, ny_lo)

# vertices: 10 (x,y) points, x2 z-planes
pts=[(x0,y2),(x0,y3),(x1,y2),(x1,y3),(x2,y2),(x2,y3),(x1,y1),(x2,y1),(x1,y0),(x2,y0)]
def vid(i,zp): return i+zp*10
verts=[]
for zp,z in enumerate((0.0,W)):
    for (x,y) in pts:
        verts.append((x,y,z))

def hexblk(a,b,c,d, div, grade):
    # a,b,c,d = (xlo,ylo)(xhi,ylo)(xhi,yhi)(xlo,yhi) point indices; z0 then zW
    return f"hex ({a} {b} {c} {d} {a+10} {b+10} {c+10} {d+10}) {div} {grade}"

blocks=[
 # B1 inlet upper  x[0,0.1] y[0.01,0.02]  pts:0,2,3,1
 hexblk(0,2,3,1, f"({nx_in} {ny_up} {nz})", "simpleGrading (1 1 1)"),
 # B2 exp upper    x[0.1,0.4] y[0.01,0.02] pts:2,4,5,3
 hexblk(2,4,5,3, f"({nx_exp} {ny_up} {nz})", "simpleGrading (1 1 1)"),
 # B3 exp lower    x[0.1,0.4] y[0,0.01]    pts:6,7,4,2   (recirc; graded to interface)
 hexblk(6,7,4,2, f"({nx_exp} {ny_lo} {nz})", f"simpleGrading (1 {R:.2f} 1)"),
 # B4 solid        x[0.1,0.4] y[-0.004,0]  pts:8,9,7,6
 hexblk(8,9,7,6, f"({nx_exp} {ny_so} {nz})", "simpleGrading (1 1 1)"),
]

def face(a,b): # a,b point indices, returns the z-extruded quad (a b b+10 a+10)
    return f"({a} {b} {b+10} {a+10})"

boundary=f"""
boundary
(
    inlet      {{ type patch; faces ( {face(0,1)} ); }}
    outlet     {{ type patch; faces ( {face(4,5)} {face(7,4)} ); }}
    topWall    {{ type wall;  faces ( {face(1,3)} {face(3,5)} ); }}
    inletFloor {{ type wall;  faces ( {face(0,2)} ); }}
    stepWall   {{ type wall;  faces ( {face(6,2)} ); }}
    heatedBase {{ type wall;  faces ( {face(8,9)} ); }}
    fluidSides {{ type wall;  faces (
        (0 2 3 1) (2 4 5 3) (6 7 4 2)
        ({0+10} {2+10} {3+10} {1+10}) ({2+10} {4+10} {5+10} {3+10}) ({6+10} {7+10} {4+10} {2+10}) ); }}
    solidSides {{ type wall;  faces (
        {face(9,7)} {face(8,6)}
        (8 9 7 6) ({8+10} {9+10} {7+10} {6+10}) ); }}
);
"""

hdr="""/*--------------------------------*- C++ -*----------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
// Backward-facing-step CHT: inlet (upper half) -> step down -> expanded heated floor over a solid slab.
scale 1;
"""
vtxt="vertices\n(\n"+"".join(f"    ( {x} {y} {z} )\n" for (x,y,z) in verts)+");\n"
btxt="blocks\n(\n"+"".join("    "+b+"\n" for b in blocks)+");\n"
out=hdr+"\n"+vtxt+"\n"+btxt+"\nedges ();\n"+boundary+"\nmergePatchPairs ();\n"
open("cht_bfs/system/blockMeshDict","w").write(out)
print(f"wrote cht_bfs/system/blockMeshDict  (B3 grading R={R:.2f}, first-cell {5e-5*1e3:.3f} mm)")
print(f"heated floor area = {(x2-x1)*W:.4f} m^2 ; step ER = {(y3-y1)/(y3-y2):.1f}")
