#!/usr/bin/env python3
"""Outlet-offset recirculation CHT: straight duct, but the OUTLET is only the upper half of the
end face; the lower half (over the heated floor) is a wall -> forces a recirculation in the
downstream floor corner -> floor hotspot. Wall-resolved (graded near y=0), solid slab below.
3 blocks: solid | fluid-lower(graded, recirc) | fluid-upper."""
import math
x0,x1=0.0,0.4
y0,y1,y2,y3=-0.004,0.0,0.010,0.020
W=0.04
nx,nz=120,16; ny_lo,ny_up,ny_so=50,25,8
def find_R(c1,L,N):
    lo,hi=1.0001,2.0
    for _ in range(200):
        k=(lo+hi)/2; c=L*(k-1)/(k**N-1)
        if c>c1: lo=k
        else: hi=k
    return k**(N-1)
R=find_R(5e-5,y2-y1,ny_lo)
pts=[(x0,y0),(x1,y0),(x0,y1),(x1,y1),(x0,y2),(x1,y2),(x0,y3),(x1,y3)]
verts=[(x,y,z) for z in (0.0,W) for (x,y) in pts]
def hx(a,b,c,d,div,g): return f"hex ({a} {b} {c} {d} {a+8} {b+8} {c+8} {d+8}) {div} {g}"
blocks=[
 hx(0,1,3,2,f"({nx} {ny_so} {nz})","simpleGrading (1 1 1)"),       # solid  y[-0.004,0]
 hx(2,3,5,4,f"({nx} {ny_lo} {nz})",f"simpleGrading (1 {R:.2f} 1)"),# fluid lower y[0,0.01] graded
 hx(4,5,7,6,f"({nx} {ny_up} {nz})","simpleGrading (1 1 1)"),       # fluid upper y[0.01,0.02]
]
def f(a,b): return f"({a} {b} {b+8} {a+8})"
boundary=f"""
boundary
(
    inlet      {{ type patch; faces ( {f(2,4)} {f(4,6)} ); }}
    outlet     {{ type patch; faces ( {f(5,7)} ); }}                 // upper half only
    outletWall {{ type wall;  faces ( {f(3,5)} ); }}                 // lower half blocked -> recirc
    topWall    {{ type wall;  faces ( {f(6,7)} ); }}
    heatedBase {{ type wall;  faces ( {f(0,1)} ); }}
    fluidSides {{ type wall;  faces (
        (2 3 5 4) (4 5 7 6) ({2+8} {3+8} {5+8} {4+8}) ({4+8} {5+8} {7+8} {6+8}) ); }}
    solidSides {{ type wall;  faces ( {f(1,3)} {f(0,2)}
        (0 1 3 2) ({0+8} {1+8} {3+8} {2+8}) ); }}
);
"""
hdr="""/*--------------------------------*- C++ -*----------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1;
"""
out=hdr+"\nvertices\n(\n"+"".join(f"    ( {x} {y} {z} )\n" for x,y,z in verts)+");\n\n"
out+="blocks\n(\n"+"".join("    "+b+"\n" for b in blocks)+");\n\nedges ();\n"+boundary+"\nmergePatchPairs ();\n"
open("cht_offset/system/blockMeshDict","w").write(out)
print(f"wrote cht_offset/system/blockMeshDict (grading R={R:.2f}); outlet=upper half, lower-floor blocked")
