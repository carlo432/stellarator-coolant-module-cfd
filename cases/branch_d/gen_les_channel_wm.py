#!/usr/bin/env python3
"""Spine L: wall-MODELED LES of a periodic plane channel, FLiBe.
Canonical building block that (a) gives the first unsteady result, (b) de-risks the LES pipeline before HPC-1,
and (c) sets up the O benchmark vs Moser-Kim-Mansour DNS. The first completed run measured Re_tau~192,
so the nearest MKM comparison is Re_tau=180.
Domain (2*pi*d, 2d, pi*d), periodic in x (streamwise) & z (spanwise), no-slip walls at y=+-d.
Driven to bulk U via fvOptions meanVelocityForce. Incompressible pimpleFoam + WALE SGS, wall functions."""
import math, os

case="les_channel_wm"
d=0.01                       # channel half-height [m]  -> full height 0.02 (matches our duct channel)
Lx, Ly, Lz = 2*math.pi*d, 2*d, math.pi*d
Ubar=1.1598                  # bulk streamwise velocity [m/s] (matches RANS/CHT)
nu=3.09278e-6                # FLiBe kinematic viscosity [m2/s]
# first shakedown measured Re_tau~192; refine x-z and average longer before claiming benchmark agreement
Nx, Ny, Nz = 32, 64, 32
gradY=12.0                   # symmetric grading toward both walls (edge/centre cell ratio per half)

for sub in ("system","constant","0"):
    os.makedirs(f"{case}/{sub}", exist_ok=True)

hdr=lambda cls,obj:("/*--------------------------------*- C++ -*----------------------------------*/\n"
 f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n")

# ---------- blockMeshDict (single block, double-graded in y) ----------
# simpleGrading y uses two-segment grading "( (0.5 0.5 gradY) (0.5 0.5 1/gradY) )" to cluster at both walls
bm=hdr("dictionary","blockMeshDict")+f"""scale 1;
vertices
(
 (0 0 0) ({Lx} 0 0) ({Lx} {Ly} 0) (0 {Ly} 0)
 (0 0 {Lz}) ({Lx} 0 {Lz}) ({Lx} {Ly} {Lz}) (0 {Ly} {Lz})
);
blocks
(
 hex (0 1 2 3 4 5 6 7) ({Nx} {Ny} {Nz})
 simpleGrading (1 ( (0.5 0.5 {gradY}) (0.5 0.5 {1.0/gradY}) ) 1)
);
edges ();
boundary
(
 bottomWall {{ type wall;   faces ( (0 1 5 4) ); }}
 topWall    {{ type wall;   faces ( (3 7 6 2) ); }}
 inlet      {{ type cyclic; neighbourPatch outlet; faces ( (0 4 7 3) ); }}
 outlet     {{ type cyclic; neighbourPatch inlet;  faces ( (1 2 6 5) ); }}
 front      {{ type cyclic; neighbourPatch back;   faces ( (4 5 6 7) ); }}
 back       {{ type cyclic; neighbourPatch front;  faces ( (0 3 2 1) ); }}
);
mergePatchPairs ();
"""
open(f"{case}/system/blockMeshDict","w").write(bm)

# ---------- transportProperties ----------
open(f"{case}/constant/transportProperties","w").write(
 hdr("dictionary","transportProperties")+f"transportModel Newtonian;\nnu {nu};\n")

# ---------- turbulenceProperties (LES / WALE) ----------
open(f"{case}/constant/turbulenceProperties","w").write(
 hdr("dictionary","turbulenceProperties")+
 "simulationType LES;\nLES\n{\n LESModel WALE;\n turbulence on;\n printCoeffs on;\n delta cubeRootVol;\n"
 " cubeRootVolCoeffs { deltaCoeff 1; }\n}\n")

# ---------- 0/U : cyclic + no-slip walls ----------
open(f"{case}/0/U","w").write(hdr("volVectorField","U")+f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform ({Ubar} 0 0);
boundaryField
{{
 "(bottomWall|topWall)" {{ type noSlip; }}
 "(inlet|outlet)" {{ type cyclic; }}
 "(front|back)"   {{ type cyclic; }}
}}
""")

# ---------- 0/p (kinematic) ----------
open(f"{case}/0/p","w").write(hdr("volScalarField","p")+"""dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 "(bottomWall|topWall)" { type zeroGradient; }
 "(inlet|outlet)" { type cyclic; }
 "(front|back)"   { type cyclic; }
}
""")

# ---------- 0/nut : SGS viscosity, Spalding wall function (wall-modeled) ----------
open(f"{case}/0/nut","w").write(hdr("volScalarField","nut")+"""dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
 "(bottomWall|topWall)" { type nutUSpaldingWallFunction; value uniform 0; }
 "(inlet|outlet)" { type cyclic; }
 "(front|back)"   { type cyclic; }
}
""")

# ---------- fvOptions : drive bulk velocity ----------
open(f"{case}/constant/fvOptions","w").write(hdr("dictionary","fvOptions")+f"""momentumSource
{{
 type meanVelocityForce;
 selectionMode all;
 fields (U);
 Ubar ({Ubar} 0 0);
}}
""")

# ---------- controlDict ----------
open(f"{case}/system/controlDict","w").write(hdr("dictionary","controlDict")+f"""application pimpleFoam;
startFrom latestTime; startTime 0; stopAt endTime;
endTime 60; deltaT 2e-4;
writeControl adjustableRunTime; writeInterval 5; purgeWrite 6;
writeFormat ascii; writePrecision 8; writeCompression off;
timeFormat general; timePrecision 6; runTimeModifiable true;
adjustTimeStep yes; maxCo 0.6; maxDeltaT 1e-3;
functions
{{
 fieldAverage1 {{ type fieldAverage; libs (fieldFunctionObjects); timeStart 20;
   fields ( U {{ mean on; prime2Mean on; base time; }} p {{ mean on; prime2Mean on; base time; }} ); }}
 yPlus {{ type yPlus; libs (fieldFunctionObjects); writeControl writeTime; }}
}}
""")

# ---------- fvSchemes (LES: 2nd-order, low-dissipation) ----------
open(f"{case}/system/fvSchemes","w").write(hdr("dictionary","fvSchemes")+"""ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
 default none;
 div(phi,U) Gauss LUST grad(U);
 div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")

# ---------- fvSolution (PIMPLE) ----------
open(f"{case}/system/fvSolution","w").write(hdr("dictionary","fvSolution")+"""solvers
{
 p { solver GAMG; smoother DICGaussSeidel; tolerance 1e-6; relTol 0.05; }
 pFinal { $p; relTol 0; }
 "(U|k)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }
}
PIMPLE { nCorrectors 2; nNonOrthogonalCorrectors 0; nOuterCorrectors 1; }
""")

print(f"wrote {case}/  (Re_b≈{Ubar*2*d/nu:.0f} on 2d; first shakedown measured Re_tau≈192 -> nearest MKM 180); "
      f"cells={Nx*Ny*Nz} ; domain ({Lx:.4f} {Ly:.4f} {Lz:.4f})")
