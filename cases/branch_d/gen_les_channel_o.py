#!/usr/bin/env python3
"""Spine O: refined wall-resolved LES of the periodic FLiBe channel to recover the log law and certify the
M harness against Moser-Kim-Mansour Re_tau~180. Same case as L but: finer streamwise/spanwise (64x64x64 ->
dx+~39, dz+~19, proper LES resolution) + longer average. ALL FOUR L-gotchas pre-fixed in the dictionaries."""
import math, os
case="les_channel_o"
d=0.01; Lx,Ly,Lz=2*math.pi*d,2*d,math.pi*d
Ubar=1.1598; nu=3.09278e-6
Nx,Ny,Nz=64,64,64        # refined x,z (was 32) -> resolve streaks; Ny unchanged (sublayer already perfect)
gradY=12.0
for sub in ("system","constant","0"): os.makedirs(f"{case}/{sub}", exist_ok=True)
hdr=lambda c,o:("/*--------------------------------*- C++ -*----------------------------------*/\n"
 f"FoamFile {{ version 2.0; format ascii; class {c}; object {o}; }}\n")

open(f"{case}/system/blockMeshDict","w").write(hdr("dictionary","blockMeshDict")+f"""scale 1;
vertices ( (0 0 0) ({Lx} 0 0) ({Lx} {Ly} 0) (0 {Ly} 0) (0 0 {Lz}) ({Lx} 0 {Lz}) ({Lx} {Ly} {Lz}) (0 {Ly} {Lz}) );
blocks ( hex (0 1 2 3 4 5 6 7) ({Nx} {Ny} {Nz}) simpleGrading (1 ( (0.5 0.5 {gradY}) (0.5 0.5 {1.0/gradY}) ) 1) );
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
""")
open(f"{case}/constant/transportProperties","w").write(hdr("dictionary","transportProperties")+f"transportModel Newtonian;\nnu {nu};\n")
open(f"{case}/constant/turbulenceProperties","w").write(hdr("dictionary","turbulenceProperties")+
 "simulationType LES;\nLES { LESModel WALE; turbulence on; printCoeffs on; delta cubeRootVol; cubeRootVolCoeffs { deltaCoeff 1; } }\n")
open(f"{case}/0/U","w").write(hdr("volVectorField","U")+f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform ({Ubar} 0 0);
boundaryField {{ "(bottomWall|topWall)" {{ type noSlip; }} "(inlet|outlet)" {{ type cyclic; }} "(front|back)" {{ type cyclic; }} }}
""")
open(f"{case}/0/p","w").write(hdr("volScalarField","p")+"""dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField { "(bottomWall|topWall)" { type zeroGradient; } "(inlet|outlet)" { type cyclic; } "(front|back)" { type cyclic; } }
""")
open(f"{case}/0/nut","w").write(hdr("volScalarField","nut")+"""dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField { "(bottomWall|topWall)" { type nutUSpaldingWallFunction; value uniform 0; } "(inlet|outlet)" { type cyclic; } "(front|back)" { type cyclic; } }
""")
open(f"{case}/constant/fvOptions","w").write(hdr("dictionary","fvOptions")+f"""momentumSource {{ type meanVelocityForce; selectionMode all; fields (U); Ubar ({Ubar} 0 0); }}
""")
# controlDict -- gotcha #1 (physical-second times), #4 (sane writeInterval, NOT 100000; fieldAverage writeTime)
open(f"{case}/system/controlDict","w").write(hdr("dictionary","controlDict")+"""application pimpleFoam;
startFrom latestTime; startTime 0; stopAt endTime;
endTime 4.0; deltaT 1e-4;          // ~74 flow-through total (tau=0.0542 s)
writeControl timeStep; writeInterval 3000; purgeWrite 3;   // periodic checkpoints incl near-end; NO per-step bloat
writeFormat ascii; writePrecision 8; writeCompression off;
timeFormat general; timePrecision 6; runTimeModifiable true;
adjustTimeStep yes; maxCo 0.6; maxDeltaT 5e-4;
functions
{
 fieldAverage1 { type fieldAverage; libs (fieldFunctionObjects); timeStart 1.5; writeControl writeTime; restartOnOutput false;
   fields ( U { mean on; prime2Mean on; base time; } p { mean on; prime2Mean on; base time; } ); }
 yPlus { type yPlus; libs (fieldFunctionObjects); writeControl writeTime; }
}
""")
open(f"{case}/system/fvSchemes","w").write(hdr("dictionary","fvSchemes")+"""ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss LUST grad(U); div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; } snGradSchemes { default corrected; }
""")
# fvSolution -- gotcha #2 (pRefCell/pRefValue) + #3 (UFinal via regex)
open(f"{case}/system/fvSolution","w").write(hdr("dictionary","fvSolution")+"""solvers
{
 p { solver GAMG; smoother DICGaussSeidel; tolerance 1e-6; relTol 0.05; }
 pFinal { $p; relTol 0; }
 "(U|k)(Final)?" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }
}
PIMPLE { nCorrectors 2; nNonOrthogonalCorrectors 0; nOuterCorrectors 1; pRefCell 0; pRefValue 0; }
""")
# setExprFields perturbation (trip turbulence)
open(f"{case}/system/setExprFieldsDict","w").write(hdr("dictionary","setExprFieldsDict")+"""expressions
(
 U { field U; dimensions [0 1 -1 0 0 0 0];
     expression #{ vector( 1.74*(1.0 - sqr(pos().y()/0.01 - 1.0)) + 0.35*sin(40.0*pos().z())*cos(3000.0*pos().y()),
                           0.35*sin(35.0*pos().x())*sin(28.0*pos().z()),
                           0.35*cos(32.0*pos().x())*sin(2500.0*pos().y()) ) #}; }
);
""")
print(f"wrote {case}/  cells={Nx*Ny*Nz} (dx+~{2*math.pi*192/Nx:.0f}, dz+~{math.pi*192/Nz:.0f} at Re_tau~192); all 4 gotchas pre-fixed")
