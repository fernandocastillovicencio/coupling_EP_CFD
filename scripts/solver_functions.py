from scripts.mesh_functions import w, w_header, w_footer
# ---------------------------------------------------------------------------- #
# -------------------------------- FOR SOLVER -------------------------------- #
# ---------------------------------------------------------------------------- #
def folder_0(casedir):
    # --------------------------------- files -------------------------------- #
    def create_U(casedir):
        with open(casedir+'/0.orig/U', 'w') as f:
            w_header(f,'volVectorField','U')
            w(f,""" 
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0.5 0 0);

boundaryField
{
    east
    {
        type            fixedValue;
        value           uniform (3 0 0);
    }

    west
    {
        type            zeroGradient;
    }
    
    ceiling
    {
        type            zeroGradient;
    }

    ".*"
    {
        type            noSlip;
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_p(casedir):
        with open(casedir+'/0.orig/p', 'w') as f:
            w_header(f,'volScalarField','p')
            w(f,""" 
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
        
    east
    {
        type            zeroGradient;
    }

    west
    {
        type            fixedValue;
        value           uniform 0;
    }
    
    ".*"
    {
        type            zeroGradient;
    }

}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_nut(casedir):
        with open(casedir+'/0.orig/nut', 'w') as f:
            w_header(f,'volScalarField','nut')
            w(f,""" 
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 1e-05;

boundaryField
{
    east
    {
        type            calculated;
        value           uniform 0;
    }

    west
    {
        type            calculated;
        value           uniform 0;
    }

    ".*"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

}

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_nuTilda(casedir):
        with open(casedir+'/0.orig/nuTilda', 'w') as f:
            w_header(f,'volScalarField','nuTilda')
            w(f,""" 
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 4e-05;

boundaryField
{
    east
    {
        type            fixedValue;
        value           uniform 0;
    }

    west
    {
        type            zeroGradient;
    }

    ".*"
    {
        type            zeroGradient;
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_k(casedir):
        with open(casedir+'/0.orig/k', 'w') as f:
            w_header(f,'volScalarField','k')
            w(f,""" 
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.375;

boundaryField
{
    east
    {
        type            fixedValue;
        value           uniform 0.375;
    }

    west
    {
        type            zeroGradient;
    }

    ".*"
    {
        type            kqRWallFunction;
        value           uniform 0.375;
    }

}


// ************************************************************************* //

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_epsilon(casedir):
        with open(casedir+'/0.orig/epsilon', 'w') as f:
            w_header(f,'volScalarField','epsilon')
            w(f,""" 
dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 14.855;

boundaryField
{
    east
    {
        type            fixedValue;
        value           uniform 14.855;
    }

    west
    {
        type            zeroGradient;
    }

    ".*"
    {
        type            epsilonWallFunction;
        value           uniform 14.855;
    }

    
}

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_omega(casedir):
        with open(casedir+'/0.orig/omega', 'w') as f:
            w_header(f,'volScalarField','omega')
            w(f,""" 
dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 440.15;

boundaryField
{
    east
    {
        type            fixedValue;
        value           $internalField;
    }

    west
    {
        type            zeroGradient;
    }

    ".*"
    {
        type            omegaWallFunction;
        value           $internalField;
    }
}

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    # --------------------------------- tasks -------------------------------- #
    create_U(casedir)
    create_p(casedir)
    # ------------------------------ turbulence ------------------------------ #
    create_nut(casedir)
    create_nuTilda(casedir)
    create_k(casedir)
    create_epsilon(casedir)
    create_omega(casedir)
    # ---------------------------------------------------------------------------- #
def folder_constant(casedir):
    # --------------------------------- files -------------------------------- #
    def create_turbulenceProperties(casedir):
        with open(casedir+'/constant/turbulenceProperties', 'w') as f:
            w_header(f,'dictionary','turbulenceProperties')
            w(f,""" 
simulationType      RAS;

RAS
{
    // Tested with kEpsilon, realizableKE, kOmega, kOmegaSST,
    // ShihQuadraticKE, LienCubicKE.
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_transportProperties(casedir):
        with open(casedir+'/constant/transportProperties', 'w') as f:
            w_header(f,'dictionary','transportProperties')
            w(f,""" 
transportModel  Newtonian;

nu              1e-05;
            """)
            w_footer(f)
    # --------------------------------- tasks -------------------------------- #
    create_turbulenceProperties(casedir)
    create_transportProperties(casedir)
# ---------------------------------------------------------------------------- #
def folder_system(casedir):
    # --------------------------------- files -------------------------------- #
    def create_fvSchemes(casedir):
        with open(casedir+'/system/fvSchemes', 'w') as f:
            w_header(f,'dictionary','fvSchemes')
            w(f,""" 
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;

    div(phi,U)      bounded Gauss linearUpwind grad(U);

    turbulence      bounded Gauss limitedLinear 1;
    div(phi,k)      $turbulence;
    div(phi,epsilon) $turbulence;
    div(phi,omega)  $turbulence;

    div(nonlinearStress) Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}


            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_fvSolution(casedir):
        with open(casedir+'/system/fvSolution', 'w') as f:
            w_header(f,'dictionary','fvSolution')
            w(f,""" 
solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    "(U|k|epsilon|omega|f|v2)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {
        p               1e-2;
        U               1e-3;
        "(k|epsilon|omega|f|v2)" 1e-3;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9; // 0.9 is more stable but 0.95 more convergent
        ".*"            0.9; // 0.9 is more stable but 0.95 more convergent
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_controlDict(casedir):
        with open(casedir+'/system/controlDict', 'w') as f:
            w_header(f,'dictionary','controlDict')
            w(f,""" 
application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1500;

deltaT          5;

writeControl    timeStep;

writeInterval   20;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_decomposeParDict(casedir):
        with open(casedir+'/system/decomposeParDict', 'w') as f:
            w_header(f,'dictionary','decomposeParDict')
            w(f,"""
numberOfSubdomains 26;
method          scotch;
            """)
            w_footer(f)
    # --------------------------------- tasks -------------------------------- #
    create_fvSchemes(casedir)
    create_fvSolution(casedir)
    create_controlDict(casedir)
    create_decomposeParDict(casedir)
# ---------------------------------------------------------------------------- #
def create_solver_files(casedir):
    folder_0(casedir)
    folder_constant(casedir)
    folder_system(casedir)
# ---------------------------------------------------------------------------- #