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

internalField   uniform (0 0 0);

boundaryField
{
    ".*"
    {
        type            noSlip;
    }
    
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_T(casedir):
        with open(casedir+'/0.orig/T', 'w') as f:
            w_header(f,'volScalarField','T')
            w(f,""" 
dimensions      [0 0 0 1 0 0 0];

internalField   uniform 273;

boundaryField
{
    inlet
    {
        type        fixedValue;
        value       uniform 270;
    }
    outlet
    {
        type        fixedValue;
        value       uniform 276;
    }
    ".*"
    {
        type        fixedValue;
        value       uniform 273;
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_p(casedir):
        with open(casedir+'/0.orig/p', 'w') as f:
            w_header(f,'volScalarField','p')
            w(f,""" 
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 1e5;

boundaryField
{   
    ".*"
    {
        type            calculated;
        value           $internalField;
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_p_rgh(casedir):
        with open(casedir+'/0.orig/p_rgh', 'w') as f:
            w_header(f,'volScalarField','p_rgh')
            w(f,""" 
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 1e5;

boundaryField
{    
    ".*"
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
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

internalField   uniform 3.75e-04;

boundaryField
{
    ".*"
    {
        type            kqRWallFunction;
        value           uniform 3.75e-04;
    }
    
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_epsilon(casedir):
        with open(casedir+'/0.orig/epsilon', 'w') as f:
            w_header(f,'volScalarField','epsilon')
            w(f,""" 
dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 4e-06;

boundaryField
{
    ".*"
    {
        type            epsilonWallFunction;
        value           uniform 4e-06;
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

internalField   uniform 0.12;

boundaryField
{
    ".*"
    {
        type            omegaWallFunction;
        value           uniform 0.12;
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

internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_alphat(casedir):
        with open(casedir+'/0.orig/alphat', 'w') as f:
            w_header(f,'volScalarField','alphat')
            w(f,""" 
dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            compressible::alphatWallFunction;
        Prt             0.85;
        value           uniform 0;
    }    
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    # --------------------------------- tasks -------------------------------- #
    create_U(casedir)
    create_T(casedir)
    create_p(casedir)
    create_p_rgh(casedir)
    # ------------------------------ turbulence ------------------------------ #
    create_k(casedir)
    create_epsilon(casedir)
    create_omega(casedir)
    create_nut(casedir)
    create_alphat(casedir)
# ---------------------------------------------------------------------------- #
def folder_constant(casedir):
    # --------------------------------- files -------------------------------- #
    def create_turbulenceProperties(casedir):
        with open(casedir+'/constant/turbulenceProperties', 'w') as f:
            w_header(f,'dictionary','turbulenceProperties')
            w(f,""" 
simulationType          RAS;

RAS
{
    RASModel            kOmegaSST;

    turbulence          on;

    printCoeffs         on;
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_thermophysicalProperties(casedir):
        with open(casedir+'/constant/thermophysicalProperties', 'w') as f:
            w_header(f,'dictionary','thermophysicalProperties')
            w(f,""" 
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       28.96;
    }
    thermodynamics
    {
        Cp              1004.4;
        Hf              0;
    }
    transport
    {
        mu              1.831e-05;
        Pr              0.705;
    }
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_g(casedir):
        with open(casedir+'/constant/g', 'w') as f:
            w_header(f,'dictionary','g')
            w(f,""" 
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);
            """)
            w_footer(f)
    # --------------------------------- tasks -------------------------------- #
    create_turbulenceProperties(casedir)
    create_thermophysicalProperties(casedir)
    create_g(casedir)
# ---------------------------------------------------------------------------- #
def folder_system(casedir):
    # --------------------------------- files -------------------------------- #
    def create_fvSchemes(casedir):
        with open(casedir+'/system/fvSchemes', 'w') as f:
            w_header(f,'dictionary','fvSchemes')
            w(f,""" 
ddtSchemes
{
    default steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;

    div(phi,U)      bounded Gauss limitedLinear 0.2;

    energy          bounded Gauss limitedLinear 0.2;
    div(phi,K)      $energy;
    div(phi,h)      $energy;

    turbulence      bounded Gauss limitedLinear 0.2;
    div(phi,k)      $turbulence;
    div(phi,epsilon) $turbulence;
    div(phi,omega) $turbulence;

    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear orthogonal;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         orthogonal;
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
    Phi
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-06;
        relTol          0.01;
    }

    p
    {
        $Phi;
    }
    p_rgh
    {
        solver           GAMG;
        tolerance        1e-7;
        relTol           0.01;
        smoother         DICGaussSeidel;
    }

    "(U|h|k|epsilon|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0.1;
    }
}

SIMPLE
{
    momentumPredictor no;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;

    residualControl
    {
        p_rgh           1e-4;
        U               1e-4;
        h               1e-4;

        // possibly check turbulence fields
        "(k|epsilon|omega)" 1e-4;
    }
}

relaxationFactors
{
    fields
    {
        rho             1.0;
        p_rgh           0.7;
    }
    equations
    {
        U               0.3;
        h               0.3;
        "(k|epsilon|omega)" 0.7;
    }
}
potentialFlow
{
    nNonOrthogonalCorrectors 5;
}

            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_fvOptions(casedir):
        with open(casedir+'/system/fvOptions', 'w') as f:
            w_header(f,'dictionary','fvOptions')
            w(f,""" 
limitT
{
    type       limitTemperature;
    min        200;
    max        400;
    selectionMode all;
}
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_controlDict(casedir):
        with open(casedir+'/system/controlDict', 'w') as f:
            w_header(f,'dictionary','controlDict')
            w(f,""" 
application     buoyantSimpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          10;

writeControl    timeStep;

writeInterval   10;

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
    create_fvOptions(casedir)
    create_controlDict(casedir)
    create_decomposeParDict(casedir)
# ---------------------------------------------------------------------------- #
def create_solver_files(casedir):
    folder_0(casedir)
    folder_constant(casedir)
    folder_system(casedir)
# ---------------------------------------------------------------------------- #