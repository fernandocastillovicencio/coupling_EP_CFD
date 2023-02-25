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

internalField   uniform (0.1 0 0);

boundaryField
{
    
    #includeEtc "caseDicts/setConstraintTypes"
    
    east
    {
        type            fixedValue;
        value           uniform (3 0 0);
    }

    west
    {
        type            pressureInletOutletVelocity;
        value           $internalField;
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
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 1e5;

boundaryField
{
        
    east
    {
        type            zeroGradient;
    }

    west
    {
        type            fixedValue;
        value           $internalField;
    }
    
    ".*"
    {
        type            zeroGradient;
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

internalField   uniform 302;

boundaryField
{
        
    east
    {
        type            fixedValue;
        value           uniform 290;
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
    def create_p_rgh(casedir):
        with open(casedir+'/0.orig/p_rgh', 'w') as f:
            w_header(f,'volScalarField','p_rgh')
            w(f,""" 
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 1e5;

boundaryField
{
        
    east
    {
        type            fixedFluxPressure;
        value           $internalField;
    }

    west
    {
        type            prghPressure;
        p               $internalField;
        value           $internalField;
    }
    
    ".*"
    {
        type            fixedFluxPressure;
        value           $internalField;
    }

}
            """)
            w_footer(f)
    # ------------------------------ turbulence ------------------------------ #
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
    def create_alphat(casedir):
        with open(casedir+'/0.orig/alphat', 'w') as f:
            w_header(f,'volScalarField','alphat')
            w(f,""" 
dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    east
    {
        type            calculated;
        value           $internalField;
    }

    west
    {
        type            calculated;
        value           $internalField;
    }

    ".*"
    {
        type            compressible::alphatJayatillekeWallFunction;
        Prt             0.85;
        value           $internalField;
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
    create_T(casedir)
    create_p_rgh(casedir)
    # ------------------------------ turbulence ------------------------------ #
    create_nut(casedir)
    create_nuTilda(casedir)
    create_alphat(casedir)
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
    equationOfState Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       28.9;
    }
    equationOfState
    {
        rho0            1;
        T0              300;
        beta            3e-03;
    }
    thermodynamics
    {
        Cp              1040;
        Hf              0;
    }
    transport
    {
        mu              1e-05;
        Pr              0.7;
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
    # ------------------------------------------------------------------------ #
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
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;

    div(phi,U)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss upwind;

    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;

    div(phi,Ekp)    bounded Gauss linear;

    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;

    div(phi,age)    bounded Gauss upwind;
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
            """)
            w_footer(f)
    # ------------------------------------------------------------------------ #
    def create_fvSolution(casedir):
        with open(casedir+'/system/fvSolution', 'w') as f:
            w_header(f,'dictionary','fvSolution')
            w(f,""" 
solvers
{
    p_rgh
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }

    "(U|h|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    age
    {
        $U;
        relTol          0.001;
    }
}

SIMPLE
{
    momentumPredictor false;
    nNonOrthogonalCorrectors 0;

    residualControl
    {
        p_rgh           1e-6;
        U               1e-6;
        h               1e-6;
        "(k|epsilon)"   1e-6;
    }
}

relaxationFactors
{
    fields
    {
        p_rgh           0.7;
    }

    equations
    {
        U               0.3;
        h               0.3;
        "(k|epsilon|R)" 0.7;
        age             1;
    }
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

endTime         60;

deltaT          0.1;

writeControl    timeStep;

writeInterval   50;

purgeWrite      10;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

functions
{
    htc
    {
        type            heatTransferCoeff;
        libs            (fieldFunctionObjects);
        field           T;
        writeControl    writeTime;
        writeInterval   1;
        htcModel        fixedReferenceTemperature;
        patches         (Wall001);
        TRef            273;
    }
}

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