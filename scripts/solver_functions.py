from scripts.mesh_functions import w, w_header, w_footer
# ---------------------------------------------------------------------------- #
# -------------------------------- FOR SOLVER -------------------------------- #
# ---------------------------------------------------------------------------- #
# def create_transProp(casedir):  
#     with open(casedir+'/constant/transportProperties', 'w') as f:
#         w_header(f,'dictionary','transportProperties')
#         w(f,""" 
# transportModel  Newtonian;

# nu              1.5e-05;
#           """)
#         # w(f, 'nu              0.01;')
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_turbProp(casedir):  
#     with open(casedir+'/constant/turbulenceProperties', 'w') as f:
#         w_header(f,'dictionary','turbulenceProperties')
#         w(f,""" 
# simulationType      RAS;

# RAS
# {
#     RASModel        kEpsilon;

#     turbulence      on;

#     printCoeffs     on;
# }
#           """)
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_0_U(casedir):
#     with open(casedir+'/0.orig/U', 'w') as f:
#         w_header(f,'volVectorField','U')
#         w(f, """ 
# Uinlet          (10 0 0);

# dimensions      [0 1 -1 0 0 0 0];

# internalField   uniform (0 0 0);

# boundaryField
# {
#     inlet
#     {
#         type            fixedValue;
#         value           uniform $Uinlet;
#     }

#     outlet
#     {
#         type            pressureInletOutletVelocity;
#         value           uniform (0 0 0);
#     }

#     wall
#     {
#         type            noSlip;
#     }
    
#     #includeEtc "caseDicts/setConstraintTypes"
# }
#           """)
#         # w(f, 'dimensions [0 1 -1 0 0 0 0];')
#         # w(f, '')
#         # w(f, 'internalField uniform (0 0 0);')
#         # w(f, '')
#         # w(f, 'boundaryField')
#         # w(f, '{')
#         # w(f, '    minX')
#         # w(f, '    {')
#         # w(f, '          type fixedValue;')
#         # w(f, '          value (1 0 0);')
#         # w(f, '    }\n')
#         # w(f, '    maxX')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    minY')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    maxY')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    minZ')
#         # w(f, '    {')
#         # w(f, '          type noSlip;')
#         # w(f, '    }\n')
#         # w(f, '    maxZ')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }')
#         # w(f, '}')
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_0_p(casedir):
#     with open(casedir+'/0.orig/p', 'w') as f:
#         w_header(f,'volScalarField','p')
#         w(f, """ 
# dimensions      [0 2 -2 0 0 0 0];

# internalField   uniform 0;

# boundaryField
# {
#     inlet
#     {
#         type            zeroGradient;
#     }

#     outlet
#     {
#         type            totalPressure;
#         p0              uniform 0;
#     }

#     wall
#     {
#         type            zeroGradient;
#     }
    
#     #includeEtc "caseDicts/setConstraintTypes"
# }

#           """)
#         # w(f, 'dimensions [0 2 -2 0 0 0 0];')
#         # w(f, '')
#         # w(f, 'internalField uniform 0;')
#         # w(f, '')
#         # w(f, 'boundaryField')
#         # w(f, '{')
#         # w(f, '    minX')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    maxX')
#         # w(f, '    {')
#         # w(f, '          type fixedValue;')
#         # w(f, '          value uniform 0;')
#         # w(f, '    }\n')
#         # w(f, '    minY')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    maxY')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    minZ')
#         # w(f, '    {')
#         # w(f, '          type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '    maxZ')
#         # w(f, '    {')
#         # w(f, '      type zeroGradient;')
#         # w(f, '    }\n')
#         # w(f, '}')
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_0_k(casedir):
#     with open(casedir+'/0.orig/k', 'w') as f:
#         w_header(f,'volScalarField','k')
#         w(f, """ 
# kInlet          1.5;   // approx k = 1.5*(I*U)^2 ; I = 0.1

# dimensions      [0 2 -2 0 0 0 0];

# internalField   uniform $kInlet;

# boundaryField
# {
#     inlet
#     {
#         type            fixedValue;
#         value           uniform $kInlet;
#     }

#     outlet
#     {
#         type            inletOutlet;
#         inletValue      uniform $kInlet;
#         value           uniform $kInlet;
#     }

#     wall
#     {
#         type            kqRWallFunction;
#         value           uniform $kInlet;
#     }
    
#     #includeEtc "caseDicts/setConstraintTypes"

# }
#         """)
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_0_epsilon(casedir):
#     with open(casedir+'/0.orig/epsilon', 'w') as f:
#         w_header(f,'volScalarField','epsilon')
#         w(f, """ 
# epsilonInlet  0.03; // Cmu^0.75 * k^1.5 / L ; L =10

# dimensions      [0 2 -3 0 0 0 0];

# internalField   uniform $epsilonInlet;

# boundaryField
# {
#     inlet
#     {
#         type            fixedValue;
#         value           uniform $epsilonInlet;
#     }

#     outlet
#     {
#         type            inletOutlet;
#         inletValue      uniform $epsilonInlet;
#         value           uniform $epsilonInlet;
#     }

#     wall
#     {
#         type            epsilonWallFunction;
#         value           uniform $epsilonInlet;
#     }
    
#     #includeEtc "caseDicts/setConstraintTypes"

# }
#         """)
#         w_footer(f)
# # ---------------------------------------------------------------------------- #
# def create_0_nut(casedir):
#     with open(casedir+'/0.orig/nut', 'w') as f:
#         w_header(f,'volScalarField','nut')
#         w(f, """ 
# dimensions      [0 2 -1 0 0 0 0];

# internalField   uniform 0;

# boundaryField
# {
#     inlet
#     {
#         type            calculated;
#         value           uniform 0;
#     }

#     outlet
#     {
#         type            calculated;
#         value           uniform 0;
#     }

#     wall
#     {
#         type            nutkWallFunction;
#         value           uniform 0;
#     }
    
#     #includeEtc "caseDicts/setConstraintTypes"

# }
#         """)
#         w_footer(f)

# ---------------------------------------------------------------------------- #
def create_0_alphat(casedir):
    with open(casedir+'/0.orig/alphat', 'w') as f:
        w_header(f,'volScalarField','alphat')
        w(f,""" 
dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    floor
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }

    ceiling
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }

    fixedWalls
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    
    Flr001
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Roof001
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Wall001
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Wall002
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Wall003
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Wall004
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
    Win001
    {
        type            compressible::alphatWallFunction;
        value           uniform 0;
    }
}
          """)
        
def create_0_epsilon(casedir):
    with open(casedir+'/0.orig/epsilon', 'w') as f:
        w_header(f,'volScalarField','epsilon')
        w(f,""" 
dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 0.01;

boundaryField
{
    floor
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }

    ceiling
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }

    fixedWalls
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    
    Flr001
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Roof001
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Wall001
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Wall002
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Wall003
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Wall004
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
    Win001
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }
}
          """)
        w_footer(f)

def create_0_G(casedir):
    with open(casedir+'/0.orig/G', 'w') as f:
        w_header(f,'volScalarField','G')
        w(f,""" 
dimensions      [1 0 -3 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    floor
    {
        type            MarshakRadiation;
        value           uniform 0;
    }

    fixedWalls
    {
        type            MarshakRadiation;
        value           uniform 0;
    }

    ceiling
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    
    Flr001
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Roof001
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Wall001
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Wall002
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Wall003
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Wall004
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
    Win001
    {
        type            MarshakRadiation;
        value           uniform 0;
    }
}
          """)
        w_footer(f)

def create_0_k(casedir):
    with open(casedir+'/0.orig/k', 'w') as f:
        w_header(f,'volScalarField','k')
        w(f, """ 
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.1;

boundaryField
{
    floor
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }

    ceiling
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }

    fixedWalls
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    
    
    Flr001
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Roof001
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Wall001
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Wall002
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Wall003
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Wall004
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
    Win001
    {
        type            kqRWallFunction;
        value           uniform 0.1;
    }
}          
          """)
        w_footer(f)

def create_0_nut(casedir):
    with open(casedir+'/0.orig/nut', 'w') as f:
        w_header(f,'volScalarField','nut')
        w(f,""" 
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    floor
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    ceiling
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    fixedWalls
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    
    Flr001
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Roof001
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Wall001
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Wall002
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Wall003
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Wall004
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    Win001
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
}          
          """)
        w_footer(f)
        
def create_0_p(casedir):
    with open(casedir+'/0.orig/p', 'w') as f:
        w_header(f,'volScalarField','p')
        w(f, """ 
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 100000;

boundaryField
{
    floor
    {
        type            calculated;
        value           $internalField;
    }

    ceiling
    {
        type            calculated;
        value           $internalField;
    }

    fixedWalls
    {
        type            calculated;
        value           $internalField;
    }

    Flr001
    {
        type            calculated;
        value           $internalField;
    }
    Roof001
    {
        type            calculated;
        value           $internalField;
    }
    Wall001
    {
        type            calculated;
        value           $internalField;
    }
    Wall002
    {
        type            calculated;
        value           $internalField;
    }
    Wall003
    {
        type            calculated;
        value           $internalField;
    }
    Wall004
    {
        type            calculated;
        value           $internalField;
    }
    Win001
    {
        type            calculated;
        value           $internalField;
    }
}          
          """)
        w_footer(f)
  
def create_0_p_rgh(casedir):
    with open(casedir+'/0.orig/p_rgh', 'w') as f:
        w_header(f,'volScalarField','p_rgh')
        w(f,""" 
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 100000;

boundaryField
{
    floor
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }

    ceiling
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }

    fixedWalls
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }

    Flr001
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Roof001
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Wall001
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Wall002
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Wall003
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Wall004
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
    Win001
    {
        type            fixedFluxPressure;
        value           uniform 100000;
    }
}
          """)

def create_0_U(casedir):
    with open(casedir+'/0.orig/U', 'w') as f:
        w_header(f,'volVectorField','U')
        w(f,""" 
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    floor
    {
        type            noSlip;
    }

    ceiling
    {
        type            noSlip;
    }

    fixedWalls
    {
        type            noSlip;
    }

    Flr001
    {
        type            noSlip;
    }
    Roof001
    {
        type            noSlip;
    }
    Wall001
    {
        type            noSlip;
    }
    Wall002
    {
        type            noSlip;
    }
    Wall003
    {
        type            noSlip;
    }
    Wall004
    {
        type            noSlip;
    }
    Win001
    {
        type            noSlip;
    }
}
          """)
        w_footer(f)
        
def create_0_T(casedir):
    with open(casedir+'/0.orig/T', 'w') as f:
        w_header(f,'volScalarField','T')
        w(f, """ 
dimensions      [0 0 0 1 0 0 0];

internalField   uniform 300;

boundaryField
{
    floor
    {
        type            fixedValue;
        value           uniform 300.0;
    }

    ceiling
    {
        type            fixedValue;
        value           uniform 300.0;
    }

    fixedWalls
    {
        type            zeroGradient;
    }

    Flr001
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Roof001
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Wall001
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Wall002
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Wall003
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Wall004
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    Win001
    {
        type            fixedValue;
        value           uniform 500.0;
    }
    
}
          """)
        w_footer(f)
   
# ---------------------------------------------------------------------------- #
def create_constant_boundaryRadiationProperties(casedir):
    with open(casedir+'/constant/boundaryRadiationProperties', 'w') as f:
        w_header(f,'dictionary','boundaryRadiationProperties')
        w(f, """ 
".*"
{
    type            lookup;
    emissivity      1.0;
    absorptivity    1.0;
}          
          """)
        w_footer(f)
        
def create_constant_g(casedir):
    with open(casedir+'/constant/g', 'w') as f:
        w_header(f,'uniformDimensionedVectorField','g')
        w(f, """ 
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);          
          """)
        w_footer(f)

def create_constant_radiationProperties(casedir):
    with open(casedir+'/constant/radiationProperties', 'w') as f:
        w_header(f,'dictionary','radiationProperties')
        w(f, """ 
radiation on;

radiationModel  P1;

// Number of flow iterations per radiation iteration
solverFreq 1;

absorptionEmissionModel constantAbsorptionEmission;

constantAbsorptionEmissionCoeffs
{
    absorptivity    absorptivity    [0 -1 0 0 0 0 0] 0.5;
    emissivity      emissivity      [0 -1 0 0 0 0 0] 0.5;
    E               E               [1 -1 -3 0 0 0 0] 0;
}

scatterModel    none;

sootModel       none;          
          """)
        w_footer(f)

def create_constant_thermophysicalProperties(casedir):
    with open(casedir+'/constant/thermophysicalProperties', 'w') as f:
        w_header(f,'dictionary','thermophysicalProperties')
        w(f, """ 
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

pRef            100000;

mixture
{
    specie
    {
        molWeight       28.9;
    }
    thermodynamics
    {
        Cp              1000;
        Hf              0;
    }
    transport
    {
        mu              1.8e-05;
        Pr              0.7;
    }
}
          """)
        w_footer(f)

def create_constant_turbulenceProperties(casedir):
    with open(casedir+'/constant/turbulenceProperties', 'w') as f:
        w_header(f,'dictionary','turbulenceProperties')
        w(f, """ 
simulationType      RAS;

RAS
{
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;
}          
          """)
        w_footer(f)

# ---------------------------------------------------------------------------- #

def create_system_controlDict(casedir):
    with open(casedir+'/system/controlDict', 'w') as f:
        w_header(f,'dictionary','controlDict')
        w(f,""" 
application     buoyantSimpleFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;
          """)
        w_footer(f)
        
def create_system_fvSchemes(casedir):
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

    energy          bounded Gauss upwind;
    div(phi,K)      $energy;
    div(phi,h)      $energy;

    turbulence      bounded Gauss upwind;
    div(phi,k)      $turbulence;
    div(phi,epsilon) $turbulence;

    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
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
          """)
        w_footer(f)
        
def create_system_fvSolution(casedir):
    with open(casedir+'/system/fvSolution', 'w') as f:
        w_header(f,'dictionary','fvSolution')
        w(f,""" 
solvers
{
    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.01;
    }

    "(U|h|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-05;
        relTol          0.1;
    }

    G
    {
        $p_rgh;
        tolerance       1e-05;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;

    residualControl
    {
        p_rgh           1e-2;
        U               1e-3;
        h               1e-3;
        G               1e-3;

        // possibly check turbulence fields
        "(k|epsilon|omega)" 1e-3;
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
        U               0.2;
        h               0.2;
        "(k|epsilon|R)" 0.5;
        G               0.7;
    }
}          
          """)
        w_footer(f)
        
# ---------------------------------------------------------------------------- #
def folder_0(casedir):
    create_0_alphat(casedir)
    create_0_epsilon(casedir)
    create_0_G(casedir)
    create_0_k(casedir)
    create_0_nut(casedir)
    create_0_p(casedir)
    create_0_p_rgh(casedir)
    create_0_U(casedir)
    create_0_T(casedir)
    
def folder_constant(casedir):
    create_constant_boundaryRadiationProperties(casedir)
    create_constant_g(casedir)
    create_constant_radiationProperties(casedir)
    create_constant_thermophysicalProperties(casedir)
    create_constant_turbulenceProperties(casedir)

def folder_system(casedir):
    create_system_controlDict(casedir)
    create_system_fvSchemes(casedir)
    create_system_fvSolution(casedir)
# ---------------------------------------------------------------------------- #
def create_solver_files(casedir):
    folder_0(casedir)
    folder_constant(casedir)
    folder_system(casedir)