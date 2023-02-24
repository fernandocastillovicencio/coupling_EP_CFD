# ----------------------------- import libraries ----------------------------- #
import os
import shutil
# ---------------------------------------------------------------------------- #
def rel(x):
    return os.path.join(os.path.dirname(__file__), x)
# ---------------------------------------------------------------------------- #
def cleanFolder(x):
    if os.path.exists(x):
        try:
            shutil.rmtree(x)
        except OSError as e:
            print("Error: %s : %s" % (x, e.strerror))
    os.mkdir(x)
# ---------------------------------------------------------------------------- #
def createDirs(casedir):
    cleanFolder(casedir+'/0.orig/')
    cleanFolder(casedir+'/constant/')
    cleanFolder(casedir+'/constant/triSurface')
    cleanFolder(casedir+'/system/')
# ---------------------------------------------------------------------------- #
def w(f, x):
    f.write(x)
    f.write('\n')
    # ------------------------------------------------------------------------ #
def w_header(f,myClass,myObject):
    w(f, '/*--------------------------------*- C++ -*----------------------------------*\\')
    w(f, '| =========                 |                                                 |')
    w(f, '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |')
    w(f, '|  \\\\    /   O peration     | Version:  v2206                                 |')
    w(f, '|   \\\\  /    A nd           | Website:  www.openfoam.com                      |')
    w(f, '|    \\\\/     M anipulation  |                                                 |')
    w(f, '\\*---------------------------------------------------------------------------*/')
    w(f, 'FoamFile')
    w(f, '{')
    w(f, '    version     2.0;')
    w(f, '    format      ascii;')
    w(f, '    class       '+myClass+';')
    w(f, '    object      '+myObject+';')
    w(f, '}')
    w(f, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //')
    w(f, '')
# ---------------------------------------------------------------------------- #
def w_footer(f):
    w(f, '\n')
    w(f, '// ************************************************************************* //')
# ---------------------------------------------------------------------------- #
# --------------------------------- BLOCKMESH -------------------------------- #
# ---------------------------------------------------------------------------- #
def create_blockMeshDict(casedir,box_limits):
    with open(casedir+'/system/blockMeshDict', 'w') as f:
        k = 5
        H = box_limits[5]-box_limits[4]
        xmin = str(box_limits[0]-k*H)
        xmax = str(box_limits[1]+k*H)
        ymin = str(box_limits[2]-k*H)
        ymax = str(box_limits[3]+k*H)
        zmin = str(box_limits[4]-k*H)
        zmax = str(box_limits[5]+k*H)
        # -------------------------------------------------------------------- #
        w_header(f,'dictionary','blockMeshDict')
        w(f, """ 
scale   1;

vertices
(
""" )
        w(f, '    ( '+xmin+'  '+ymin+'  '+zmin+' )')
        w(f, '    ( '+xmax+'  '+ymin+'  '+zmin+' )')
        w(f, '    ( '+xmax+'  '+ymax+'  '+zmin+' )')
        w(f, '    ( '+xmin+'  '+ymax+'  '+zmin+' )')
        w(f, '    ( '+xmin+'  '+ymin+'  '+zmax+' )')
        w(f, '    ( '+xmax+'  '+ymin+'  '+zmax+' )')
        w(f, '    ( '+xmax+'  '+ymax+'  '+zmax+' )')
        w(f, '    ( '+xmin+'  '+ymax+'  '+zmax+' )')
        w(f, """ 
);

edges
(
);

blocks
(
        """ ) 
        # ------------------------------ grading ----------------------------- #
        mesh_division_size=25
        ds = (box_limits[5]-box_limits[4])/mesh_division_size
        xcells = str(round( (box_limits[1]-box_limits[0])/ds ))
        ycells = str(round( (box_limits[3]-box_limits[2])/ds ))
        zcells = str(round( (box_limits[5]-box_limits[4])/ds ))
        w( f, 'hex (0 1 2 3 4 5 6 7) ('+xcells+' '+ycells+' '+zcells+') simpleGrading (1 1 1)')
        # ----------------------------- next part ---------------------------- #
        w(f, """
);

boundary
(

    inlet
    {
        type wall;
        faces
        (
            (4 7 3 0)
        );
    }
    
    outlet
    {
        type wall;
        faces
        (
            (6 5 1 2)
        );
    }
    
    walls
    {
        type wall;
        faces
        (
            (0 1 5 4)
            (2 3 7 6)
            (4 5 6 7)
            (3 2 1 0)
        );
    }

);

mergePatchPairs
(
);
          """)
        w_footer(f)
# ---------------------------------------------------------------------------- #
def create_ctlDict(casedir):
    with open(casedir+'/system/controlDict', 'w') as f:
        w_header(f,'dictionary','controlDict')
        w(f, """ 
application     icoFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         0.5;

deltaT          0.005;

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
# ---------------------------------------------------------------------------- #
# ------------------------------ SNAPPY HEX MESH ----------------------------- #
# ---------------------------------------------------------------------------- #
def create_fvSch(casedir):
   with open(casedir+'/system/fvSchemes', 'w') as f:
        w_header(f,'dictionary', 'fvSchemes')
        w(f, """ 
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
    grad(p)         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linear;
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
# ---------------------------------------------------------------------------- #
def create_fvSol(casedir):
   with open(casedir+'/system/fvSolution', 'w') as f:
       w_header(f,'dictionary','fvSolution')
       w(f, """ 
solvers
{
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
}

PISO
{
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
}

         """)
       w_footer(f)
# ---------------------------------------------------------------------------- #
def get_geometry_files(stl_files):
    for file in stl_files:
        src = rel('../files/geometry/'+file+'.stl')
        dst = rel('../files/foamCase/constant/triSurface/'+file+'.stl')
        shutil.copyfile(src, dst)
# ---------------------------------------------------------------------------- #
def create_surfaceFeatureExtractDict(casedir):
    with open(casedir+'/system/surfaceFeatureExtractDict', 'w') as f:
        w_header(f,'dictionary','surfaceFeatureExtractDict')
        w(f, 'body.stl')
        w(f, '{')
        w(f, '      extractionMethod        extractFromSurface;')
        w(f, '      extractFromSurfaceCoeffs')
        w(f, '      {')
        w(f, '              includedAngle   150;')
        w(f, '      }')
        w(f, '      subsetFeatures')
        w(f, '      {')
        w(f, '          nonManifoldEdges    yes;')
        w(f, '          openEdges           yes;')
        w(f, '      }')
        w(f, '      writeObj                yes;')
        w(f, '}')
        w_footer(f)
# ---------------------------------------------------------------------------- #
def box_geometry(f,H,box_limits,k,i):
    w(f, '      refinementBox'+i)
    w(f, '      {')
    w(f, '              type box;')
    xmin = [ box_limits[0]-k*H , box_limits[2]-k*H , box_limits[4]-k*H ]
    x = str(xmin[0]) ; y = str(xmin[1]) ; z = str(xmin[2]);
    w(f, '              min ( '+x+' '+y+' '+z+' );')
    xmax = [ box_limits[1]+k*H , box_limits[3]+k*H , box_limits[5]+k*H ]
    x = str(xmax[0]) ; y = str(xmax[1]) ; z = str(xmax[2]);
    w(f, '              max ( '+x+' '+y+' '+z+' );')
    w(f, '      }')
# ---------------------------------------------------------------------------- #
def box_refinement(f,n):
        w(f, '        refinementBox'+str(n))
        w(f, '        {')
        w(f, '            mode inside;')
        w(f, '            levels (('+str(n)+' '+str(n)+'));')
        w(f, '        }')
# ---------------------------------------------------------------------------- #
def ground_region(f,r):
    w(f,'   ground')
    w(f,'   {')
    w(f,'           type searchableDisk;')
    w(f,'           origin (0 0 0);')
    w(f,'           normal (0 0 1);')
    w(f,'           radius '+str(r)+';')
    w(f,'   }')
# ---------------------------------------------------------------------------- #
def ground_refinement(f):
    w(f, '      ground')
    w(f, '      {')
    w(f, '             mode distance;')
    w(f, '             levels ((2 2));')
    w(f, '      }')
# ---------------------------------------------------------------------------- #
def ground_layers(f):
    w(f, '              ground')
    w(f, '              {')
    w(f, '                      nSurfaceLayers 10;')
    w(f, '              }')
# ---------------------------------------------------------------------------- #
def stl_geometry(f,file):
    w(f, '      '+file+'.stl')
    w(f, '      {')
    w(f, '              type triSurfaceMesh;')
    w(f, '              name '+file+'; ')
    w(f, '      }')
# ---------------------------------------------------------------------------- #
def stl_refinement_surfaces(f,file):
    w(f,'        '+file)
    w(f,'        {')
    w(f,'               level (3 3);')
    w(f,'        }')
# ---------------------------------------------------------------------------- #
def stl_refinement_layers(f,file):
    w(f, '              '+file)
    w(f, '              {')
    w(f, '                      nSurfaceLayers 10;')
    w(f, '              }')
# ---------------------------------------------------------------------------- #
def create_snappyHexMeshDict(casedir,box_limits,stl_files):
    with open(casedir+'/system/snappyHexMeshDict', 'w') as f:
        # ----------------------------- variables ---------------------------- #
        factor = 5
        H = box_limits[5]-box_limits[4]
        # -------------------------------------------------------------------- #
        n = 1 ; r=1/4
        box_scales = [1]*n
        for i in range(n):
            box_scales[i]*=r**(i)
        # -------------------------------------------------------------------- #
        w_header(f,'dictionary','snappyHexMeshDict')
        # ----------------------- initial configuration ---------------------- #
        w(f,'castellatedMesh true;')
        w(f,'snap            false;')
        w(f,'addLayers       true;')
        w(f, '\n')
        w(f, 'geometry')
        w(f, '{')
        # -------------------------- read stl files -------------------------- #
        for file in stl_files:
            stl_geometry(f,file)
        # --------------------------- insert ground -------------------------- #
        r = factor*max(box_limits[1]-box_limits[0],box_limits[3]-box_limits[2],\
            box_limits[5]-box_limits[4])
        ground_region(f,r)
        w(f, '};')
        w(f, '')
        w(f, '')
        w(f, '')
        # -------------------------------------------------------------------- #
        # --------------------- castellated mesh controls -------------------- #
        w(f, 'castellatedMeshControls')
        w(f, '{')
        w(f, '    nCellsBetweenLevels 2;')
        w(f, '    allowFreeStandingZoneFaces true;')
        w(f, '    maxGlobalCells 2e8;')
        w(f, '    maxLoadUnbalance 0.25;')
        w(f, '    maxLocalCells 1e8;')
        w(f, '    minRefinementCells 1;')
        w(f, '    planarAngle 30;')
        w(f, '    resolveFeatureAngle 30;')
        w(f, '')
        w(f, '    features')
        w(f, '    (')
        # w(f, '        { file \"body.eMesh\"; level 0; }')
        w(f, '    );')
        w(f, '    refinementSurfaces')
        w(f, '    {')
        for file in stl_files:
            stl_refinement_surfaces(f,file)
        w(f, '    }')
        w(f, '    refinementRegions')
        w(f, '    {')
        # ground_refinement(f)
        w(f, '    }')
        locx = box_limits[1]+factor/2*H
        locy = box_limits[3]+factor/2*H
        locz = box_limits[5]+factor/2*H
        w(f, '    locationInMesh ('+\
            str(locx)+' '+str(locy)+' '+str(locz)+'); ')
        w(f, '}')
        w(f, '')
        # -------------------------------------------------------------------- #
        # --------------------------- snap controls -------------------------- #
        w(f, 'snapControls')
        w(f, '{')
        w(f, '      explicitFeatureSnap true;')
        w(f, '      implicitFeatureSnap false;')
        w(f, '      multiRegionFeatureSnap false;')
        w(f, '      nFeatureSnapIter 10;')
        w(f, '      nRelaxIter 5;')
        w(f, '      nSmoothPatch 3;')
        w(f, '      nSolveIter 100;') #300
        w(f, '      tolerance 2.0;')
        w(f, '}')
        w(f, '')
        # -------------------------------------------------------------------- #
        # -------------------------- layers control -------------------------- #
        w(f, '// Settings for the layer addition.')
        w(f, 'addLayersControls')
        w(f, '{')
        # ------------------------- basic parameters ------------------------- #
        w(f, '        relativeSizes true;')
        w(f, '        expansionRatio 1.25;')
        w(f, '        featureAngle 130;')
        # w(f, '        featureAngle              270;	')
        # w(f, '        finalLayerThickness       0.5;')
        # w(f, '        firstLayerThickness       0.1;')
        w(f, '        minThickness 0.05;')
        w(f, '        nGrow 0;')
        # w(f, '        nSurfaceLayers 			3;')
        w(f, '        thickness 1;')
        # ----------------------------- advanced 1---------------------------- #
        w(f, '        maxFaceThicknessRatio 0.25;') # Stop layer growth on highly warped cells
        # ------------------- advanced: patch displacement ------------------- #
        w(f, '        nSmoothSurfaceNormals 0;')
        w(f, '        nSmoothThickness 10;')
        # ------------------ advanced: medial axis analysis ------------------ #
        w(f, '        maxThicknessToMedialRatio 0.3;') #
        w(f, '        minMedialAxisAngle 120;')
        w(f, '        minMedianAxisAngle 120;')
        # w(f, '        nMedialAxisIter             10;')
        # w(f, '        nSmoothDisplacement 		    90;')
        # w(f, '        nSmoothNormals 				15;')
        w(f, '        nSmoothNormals 30;')
        # -------------------------- mesh shrinking -------------------------- #
        w(f, '        nLayerIter 50;')
        w(f, '        nRelaxedIter 20;')
        w(f, '        nRelaxIter 5;') #5 or 25
        w(f, '        slipFeatureAngle 30;') #130 or 30
        w(f, '        nBufferCellsNoExtrude 0;') 
        
        # w(f, '        additionalReporting         true;')
        # ------------------------------ unknown ----------------------------- #
        # w(f, '        nBufferCellsNoExtrude       0;')
        # ------------------------------ layers ------------------------------ #
        w(f, '        layers')
        w(f, '        {')
        for file in stl_files:
            stl_refinement_layers(f,file)
        # ground_layers(f)
        w(f, '        }')
        w(f, '}')
        w(f, '')
        # ----------------------- mesh quality controls ---------------------- #
        w(f, 'meshQualityControls')
        w(f, '{')
        w(f, '      #include \"meshQualityDict\"')
        w(f, '      relaxed')
        w(f, '      {')
        w(f, '              maxNonOrtho 80;')
        w(f, '      }')
        w(f, '      nSmoothScale 4;')
        w(f, '      errorReduction 0.75;')
        w(f, '}')
        w(f, 'debug 0;')
        w(f, 'mergeTolerance 1e-6;')
        w_footer
# ---------------------------------------------------------------------------- #
def create_meshQualityDict(casedir):
    with open(casedir+'/system/meshQualityDict', 'w') as f:
        w_header(f, 'dictionary', 'meshQualityDict')
        # w(f, 'maxNonOrtho 80;')
        w(f, 'maxNonOrtho 70;')
        w(f, 'maxBoundarySkewness 4;	//original 20')
        w(f, 'maxInternalSkewness 4;')
        w(f, 'maxConcave 70;')
        w(f, 'minVol 1e-13;')
        w(f, 'minTetQuality -1e30; 	//for best layer insertion	')
        w(f, 'minArea -1;')
        w(f, 'minTwist 0.02;')
        w(f, 'minDeterminant 0.001;')
        w(f, 'minFaceWeight 0.05;')
        w(f, 'minVolRatio 0.01;')
        w(f, 'minTriangleTwist -1;')
        w_footer(f)
# ---------------------------------------------------------------------------- #
def create_decomposeParDict(casedir):
    with open(casedir+'/system/decomposeParDict', 'w') as f:
        w_header(f,'dictionary','decomposeParDict')
        w(f, 'numberOfSubdomains 26;')
        w(f, 'method          scotch;')
        w_footer(f)
# # ---------------------------------------------------------------------------- #
# # -------------------------------- FOR SOLVER -------------------------------- #
# # ---------------------------------------------------------------------------- #
def createDirsAndFiles(casedir, box_limits, stl_files):
    # ----------------------------- for blockMesh ---------------------------- #
    createDirs(casedir)
    create_blockMeshDict(casedir,box_limits)
    create_ctlDict(casedir)
    # --------------------------- for snappyHexMesh -------------------------- #
    create_fvSch(casedir)
    create_fvSol(casedir)
    get_geometry_files(stl_files)
    # create_surfaceFeatureExtractDict(casedir)
    create_snappyHexMeshDict(casedir,box_limits,stl_files)
    create_meshQualityDict(casedir)
    create_decomposeParDict(casedir)