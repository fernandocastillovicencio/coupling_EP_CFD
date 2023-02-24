import os
from scripts.mesh_functions import rel
from scripts.solver_functions import create_solver_files

def mkSolver():    
    # -------------------------- Create solver files ------------------------- #
    casedir = rel('../files/foamCase')
    create_solver_files(casedir)
    
     
    # --------------------------- Prepare directory -------------------------- #
    prepare_commands="""
    # cd files/foamCase
    rm -rf 0 > /dev/null 2>&1
    # cp -r 0.orig 0
    # """
    os.system(prepare_commands)
    # ------------------------------ Run Solver ------------------------------ #
    # ------------------------------------------------------------------------ #
    single_solver_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    cp -r 0.orig 0
    # rm -rf processor*    
    potentialFoam
    simpleFoam
    paraFoam
    """
    # ------------------------------------------------------------------------ #
    multi_solver_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    cp -r 0.orig 0
    rm -rf processor*
    potentialFoam
    decomposePar -force
    mpirun --use-hwthread-cpus -np 26 simpleFoam -parallel > log.foamRun
    reconstructPar
    touch foamCase.foam
    """
    # ------------------------------------------------------------------------ #
    os.system(multi_solver_commands)
    # ------------------------------------------------------------------------ #
    # ------------------------------------------------------------------------ #
    # ---------------------------- Clean directory --------------------------- #
    cleanCmd="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    # foamCleanTutorials
    # foamCleanPolyMesh
    rm -rf constant/extendedFeatureEdgeMesh 1* 2* 3* 4* 5* 6* 7* 8* 9* processor* 0
    """
    os.system(cleanCmd)
    