import os
from scripts.mesh_functions import rel
from scripts.solver_functions import create_solver_files

def mkSolver():    
    # -------------------------- Create solver files ------------------------- #
    casedir = rel('../files/foamCase')
    create_solver_files(casedir)
    
     
    # --------------------------- Prepare directory -------------------------- #
    prepare_commands="""
    cd files/foamCase
    rm -rf [1-9]* 0.[0-9]* processor* 0 h* postProcessing
    rm -rf 0 > /dev/null 2>&1
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
    buoyantSimpleFoam > log.foamRun
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
    mpirun --use-hwthread-cpus -np 26 buoyantSimpleFoam -parallel > log.foamRun
    reconstructPar > log.reconstructPar
    touch foamCase.foam
    rm -rf processor*
    # paraFoam
    """
    # ------------------------------------------------------------------------ #
    os.system(multi_solver_commands)
    # ------------------------------------------------------------------------ #
    # ------------------------------------------------------------------------ #
    # ---------------------------- Post Processing --------------------------- #
    postprocessing_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    postProcess -func 'patchAverage(name=Wall001,heatTransferCoeff(T))' -latestTime
    postProcess -func 'patchAverage(name=Wall002,heatTransferCoeff(T))' -latestTime
    postProcess -func 'patchAverage(name=Wall003,heatTransferCoeff(T))' -latestTime
    postProcess -func 'patchAverage(name=Wall004,heatTransferCoeff(T))' -latestTime
    postProcess -func 'patchAverage(name=Roof001,heatTransferCoeff(T))' -latestTime
    postProcess -func 'patchAverage(name=Win001,heatTransferCoeff(T))' -latestTime
    paraFoam
    """
    # ------------------------------------------------------------------------ #
    os.system(postprocessing_commands)
    # ------------------------------------------------------------------------ #
    # ------------------------------------------------------------------------ #
    # ---------------------------- Clean directory --------------------------- #
    cleanCmd="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    # foamCleanTutorials
    # foamCleanPolyMesh
    rm -rf constant/extendedFeatureEdgeMesh [1-9]* 0.[0-9]*  processor* 0
    """
    # os.system(cleanCmd)
    