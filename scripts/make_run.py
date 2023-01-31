import os
def mkRun():
    # --------------------------- Prepare directory -------------------------- #
    prepare_commands="""
    # cd files/foamCase
    rm -rf 0 > /dev/null 2>&1
    # cp -r 0.orig 0
    # """
    # os.system(prepare_commands)
    # ------------------------------ Run Solver ------------------------------ #
    solver_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    cp -r 0.orig 0
    rm -rf processor*
    # simpleFoam > log.foamRun
    decomposePar -force
    mpirun --use-hwthread-cpus -np 26 simpleFoam -parallel > log.foamRun
    reconstructPar
    # paraFoam
    """
    os.system(solver_commands)
    # ------------------------------------------------------------------------ #
    # ------------------------------------------------------------------------ #
    # ---------------------------- Clean directory --------------------------- #
    cleanCmd="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    foamCleanTutorials
    foamCleanPolyMesh
    rm -rf constant/extendedFeatureEdgeMesh
    """
    os.system(cleanCmd)
    