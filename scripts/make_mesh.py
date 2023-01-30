def mkMesh():
    # ----------------------- import standard libraries ---------------------- #
    import os
    import pickle
    # ------------------------ import custom libraries ----------------------- #
    from scripts.mesh_functions import rel,createDirsAndFiles
    # ---------------------- changing working directory ---------------------- #
    # os.chdir(rel('../files/foamCase'))
    casedir = rel('../files/foamCase')
    # ------------------------ get bounding box limits ----------------------- #
    with open (rel('../files/geometry/box_limits.pickle'), 'rb') as fp:
        box_limits = pickle.load(fp)
    # ------------------------------------------------------------------------ #
    # ----- CREATE NECESSARY FILES AND FOLDERS FOR MESHING AND PROCESSING ---- #
    os.system('cd files/foamCase && rm -rf *')
    createDirsAndFiles(casedir,box_limits)
    # ------------------------------ PREPARE DIR ----------------------------- #
    prepare_commands=""" 
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    rm -rf processor*
    foamCleanTutorials
    foamCleanPolyMesh
    """
    os.system(prepare_commands)
    # print('PASSO1')
    # ------------------------------- BLOCKMESH ------------------------------ #
    blockMesh_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    blockMesh
    # paraFoam
    """
    os.system(blockMesh_commands)
    # ----------------------------- SNAPPYHEXMESH ---------------------------- #
    multi_snappy_commands= """
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    surfaceFeatureExtract -noFunctionObjects
    decomposePar
    mpirun --use-hwthread-cpus -np 24 snappyHexMesh -parallel -noFunctionObjects -overwrite
    reconstructParMesh -constant
    mpirun --use-hwthread-cpus -np 24 checkMesh -parallel -latestTime
    paraFoam
    foamCleanTutorials
    foamCleanPolyMesh
    """
    single_snappy_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    surfaceFeatureExtract -noFunctionObjects
    snappyHexMesh -noFunctionObjects -overwrite
    paraFoam
    foamCleanTutorials && foamCleanPolyMesh && ls && rm -rf constant/extendedFeatureEdgeMesh
    """
    os.system(multi_snappy_commands)

    