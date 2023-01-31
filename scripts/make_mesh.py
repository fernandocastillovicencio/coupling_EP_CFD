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
    # ----------------- get the list of the geomtry stl files ---------------- #
    with open (rel('../files/geometry/stl_files.pickle'), 'rb') as fp:
        stl_files = pickle.load(fp)
    # ------------------------------------------------------------------------ #
    # ----- CREATE NECESSARY FILES AND FOLDERS FOR MESHING AND PROCESSING ---- #
    os.system('cd files/foamCase && rm -rf *')
    createDirsAndFiles(casedir,box_limits,stl_files)
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
    # surfaceFeatureExtract -noFunctionObjects
    decomposePar
    mpirun --use-hwthread-cpus -np 26 snappyHexMesh -parallel -noFunctionObjects -overwrite > log.snappyHexMesh
    reconstructParMesh -constant
    mpirun --use-hwthread-cpus -np 26 checkMesh -parallel -latestTime > log.checkMesh
    paraFoam
    """
    single_snappy_commands="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    # surfaceFeatureExtract -noFunctionObjects
    snappyHexMesh -noFunctionObjects -overwrite > log.snappyHexMesh
    # checkMesh
    paraFoam
    """
    os.system(multi_snappy_commands)
    # ------------------------------------------------------------------------ #
    cleanCmd="""
    . /usr/lib/openfoam/openfoam2206/etc/bashrc
    cd files/foamCase
    foamCleanTutorials
    foamCleanPolyMesh
    rm -rf constant/extendedFeatureEdgeMesh
    """
    os.system(cleanCmd)
    