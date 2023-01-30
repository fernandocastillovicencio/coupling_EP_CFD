#!/bin/bash

# --------------------------------- GEOMETRY --------------------------------- #
rm -rf files/geometry/*
mkdir -p files/geometry/components
salome -t scripts/main.py

# ----------------------------------- MESH ----------------------------------- #
# ------------------------------ prepare folder ------------------------------ #
# # workon energyenv
# rm -rf files/geometry/*
# mkdir -p files/foamCase
# cd files/foamCase
# # --------------------------- create files and dirs -------------------------- #
# python3 scripts/make_mesh.py
# ------------------------------ clean foam case ----------------------------- #
# rm -rf 0 > /dev/null 2>&1
# rm *.obj > /dev/null 2>&1
# rm -rf processor*
# rm -rf sol* figs* */*/*.eMesh constant/extended* geometry
# foamCleanTutorials
# foamCleanPolyMesh
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
cd ../..
