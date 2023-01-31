# ----------------------------- import libraries ----------------------------- #
# import numpy as np
# import pickle
# from eppy.function_helpers import getcoords
# -------------------------- import local libraries -------------------------- #
# from scripts.read_case import get_idf_file
# from make_geometry import mkGeometry
# -------------------------- import salome libraries ------------------------- #
# import salome
# salome.salome_init()
# from salome.geom import geomBuilder
# geompy = geomBuilder.New()
# ------------------------------- get IDF file ------------------------------- #
# idf = get_idf_file()
# print(type(idf))
# # ------------------------------- Make Geometry ------------------------------ #
# mkGeometry(idf1)
# # ------------------------------- Make Geometry ------------------------------ 
# def mkGeometry(idf):
#     print(idf.idfobjects['BuildingSurface:Detailed'])
            

import os
def mkGeometry():
    salome_command = 'salome -t scripts/geometry_functions.py'
    # salome_command = 'ls'
    os.system(salome_command)