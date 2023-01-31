# ----------------------------- import libraries ----------------------------- #
import numpy as np
import pickle
from eppy.function_helpers import getcoords
from read_case import get_idf_file
import os
# os.chdir('/home/fernando/workspace/research/phd/webapp/wip/coupling_EP_CFD/')
# ---------------------------------------------------------------------------- #
cmd = """ 
rm -rf files/geometry/*
mkdir -p files/geometry
"""
os.system(cmd)
# -------------------------- import salome libraries ------------------------- #
import salome
salome.salome_init()
from salome.geom import geomBuilder
geompy = geomBuilder.New()
# ---------------------------------------------------------------------------- #
# --------------------------------- functions -------------------------------- #
vxname = lambda i: 'vertex'+str(i).zfill(2)
edname = lambda i: 'edge'+str(i).zfill(2)
# ---------------------------------------------------------------------------- #
geom_dir = 'files/geometry/'
# --------------------------------- boxlimits -------------------------------- #
def calc_box_limits(x,box_limits):
    for i in range(len(x)):
        box_limits[0] = min( x[i][0] , box_limits[0] ) #xmin
        box_limits[1] = max( x[i][0] , box_limits[1] ) #xmax 
        box_limits[2] = min( x[i][1] , box_limits[2] ) #ymin
        box_limits[3] = max( x[i][1] , box_limits[3] ) #ymax
        box_limits[4] = min( x[i][2] , box_limits[4] ) #zmin
        box_limits[5] = max( x[i][2] , box_limits[5] ) #zmax
    return box_limits
# ----------------------------- Geometry function ---------------------------- #
def mkGeometry(idf):
    box_limits=np.zeros(6)
    stl_files = []
    
    for surface in idf.idfobjects['BuildingSurface:Detailed']:
        # ------------------------------ labels ------------------------------ #
        surface_name = surface['Name'][6:].replace(' ','_')
        # ---------------------- get vertex coordinates ---------------------- #
        x = getcoords(surface)
        # ---------------------------- box limits ---------------------------- #
        box_limits = calc_box_limits(x,box_limits)
        # ----------------------------- VERTICES ----------------------------- #
        for i in range(len(x)):
            # -------------------------- create vertices ------------------------- #
            cmd = vxname(i)+' = geompy.MakeVertex('+str(x[i][0])+','+str(x[i][1])+\
                ' , '+str(x[i][2])+' )'
            exec(cmd)
        # ------------------------------- EDGES ------------------------------ #
        for i in range(len(x)):
            # -------------- select the vertices to make an edge ------------- #
            v1 = vxname(i);
            v2 = vxname(i+1) if i<3 else vxname(0)
            # ------------------------- create Edges ------------------------- #
            cmd = edname(i)+'= geompy.MakeEdge('+v1+', '+v2+')'
            exec(cmd)
        # ----------------------------- SURFACES ----------------------------- #
        cmd = surface_name+' = geompy.MakeFaceWires(['
        for i in range(len(x)): 
            cmd+=edname(i)+','
        cmd+='] , 1 )'
        exec(cmd)
        # cmd = "geompy.addToStudy("+surface_name+",\'"+surface_name+"\')"
        # exec(cmd)
        cmd = "geompy.ExportSTL("+surface_name+", \'"+geom_dir+surface_name+\
            ".stl\'"+", True, 0.001, True)"   
        exec(cmd)
        stl_files.append(surface_name)
    # ---------------------------- save box limits --------------------------- #
    with open(geom_dir+'box_limits.pickle', 'wb') as f:
        pickle.dump(box_limits, f)
    # ------------------------------------------------------------------------ #
    # ----------------------------- FENESTRATION ----------------------------- #
    for fenestration in idf.idfobjects['FenestrationSurface:Detailed']:
        fenestration_name = fenestration['Name'].split(':')[2]
        building_surf_name = fenestration['Name'].split(':')[1]
        # -------------------------- get coordinates ------------------------- #
        x = getcoords(fenestration)
        # ------------------------------------------------------------------------ #
        # ----------------------------- VERTICES ----------------------------- #
        for i in range(len(x)):
            # -------------------------- create vertices ------------------------- #
            cmd = vxname(i)+' = geompy.MakeVertex('+str(x[i][0])+','+str(x[i][1])+\
                ' , '+str(x[i][2])+' )'
            exec(cmd)
        # ------------------------------- EDGES ------------------------------ #
        for i in range(len(x)):
            # -------------- select the vertices to make an edge ------------- #
            v1 = vxname(i);
            v2 = vxname(i+1) if i<3 else vxname(0)
            # ------------------------- create Edges ------------------------- #
            cmd = edname(i)+'= geompy.MakeEdge('+v1+', '+v2+')'
            exec(cmd)
        # ----------------------------- SURFACES ----------------------------- #
        cmd = fenestration_name+' = geompy.MakeFaceWires(['
        for i in range(len(x)): 
            cmd+=edname(i)+','
        cmd+='] , 1 )'
        exec(cmd)
        # ------------------------------------------------------------ #
        # cmd = "geompy.addToStudy("+fenestration_name+",\'"+fenestration_name+"\')"
        # exec(cmd)
        cmd = "geompy.ExportSTL("+fenestration_name+", \'"+geom_dir+\
            fenestration_name+".stl\'"+", True, 0.001, True)"   
        exec(cmd)
        stl_files.append(fenestration_name)
        # -------------------------------------------------------------------- #
        # ------------- Combine Fenestration and BuildingSurfaces ------------ #
        cmd = building_surf_name+' = geompy.MakeCutList('+building_surf_name+\
                ', ['+fenestration_name+'], True)'
        exec(cmd)
        # cmd = "geompy.addToStudy("+building_surf_name+",\'"+building_surf_name+"\')"
        # exec(cmd)
        cmd = "geompy.ExportSTL("+building_surf_name+", \'"+geom_dir+\
            building_surf_name+".stl\'"+", True, 0.001, True)"   
        exec(cmd)
        # ----------------------- Save list of surfaces ---------------------- #
        # ---------------------------- save box limits --------------------------- #
    with open(geom_dir+'stl_files.pickle', 'wb') as f:
        pickle.dump(stl_files, f)
    #    # -------------------------------------------------------------------- #
# --------------------------------- READ FILE -------------------------------- #
idf = get_idf_file()
# ---------------------------------------------------------------------------- #
# --------------------------------- GEOMETRY --------------------------------- #
mkGeometry(idf)