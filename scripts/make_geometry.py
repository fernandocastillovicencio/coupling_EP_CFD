# ----------------------------- import libraries ----------------------------- #
import numpy as np
from eppy.function_helpers import getcoords
# ---------------------------------------------------------------------------- #

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
# ----------------------------- Geometry function ---------------------------- #
def mkGeometry(idf):
    stl_files = []
    
    for surface in idf.idfobjects['BuildingSurface:Detailed']:
        # ------------------------------ labels ------------------------------ #
        surface_name = surface['Name'][6:].replace(' ','_')
        # ---------------------- get vertex coordinates ---------------------- #
        x = getcoords(surface)
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
            # cmd = "geompy.addToStudy( "+edname(i)+", \'"+edname(i)+"\' )"
            # exec(cmd)
        # ----------------------------- SURFACES ----------------------------- #
        cmd = surface_name+' = geompy.MakeFaceWires(['
        for i in range(len(x)): 
            cmd+=edname(i)+','
        cmd+='] , 1 )'
        exec(cmd)
        
        cmd = "geompy.addToStudy("+surface_name+",\'"+surface_name+"\')"
        exec(cmd)
        
        
        cmd = "geompy.ExportSTL("+surface_name+", \'files/"+surface_name+".stl\'"+\
            ", True, 0.001, True)"   
        exec(cmd)
        stl_files.append(surface_name)
    
    # ------------------------------------------------------------------------ #
    # ----------------------------- FENESTRATION ----------------------------- #
    for fenestration in idf.idfobjects['FenestrationSurface:Detailed']:
        fenestration_name = fenestration['Name'].split(':')[2]
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
        cmd = "geompy.ExportSTL("+fenestration_name+", \'files/"+fenestration_name+".stl\'"+\
            ", True, 0.001, True)"   
        exec(cmd)
        stl_files.append(fenestration_name)
        # -------------------------------------------------------------------- #
        # -------------------------- Merge STL files ------------------------- #
        with open("files/geometry.stl", "w") as f:
            for surface in stl_files:
                # print(surface)
                i = 0
                for line in open('files/'+surface+'.stl', 'r').readlines():
                    i+=1
                    if 'solid' in str(line):
                        line =  line.strip()+' '+surface+'\n'
                    f.write(line)
