import os
from scripts.make_geometry import mkGeometry
from scripts.make_mesh import mkMesh
from scripts.make_run import mkRun
# --------------------------------- PREAMBLE --------------------------------- #
os.chdir('/home/fernando/workspace/research/phd/webapp/wip/coupling_EP_CFD')
# --------------------------------- GEOMETRY --------------------------------- #
mkGeometry()
# ----------------------------------- MESH ----------------------------------- #
mkMesh()
# ------------------------------------ RUN ----------------------------------- #
mkRun()