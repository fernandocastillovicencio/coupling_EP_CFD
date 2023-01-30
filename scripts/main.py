# ------------------------- import standard libraries ------------------------ #
import sys
# ------------------------------ add local path ------------------------------ #
import os
def set_proj_dir():
    os.chdir('/home/fernando/workspace/research/phd/webapp/wip/coupling_EP_CFD')
sys.path.append('scripts')
# -------------------------- import local libraries -------------------------- #
from read_case import get_idf_file
from make_geometry import mkGeometry
# ------------------------------- get IDF file ------------------------------- #
set_proj_dir()
idf1 = get_idf_file()
# ------------------------------- Make Geometry ------------------------------ #
mkGeometry(idf1)
# ------------------------------- Make Geometry ------------------------------ 