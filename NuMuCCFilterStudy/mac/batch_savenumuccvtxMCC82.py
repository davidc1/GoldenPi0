#
# Example PyROOT script to run analysis module, ana_base.
# The usage is same for inherited analysis class instance.
#

# Load libraries
import os, ROOT, sys
from ROOT import *
from larlite import larlite as fmwk


INPATH = '/pnfs/uboone/scratch/users/davidc1/MCC8_2_BNB_5E19_Sel2Lite_test_out/v06_45_00/'

OUTPATH = '/uboone/data/users/davidc1/mcc82/5e19/'

dirs = os.listdir(INPATH)

ctr = 0

for directory in dirs:

    DIRPATH = INPATH + directory

    if (os.path.isdir(DIRPATH) == False):
        continue

    cmd = 'python savesel2vtxtrkMCC82.py %s/larlite_*.root larlite_%04i.root'%(DIRPATH,ctr)
    os.system(cmd)

    print cmd
    ctr += 1

        
sys.exit(0);
