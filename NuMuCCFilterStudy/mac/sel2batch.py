import sys, os

import ROOT
from larlite import larlite as fmwk

INPATH  = '/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/MC/sel2_MCC8_BNB_COSMIC_FULL_LARLITE_v04/'
OUTPATH = INPATH

ctr = 1

while (ctr < 961):

    ll_files = []

    fname = '%s/larlite_mcinfo_%04i.root'%(INPATH,ctr)
    if (os.path.isfile(fname) == True):
        ll_files.append(fname)

    fname = '%s/larlite_reco2d_%04i.root'%(INPATH,ctr)
    if (os.path.isfile(fname) == True):
        ll_files.append(fname)

    fname = '%s/larlite_pandora_%04i.root'%(INPATH,ctr)
    if (os.path.isfile(fname) == True):
        ll_files.append(fname)

    if (len(ll_files) != 3): 
        ctr += 1
        continue

    cmd = 'python savesel2vtxtrk.py %s %s %s'%(ll_files[0],ll_files[1],ll_files[2])
    print cmd
    os.system(cmd)
    cmd = 'mv larlite_numucc.root %s/larlite_numucc_%04i.root'%(OUTPATH,ctr)
    print cmd
    os.system(cmd)

    ctr += 1

sys.exit()

