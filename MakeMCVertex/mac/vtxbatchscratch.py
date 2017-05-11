import sys, os

import ROOT
from larlite import larlite as fmwk

INPATH  = '/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/MC/sel2_MCC8_BNB_COSMIC_FULL_LARLITE_v01/'
OUTPATH = INPATH #'/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/MC/sel2_MCC8_BNB_COSMIC_FULL_LARLITE_compressed/'

ctr = 1

while (ctr < 348):

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

    cmd = 'python makemcvertex.py %s %s %s'%(ll_files[0],ll_files[1],ll_files[2])
    print cmd
    os.system(cmd)
    cmd = 'mv larlite_vertex.root %s/larlite_vertex_%04i.root'%(OUTPATH,ctr)
    print cmd
    os.system(cmd)

    ctr += 1

    #if (ctr >= 1):
    #    break
    

sys.exit()

