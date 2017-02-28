import sys, os

for x in xrange(336):
    
    cmd = 'python dataReleaseMaker.py /pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/HandScannedPi0/selection2/larlite/larlite_wire_%04i.root /pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/HandScannedPi0/selection2/larlite/larlite_neutrino_%04i.root'%(x,x)
    print cmd
    os.system(cmd)
