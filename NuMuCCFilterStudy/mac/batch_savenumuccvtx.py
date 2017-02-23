#
# Example PyROOT script to run analysis module, ana_base.
# The usage is same for inherited analysis class instance.
#

# Load libraries
import os, ROOT, sys
from ROOT import *
from larlite import larlite as fmwk


minfile = int(sys.argv[-2])
maxfile = int(sys.argv[-1])


for n in xrange(minfile,maxfile):

    # Create ana_processor instance
    my_proc=fmwk.ana_processor()

    # Specify IO mode
    my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

    my_proc.add_input_file('/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/NuMuCC_Sel2_reco1and2_larlite/larlite_reco2d_%04i.root'%n)

    my_proc.set_ana_output_file("ana.root")
    my_proc.set_output_file("/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/NuMuCC_Sel2_reco1and2_larlite//larlite_neutrino_%04i.root"%n)
    
    ana = fmwk.SaveNumuCCVtx()
    
    fin = open('/pnfs/uboone/persistent/users/oscillations_group/GoldenPi0/NuMuCC_Sel2/numucc_sel2_5E19_info.txt')
        
    for line in fin:
        words = line.split()
        if (len(words) != 12): continue
        run   = int(words[0])
        event = int(words[2])
        x     = float(words[3])
        y     = float(words[4])
        z     = float(words[5])
        ana.addVertex(run,event,x,y,z)

    my_proc.add_process(ana)

    my_proc.set_data_to_write(fmwk.data.kVertex,'numuCC_vertex')

    my_proc.run()

sys.exit(0);
