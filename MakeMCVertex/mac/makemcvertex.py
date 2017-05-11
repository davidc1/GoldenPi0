import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


import ROOT
from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("ana.root");

# Specify data output root file name
my_proc.set_output_file("larlite_vertex.root")

mcvtx = fmwk.MakeMCVertex()
mcvtx.SetXOffset(0.7)
my_proc.add_process(mcvtx)

my_proc.set_data_to_write(fmwk.data.kMCShower,    "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kMCTrack,     "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kMCTruth,     "generator"     )
my_proc.set_data_to_write(fmwk.data.kVertex,      "mcvertex"      )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"       )

#my_proc.set_data_to_write(fmwk.data.kVertex,'mcvertex')

print
print  "Finished configuring ana_processor. Start event loop!"
print

ctr = 0

my_proc.enable_filter(True)

my_proc.run()

sys.exit()

