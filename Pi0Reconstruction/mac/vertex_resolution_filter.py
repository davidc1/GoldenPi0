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
my_proc.set_ana_output_file("vertex_resolution_ana.root");
my_proc.set_output_file("vertex_resolution_out.root");

vtx = fmwk.VertexResolution()
vtx.setResolutionFilter(4.)
my_proc.add_process(vtx)

my_proc.enable_filter(True)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

