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
my_proc.enable_filter(True)

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("ana_out.root");
my_proc.set_output_file("out.root");

# Specify data output root file name
#my_proc.set_output_file("rawdigit_hits.root")

search = fmwk.SearchPFPartHierarchy()
search.SetVerbose(False)
search.SetFilterShowers(True,2)
search.SetFilterTracks(False)
my_proc.add_process(search)

#my_proc.set_data_to_write(fmwk.data.kHit,'rawhit')
#my_proc.set_data_to_write(fmwk.data.kRawDigit,'daq')

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

