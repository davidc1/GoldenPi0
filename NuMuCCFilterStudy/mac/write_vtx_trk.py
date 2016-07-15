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
my_proc.set_ana_output_file("calib_ana.root");

# Specify data output root file name
my_proc.set_output_file("larlite_neutrino.root")

search = fmwk.SearchPFPartHierarchy()
my_proc.add_process(search)

my_proc.set_data_to_write(fmwk.data.kVertex,'numuCC_vertex')
my_proc.set_data_to_write(fmwk.data.kTrack,'numuCC_track')

print
print  "Finished configuring ana_processor. Start event loop!"
print

ctr = 0

while my_proc.process_event():

    print 'Counter : ',ctr

    usrinput = raw_input("Hit Enter: next evt  || int: go to event number ||  q: exit viewer\n")
    if ( usrinput == "q" ):
        sys.exit(0)

    ctr += 1

sys.exit()

