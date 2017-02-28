

import ROOT
import larlite
from ROOT import evd, larlite, larutil
from evdmanager import geometry

# Needed for the process to read data and filter noise:

# Needed for array manipulation and color mapping:
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

import sys
import time

import numpy

import math

import Image
import ImageFont
import ImageDraw
import matplotlib
from matplotlib import pyplot as plt

# This defines the official argoneut color map
# in a way usable by matplotlib colormaps
collection_cdict = {'red':
                    ((0.0, 030./255, 030./255),
                     (0.2, 000./255, 000./255),
                     (0.3, 030./255, 030./255),
                     (1.0, 255./255, 255./255)),
                    'green':
                    ((0.0, 030./255, 030./255),
                     (0.2, 000./255, 000./255),
                     (0.3, 255./255, 255./255),
                     (1.0, 000./255, 000./255)),
                    'blue':
                    ((0.0, 255./255, 255./255),
                     (1.0, 000./255, 000./255))
                    }

induction_cdict = {'red':
                    ((0.0, 030./255, 030./255),
                     (0.32, 000./255, 000./255),
                     (0.8, 000./255, 000./255),
                     (1.0, 255./255, 255./255)),
                    'green':
                    ((0.0, 030./255, 030./255),
                     (0.32, 255./255, 255./255),
                     (0.8, 255./255, 255./255),
                     (1.0, 000./255, 000./255)),
                    'blue':
                    ((0.0, 255./255, 255./255),
                     (0.32, 255./255, 255./255),
                     (0.8, 000./255, 000./255),
                     (1.0, 000./255, 000./255))
                    }


# This takes the data and maps it on to an RGBalpha array
# Create a color map using the official color scheme
coll_map = LinearSegmentedColormap('collection_map', collection_cdict)
ind_map = LinearSegmentedColormap('induction_map', induction_cdict)
# register the color map:
plt.register_cmap(cmap=coll_map)
plt.register_cmap(cmap=ind_map)

# this is a utility that maps arrays into rgb arrays
collectionMapper = cm.ScalarMappable()
collectionMapper.set_cmap('collection_map')
# Set the levels of the color mapping:
collectionMapper.set_clim(vmin=-25, vmax=100)

inductionMapper = cm.ScalarMappable()
inductionMapper.set_cmap('induction_map')
# Set the levels of the color mapping:
inductionMapper.set_clim(vmin=-25, vmax=100)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def getWireRange(vtx):

    geomH = larutil.GeometryHelper.GetME()
    
    vtx_w = 200
    vtx_t = 200

    buffer_w = int(70 / 0.3)
    buffer_t = int(70 / 0.057) 
    
    vtxXYZ = ROOT.TVector3( vtx.X(), vtx.Y(), vtx.Z() )
    vtxWT  = geomH.Point_3Dto2D(vtxXYZ, 2)
    vtx_w = int(vtxWT.w / 0.3)
    vtx_t = int(vtxWT.t / 0.057) + 800
    #print 'vertex being projected to wire : %i'%(vtx_w)
    #print 'vertex being projected to time : %i'%(vtx_t)

    wire_min = vtx_w - buffer_w
    if (wire_min < 0):
        wire_min = 0
    wire_max = vtx_w + buffer_w
    if (wire_max > 3456):
        wire_max = 3456
    time_min = vtx_t - buffer_t
    if (time_min < 0):
        time_min = 0
    time_max = vtx_t + buffer_t
    if (time_max > 6000):
        time_max = 6000

    return (wire_min, wire_max), (time_min, time_max)
    

def drawEvent(wire_range, time_range, coll_image, output_folder, run, event):

    #print bcolors.OKBLUE, "Initializing data ...", bcolors.ENDC

    w_width = wire_range[1]-wire_range[0]
    t_width = time_range[1]-time_range[0]
    t_cm = t_width * 0.057
    print 'wire width: ',w_width
    print 'time width: ',t_width
    coll_image = coll_image[wire_range[0] : wire_range[1] , time_range[0] : time_range[1] ]

    coll_image = numpy.flipud(numpy.transpose(coll_image))

    coll_image = numpy.repeat( coll_image, 6, axis=1)

    # This function maps the data onto the color scheme and makes and RGB array
    #print 'map collection-plane array'
    coll_array = collectionMapper.to_rgba(coll_image)
    coll_array = numpy.uint8(coll_array*255)

    collfileName = "Run_%04i_Event_%05i_collection.png"%(run,event)
    #print collfileName

    coll_image = Image.fromarray(coll_array)
    draw = ImageDraw.Draw(coll_image)

    run_info = 'RUN %i EVENT %i'%(run,event)
    fontsize = 60
    font = ImageFont.truetype('courier10.ttf', 60)
    draw.text((100, t_width  - fontsize - 100),run_info,(255,255,255), font=font)
    draw.text((t_width - 300, t_width  - fontsize - 100),'%i cm'%(int(t_cm*(450./t_width) ) ) ,(255,255,255), font=font)
    draw.line( [ (t_width - 430, t_width - fontsize - 170), (t_width+20, t_width - fontsize - 170) ], "white", width=10 )

    coll_image.save(output_folder + collfileName)


def main():


    larutil.LArUtilManager.Reconfigure(larlite.geo.kMicroBooNE)

    drawWire = evd.DrawWire()
    #print 'initize image processing'
    drawWire.initialize()
    
    ana_processor = larlite.ana_processor()
    ana_processor.set_io_mode(larlite.storage_manager.kREAD)
    ana_processor.add_process(drawWire)
    ana_processor.add_input_file( sys.argv[1] )


    manager = larlite.storage_manager()
    manager.reset()
    manager.add_in_filename( sys.argv[2] )
    manager.set_io_mode(larlite.storage_manager.kREAD)
    manager.open()

    folder = "/uboone/data/users/davidc1/numuCCFilterOutput/selection2/handscanned/images/"

    froot = ROOT.TFile(sys.argv[1])
    tree  = froot.Get("larlite_id_tree")
    nentries = tree.GetEntries()
    froot.Close()
    
    i = 0
    for entry in xrange(nentries):
        #print 'processing event %i'%entry

        ana_processor.process_event( entry )
        manager.go_to( entry )

        # grab the vertex:
        vertex_v = manager.get_data(larlite.data.kVertex,"numuCC_vertex")
        if (len(vertex_v) != 1):
            print 'Number of vertices != 1 : continue...'
            continue
        vtx = vertex_v[0]
        wire_range, time_range = getWireRange(vtx)

        run = manager.run_id()
        event = manager.event_id()

        image = drawWire.getArrayByPlane(2)
        
        drawEvent(wire_range, time_range, image, folder, run, event)


if __name__ == '__main__':

    main()
