#ifndef LARLITE_VERTEXEXAMPLE_CXX
#define LARLITE_VERTEXEXAMPLE_CXX

#include "VertexExample.h"
// without this line the module will not know about vertices
// and therfore we won't be able to read vertex info from
// the input file
#include "DataFormat/vertex.h" 
// include utility functions from GeoHelper
#include "LArUtil/GeometryHelper.h"
// include LArProperties utilities
#include "LArUtil/LArProperties.h"

namespace larlite {

  VertexExample::VertexExample()
  {
    // setting the producer name for the vertex we want to load
    // this can be overwritten by the SetVertexProducer function
    // at run time.
    _vertex_producer = "numuCC_vertex";
    _name="VertexExample";
    _fout=0;
  }

  bool VertexExample::initialize() {

    return true;
  }
  
  bool VertexExample::analyze(storage_manager* storage) {

    // get the drift velocity
    double efield   = larutil::LArProperties::GetME()->Efield(); // kV/cm
    double temp     = larutil::LArProperties::GetME()->Temperature(); // Kelvin
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(efield,temp); // [cm/us]

    std::cout << "Drift velocity is set to : " << driftVel << " [cm/us]" << std::endl;
  
    // load GeometryHelper utility
    auto geomHelper = ::larutil::GeometryHelper::GetME();

    // grab the vertex data-product
    auto *event_vertex = storage->get_data<larlite::event_vertex>(_vertex_producer);

    // make sure there are vertices in this event!
    // if not exit and send an error message
    if (!event_vertex or (event_vertex->size() == 0) ){
      std::cout << "No vertex found in this file. Quit this event." << std::endl;
      return false;
    }

    std::cout << "In this event we found " << event_vertex->size() << " vertices. Take a look at each one..." << std::endl;

    // loop through the vertices and print out some basic information
    for (size_t n=0; n < event_vertex->size(); n++){

      // grab this specific vertex
      auto const& vertex = event_vertex->at(n);

      // what information is available in the vertex data-product?
      // check it out in this file [path is relative to the top LArLite directory]
      // core/DataFormat/vertex.h

      // grab vertex X,Y,Z information
      // create a c++ array of length 3 where to store the XYZ coordinates
      double * xyz = new double[3];
      // this fuction assigned to the variable 'xyz' the vertex information
      vertex.XYZ(xyz);
      
      std::cout << "Reco'd vertex @ [X, Y, Z] -> [" 
		<< xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]" << std::endl;
      
      // project this point onto the various planes
      for (int plane = 0; plane < 3; plane++){
	
	// geoHelper is a set of utility functions that help with geometric stuff..
	auto const& projection2D = geomHelper->Point_3Dto2D(xyz, plane);
	//projection2D is the projection in cm/cm space.
	std::cout << "\tOn Plane " << plane << " this corresponds to [W,T] -> ["
		  << projection2D.w << ", " << projection2D.t << "] in cm." << std::endl;
	// get the same 2D projection information but in wire/time coordinates
	int wire = (int) ( projection2D.w / geomHelper->WireToCm() );
	int time = (int) ( projection2D.t / geomHelper->TimeToCm() );
	std::cout << "\t in wire/time coordinate this corresponds to [W,T] -> [" 
		  << wire << ", " << time << "]" << std::endl;

      }// for all planes

    }// for all vertices in the event
	

  
    return true;
  }

  bool VertexExample::finalize() {

  
    return true;
  }

}
#endif
