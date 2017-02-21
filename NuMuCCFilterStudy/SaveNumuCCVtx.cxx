#ifndef LARLITE_SAVENUMUCCVTX_CXX
#define LARLITE_SAVENUMUCCVTX_CXX

#include "SaveNumuCCVtx.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool SaveNumuCCVtx::initialize() {

    //event_file.open("event_file.txt");

    return true;
  }
  
  bool SaveNumuCCVtx::analyze(storage_manager* storage) {

    auto run = storage->run_id();
    auto event = storage->event_id();
    auto subrun = storage->subrun_id();

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");

    //if (ev_vtx->size() != 1) return true;
    //auto const& vtx = ev_vtx->at(0);

    
    std::cout << "RUN : " << run << " EVENT : " << event << std::endl;

    storage->set_id(run, subrun, event);

    std::pair<int,int> evinfo;
    evinfo = std::make_pair(run,event);

    double xyz[3];

    // loop through all entries in map
    for (auto const& element : _vtxmap)
      {
	
	if ( element.first != evinfo) continue;
	
	auto const& xyzloc = element.second;
	
	xyz[0] = xyzloc.at(0);
	xyz[1] = xyzloc.at(1);
	xyz[2] = xyzloc.at(2);

	vertex thisvtx(xyz);

	ev_vtx->emplace_back( thisvtx );

	return true;
	
      }

    /*
    event_file << run << " " << subrun << " " << event << " "
	       << vtx.X() << " " << vtx.Y() << " " << vtx.Z()
	       << "\n";
    */

    std::cout << "could not find vertex for event" << std::endl;

    return true;
  }

  bool SaveNumuCCVtx::finalize() {

    //event_file.close();

    return true;
  }

  void SaveNumuCCVtx::addVertex(const int& run, const int& event, const double& x, const double& y, const double& z) {

    std::vector<double> vtx = {x,y,z};

    std::pair<int,int> evinfo;
    evinfo = std::make_pair(run,event);

    _vtxmap[evinfo] = vtx;
    
    return;

  }

}
#endif
