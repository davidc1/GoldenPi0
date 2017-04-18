#ifndef LARLITE_MAKEMCVERTEX_CXX
#define LARLITE_MAKEMCVERTEX_CXX

#include "MakeMCVertex.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"

#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool MakeMCVertex::initialize() {

    _SCE = new larutil::SpaceChargeMicroBooNE();

    _time2cm = larutil::GeometryHelper::GetME()->TimeToCm();;

    return true;
  }
  
  bool MakeMCVertex::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator");
    if(!ev_mctruth || !ev_mctruth->size() ) return false;
    
    auto nu = ev_mctruth->at(0).GetNeutrino();
    auto parts = ev_mctruth->at(0).GetParticles();
    
    auto ev_vtx = storage->get_data<event_vertex>("mcvertex");
    
    storage->set_id(storage->run_id(), storage->subrun_id(), storage->event_id());
    
    ev_vtx->reserve(1);
    double xyz[3] = {0.};

    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();
    auto tvtx = traj.at(traj.size()-1).T(); // ns
    // time in tick :
    auto vtxtick = (tvtx / 1000.) * 2.;
    // time in cm :
    auto vtxtimecm = vtxtick * _time2cm; 

    //std::cout << "neutrino time @ " << traj.at(traj.size()-1).T() << std::endl;
    //std::cout << "correction in cm " << vtxtimecm << std::endl;
    //std::cout << std::endl;

    // get spacecharge correction
    auto sce_corr = _SCE->GetPosOffsets(xvtx,yvtx,zvtx);
    
    xyz[0] = xvtx + _offset - sce_corr.at(0) + vtxtimecm;
    xyz[1] = yvtx + sce_corr.at(1);
    xyz[2] = zvtx + sce_corr.at(2);

    if ( (xyz[0] < 0) || (xyz[0] > 256) || (xyz[1] < -116) || (xyz[1] > 116) || (xyz[2] < 0) || (xyz[2] > 1036) )
      return false;

    //std::cout << "created vertex @ [ " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "] " << std::endl;
    
    vertex new_vtx(xyz) ;
    ev_vtx->push_back(new_vtx);

    return true;
  }

  bool MakeMCVertex::finalize() {

    return true;
  }

}
#endif
