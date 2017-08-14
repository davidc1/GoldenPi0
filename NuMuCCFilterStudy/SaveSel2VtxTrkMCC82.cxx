#ifndef LARLITE_SAVESEL2VTXTRKMCC82_CXX
#define LARLITE_SAVESEL2VTXTRKMCC82_CXX

#include "SaveSel2VtxTrkMCC82.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"

namespace larlite {

  SaveSel2VtxTrkMCC82::SaveSel2VtxTrkMCC82()
  {
    _name="SaveSel2VtxTrkMCC82";
    _fout=0;
    _verbose = true;
  }

  bool SaveSel2VtxTrkMCC82::initialize() {

    return true;
  }
  
  bool SaveSel2VtxTrkMCC82::analyze(storage_manager* storage) {

    // keep track of number of showers and tracks found
    int n_showers = 0;
    int n_tracks  = 0;

    // print event information
    if (_verbose){
      std::cout << "Run   " << storage->run_id() << std::endl;
      std::cout << "Event " << storage->event_id() << std::endl;
    }
    
    // get a handle to the association
    auto ev_ass = storage->get_data<larlite::event_ass>("NuMuCCSelectionII");

    // get track associations
    auto const& ass_trk_v = ev_ass->association(3);
    // and vertex associations
    auto const& ass_vtx_v = ev_ass->association(1);

    //std::cout << "ass trk size : " << ass_trk_v.size() << std::endl;
    //std::cout << "ass vtx size : " << ass_vtx_v.size() << std::endl;

    // get tracks and vertices
    auto ev_trk = storage->get_data<larlite::event_track> ("pandoraNu");
    auto ev_vtx = storage->get_data<larlite::event_vertex>("pandoraNu");

    // get tracks and vertices
    auto ev_trk2 = storage->get_data<larlite::event_track> ("NuMuCCSelectionII");
    auto ev_vtx2 = storage->get_data<larlite::event_vertex>("NuMuCCSelectionII");

    // also, create track & vertex associated with the neutrino interaction
    auto ev_nu_trk = storage->get_data<larlite::event_track>("numuCC_track");
    auto ev_nu_vtx = storage->get_data<larlite::event_vertex>("numuCC_vertex");

    storage->set_id( ev_ass->run(), ev_ass->subrun(), ev_ass->event_id() );

    
    // are there tracks? are there vertices?
    if (!ev_trk or (ev_trk->size() == 0)){
      std::cout << "No track! exit" << std::endl;
      return false;
    }
    if (!ev_vtx or (ev_vtx->size() == 0)){
      std::cout << "No vertex! exit" << std::endl;
      return false;
    }


    // grab associations by product id
    auto const& ass_trk_vtx_v = ev_ass->association(ev_trk2->id(),ev_vtx->id());
    //std::cout << "trk -> vtx ass returns " << ass_trk_vtx_v.size() << " associations" << std::endl;
    for (size_t kk=0; kk < ass_trk_vtx_v.size(); kk++) {
      auto const& ass_trk_vtx = ass_trk_vtx_v.at(kk);
      //std::cout << "\t entry : " << kk << " has " << ass_trk_vtx.size() << " values" << std::endl;
      if (ass_trk_vtx.size() == 1) {
	//std::cout << "\t\t associated w/ entry " << ass_trk_vtx.at(0) << std::endl;
	ev_nu_vtx->emplace_back( ev_vtx->at( ass_trk_vtx[0] ) );
      }
    }
    
    auto const& ass_vtx_trk_v = ev_ass->association(ev_vtx2->id(),ev_trk->id());
    //std::cout << "vtx -> trk ass returns " << ass_vtx_trk_v.size() << " associations" << std::endl;
    for (size_t kk=0; kk < ass_vtx_trk_v.size(); kk++) {
      auto const& ass_vtx_trk = ass_vtx_trk_v.at(kk);
      //std::cout << "\t entry : " << kk << " has " << ass_vtx_trk.size() << " values" << std::endl;
      if (ass_vtx_trk.size() == 1) {
	//std::cout << "\t\t associated w/ entry " << ass_vtx_trk.at(0) << std::endl;
	ev_nu_trk->emplace_back( ev_trk->at( ass_vtx_trk[0] ) );
      }
    }


    //std::cout << "number of tracks   : " << ev_trk->size() << std::endl;
    //std::cout << "number of vertices : " << ev_vtx->size() << std::endl;

    return true;
  }
  
  bool SaveSel2VtxTrkMCC82::finalize() {
    
    return true;
  }

}
#endif
