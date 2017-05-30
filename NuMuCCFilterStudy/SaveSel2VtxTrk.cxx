#ifndef LARLITE_SAVESEL2VTXTRK_CXX
#define LARLITE_SAVESEL2VTXTRK_CXX

#include "SaveSel2VtxTrk.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"

namespace larlite {

  SaveSel2VtxTrk::SaveSel2VtxTrk()
  {
    _name="SaveSel2VtxTrk";
    _fout=0;
    _verbose = true;
  }

  bool SaveSel2VtxTrk::initialize() {

    return true;
  }
  
  bool SaveSel2VtxTrk::analyze(storage_manager* storage) {

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
    
    // get the association keys
    auto const& ass_keys = ev_ass->association_keys();
    
    larlite::AssSet_t ass_trk_vtx_v;
    larlite::event_track *ev_trk = nullptr;
    ass_trk_vtx_v = storage->find_one_ass( ass_keys[0].second, ev_trk, ev_ass->name() );

    larlite::AssSet_t ass_vtx_trk_v;
    larlite::event_vertex *ev_vtx = nullptr;
    ass_vtx_trk_v = storage->find_one_ass( ass_keys[0].first, ev_vtx, ev_ass->name() );

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


    if (_verbose){
      std::cout << "Associations between vtx and trk : " << ass_vtx_trk_v.size() << std::endl;
      std::cout << "Associations between trk and vtx : " << ass_trk_vtx_v.size() << std::endl;
    }

    // find the track and vertex associated to the neutrino
    for (size_t i=0; i < ass_vtx_trk_v.size(); i++){
      
      if (ass_vtx_trk_v[i].size() == 0){
	//std::cout << "vtx->trk association is empty..." << std::endl;
	continue;
      }
      if (_verbose){
	std::cout << "trk " << i << " associated to vtx " << ass_vtx_trk_v[i][0] << std::endl;
	std::cout << ev_trk->size() << " tracks present.." << std::endl;
	std::cout << ev_vtx->size() << " vertices present.." << std::endl;
      }
      auto const& nutrk = ev_trk->at(i);
      auto const& nuvtx = ev_vtx->at( ass_vtx_trk_v[i][0] );

      ev_nu_trk->emplace_back( nutrk );
      ev_nu_vtx->emplace_back( nuvtx );

    }

    return true;
  }
  
  bool SaveSel2VtxTrk::finalize() {
    
    return true;
  }

}
#endif
