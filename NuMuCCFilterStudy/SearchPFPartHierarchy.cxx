#ifndef LARLITE_SEARCHPFPARTHIERARCHY_CXX
#define LARLITE_SEARCHPFPARTHIERARCHY_CXX

#include "SearchPFPartHierarchy.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"

namespace larlite {

  SearchPFPartHierarchy::SearchPFPartHierarchy()
  {
    _name="SearchPFPartHierarchy";
    _fout=0;
    _verbose = true;
    _n_tracks = 0;
    _n_showers = 0;
    _n_daughters = 0;
    _filter_showers = false;
    _filter_tracks  = false;
    _min_tracks     = false;
    _min_showers    = false;
    _min_daughters  = false;
  }

  bool SearchPFPartHierarchy::initialize() {

    return true;
  }
  
  bool SearchPFPartHierarchy::analyze(storage_manager* storage) {

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

    // are there tracks? are there vertices?
    if (!ev_trk or (ev_trk->size() == 0)){
      std::cout << "No track! exit" << std::endl;
      return false;
    }
    if (!ev_vtx or (ev_vtx->size() == 0)){
      std::cout << "No vertex! exit" << std::endl;
      return false;
    }

    /*
    // grab PFParticles associated with these tracks
    larlite::AssSet_t ass_trk_pfpart_v;
    larlite::event_pfpart *ev_pfpart = nullptr;
    ass_trk_pfpart_v = storage->find_one_ass( ev_trk->id(), ev_pfpart, ev_trk->name() );

    if (!ev_pfpart or (ev_pfpart->size() == 0)){
      std::cout << "No pfpart! exit" << std::endl;
      return false;
    }
    */

    storage->set_id( ev_ass->run(), ev_ass->subrun(), ev_ass->event_id() );

    if (_verbose){
      std::cout << "Associations between vtx and trk : " << ass_vtx_trk_v.size() << std::endl;
      std::cout << "Associations between trk and vtx : " << ass_trk_vtx_v.size() << std::endl;
    }

    // find the track and vertex associated to the neutrino
    for (size_t i=0; i < ass_vtx_trk_v.size(); i++){
      
      if (ass_vtx_trk_v[i].size() == 0){
	std::cout << "vtx->trk association is empty..." << std::endl;
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

    // find the track and vertex associated to the neutrino
    for (size_t i=0; i < ass_trk_vtx_v.size(); i++){
      
      if (ass_trk_vtx_v[i].size() == 0){
	std::cout << "trk->vtx association is empty..." << std::endl;
	continue;
      }
      if (_verbose){
	std::cout << "vtx " << i << " associated to trk " << ass_trk_vtx_v[i][0] << std::endl;
	std::cout << ev_trk->size() << " tracks present.." << std::endl;
	std::cout << ev_vtx->size() << " vertices present.." << std::endl;
      }
      /*
      auto const& nutrk = ev_trk->at(i);
      auto const& nuvtx = ev_vtx->at( ass_vtx_trk_v[i][0] );

      ev_nu_trk->emplace_back( nutrk );
      ev_nu_vtx->emplace_back( nuvtx );
      */

    }
      /*
      
      // grab the PFParticle associated with this muon
      if (ass_trk_pfpart_v.size() <= i)
	return false;
      if (_verbose)
	std::cout << "PFParts associated with muon : " <<  ass_trk_pfpart_v[i].size() << std::endl;
      auto pfpart_idx = ass_trk_pfpart_v[i][0];
      auto muon = ev_pfpart->at(pfpart_idx);

      if (_verbose)
	std::cout << "Muon PFPart info :" << std::endl
		  << "\tPDG code   : " << muon.PdgCode() << std::endl
		  << "\tDaughters? : " << muon.NumDaughters() << std::endl
		  << "\t Parent?   : " << muon.Parent() << std::endl;
      
      // grab parent
      if (muon.Parent() >= ev_pfpart->size()){
	if (_verbose)
	  std::cout << "Muon parent not here..." << std::endl;
	return false;
      }
      
      auto neutrino = ev_pfpart->at( muon.Parent() );

      if (_verbose)
	std::cout << "Neutrino PFPart info :" << std::endl
		  << "\tPDG code   : " << neutrino.PdgCode() << std::endl
		  << "\tDaughters? : " << neutrino.NumDaughters() << std::endl
		  << "\t Parent?   : " << neutrino.Parent() << std::endl;
      
      // print neutrino daughters
      for (auto daughter_idx : neutrino.Daughters() ){
	auto daughter = ev_pfpart->at(daughter_idx);
	if (_verbose)
	  std::cout << "daughter PFPart info :" << std::endl
		    << "\tPDG code   : " << daughter.PdgCode() << std::endl
		    << "\tDaughters? : " << daughter.NumDaughters() << std::endl
		    << "\t Parent?   : " << daughter.Parent() << std::endl;
	if (daughter.PdgCode() == 11)
	  n_showers += 1;
	if (daughter.PdgCode() == 13)
	  n_tracks += 1;
      }
      
      // found the muon -> so exit track loop...

      if (_verbose)
	std::cout << std::endl << std::endl << std::endl;
      
      break;
    }

    if (_min_daughters){
      if ( (n_showers + n_tracks) >= _n_daughters)
	return true;
      else
	return false;
    }// if we want to set a minimum number of daughters

    if (!_min_showers){
      if (_filter_showers and (n_showers != _n_showers) )
	return false;
    }
    else{
      if (_filter_showers and (n_showers < _n_showers) )
	return false;
    }
    if (!_min_tracks){
      if (_filter_tracks and (n_tracks != _n_tracks) )
	return false;
    }
    else{
      if (_filter_tracks and (n_tracks < _n_tracks) )
	return false;
    }
    
      */
    
    return true;
  }
  
  bool SearchPFPartHierarchy::finalize() {
    
    return true;
  }

}
#endif
