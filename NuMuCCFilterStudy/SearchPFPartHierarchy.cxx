#ifndef LARLITE_SEARCHPFPARTHIERARCHY_CXX
#define LARLITE_SEARCHPFPARTHIERARCHY_CXX

#include "SearchPFPartHierarchy.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"

namespace larlite {

  bool SearchPFPartHierarchy::initialize() {

    return true;
  }
  
  bool SearchPFPartHierarchy::analyze(storage_manager* storage) {

    // print event information
    std::cout << "Run   " << storage->run_id() << std::endl;
    std::cout << "Event " << storage->event_id() << std::endl;
    
    // get a handle to the association
    auto ev_ass = storage->get_data<larlite::event_ass>("NuMuCCInclusive");
    
    // get the association keys
    auto const& ass_keys = ev_ass->association_keys();
    
    larlite::AssSet_t ass_trk_vtx_v;
    larlite::event_track *ev_trk = nullptr;
    ass_trk_vtx_v = storage->find_one_ass( ass_keys[0].second, ev_trk, ev_ass->name() );

    larlite::AssSet_t ass_vtx_trk_v;
    larlite::event_track *ev_vtx = nullptr;
    ass_vtx_trk_v = storage->find_one_ass( ass_keys[0].second, ev_vtx, ev_ass->name() );

    // are there tracks? are there vertices?
    if (!ev_trk or (ev_trk->size() == 0)){
      std::cout << "No track! exit" << std::endl;
      return false;
    }
    if (!ev_vtx or (ev_vtx->size() == 0)){
      std::cout << "No vertex! exit" << std::endl;
      return false;
    }

    // grab PFParticles associated with these tracks
    larlite::AssSet_t ass_trk_pfpart_v;
    larlite::event_pfpart *ev_pfpart = nullptr;
    ass_trk_pfpart_v = storage->find_one_ass( ev_trk->id(), ev_pfpart, ev_trk->name() );

    if (!ev_pfpart or (ev_pfpart->size() == 0)){
      std::cout << "No pfpart! exit" << std::endl;
      return false;
    }

    std::cout << "Associations between vtx and track : " << ass_vtx_trk_v.size() << std::endl;

    // find the track and vertex associated to the neutrino
    for (size_t i=0; i < ass_vtx_trk_v.size(); i++){
      
      if (ass_vtx_trk_v[i].size() == 0){
	std::cout << "vtx->trk association is empty..." << std::endl;
	continue;
      }
      
      std::cout << "trk " << i << " associated to vtx " << ass_vtx_trk_v[i][0] << std::endl;
      auto const& nutrk = ev_trk->at(i);
      auto const& nuvtx = ev_vtx->at( ass_vtx_trk_v[i][0] );
      
      // grab the PFParticle associated with this muon
      std::cout << "PFParts associated with muon : " <<  ass_trk_pfpart_v[i].size() << std::endl;
      auto pfpart_idx = ass_trk_pfpart_v[i][0];
      auto muon = ev_pfpart->at(pfpart_idx);
      
      std::cout << "Muon PFPart info :" << std::endl
		<< "\tPDG code   : " << muon.PdgCode() << std::endl
		<< "\tDaughters? : " << muon.NumDaughters() << std::endl
		<< "\t Parent?   : " << muon.Parent() << std::endl;
      
      // grab parent
      if (muon.Parent() >= ev_pfpart->size()){
	std::cout << "Muon parent not here..." << std::endl;
	return false;
      }
      
      auto neutrino = ev_pfpart->at( muon.Parent() );
      
      std::cout << "Neutrino PFPart info :" << std::endl
		<< "\tPDG code   : " << neutrino.PdgCode() << std::endl
		<< "\tDaughters? : " << neutrino.NumDaughters() << std::endl
		<< "\t Parent?   : " << neutrino.Parent() << std::endl;
      
      // print neutrino daughters
      for (auto daughter_idx : neutrino.Daughters() ){
	auto daughter = ev_pfpart->at(daughter_idx);
	std::cout << "daughter PFPart info :" << std::endl
		  << "\tPDG code   : " << daughter.PdgCode() << std::endl
		  << "\tDaughters? : " << daughter.NumDaughters() << std::endl
		  << "\t Parent?   : " << daughter.Parent() << std::endl;
      }
      
      // found the muon -> so exit track loop...

      std::cout << std::endl << std::endl << std::endl;
      
      break;
    }
    
    return true;
  }
  
  bool SearchPFPartHierarchy::finalize() {
    
    return true;
  }

}
#endif
