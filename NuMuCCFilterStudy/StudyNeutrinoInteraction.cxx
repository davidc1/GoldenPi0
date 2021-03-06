#ifndef LARLITE_STUDYNEUTRINOINTERACTION_CXX
#define LARLITE_STUDYNEUTRINOINTERACTION_CXX

#include "StudyNeutrinoInteraction.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"


namespace larlite {

  StudyNeutrinoInteraction::StudyNeutrinoInteraction()
    : _tree(nullptr)
  {
    _name="StudyNeutrinoInteraction";
    _fout=0;
    _verbose = true;
    _n_tracks = 0;
    _n_showers = 0;
    _filter_showers = false;
    _filter_tracks  = false;
    _min_tracks     = false;
    _min_showers    = false;
  }

  bool StudyNeutrinoInteraction::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_muon_length",&_muon_length,"muon_length/D");

    return true;
  }
  
  bool StudyNeutrinoInteraction::analyze(storage_manager* storage) {

    // keep track of number of showers and tracks found
    int n_showers = 0;
    int n_tracks  = 0;
    // keep track of the reconstructed tracks produced in the neutrino
    // interaction. This vector contains all reco tracks with parent
    // the reconstructed neutrino
    std::vector<larlite::track> nu_trk_v;

    // print event information
    if (_verbose){
      std::cout << "Run   " << storage->run_id() << std::endl;
      std::cout << "Event " << storage->event_id() << std::endl;
    }
    
    // get a handle to the association
    auto ev_ass = storage->get_data<larlite::event_ass>("NuMuCCInclusive");
    
    // get the association keys
    auto const& ass_keys = ev_ass->association_keys();
    
    larlite::AssSet_t ass_trk_vtx_v;
    larlite::event_track *ev_trk = nullptr;
    ass_trk_vtx_v = storage->find_one_ass( ass_keys[0].second, ev_trk, ev_ass->name() );

    larlite::AssSet_t ass_vtx_trk_v;
    larlite::event_vertex *ev_vtx = nullptr;
    ass_vtx_trk_v = storage->find_one_ass( ass_keys[0].first, ev_vtx, ev_ass->name() );

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
    bool pfpart = true; // have the PFParticles been found?
    ass_trk_pfpart_v = storage->find_one_ass( ev_trk->id(), ev_pfpart, ev_trk->name() );

    if (!ev_pfpart or (ev_pfpart->size() == 0)){
      std::cout << "No pfpart! exit" << std::endl;
      pfpart = false;
      return false;
    }

    // and now grab tracks associated to the same PFParts
    larlite::AssSet_t ass_pfpart_trk_v;
    larlite::event_track *ev_trk_2 = nullptr;
    bool tracks = true; // have the tracks been found?
    if (pfpart){
      ass_pfpart_trk_v = storage->find_one_ass( ev_pfpart->id(), ev_trk_2, ev_pfpart->name() );
      
      if (!ev_trk_2 or (ev_trk_2->size() == 0)){
	std::cout << "No track associated to PFPart! exit" << std::endl;
	tracks = false;
      }
    }

    if (_verbose)
      std::cout << "Associations between vtx and track : " << ass_vtx_trk_v.size() << std::endl;

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

      // get the muon length
      _muon_length = nutrk.Length();
      _tree->Fill();
      
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

	if (daughter.PdgCode() == 13){
	  n_tracks += 1;
	  if (tracks == false) // don't search for reco tracks if they have not been found.
	    continue;
	  auto trk_idx = ass_pfpart_trk_v[daughter_idx];
	  if (trk_idx.size() != 1){
	    std::cout << "no associated track...skip" << std::endl;
	    continue;
	  }
	  auto trk = ev_trk_2->at(trk_idx[0]);
	  std::cout << "Added daughter to neutrino tracks" << std::endl;
	  // add this track to the list of track-outputs of the neutrino interaction
	  nu_trk_v.push_back( trk );
	}// if tracks
      }
      
      // found the muon -> so exit track loop...
      if (_verbose)
	std::cout << std::endl << std::endl << std::endl;
      
      break;
    }

    return true;
  }
  
  bool StudyNeutrinoInteraction::finalize() {

    _tree->Write();
    
    return true;
  }

}
#endif
