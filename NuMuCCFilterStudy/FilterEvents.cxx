#ifndef LARLITE_FILTEREVENTS_CXX
#define LARLITE_FILTEREVENTS_CXX

#include "FilterEvents.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"
#include "LArUtil/GeometryHelper.h"

namespace larlite {

  FilterEvents::FilterEvents()
  {
    _name="FilterEvents";
    _fout=0;
    _verbose = true;
    _n_tracks = 0;
    _n_showers = 0;
    _filter_showers = false;
    _filter_tracks  = false;
    _min_tracks     = false;
    _min_showers    = false;
    _max_dist = 10;
  }

  bool FilterEvents::initialize() {

    return true;
  }
  
  bool FilterEvents::analyze(storage_manager* storage) {

    auto geoHelper = larutil::GeometryHelper::GetME();

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

    larlite::vertex nuvtx;

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
      nuvtx = ev_vtx->at( ass_vtx_trk_v[i][0] );
      
      // grab the PFParticle associated with this muon
      if (ass_trk_pfpart_v.size() <= i)
	return false;
      if (_verbose)
	std::cout << "PFParts associated with muon : " <<  ass_trk_pfpart_v[i].size() << std::endl;
      auto pfpart_idx = ass_trk_pfpart_v[i][0];
      auto muon = ev_pfpart->at(pfpart_idx);

      // grab parent
      if (muon.Parent() >= ev_pfpart->size()){
	if (_verbose)
	  std::cout << "Muon parent not here..." << std::endl;
	return false;
      }
      
      auto neutrino = ev_pfpart->at( muon.Parent() );

      // print neutrino daughters
      for (auto daughter_idx : neutrino.Daughters() ){
	auto daughter = ev_pfpart->at(daughter_idx);

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
    }
    
    // if there is a single track -> keep the event
    if (nu_trk_v.size() == 1)
      return true;
    
    // if more than 1 track -> require that at least 2 have a "common origin"
    // grab the vertex coordinates on the Y plane
    larutil::Point2D vertex2D;
    double * xyz = new double[3];
    nuvtx.XYZ(xyz);
    vertex2D = geoHelper -> Point_3Dto2D(xyz, 2);
    
    // number of tracks close to the vertex
    int numClose = 0;
    
    // lopp through all trakcs -> at least 2 (muon + other) have to lie within some distance
    // of the vertex
    for (auto const& trk : nu_trk_v){
      
      auto start = geoHelper->Point_3Dto2D( trk.LocationAtPoint(0), 2 );
      auto end   = geoHelper->Point_3Dto2D( trk.LocationAtPoint( trk.NumberTrajectoryPoints() - 1 ), 2 );
      double startDist = sqrt( (start.w - vertex2D.w) * (start.w - vertex2D.w) +
			       (start.t - vertex2D.t) * (start.t - vertex2D.t) );
      double endDist   = sqrt( (end.w - vertex2D.w) * (end.w - vertex2D.w) +
			       (end.t - vertex2D.t) * (end.t - vertex2D.t) );
      
      if ( (startDist < _max_dist) or (endDist < _max_dist) )
	numClose += 1;
      
    }// for all tracks
    
    // if the number of close tracks is less than 2 -> return false
    if (numClose < 2)
      return false;
    

    return true;
  }
  
  bool FilterEvents::finalize() {
    
    return true;
  }

}
#endif
