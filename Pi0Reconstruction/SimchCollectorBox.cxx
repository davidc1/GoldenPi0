#ifndef LARLITE_SIMCHCOLLECTORBOX_CXX
#define LARLITE_SIMCHCOLLECTORBOX_CXX

#include "SimchCollectorBox.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include "DataFormat/simch.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool SimchCollectorBox::initialize() {

    _ide_tree = new TTree("_ide_tree","IDE TREE");
    _ide_tree->Branch("_x",&_x,"x/D");
    _ide_tree->Branch("_y",&_y,"y/D");
    _ide_tree->Branch("_z",&_z,"z/D");
    _ide_tree->Branch("_e",&_e,"e/D");

    _hit_tree = new TTree("_hit_tree","HIT TREE");
    _hit_tree->Branch("_w",&_w,"w/D");
    _hit_tree->Branch("_t",&_t,"t/D");
    _hit_tree->Branch("_adc",&_adc,"adc/D");

    return true;
  }
  
  bool SimchCollectorBox::analyze(storage_manager* storage) {

    std::cout << "Event : " << storage->event_id() << std::endl;
    
    auto *ev_mctruth  = storage->get_data<event_mctruth>("generator");
    auto *ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto *ev_vertex   = storage->get_data<event_vertex>("numuCC_vertex");
    auto *ev_shower   = storage->get_data<event_shower>("showerreco");
    auto *ev_simch    = storage->get_data<event_simch>("largeant");
    auto *ev_hit      = storage->get_data<event_hit>("gaushit");
    
    if ( (!ev_mctruth) or (ev_mctruth->size() == 0) )
      return false; 
    
    // get all MCParticles
    auto mctruth = ev_mctruth->at(0);
    
    auto part_v = mctruth.GetParticles();
    
    int n_pi0 = 0;
    
    size_t pi0_idx = 0;
    int    pi0_trkid = 0;

    double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
    
    for (size_t i=0; i < part_v.size(); i++){
      auto const& part = part_v[i];
      //std::cout << "found " << part.PdgCode() << std::endl;
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	//std::cout << "found " << part.PdgCode() << "\tw/ energy : " << part.Trajectory().at(0).E() <<  std::endl;
	n_pi0    += 1;
	pi0_idx   = i;
	pi0_trkid = part.TrackId();
	_mc_vtx_x = part.Trajectory().at( 0 ).X();
	_mc_vtx_y = part.Trajectory().at( 0 ).Y();
	_mc_vtx_z = part.Trajectory().at( 0 ).Z();
      }
    }
    
    if (n_pi0 != 1)
      return false;
    
    // is the vertex in the TPC?
    auto const& vtx = part_v.at(pi0_idx).Trajectory().at(0);
    if ( (vtx.X() < 0) or (vtx.X() > 256) or (vtx.Y() < -116) or (vtx.Y() > 116) or (vtx.Z() < 0) or (vtx.Z() > 1036) ){
      std::cout << "Vertex @ x : " << vtx.X() << "\t y : " << vtx.Y() << "\t z : " << vtx.Z() << std::endl;
      return false;
    }
    
    size_t pi0_ID_1 = 0;
    size_t pi0_ID_2 = 0;
    size_t idx_1 = 0;
    size_t idx_2 = 0;
    size_t n_found = 0;
    std::cout << "there are " << ev_mcshower->size() << " MCShowers" << std::endl;
    for (size_t i=0; i < ev_mcshower->size(); i++){
      auto const& mcs = ev_mcshower->at(i);
      //std::cout << "MCShower energy : " << mcs.DetProfile().E() << std::endl;
      // distance from vertex
      double x = mcs.Start().X();
      double y = mcs.Start().Y();
      double z = mcs.Start().Z();
      double d = sqrt( ( (_mc_vtx_x - x) * (_mc_vtx_x - x) ) +
		       ( (_mc_vtx_y - y) * (_mc_vtx_y - y) ) +
		       ( (_mc_vtx_z - z) * (_mc_vtx_z - z) ) );
      if ( d < 0.01 ){
	if (n_found == 0){
	  pi0_ID_1 = mcs.MotherTrackID();
	  idx_1 = i;
	  n_found += 1;
	}
	else if (n_found == 1){
	  pi0_ID_2 = mcs.MotherTrackID();
	  idx_2 = i;
	  n_found += 1;
	}
	else
	  n_found += 1;
      }// if mother is a Pi0
    }// for all MCShowers

    std::cout << "Found " << n_found << " Gammas" << std::endl;
    std::cout << "Pi0 ID  : " << pi0_trkid << std::endl;

    if ( n_found != 2 ) return false;
    
    auto const& shr1 = ev_mcshower->at(idx_1);
    auto const& shr2 = ev_mcshower->at(idx_2);

    std::cout << "shr1 E : " << shr1.Start().E() << std::endl;
    std::cout << "shr2 E : " << shr2.Start().E() << std::endl;

    //std::cout << "shr1 : TrackID: " << shr1.TrackID() << std::endl;
    //std::cout << "shr2 : TrackID: " << shr2.TrackID() << std::endl;

    std::vector<unsigned int> shr1_trackid_v_all;
    shr1_trackid_v_all.push_back(shr1.TrackID());
    std::vector<unsigned int> shr2_trackid_v_all;
    shr2_trackid_v_all.push_back(shr2.TrackID());
    
    auto const& shr1_trackid_v = shr1.DaughterTrackID();
    auto const& shr2_trackid_v = shr2.DaughterTrackID();

    for (auto const& id : shr1_trackid_v)
      shr1_trackid_v_all.push_back(id);
    for (auto const& id : shr2_trackid_v)
      shr2_trackid_v_all.push_back(id);

    //std::cout << "number of particles for shr1: " << shr1_trackid_v_all.size() << std::endl;
    //std::cout << "number of particles for shr2: " << shr2_trackid_v_all.size() << std::endl;

    double ide_e_sum_1 = 0;
    double ide_e_sum_2 = 0;
    
    double ide_e_shr1, ide_e_shr2;

    for (size_t l=0; l < ev_simch->size(); l++){
      
      auto const& simch = ev_simch->at(l);
      
      if (simch.Channel() < 4800)
	continue;
      
      auto const& all_ide = simch.TrackIDsAndEnergies(0,19600);
      
      for (size_t j=0; j < all_ide.size(); j++){
	
	auto const& ide = all_ide[j];

	_x = ide.x;
	_y = ide.y;
	_z = ide.z;
	_e = ide.energy;

	if ( (_z > 534) && (_z < 634) && (_x > 0) && (_x < 100) ){
	  
	  _ide_tree->Fill();
	  
	  if ( (_x > 5) && (_x < 45) && (_z > 540) && (_z < 575) ){
	    ide_e_sum_1 += _e;
	  }
	  if ( (_x > 55) && (_x < 90) && (_z > 540) && (_z < 575) ){
	    ide_e_sum_2 += _e;
	  }
	}

	/*

	if (std::find(shr1_trackid_v_all.begin(),
		      shr1_trackid_v_all.end(),
		      ide.trackID) != shr1_trackid_v_all.end())
	  ide_e_shr1 += ide.energy;
	
	if (std::find(shr2_trackid_v_all.begin(),
		      shr2_trackid_v_all.end(),
		      ide.trackID) != shr2_trackid_v_all.end())
	  ide_e_shr2 += ide.energy;
	*/
	
      }// for all IDes
      
    }// for all simchannels

    std::cout << "IDE sum 1 : " << ide_e_sum_1 << std::endl;
    std::cout << "IDE sum 2 : " << ide_e_sum_2 << std::endl;

    double adc_e_sum = 0;

    for (size_t h=0; h < ev_hit->size(); h++){

      auto const& hit = ev_hit->at(h);

      if (hit.WireID().Plane != 2) continue;
      
      _w = hit.WireID().Wire * 0.3;
      _t = hit.PeakTime() * 0.05;
      _adc = hit.Integral();
      _hit_tree->Fill();

      adc_e_sum += _adc;

    }

    std::cout << "ADC shr 1 : " << adc_e_sum << std::endl;
      
    
    return true;
  }

  bool SimchCollectorBox::finalize() {

    _ide_tree->Write();
    _hit_tree->Write();

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
