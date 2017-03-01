#ifndef LARLITE_BNBPI0PROPERTIES_CXX
#define LARLITE_BNBPI0PROPERTIES_CXX

#include "BNBPi0Properties.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"

namespace larlite {

  BNBPi0Properties::BNBPi0Properties()
    : _tree(nullptr)
  {
    _name="BNBPi0Properties";
    _fout=0;
  }

  bool BNBPi0Properties::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","shower study tree");

    _tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
    
    _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
    _tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");

    _tree->Branch("_event",&_event,"event/I");

    // vertex info
    _tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
    _tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
    _tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");

    // 1st shower info
    _tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
    _tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
    _tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
    _tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
    _tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
    _tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
    _tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");

    // 2nd shower info
    _tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
    _tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
    _tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
    _tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
    _tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
    _tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
    _tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");

    return true;
  }
  
  bool BNBPi0Properties::analyze(storage_manager* storage) {

    auto *ev_mctruth  = storage->get_data<event_mctruth>("generator");
    auto *ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto *ev_shower   = storage->get_data<event_shower>("showerreco");

    if ( (!ev_mctruth) or (ev_mctruth->size() == 0) )
      return false; 
    
    // get all MCParticles
    auto mctruth = ev_mctruth->at(0);
    
    _nu_e = mctruth.GetNeutrino().Nu().Trajectory().at(0).E();
    
    auto part_v = mctruth.GetParticles();
    
    int n_pi0 = 0;
    
    size_t pi0_idx = 0;
    int    pi0_trkid = 0;
    
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
	_pi0_e    = part.Trajectory().at( 0 ).E();
      }
    }
    
    size_t pi0_ID_1 = 0;
    size_t pi0_ID_2 = 0;
    size_t idx_1 = 0;
    size_t idx_2 = 0;
    size_t n_found = 0;
    std::cout << "there are " << ev_mcshower->size() << " MCShowers" << std::endl;
    for (size_t i=0; i < ev_mcshower->size(); i++){
      auto const& mcs = ev_mcshower->at(i);
      //std::cout << "MCShower energy : " << mcs.Start().E() << std::endl;
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
    
    if ( (n_found == 2) ){//and (pi0_ID_1 == pi0_trkid) and (pi0_ID_2 == pi0_trkid ) ){
      
      auto const& shr1 = ev_mcshower->at(idx_1);
      auto const& shr2 = ev_mcshower->at(idx_2);
      
      //_pi0_e  = shr1.MotherEnd().E();
      
      _mc_shr1_e  = shr1.Start().E();
      std::cout << "Start E1 : " <<  shr1.Start().E() << std::endl;
      _mc_shr1_x  = shr1.Start().X();
      _mc_shr1_y  = shr1.Start().Y();
      _mc_shr1_z  = shr1.Start().Z();
      double mom1 = sqrt ( ( shr1.Start().Px() * shr1.Start().Px() ) +
			   ( shr1.Start().Py() * shr1.Start().Py() ) +
			   ( shr1.Start().Pz() * shr1.Start().Pz() ) );
      
      _mc_shr1_px = shr1.Start().Px() / mom1;
      _mc_shr1_py = shr1.Start().Py() / mom1;
      _mc_shr1_pz = shr1.Start().Pz() / mom1;
      
      _mc_shr2_e  = shr2.Start().E();
      std::cout << "Start E2 : " <<  shr2.Start().E() << std::endl;
      _mc_shr2_x  = shr2.Start().X();
      _mc_shr2_y  = shr2.Start().Y();
      _mc_shr2_z  = shr2.Start().Z();
      double mom2 = sqrt ( ( shr2.Start().Px() * shr2.Start().Px() ) +
			   ( shr2.Start().Py() * shr2.Start().Py() ) +
			   ( shr2.Start().Pz() * shr2.Start().Pz() ) );
      
      _mc_shr2_px = shr2.Start().Px() / mom2;
      _mc_shr2_py = shr2.Start().Py() / mom2;
      _mc_shr2_pz = shr2.Start().Pz() / mom2;
      
      } // if we found 2 mc showers that come from the beam pi0

    _tree->Fill();
  
    return true;
  }

  bool BNBPi0Properties::finalize() {

    _tree->Write();
  
    return true;
  }

}
#endif
