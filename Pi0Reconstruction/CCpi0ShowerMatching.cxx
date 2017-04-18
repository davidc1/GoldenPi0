#ifndef LARLITE_CCPI0SHOWERMATCHING_CXX
#define LARLITE_CCPI0SHOWERMATCHING_CXX

#include "CCpi0ShowerMatching.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool CCpi0ShowerMatching::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","shower study tree");

    _tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
    
    _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
    _tree->Branch("_pi_e",&_pi0_e,"pi0_e/D");

    _tree->Branch("_event",&_event,"event/I");

    // vertex info
    _tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
    _tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
    _tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
    _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
    _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
    _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");

    // 1st shower info
    _tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
    _tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
    _tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
    _tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
    _tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
    _tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
    _tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
    _tree->Branch("_rc_shr1_x",&_rc_shr1_x,"rc_shr1_x/D");
    _tree->Branch("_rc_shr1_y",&_rc_shr1_y,"rc_shr1_y/D");
    _tree->Branch("_rc_shr1_z",&_rc_shr1_z,"rc_shr1_z/D");
    _tree->Branch("_rc_shr1_e",&_rc_shr1_e,"rc_shr1_e/D");
    _tree->Branch("_rc_shr1_px",&_rc_shr1_px,"rc_shr1_px/D");
    _tree->Branch("_rc_shr1_py",&_rc_shr1_py,"rc_shr1_py/D");
    _tree->Branch("_rc_shr1_pz",&_rc_shr1_pz,"rc_shr1_pz/D");

    // 2nd shower info
    _tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
    _tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
    _tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
    _tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
    _tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
    _tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
    _tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");
    _tree->Branch("_rc_shr2_x",&_rc_shr2_x,"rc_shr2_x/D");
    _tree->Branch("_rc_shr2_y",&_rc_shr2_y,"rc_shr2_y/D");
    _tree->Branch("_rc_shr2_z",&_rc_shr2_z,"rc_shr2_z/D");
    _tree->Branch("_rc_shr2_e",&_rc_shr2_e,"rc_shr2_e/D");
    _tree->Branch("_rc_shr2_px",&_rc_shr2_px,"rc_shr2_px/D");
    _tree->Branch("_rc_shr2_py",&_rc_shr2_py,"rc_shr2_py/D");
    _tree->Branch("_rc_shr2_pz",&_rc_shr2_pz,"rc_shr2_pz/D");

    // shower correlations
    _tree->Branch("_dot1",&_dot1,"dot1/D");
    _tree->Branch("_dot2",&_dot2,"dot2/D");

    _tree->Branch("_rc_oangle",&_rc_oangle,"rc_oangle/D");
    _tree->Branch("_mc_oangle",&_mc_oangle,"mc_oangle/D");
    _tree->Branch("_rc_mass"  ,&_rc_mass  ,"_rc_mass/D" );


    return true;
  }
  
  bool CCpi0ShowerMatching::analyze(storage_manager* storage) {

    Reset();

    _event = storage->event_id();

    // start with mc info

    auto *ev_mctruth  = storage->get_data<event_mctruth>("generator");
    auto *ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    //auto *ev_vertex   = storage->get_data<event_vertex>("numuCC_vertex");
    auto *ev_vertex   = storage->get_data<event_vertex>("mcvertex");
    auto *ev_shower   = storage->get_data<event_shower>("showerreco");

    // skip if less than 2 showers reconstructed
    if ( !ev_shower ) {
      print(larlite::msg::kERROR,__FUNCTION__," missing reco information.");
      return false;
    }
      
    if ( (!ev_mctruth) or (ev_mctruth->size() == 0) ) {
      print(larlite::msg::kERROR,__FUNCTION__," missing truth information.");
      return false;
    }
    
    // get all MCParticles
    auto mctruth = ev_mctruth->at(0);
    
    _nu_e = mctruth.GetNeutrino().Nu().Trajectory().at(0).E();
    auto part_v = mctruth.GetParticles();
    
    int    n_pi0 = 0;
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

    if ( (n_found != 2) ) return false;
      
    auto const& mcshr1 = ev_mcshower->at(idx_1);
    auto const& mcshr2 = ev_mcshower->at(idx_2);
    
    _pi0_e  = mcshr1.MotherEnd().E();
    
    _mc_shr1_e  = mcshr1.Start().E();
    _mc_shr1_x  = mcshr1.Start().X();
    _mc_shr1_y  = mcshr1.Start().Y();
    _mc_shr1_z  = mcshr1.Start().Z();
    double mom1 = sqrt ( ( mcshr1.Start().Px() * mcshr1.Start().Px() ) +
			 ( mcshr1.Start().Py() * mcshr1.Start().Py() ) +
			 ( mcshr1.Start().Pz() * mcshr1.Start().Pz() ) );
    
    _mc_shr1_px = mcshr1.Start().Px() / mom1;
    _mc_shr1_py = mcshr1.Start().Py() / mom1;
    _mc_shr1_pz = mcshr1.Start().Pz() / mom1;
    
    _mc_shr2_e  = mcshr2.Start().E();
    _mc_shr2_x  = mcshr2.Start().X();
    _mc_shr2_y  = mcshr2.Start().Y();
    _mc_shr2_z  = mcshr2.Start().Z();
    double mom2 = sqrt ( ( mcshr2.Start().Px() * mcshr2.Start().Px() ) +
			 ( mcshr2.Start().Py() * mcshr2.Start().Py() ) +
			 ( mcshr2.Start().Pz() * mcshr2.Start().Pz() ) );
    
    _mc_shr2_px = mcshr2.Start().Px() / mom2;
    _mc_shr2_py = mcshr2.Start().Py() / mom2;
    _mc_shr2_pz = mcshr2.Start().Pz() / mom2;

    _mc_oangle  = mcshr1.Start().Momentum().Vect().Dot( mcshr1.Start().Momentum().Vect() );
    _mc_oangle /= mcshr1.Start().Momentum().Vect().Mag();
    _mc_oangle /= mcshr2.Start().Momentum().Vect().Mag();

    std::vector<larlite::mcshower> pi0_mcshower_v = {mcshr1,mcshr2};

    _n_reco_showers = ev_shower->size();
    
    if (ev_shower->size() < 2) {
      print(larlite::msg::kERROR,__FUNCTION__," less than 2 reco showers. Skip.");
      _tree->Fill();
      return false;
    }
    
    // moving on to reconstruction
    if (ev_vertex->size() == 1){
      auto const& vtx = ev_vertex->at(0);
      _rc_vtx_x = vtx.X();
      _rc_vtx_y = vtx.Y();
      _rc_vtx_z = vtx.Z();
    }

    std::vector<larlite::shower> reco_shower_v;
    for (size_t i=0; i < ev_shower->size(); i++) 
      reco_shower_v.push_back( ev_shower->at(i) );
    
    // MC <-> RC matching
    auto MCRCmatch = Match(pi0_mcshower_v, reco_shower_v);
    
    auto const& rcshr1 = ev_shower->at( MCRCmatch.first );
    
    _rc_shr1_e = rcshr1.Energy();
    _rc_shr1_x = rcshr1.ShowerStart().X();
    _rc_shr1_y = rcshr1.ShowerStart().Y();
    _rc_shr1_z = rcshr1.ShowerStart().Z();
    mom1 = sqrt( ( rcshr1.Direction().X() * rcshr1.Direction().X() ) +
		 ( rcshr1.Direction().Y() * rcshr1.Direction().Y() ) +
		 ( rcshr1.Direction().Z() * rcshr1.Direction().Z() ) );
    _rc_shr1_px = rcshr1.Direction().X() / mom1;
    _rc_shr1_py = rcshr1.Direction().Y() / mom1;
    _rc_shr1_pz = rcshr1.Direction().Z() / mom1;


    _dot1  = rcshr1.Direction().Dot( mcshr1.Start().Momentum().Vect() );
    _dot1 /= mcshr1.Start().Momentum().Vect().Mag();
    _dot1 /= rcshr1.Direction().Mag();
    
    
    auto const& rcshr2 = ev_shower->at( MCRCmatch.second );
    
    _rc_shr2_e = rcshr2.Energy();
    _rc_shr2_x = rcshr2.ShowerStart().X();
    _rc_shr2_y = rcshr2.ShowerStart().Y();
    _rc_shr2_z = rcshr2.ShowerStart().Z();
    mom2 = sqrt( ( rcshr2.Direction().X() * rcshr2.Direction().X() ) +
		 ( rcshr2.Direction().Y() * rcshr2.Direction().Y() ) +
		 ( rcshr2.Direction().Z() * rcshr2.Direction().Z() ) );
    _rc_shr2_px = rcshr2.Direction().X() / mom2;
    _rc_shr2_py = rcshr2.Direction().Y() / mom2;
    _rc_shr2_pz = rcshr2.Direction().Z() / mom2;

    _dot2  = rcshr2.Direction().Dot( mcshr2.Start().Momentum().Vect() );
    _dot2 /= mcshr2.Start().Momentum().Vect().Mag();
    _dot2 /= rcshr2.Direction().Mag();

    _rc_oangle  = rcshr1.Direction().Dot( rcshr2.Direction() );
    _rc_oangle /= rcshr1.Direction().Mag();
    _rc_oangle /= rcshr2.Direction().Mag();

    _rc_mass = sqrt( 2 * _rc_shr1_e * _rc_shr2_e * ( 1 - _rc_oangle ) );
    
    _tree->Fill();

    return true;
  }
  
  bool CCpi0ShowerMatching::finalize() {

    if (_fout) _fout->cd();
    _tree->Write();
    
    return true;
  }

  void CCpi0ShowerMatching::Reset() {

    _n_reco_showers = 0;
    _nu_e = 0;
    _pi0_e = 0;

    _mc_vtx_x= _mc_vtx_y= _mc_vtx_z= 0;
    _rc_vtx_x= _rc_vtx_y= _rc_vtx_z= 0;
    
    _mc_shr1_x=  _mc_shr1_y=  _mc_shr1_z= 0;
    _mc_shr1_px= _mc_shr1_py= _mc_shr1_pz= 0;
    _mc_shr1_e= 0;
    _rc_shr1_x=  _rc_shr1_y=  _rc_shr1_z= 0;
    _rc_shr1_px= _rc_shr1_py= _rc_shr1_pz= 0;
    _rc_shr1_e= 0;
    
    _mc_shr2_x=  _mc_shr2_y=  _mc_shr2_z= 0;
    _mc_shr2_px= _mc_shr2_py= _mc_shr2_pz= 0;
    _mc_shr2_e= 0;
    _rc_shr2_x=  _rc_shr2_y=  _rc_shr2_z= 0;
    _rc_shr2_px= _rc_shr2_py= _rc_shr2_pz= 0;
    _rc_shr2_e= 0;
    
    return;
  }

  std::pair<int,int> CCpi0ShowerMatching::Match(const std::vector<larlite::mcshower>& mcs_v,
						const std::vector<larlite::shower>&   shr_v) {

    // STEP 1
    // sort reco showers by energy
    std::vector<double> shr_energy_v;

    for (auto const& shr : shr_v)
      shr_energy_v.push_back( shr.Energy() );


    std::reverse(shr_energy_v.begin(), shr_energy_v.end());

    std::vector<size_t> sorted_idx_v;

    for (auto const& E : shr_energy_v) {
      for (size_t idx=0; idx < shr_v.size(); idx++){
	if (E == shr_v.at(idx).Energy())
	  sorted_idx_v.push_back( idx );
      }// for all reco showers
    }// for all energy values

      
    // if sorted idx list size different
    // than shower vector -> error!

    if (sorted_idx_v.size() != shr_v.size() )
      print(larlite::msg::kERROR,__FUNCTION__," did not sort all showers successfully...");

    // STEP 2
    // now match to true showers
    std::vector<size_t> matched_indices;

    // find best matching reco shower for each MC shower
    for (auto const& mcshr : mcs_v) {

      double dotmax = -1.;
      size_t idxmax = 0.;
      
      // loop through reco indices
      for (auto const& idx : sorted_idx_v) {

	// has this index already been used?
	if (std::find(matched_indices.begin(), matched_indices.end(), idx) != matched_indices.end() ) continue;
	
	// grab reco shower
	auto const& rcshr = shr_v.at(idx);
	
	double dot = rcshr.Direction().Dot( mcshr.Start().Momentum().Vect() );
	dot /= mcshr.Start().Momentum().Vect().Mag();
	dot /= rcshr.Direction().Mag();
	
	if (dot > dotmax) { dotmax = dot; idxmax = idx; }
	
      }// for all sorted indices

      std::cout << "dot product max = " << dotmax << std::endl;
      
      matched_indices.push_back( idxmax );

    }// for all true showers
    
    // did we not find two matching showers?
    if (matched_indices.size() != 2)
      print(larlite::msg::kERROR,__FUNCTION__," did not find two reco matches!");
    
    return std::make_pair( matched_indices[0], matched_indices[1] );
    
  }// and of function
  
}
#endif
