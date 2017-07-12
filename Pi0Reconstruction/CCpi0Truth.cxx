#ifndef LARLITE_CCPI0TRUTH_CXX
#define LARLITE_CCPI0TRUTH_CXX

#include "CCpi0Truth.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool CCpi0Truth::initialize() {

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
    _tree->Branch("_mc_shr1_edep",&_mc_shr1_edep,"mc_shr1_edep/D");
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
    _tree->Branch("_mc_shr2_edep",&_mc_shr2_edep,"mc_shr2_edep/D");
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

    _tree->Branch("_rcradlen1",&_rcradlen1,"rcradlen1/D");
    _tree->Branch("_rcradlen2",&_rcradlen2,"rcradlen2/D");
    _tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
    _tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

    // shower correlations
    _tree->Branch("_dot1",&_dot1,"dot1/D");
    _tree->Branch("_dot2",&_dot2,"dot2/D");
    _tree->Branch("_strt1",&_strt1,"strt1/D");
    _tree->Branch("_strt2",&_strt2,"strt2/D");

    _tree->Branch("_rc_oangle",&_rc_oangle,"rc_oangle/D");
    _tree->Branch("_mc_oangle",&_mc_oangle,"mc_oangle/D");
    _tree->Branch("_mc_mass"  ,&_mc_mass  ,"_mc_mass/D" );
    _tree->Branch("_mc_mass_edep"  ,&_mc_mass_edep  ,"_mc_mass_edep/D" );
    _tree->Branch("_rc_mass"  ,&_rc_mass  ,"_rc_mass/D" );

    _tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");
    
    
    return true;
  }
  
  bool CCpi0Truth::analyze(storage_manager* storage) {

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

    _n_reco_showers = ev_shower->size();
    
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

   MinDWall();
    
   if (n_pi0 != 1){
     _tree->Fill();
      return false;
   }
    
    // is the vertex in the TPC?
    auto const& vtx = part_v.at(pi0_idx).Trajectory().at(0);
    if ( (vtx.X() < 0) or (vtx.X() > 256) or (vtx.Y() < -116) or (vtx.Y() > 116) or (vtx.Z() < 0) or (vtx.Z() > 1036) ){
      std::cout << "Vertex @ x : " << vtx.X() << "\t y : " << vtx.Y() << "\t z : " << vtx.Z() << std::endl;
      _tree->Fill();
      return false;
    }
    
    size_t pi0_ID_1 = 0;
    size_t pi0_ID_2 = 0;
    size_t idx_1 = 0;
    size_t idx_2 = 0;
    size_t n_found = 0;
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

    size_t idxLARGE = idx_1;
    size_t idxSMALL = idx_2;

    if (ev_mcshower->at(idx_1).Start().E() < ev_mcshower->at(idx_2).Start().E() )
      { idxLARGE = idx_2; idxSMALL = idx_1; }
    
    auto const& mcshr1 = ev_mcshower->at(idxLARGE);
    auto const& mcshr2 = ev_mcshower->at(idxSMALL);
    
    _pi0_e  = mcshr1.MotherEnd().E();
    
    _mc_shr1_e  = mcshr1.Start().E();
    _mc_shr1_edep  = mcshr1.DetProfile().E();
    _mc_shr1_x  = mcshr1.DetProfile().X();
    _mc_shr1_y  = mcshr1.DetProfile().Y();
    _mc_shr1_z  = mcshr1.DetProfile().Z();
    double mom1 = sqrt ( ( mcshr1.Start().Px() * mcshr1.Start().Px() ) +
			 ( mcshr1.Start().Py() * mcshr1.Start().Py() ) +
			 ( mcshr1.Start().Pz() * mcshr1.Start().Pz() ) );
    
    _mc_shr1_px = mcshr1.Start().Px() / mom1;
    _mc_shr1_py = mcshr1.Start().Py() / mom1;
    _mc_shr1_pz = mcshr1.Start().Pz() / mom1;
    
    _mc_shr2_e  = mcshr2.Start().E();
    _mc_shr2_edep  = mcshr2.DetProfile().E();
    _mc_shr2_x  = mcshr2.DetProfile().X();
    _mc_shr2_y  = mcshr2.DetProfile().Y();
    _mc_shr2_z  = mcshr2.DetProfile().Z();
    double mom2 = sqrt ( ( mcshr2.Start().Px() * mcshr2.Start().Px() ) +
			 ( mcshr2.Start().Py() * mcshr2.Start().Py() ) +
			 ( mcshr2.Start().Pz() * mcshr2.Start().Pz() ) );
    
    _mc_shr2_px = mcshr2.Start().Px() / mom2;
    _mc_shr2_py = mcshr2.Start().Py() / mom2;
    _mc_shr2_pz = mcshr2.Start().Pz() / mom2;

    _mc_oangle  = mcshr1.Start().Momentum().Vect().Dot( mcshr2.Start().Momentum().Vect() );
    _mc_oangle /= mcshr1.Start().Momentum().Vect().Mag();
    _mc_oangle /= mcshr2.Start().Momentum().Vect().Mag();

    _mcradlen1 = sqrt( ( (_mc_shr1_x - _mc_vtx_x) * (_mc_shr1_x - _mc_vtx_x) ) +
		       ( (_mc_shr1_y - _mc_vtx_y) * (_mc_shr1_y - _mc_vtx_y) ) +
		       ( (_mc_shr1_z - _mc_vtx_z) * (_mc_shr1_z - _mc_vtx_z) ) );

    _mcradlen2 = sqrt( ( (_mc_shr2_x - _mc_vtx_x) * (_mc_shr2_x - _mc_vtx_x) ) +
		       ( (_mc_shr2_y - _mc_vtx_y) * (_mc_shr2_y - _mc_vtx_y) ) +
		       ( (_mc_shr2_z - _mc_vtx_z) * (_mc_shr2_z - _mc_vtx_z) ) );

    std::vector<larlite::mcshower> pi0_mcshower_v = {mcshr1,mcshr2};

    _mc_mass = sqrt( 2 * _mc_shr1_e * _mc_shr2_e * ( 1 - _mc_oangle ) );
    _mc_mass_edep = sqrt( 2 * _mc_shr1_edep * _mc_shr2_edep * ( 1 - _mc_oangle ) );
    
    _tree->Fill();

    return true;
  }
  
  bool CCpi0Truth::finalize() {

    if (_fout) _fout->cd();
    _tree->Write();
    
    return true;
  }

  void CCpi0Truth::Reset() {

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

  std::pair<int,int> CCpi0Truth::Match(const std::vector<larlite::mcshower>& mcs_v,
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
	
	double dot = rcshr.Direction().Dot( mcshr.DetProfile().Momentum().Vect() );
	dot /= mcshr.DetProfile().Momentum().Vect().Mag();
	dot /= rcshr.Direction().Mag();
	
	if (dot > dotmax) { dotmax = dot; idxmax = idx; }
	
      }// for all sorted indices

      matched_indices.push_back( idxmax );

    }// for all true showers
    
    // did we not find two matching showers?
    if (matched_indices.size() != 2)
      print(larlite::msg::kERROR,__FUNCTION__," did not find two reco matches!");
    
    return std::make_pair( matched_indices[0], matched_indices[1] );
    
  }// and of function

  void CCpi0Truth::MinDWall() {

    _dwallmin = 1000;

    if ( (_mc_vtx_x > 0) && (_mc_vtx_x < _dwallmin) )
      _dwallmin = _mc_vtx_x;

    if ( (_mc_vtx_x < 256) && ( (256-_mc_vtx_x) < _dwallmin) )
      _dwallmin = (256.-_mc_vtx_x);

    if ( (_mc_vtx_z > 0) && (_mc_vtx_z < _dwallmin) )
      _dwallmin = _mc_vtx_z;

    if ( (_mc_vtx_z < 1030) && ( (1030-_mc_vtx_z) < _dwallmin) )
      _dwallmin = (1030.-_mc_vtx_z);

    if ( (_mc_vtx_y > -116) && ( (_mc_vtx_y+116) < _dwallmin) )
      _dwallmin = (_mc_vtx_y+116);

    if ( (_mc_vtx_y < 116) && ( (116-_mc_vtx_y) < _dwallmin) )
      _dwallmin = (116.-_mc_vtx_y);
  }
  
}
#endif
