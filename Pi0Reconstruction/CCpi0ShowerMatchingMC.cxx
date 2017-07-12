#ifndef LARLITE_CCPI0SHOWERMATCHINGMC_CXX
#define LARLITE_CCPI0SHOWERMATCHINGMC_CXX

#include "CCpi0ShowerMatchingMC.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool CCpi0ShowerMatchingMC::initialize() {

    _w2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _t2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    if (_tree) delete _tree;
    _tree = new TTree("_tree","shower study tree");

    _tree->Branch("_n_reco_showers",&_n_reco_showers,"n_reco_showers/I");
    _tree->Branch("_event",&_event,"event/I");

    // vertex info
    _tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
    _tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
    _tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");
    _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D");
    _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D");
    _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D");

    // mc shower info
    _tree->Branch("_mc_shr1_x",&_mc_shr1_x,"mc_shr1_x/D");
    _tree->Branch("_mc_shr1_y",&_mc_shr1_y,"mc_shr1_y/D");
    _tree->Branch("_mc_shr1_z",&_mc_shr1_z,"mc_shr1_z/D");
    _tree->Branch("_mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/D");
    _tree->Branch("_mc_shr1_px",&_mc_shr1_px,"mc_shr1_px/D");
    _tree->Branch("_mc_shr1_py",&_mc_shr1_py,"mc_shr1_py/D");
    _tree->Branch("_mc_shr1_pz",&_mc_shr1_pz,"mc_shr1_pz/D");
    _tree->Branch("_mc_shr2_x",&_mc_shr2_x,"mc_shr2_x/D");
    _tree->Branch("_mc_shr2_y",&_mc_shr2_y,"mc_shr2_y/D");
    _tree->Branch("_mc_shr2_z",&_mc_shr2_z,"mc_shr2_z/D");
    _tree->Branch("_mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/D");
    _tree->Branch("_mc_shr2_px",&_mc_shr2_px,"mc_shr2_px/D");
    _tree->Branch("_mc_shr2_py",&_mc_shr2_py,"mc_shr2_py/D");
    _tree->Branch("_mc_shr2_pz",&_mc_shr2_pz,"mc_shr2_pz/D");

    // reco shower info
    _tree->Branch("_mc_shr_x",&_mc_shr_x,"mc_shr_x/D");
    _tree->Branch("_mc_shr_y",&_mc_shr_y,"mc_shr_y/D");
    _tree->Branch("_mc_shr_z",&_mc_shr_z,"mc_shr_z/D");
    _tree->Branch("_mc_shr_e",&_mc_shr_e,"mc_shr_e/D");
    _tree->Branch("_mc_shr_px",&_mc_shr_px,"mc_shr_px/D");
    _tree->Branch("_mc_shr_py",&_mc_shr_py,"mc_shr_py/D");
    _tree->Branch("_mc_shr_pz",&_mc_shr_pz,"mc_shr_pz/D");

    _tree->Branch("_mcradlen",&_mcradlen,"mcradlen/D");
    _tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
    _tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

    // reco cluster info
    _tree->Branch("_ip",&_ip,"ip/D");
    _tree->Branch("_lin",&_lin,"lin/D");
    _tree->Branch("_ssv",&_ssv,"ssv/D");
    _tree->Branch("_slope",&_slope,"slope/D");

    // MC -> RC shower comparisons
    _tree->Branch("_dot",&_dot,"dot/D");
    _tree->Branch("_strt",&_strt,"strt/D");
    _tree->Branch("_erc",&_erc,"erc/D");

    // pi0 related MC information
    _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
    _tree->Branch("_n_trk",&_n_trk,"n_trk/I");
    _tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");
    _tree->Branch("_mc_oangle",&_mc_oangle,"mc_oangle/D");
    _tree->Branch("_mc_mass"  ,&_mc_mass  ,"_mc_mass/D" );
    _tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");
    
    
    return true;
  }
  
  bool CCpi0ShowerMatchingMC::analyze(storage_manager* storage) {

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

    if (loadVertex(ev_vertex) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
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
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	n_pi0    += 1;
	pi0_idx   = i;
	pi0_trkid = part.TrackId();
	_mc_vtx_x = part.Trajectory().at( 0 ).X();
	_mc_vtx_y = part.Trajectory().at( 0 ).Y();
	_mc_vtx_z = part.Trajectory().at( 0 ).Z();
      }
      if ( ((part.PdgCode() == 13) || (part.PdgCode() == 2212) || (fabs(part.PdgCode()) == 211)) and
	   (part.StatusCode() == 1) and
	   (part.Trajectory().at(0).Momentum().Vect().Mag() > 0.3) )
	_n_trk += 1;
    }

   MinDWall();
    
    if (n_pi0 != 1)
      return false;
    
    // is the vertex in the TPC?
    auto const& vtx = part_v.at(pi0_idx).Trajectory().at(0);
    if ( (vtx.X() < 0) or (vtx.X() > 256) or (vtx.Y() < -116) or (vtx.Y() > 116) or (vtx.Z() < 0) or (vtx.Z() > 1036) ){
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
    
    _mc_shr1_e  = mcshr1.DetProfile().E();
    _mc_shr1_x  = mcshr1.DetProfile().X();
    _mc_shr1_y  = mcshr1.DetProfile().Y();
    _mc_shr1_z  = mcshr1.DetProfile().Z();
    double mom1 = sqrt ( ( mcshr1.DetProfile().Px() * mcshr1.DetProfile().Px() ) +
			 ( mcshr1.DetProfile().Py() * mcshr1.DetProfile().Py() ) +
			 ( mcshr1.DetProfile().Pz() * mcshr1.DetProfile().Pz() ) );
    
    _mc_shr1_px = mcshr1.DetProfile().Px() / mom1;
    _mc_shr1_py = mcshr1.DetProfile().Py() / mom1;
    _mc_shr1_pz = mcshr1.DetProfile().Pz() / mom1;
    
    _mc_shr2_e  = mcshr2.DetProfile().E();
    _mc_shr2_x  = mcshr2.DetProfile().X();
    _mc_shr2_y  = mcshr2.DetProfile().Y();
    _mc_shr2_z  = mcshr2.DetProfile().Z();
    double mom2 = sqrt ( ( mcshr2.DetProfile().Px() * mcshr2.DetProfile().Px() ) +
			 ( mcshr2.DetProfile().Py() * mcshr2.DetProfile().Py() ) +
			 ( mcshr2.DetProfile().Pz() * mcshr2.DetProfile().Pz() ) );
    
    _mc_shr2_px = mcshr2.DetProfile().Px() / mom2;
    _mc_shr2_py = mcshr2.DetProfile().Py() / mom2;
    _mc_shr2_pz = mcshr2.DetProfile().Pz() / mom2;

    _mc_oangle  = mcshr1.DetProfile().Momentum().Vect().Dot( mcshr2.DetProfile().Momentum().Vect() );
    _mc_oangle /= mcshr1.DetProfile().Momentum().Vect().Mag();
    _mc_oangle /= mcshr2.DetProfile().Momentum().Vect().Mag();

    _mcradlen1 = sqrt( ( (_mc_shr1_x - _mc_vtx_x) * (_mc_shr1_x - _mc_vtx_x) ) +
		       ( (_mc_shr1_y - _mc_vtx_y) * (_mc_shr1_y - _mc_vtx_y) ) +
		       ( (_mc_shr1_z - _mc_vtx_z) * (_mc_shr1_z - _mc_vtx_z) ) );

    _mcradlen2 = sqrt( ( (_mc_shr2_x - _mc_vtx_x) * (_mc_shr2_x - _mc_vtx_x) ) +
		       ( (_mc_shr2_y - _mc_vtx_y) * (_mc_shr2_y - _mc_vtx_y) ) +
		       ( (_mc_shr2_z - _mc_vtx_z) * (_mc_shr2_z - _mc_vtx_z) ) );

    std::vector<larlite::mcshower> pi0_mcshower_v = {mcshr1,mcshr2};

    _mc_mass = sqrt( 2 * _mc_shr1_e * _mc_shr2_e * ( 1 - _mc_oangle ) );
    
    _n_reco_showers = ev_shower->size();
    
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

    //std::cout << "There are " << reco_shower_v.size() << " reco showers in the event" << std::endl;
    
    // MC <-> RC matching
    auto MCRCmatch_v = Match(pi0_mcshower_v, reco_shower_v);

    for (size_t mcidx = 0; mcidx < MCRCmatch_v.size(); mcidx++) {

      auto const& mcshr = pi0_mcshower_v.at( mcidx );

      auto const& rcidx = MCRCmatch_v[mcidx];
 
      _mc_shr_e = mcshr.DetProfile().E();
      _mc_shr_x = mcshr.DetProfile().X();
      _mc_shr_y = mcshr.DetProfile().Y();
      _mc_shr_z = mcshr.DetProfile().Z();

      if (_mc_shr_e <= 10) continue;

      double momentum = sqrt ( ( mcshr.DetProfile().Px() * mcshr.DetProfile().Px() ) +
			       ( mcshr.DetProfile().Py() * mcshr.DetProfile().Py() ) +
			       ( mcshr.DetProfile().Pz() * mcshr.DetProfile().Pz() ) );
      
      _mc_shr_px = mcshr.DetProfile().Px() / momentum;
      _mc_shr_py = mcshr.DetProfile().Py() / momentum;
      _mc_shr_pz = mcshr.DetProfile().Pz() / momentum;

      _mcradlen = sqrt( ( (_mc_shr_x - _mc_vtx_x) * (_mc_shr_x - _mc_vtx_x) ) +
			( (_mc_shr_y - _mc_vtx_y) * (_mc_shr_y - _mc_vtx_y) ) +
			( (_mc_shr_z - _mc_vtx_z) * (_mc_shr_z - _mc_vtx_z) ) );

      if (rcidx == -1) {
	//std::cout << "\t\t MCshower did not find a match..." << std::endl;
	_dot  = -1;
	_erc  = -1;
	_strt = -1;
	_tree->Fill();
	continue;
      }

      auto const& rcshr = reco_shower_v[ rcidx ];
      
      _dot  = rcshr.Direction().Dot( mcshr.DetProfile().Momentum().Vect() );
      _dot /= mcshr.DetProfile().Momentum().Vect().Mag();
      _dot /= rcshr.Direction().Mag();

      _erc = rcshr.Energy();

      _strt = (rcshr.ShowerStart() - mcshr.DetProfile().Position().Vect()).Mag();     

      _tree->Fill();


    }// loop through RC showers

    return true;
  }
  
  bool CCpi0ShowerMatchingMC::finalize() {

    if (_fout) _fout->cd();
    _tree->Write();
    
    return true;
  }

  void CCpi0ShowerMatchingMC::Reset() {

    _n_trk = 0;
    _n_reco_showers = 0;
    _nu_e = 0;
    _pi0_e = 0;

    _mc_vtx_x= _mc_vtx_y= _mc_vtx_z= 0;
    _rc_vtx_x= _rc_vtx_y= _rc_vtx_z= 0;
    
    _mc_shr1_x=  _mc_shr1_y=  _mc_shr1_z= 0;
    _mc_shr1_px= _mc_shr1_py= _mc_shr1_pz= 0;
    _mc_shr1_e= 0;
    
    _mc_shr2_x=  _mc_shr2_y=  _mc_shr2_z= 0;
    _mc_shr2_px= _mc_shr2_py= _mc_shr2_pz= 0;
    _mc_shr2_e= 0;

    _mc_shr_x=  _mc_shr_y=  _mc_shr_z= 0;
    _mc_shr_px= _mc_shr_py= _mc_shr_pz= 0;
    _mc_shr_e= 0;
    
    return;
  }

  std::vector<int> CCpi0ShowerMatchingMC::Match(const std::vector<larlite::mcshower>& mcs_v,
						const std::vector<larlite::shower>&   shr_v) {

    // now match to true showers
    std::vector<int> matched_indices;

    // find best matching RC shower for each MC shower
    for (auto const& mcshr : mcs_v) {

      double dotmax = -2.;
      int    idxmax = -1;
      
      // loop through mc indices
      for (size_t i=0; i < shr_v.size(); i++) {

	auto const& rcshr = shr_v[i];

	// in this module we want all reco'd showers to be matched to the
	// most compatible true shower. Don't skip already matched mcshowers
	
	double dot = rcshr.Direction().Dot( mcshr.DetProfile().Momentum().Vect() );
	dot /= mcshr.DetProfile().Momentum().Vect().Mag();
	dot /= rcshr.Direction().Mag();
	
	if (dot > dotmax) { dotmax = dot; idxmax = i; }
	
      }// for all sorted indices
      
      matched_indices.push_back( idxmax );

    }// for all true showers
    
    return matched_indices;
    
  }// and of function

  void CCpi0ShowerMatchingMC::MinDWall() {

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


  bool CCpi0ShowerMatchingMC::loadVertex(event_vertex* ev_vtx) {
    
    if (ev_vtx->size() != 1) return false;
    
    // get vertex position on each plane
    if ( (ev_vtx->size() == 1) ){
      auto const& vtx = ev_vtx->at(0);

      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      
      auto geoH = larutil::GeometryHelper::GetME();
      auto geom = larutil::Geometry::GetME();

      for (size_t pl = 0; pl < 3; pl++){
	double *origin;
	origin = new double[3];
	geom->PlaneOriginVtx(pl,origin);
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _t2cm - origin[0];
      }
    }

    return true;
  }
  
}
#endif
