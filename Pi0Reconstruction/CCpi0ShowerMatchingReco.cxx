#ifndef LARLITE_CCPI0SHOWERMATCHINGRECO_CXX
#define LARLITE_CCPI0SHOWERMATCHINGRECO_CXX

#include "CCpi0ShowerMatchingReco.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool CCpi0ShowerMatchingReco::initialize() {

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
    _tree->Branch("_rc_shr_x",&_rc_shr_x,"rc_shr_x/D");
    _tree->Branch("_rc_shr_y",&_rc_shr_y,"rc_shr_y/D");
    _tree->Branch("_rc_shr_z",&_rc_shr_z,"rc_shr_z/D");
    _tree->Branch("_rc_shr_e",&_rc_shr_e,"rc_shr_e/D");
    _tree->Branch("_rc_shr_px",&_rc_shr_px,"rc_shr_px/D");
    _tree->Branch("_rc_shr_py",&_rc_shr_py,"rc_shr_py/D");
    _tree->Branch("_rc_shr_pz",&_rc_shr_pz,"rc_shr_pz/D");

    _tree->Branch("_rcradlen",&_rcradlen,"rcradlen/D");
    _tree->Branch("_mcradlen1",&_mcradlen1,"mcradlen1/D");
    _tree->Branch("_mcradlen2",&_mcradlen2,"mcradlen2/D");

    // reco cluster info
    _tree->Branch("_ip",&_ip,"ip/D");
    _tree->Branch("_lin",&_lin,"lin/D");
    _tree->Branch("_ssv",&_ssv,"ssv/D");
    _tree->Branch("_slope",&_slope,"slope/D");
    _tree->Branch("_slopedirangle",&_slopedirangle,"slopedirangle/D");

    // MC -> RC shower comparisons
    _tree->Branch("_dot",&_dot,"dot/D");
    _tree->Branch("_strt",&_strt,"strt/D");
    _tree->Branch("_emc",&_emc,"emc/D");

    // pi0 related MC information
    _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
    _tree->Branch("_pi_e",&_pi0_e,"pi0_e/D");
    _tree->Branch("_mc_oangle",&_mc_oangle,"mc_oangle/D");
    _tree->Branch("_mc_mass"  ,&_mc_mass  ,"_mc_mass/D" );
    _tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");
    
    
    return true;
  }
  
  bool CCpi0ShowerMatchingReco::analyze(storage_manager* storage) {

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


    if (ev_shower->size() == 0) return true;
    
    // grab clusters associated with the reconstructed shower
    auto reco_shr_cluster_v = GetRecoShowerClusters(storage, ev_shower);
    
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
    
    // MC <-> RC matching
    auto MCRCmatch_v = Match(pi0_mcshower_v, reco_shower_v);

    for (size_t rcidx = 0; rcidx < MCRCmatch_v.size(); rcidx++) {

      auto const& mcidx = MCRCmatch_v[rcidx];

      auto const& rcshr = ev_shower->at( rcidx );
      auto const& mcshr = pi0_mcshower_v.at( mcidx );

      // get clusters associated to this shower
      auto reco_clusters = reco_shr_cluster_v.at( rcidx );
      // fill cluster linearity info
      std::vector<double> hit_w_v, hit_t_v;
      //hit_w_v.resize( reco_clusters.at(2).size() );
      //hit_t_v.resize( reco_clusters.at(2).size() );
      for (auto const& hit : reco_clusters.at(2)) {
	hit_w_v.push_back( hit.WireID().Wire * _w2cm - _vtx_w_cm[2] );
	hit_t_v.push_back( hit.PeakTime()    * _t2cm - _vtx_t_cm[2] );
      }

      ::twodimtools::Linearity clusterlin(hit_w_v,hit_t_v);
      _ip = clusterlin.IP(0.,0.);
      _ssv = clusterlin._summed_square_variance;
      _lin = clusterlin._local_lin_truncated_avg;
      _slope = clusterlin._slope;

      double slope3D = rcshr.Direction().X()/rcshr.Direction().Z();
      slope3D       /= sqrt( ( rcshr.Direction().X() * rcshr.Direction().X() ) +
			     ( rcshr.Direction().Z() * rcshr.Direction().Z() ) );
      
      _slopedirangle = atan(slope3D-_slope)/(1+_slope*slope3D);
      
      _rc_shr_e = rcshr.Energy();
      _rc_shr_x = rcshr.ShowerStart().X();
      _rc_shr_y = rcshr.ShowerStart().Y();
      _rc_shr_z = rcshr.ShowerStart().Z();

      double mom = sqrt( ( rcshr.Direction().X() * rcshr.Direction().X() ) +
			 ( rcshr.Direction().Y() * rcshr.Direction().Y() ) +
			 ( rcshr.Direction().Z() * rcshr.Direction().Z() ) );
      
      _rc_shr_px = rcshr.Direction().X() / mom;
      _rc_shr_py = rcshr.Direction().Y() / mom;
      _rc_shr_pz = rcshr.Direction().Z() / mom;

      _rcradlen = sqrt( ( (_rc_shr_x - _rc_vtx_x) * (_rc_shr_x - _rc_vtx_x) ) +
			( (_rc_shr_y - _rc_vtx_y) * (_rc_shr_y - _rc_vtx_y) ) +
			( (_rc_shr_z - _rc_vtx_z) * (_rc_shr_z - _rc_vtx_z) ) );
      
      
      _dot  = rcshr.Direction().Dot( mcshr.DetProfile().Momentum().Vect() );
      _dot /= mcshr.DetProfile().Momentum().Vect().Mag();
      _dot /= rcshr.Direction().Mag();

      _emc = mcshr.DetProfile().E();

      _strt = (rcshr.ShowerStart() - mcshr.DetProfile().Position().Vect()).Mag();
      
      _tree->Fill();


    }// loop through RC showers

    return true;
  }
  
  bool CCpi0ShowerMatchingReco::finalize() {

    if (_fout) _fout->cd();
    _tree->Write();
    
    return true;
  }

  void CCpi0ShowerMatchingReco::Reset() {

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

    _rc_shr_x=  _rc_shr_y=  _rc_shr_z= 0;
    _rc_shr_px= _rc_shr_py= _rc_shr_pz= 0;
    _rc_shr_e= 0;
    
    return;
  }

  std::vector<int> CCpi0ShowerMatchingReco::Match(const std::vector<larlite::mcshower>& mcs_v,
						  const std::vector<larlite::shower>&   shr_v) {

    // now match to true showers
    std::vector<int> matched_indices;

    // find best matching MC shower for each reco shower
    for (auto const& rcshr : shr_v) {

      double dotmax = -1.;
      size_t idxmax = 0.;
      
      // loop through mc indices
      for (size_t i=0; i < mcs_v.size(); i++) {

	auto const& mcshr = mcs_v[i];

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

  void CCpi0ShowerMatchingReco::MinDWall() {

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


  std::vector< std::vector< std::vector< larlite::hit > > > CCpi0ShowerMatchingReco::GetRecoShowerClusters(larlite::storage_manager* storage, larlite::event_shower* ev_shower) {


    larlite::event_pfpart*  ev_pfpart;
    larlite::event_cluster* ev_cluster;
    larlite::event_hit*     ev_hit;

    auto ass_shr_pfp_v = storage->find_one_ass( ev_shower->id() , ev_pfpart , ev_shower->name()  );
    auto ass_pfp_cls_v = storage->find_one_ass( ev_pfpart->id() , ev_cluster, ev_pfpart->name()  );
    auto ass_cls_hit_v = storage->find_one_ass( ev_cluster->id(), ev_hit    , ev_cluster->name() );

    std::vector< std::vector< std::vector< larlite::hit > > > shr_v_hits;

    // for every PFParticle associated with each shower
    for (auto const& ass_shr_pfp : ass_shr_pfp_v) {

      // new vector for this shower (one entry per plane)
      std::vector< std::vector< larlite::hit> > shr_hits;
      shr_hits.resize(3);

      // for every list of clusters associated with the PFParticle
      for (auto const& ass_pfp_clus : ass_shr_pfp) {

	for (auto const& clus_idx : ass_pfp_cls_v[ass_pfp_clus]) {
	  
	  // new cluster -> create vector of hits for this cluster
	  std::vector< larlite::hit > hit_v;
	  
	  // grab the hit indices associated with this cluster
	  auto const& hit_idx_v = ass_cls_hit_v.at( clus_idx );
	  
	  hit_v.reserve(hit_idx_v.size());
	  
	  for  (auto const& hit_idx : hit_idx_v) {
	    
	    hit_v.push_back( ev_hit->at(hit_idx) );
	    
	  }// for all hit indices

	  shr_hits.at( ev_hit->at(hit_idx_v[0]).WireID().Plane ) =  hit_v ;

	}// for all clusters associated to the PFPart
	
      }// for all PFParts

      shr_v_hits.push_back( shr_hits );
    }// for all showers

    return shr_v_hits;
  }

  bool CCpi0ShowerMatchingReco::loadVertex(event_vertex* ev_vtx) {
    
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
