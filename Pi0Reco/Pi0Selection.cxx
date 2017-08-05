#ifndef LARLITE_PI0SELECTION_CXX
#define LARLITE_PI0SELECTION_CXX

#include "Pi0Selection.h"

namespace larlite {

  Pi0Selection::Pi0Selection()
    : _tree(nullptr)
    , _pi0_tree(nullptr)
    , _shower_tree(nullptr)
  { 
    _name="Pi0Selection";
    _fout=0;
    _edepmin = 0;
    _containmentcorrection = false;
  }

  bool Pi0Selection::initialize() {

    _ctr = -1;

    _TPC = ::geoalgo::AABox_t(0,-116,0,256,116,1036);

    if (_pi0_tree) delete _pi0_tree;
    _pi0_tree = new TTree("_pi0_tree","pi0 tree");
    _pi0_tree->Branch("ip",&_ip,"ip/D");
    _pi0_tree->Branch("el",&_el,"el/D");
    _pi0_tree->Branch("eh",&_eh,"eh/D");
    _pi0_tree->Branch("rl",&_rl,"rl/D");
    _pi0_tree->Branch("rh",&_rh,"rh/D");
    _pi0_tree->Branch("angle",&_angle,"angle/D");
    _pi0_tree->Branch("mass",&_mass,"mass/D");
    _pi0_tree->Branch("run",&_run,"run/I");
    _pi0_tree->Branch("sub",&_sub,"sub/I");
    _pi0_tree->Branch("evt",&_evt,"evt/I");
    _pi0_tree->Branch("ctr",&_ctr,"ctr/I");
    _pi0_tree->Branch("npi0",&_npi0,"npi0/I");
    _pi0_tree->Branch("ngamma",&_ngamma,"ngamma/I");
    _pi0_tree->Branch("ipvtx",&_ipvtx,"ipvtx/D");
    
    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("ip",&_ip,"ip/D");
    _tree->Branch("el",&_el,"el/D");
    _tree->Branch("eh",&_eh,"eh/D");
    _tree->Branch("rl",&_rl,"rl/D");
    _tree->Branch("rh",&_rh,"rh/D");
    _tree->Branch("angle",&_angle,"angle/D");
    _tree->Branch("mass",&_mass,"mass/D");

    // reconstructed positions
    _tree->Branch("rc0x",&_rc0x,"rc0x/D");
    _tree->Branch("rc0y",&_rc0y,"rc0y/D");
    _tree->Branch("rc0z",&_rc0z,"rc0z/D");
    _tree->Branch("rc1x",&_rc1x,"rc1x/D");
    _tree->Branch("rc1y",&_rc1y,"rc1y/D");
    _tree->Branch("rc1z",&_rc1z,"rc1z/D");

    // true shower positions
    _tree->Branch("mc0x",&_mc0x,"mc0x/D");
    _tree->Branch("mc0y",&_mc0y,"mc0y/D");
    _tree->Branch("mc0z",&_mc0z,"mc0z/D");
    _tree->Branch("mc1x",&_mc1x,"mc1x/D");
    _tree->Branch("mc1y",&_mc1y,"mc1y/D");
    _tree->Branch("mc1z",&_mc1z,"mc1z/D");

    // vertex (true and reco)
    _tree->Branch("mcvtxx",&_mcvtxx,"mcvtxx/D");
    _tree->Branch("mcvtxy",&_mcvtxy,"mcvtxy/D");
    _tree->Branch("mcvtxz",&_mcvtxz,"mcvtxz/D");
    _tree->Branch("rcvtxx",&_rcvtxx,"rcvtxx/D");
    _tree->Branch("rcvtxy",&_rcvtxy,"rcvtxy/D");
    _tree->Branch("rcvtxz",&_rcvtxz,"rcvtxz/D");

    // vertex -> IP midpoint separation [reconstructed quantities only]
    _tree->Branch("ipvtx",&_ipvtx,"ipvtx/D");

    // rc -> mc separation
    _tree->Branch("d0",&_d0,"d0/D");
    _tree->Branch("d1",&_d1,"d1/D");

    // rc -> mc angle
    _tree->Branch("angle0",&_angle0,"angle0/D");
    _tree->Branch("angle1",&_angle1,"angle1/D");

    // reco & true energies
    _tree->Branch("rce0",&_rce0,"rce0/D");
    _tree->Branch("rce1",&_rce1,"rce1/D");
    _tree->Branch("mce0",&_mce0,"mce0/D");
    _tree->Branch("mce1",&_mce1,"mce1/D");
    _tree->Branch("mcedep0",&_mcedep0,"mcedep0/D");
    _tree->Branch("mcedep1",&_mcedep1,"mcedep1/D");

    // true pi0 information
    _tree->Branch("pi0px",&_pi0px,"pi0px/D");
    _tree->Branch("pi0py",&_pi0py,"pi0py/D");
    _tree->Branch("pi0pz",&_pi0pz,"pi0pz/D");
    _tree->Branch("pi0e", &_pi0e, "pi0e/D" );
    _tree->Branch("pi0a", &_pi0a, "pi0a/D" );
    _tree->Branch("npi0",&_npi0,"npi0/I");
    _tree->Branch("ngamma",&_ngamma,"ngamma/I");

    _tree->Branch("nrecoshr",&_nrecoshr,"nrecoshr/I");
    _tree->Branch("nrecoshrcut",&_nrecoshrcut,"nrecoshrcut/I");

    // event-wise information
    _tree->Branch("run",&_run,"run/I");
    _tree->Branch("sub",&_sub,"sub/I");
    _tree->Branch("evt",&_evt,"evt/I");
    _tree->Branch("ctr",&_ctr,"ctr/I");

    /// 
    /// Tree to hold shower-by-shower information
    /// for showers from a reconstructed pi0
    ///
    if (_shower_tree) delete _shower_tree;
    _shower_tree = new TTree("_shower_tree","shower tree");
    // rc -> mc angle
    _shower_tree->Branch("anglediff",&_anglediff,"anglediff/D");
    // reco & true energies
    _shower_tree->Branch("rce",&_rce,"rce/D");
    _shower_tree->Branch("mce",&_mce,"mce/D");
    _shower_tree->Branch("mcedep",&_mcedep,"mcedep/D");
    // shower position
    _shower_tree->Branch("rcx",&_rcx,"rcx/D");
    _shower_tree->Branch("rcy",&_rcy,"rcy/D");
    _shower_tree->Branch("rcz",&_rcz,"rcz/D");
    /// shower momentum
    _shower_tree->Branch("rcpx",&_rcpx,"rcpx/D");
    _shower_tree->Branch("rcpy",&_rcpy,"rcpy/D");
    _shower_tree->Branch("rcpz",&_rcpz,"rcpz/D");
    // shower distance to wall
    _shower_tree->Branch("dwall",&_dwall,"dwall/D");
    // neutrino vertex
    _shower_tree->Branch("rcvtxx",&_rcvtxx,"rcvtxx/D");
    _shower_tree->Branch("rcvtxy",&_rcvtxy,"rcvtxy/D");
    _shower_tree->Branch("rcvtxz",&_rcvtxz,"rcvtxz/D");
    
    
    return true;
  }
  
  bool Pi0Selection::analyze(storage_manager* storage) {

    _npi0     = 0;
    _ngamma   = 0;
    _mass     = -1;

    _ctr += 1;

    std::vector<larlite::mcshower> pi0_gamma_v;

    auto ev_shr = storage->get_data<event_shower>("showerreco");
    auto ev_trk = storage->get_data<event_track> ("pi0reco");
    auto ev_vtx = storage->get_data<event_vertex>(_vtx_producer);
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctruth>("generator");

    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    _evt = storage->event_id();
    _sub = storage->subrun_id();
    _run = storage->run_id();

    _nrecoshr = ev_shr->size();

    if ( !ev_vtx || (ev_vtx->size() != 1) ) {
      print(larlite::msg::kWARNING,__FUNCTION__,"number of vertices != 1 -> skip event");
      return true;
    }

    // do we want to do MC matching?
    if (_mcmatch) {

      if ( !ev_mcs->size() || !ev_mct->size() ) {
	print(larlite::msg::kWARNING,__FUNCTION__,"No truth info available and mcmatching on -> quit.");
	return true;
      }

      // grab true vertex
      auto mctruth = ev_mct->at(0);
      auto part_v  = mctruth.GetParticles();

      for (auto const& p : part_v) {
	if ( (p.PdgCode() == 111) and (p.StatusCode() == 1) ) {
	  _npi0 += 1;
	  auto const& pvtx = p.Trajectory().at(0);
	  _mcvtx = ::geoalgo::Point_t( pvtx.X(), pvtx.Y(), pvtx.Z() );
	  _mcvtxx = _mcvtx[0];
	  _mcvtxy = _mcvtx[1];
	  _mcvtxz = _mcvtx[2];
	  _pi0e   = pvtx.E();
	}
	if ( (p.StatusCode() == 1) && ( (p.PdgCode() == 22) || (p.PdgCode() == 11) ) )
	  _ngamma += 1;
      }// for all particles

      // get ready to store the showers from the pi0
      for (auto const& mcs : *ev_mcs) {

	if ( mcs.AncestorPdgCode() != 111 ) continue;

	::geoalgo::Point_t shrvtx(mcs.Start().X(), mcs.Start().Y(), mcs.Start().Z() );
	if ( shrvtx.Dist( _mcvtx ) < 0.01 )
	  
	  pi0_gamma_v.push_back( mcs );

      }// for all MCShowers

      // do we have at least 2 truth-showers with minimum edep?
      if (_edepmin > 0) {
	if (pi0_gamma_v.size() < 2) return true;
	if ( (pi0_gamma_v[0].DetProfile().E() < _edepmin) || (pi0_gamma_v[1].DetProfile().E() < _edepmin) )
	  return true;
      }// if we should check for a minimum edep

    }// if we should do mc-matching

    auto vtx = ev_vtx->at(0);
    _rcvtx  = ::geoalgo::Point_t(vtx.X(),vtx.Y(),vtx.Z());
    _rcvtxx = _rcvtx[0];
    _rcvtxy = _rcvtx[1];
    _rcvtxz = _rcvtx[2];

    // first, filter showers and identify subset to be used for pi0 selection
    auto selected_shower_idx_v = FilterShowers(ev_shr);

    _nrecoshrcut = selected_shower_idx_v.size();

    // keep track of which pair leads to the best match, given the IP value
    double bestIP    = 40000.0;
    double bestAngle = -1;
    double bestMass  = -1;
    int    bestPair  = -1;
    _mass = -1;

    // if we have at least 2 shower candidates:
    if (selected_shower_idx_v.size() >= 2) {
    
      // get shower-pair combinatorics
      auto shr_pairs_idx_v = Combinatorics( selected_shower_idx_v );
      
      for (size_t pidx = 0; pidx < shr_pairs_idx_v.size(); pidx++ ) {
	
	auto const& shr_pair = shr_pairs_idx_v.at(pidx);
	
	auto const& shr1 = ev_shr->at(shr_pair.first);
	auto const& shr2 = ev_shr->at(shr_pair.second);
	
	// transform into half-line objects
	auto const& srt1 = shr1.ShowerStart();
	auto const& dir1 = shr1.Direction();
	auto const& HL1  = geoalgo::HalfLine( srt1.X(), srt1.Y(), srt1.Z(),
					      -dir1.X(),-dir1.Y(),-dir1.Z());
	
	auto const& srt2 = shr2.ShowerStart();
	auto const& dir2 = shr2.Direction();
	auto const& HL2  = geoalgo::HalfLine( srt2.X(), srt2.Y(), srt2.Z(),
					      -dir2.X(),-dir2.Y(),-dir2.Z());
	
	if ( shr1.Energy() < shr2.Energy() ) {
	  _el = shr1.Energy();
	  _eh = shr2.Energy();
	  _rl = HL1.Start().Dist( _rcvtx );
	  _rh = HL2.Start().Dist( _rcvtx );
	}
	else {
	  _el = shr2.Energy();
	  _eh = shr1.Energy();
	  _rl = HL2.Start().Dist( _rcvtx );
	  _rh = HL1.Start().Dist( _rcvtx );
	}
	
	::geoalgo::Point_t pt1, pt2;
	_ip = sqrt( _geoAlgo.SqDist(HL1,HL2,pt1,pt2) );
	_angle = HL1.Dir().Angle( HL2.Dir() );
	
	_ipvtx = _rcvtx.Dist( ((pt1+pt2)/2.) );

	// save mass value regardless
	_mass  = sqrt( 2 * shr1.Energy() * shr2.Energy() * ( 1 - cos(_angle) ) );
	
	// does it pass out cuts?
	if ( ( _ip < _ipmax ) && ( _angle > _anglemin * (3.14/180.) ) ) {
	  // is it the best match?
	  if (_ip < bestIP)
	    { 
	      bestIP    = _ip; 
	      bestPair  = pidx; 
	      bestAngle = _angle;
	      bestMass  = _mass;
	    }
	}
	
	_pi0_tree->Fill();
	
      }// for all shower-pairs

      
      _mass  = bestMass;
      _angle = bestAngle * 180. / 3.14;
      _ip    = bestIP;

      if (bestPair != -1) {
	std::vector<larlite::shower> rcshr_v = { ev_shr->at( shr_pairs_idx_v[bestPair].first) , ev_shr->at( shr_pairs_idx_v[bestPair].second) };
	// we have all it takes to fill shower-level comparisons, let's start

	auto shr0 = rcshr_v.at( 0 );
	auto shr1 = rcshr_v.at( 1 );

	_rce0    = shr0.Energy();
	_rce1    = shr1.Energy();

	if (_containmentcorrection) {
	  _rce0 *= ContainmentCorr(
	  

	::geoalgo::Point_t rc0strt(shr0.ShowerStart().X(), shr0.ShowerStart().Y(), shr0.ShowerStart().Z() );
	::geoalgo::Point_t rc1strt(shr1.ShowerStart().X(), shr1.ShowerStart().Y(), shr1.ShowerStart().Z() );

	::geoalgo::Point_t rc0dir(shr0.Direction().X(), shr0.Direction().Y(), shr0.Direction().Z() );
	::geoalgo::Point_t rc1dir(shr1.Direction().X(), shr1.Direction().Y(), shr1.Direction().Z() );

	rc0dir.Normalize();
	rc1dir.Normalize();

	auto pi0mom = (rc0dir * _rce0) + (rc1dir * _rce1);
	pi0mom.Normalize();
	_pi0px = pi0mom[0];
	_pi0py = pi0mom[1];
	_pi0pz = pi0mom[2];

	_rc0x = rc0strt[0];
	_rc0y = rc0strt[1];
	_rc0z = rc0strt[2];
	
	_rc1x = rc1strt[0];
	_rc1y = rc1strt[1];
	_rc1z = rc1strt[2];
	
	// did we find a pi0 in this event? if so do RC <-> MC matching and fill additional TTree information
	if ( pi0_gamma_v.size() == 2 ) {
	  
	  auto MCRCmatch = MCRCMatch( rcshr_v, pi0_gamma_v);

	  //auto shr0 = rcshr_v.at( MCRCmatch[0].first );
	  //auto shr1 = rcshr_v.at( MCRCmatch[1].first );
	  auto shr0 = rcshr_v.at( MCRCmatch[0].first );
	  auto shr1 = rcshr_v.at( MCRCmatch[1].first );
	  auto mcs0 = pi0_gamma_v.at( MCRCmatch[0].second );
	  auto mcs1 = pi0_gamma_v.at( MCRCmatch[1].second );
	  _mce0    = mcs0.Start().E();
	  _mce1    = mcs1.Start().E();
	  _mcedep0 = mcs0.DetProfile().E();
	  _mcedep1 = mcs1.DetProfile().E();
	  _pi0a    = mcs0.Start().Momentum().Vect().Angle( mcs1.Start().Momentum().Vect() );
	  
	  ::geoalgo::Point_t mc0strt(mcs0.DetProfile().X(), mcs0.DetProfile().Y(), mcs0.DetProfile().Z() );
	  ::geoalgo::Point_t mc1strt(mcs1.DetProfile().X(), mcs1.DetProfile().Y(), mcs1.DetProfile().Z() );
	  
	  ::geoalgo::Point_t mc0dir(mcs0.DetProfile().Px(), mcs0.DetProfile().Py(), mcs0.DetProfile().Pz() );
	  ::geoalgo::Point_t mc1dir(mcs1.DetProfile().Px(), mcs1.DetProfile().Py(), mcs1.DetProfile().Pz() );
	  
	  mc0dir.Normalize();
	  mc1dir.Normalize();
	
	  _angle0 = rc0dir.Angle( mc0dir );
	  _angle1 = rc1dir.Angle( mc1dir );
	  
	  _mc0x = mc0strt[0];
	  _mc0y = mc0strt[1];
	  _mc0z = mc0strt[2];
	  
	  _mc1x = mc1strt[0];
	  _mc1y = mc1strt[1];
	  _mc1z = mc1strt[2];
	  
	  _d0 = rc0strt.Dist( mc0strt );
	  _d1 = rc1strt.Dist( mc1strt );
	  
	  _ip    = bestIP;	
	  _angle = 180. * rc0dir.Angle( rc1dir ) / 3.14;
	  _mass  = sqrt( 2 * _rce0 * _rce1 * ( 1 - cos( rc0dir.Angle( rc1dir ) ) ) );
	  
	  // fill shower-by-shower tree
	  _rce = _rce0;
	  _mce = _mce0;
	  _mcedep = _mcedep0;
	  _anglediff = _angle0;
	  _rcx = _rc0x;
	  _rcy = _rc0y;
	  _rcz = _rc0z;
	  _rcpx = rc0dir[0];
	  _rcpy = rc0dir[1];
	  _rcpz = rc0dir[2];
	  auto HL  = geoalgo::HalfLine( _rcx, _rcy, _rcz, _rcpx, _rcpy, _rcpz);
	  auto pts = _geoAlgo.Intersection(HL,_TPC);
	  _dwall = 1036.;
	  if (pts.size() == 1)
	    _dwall = pts[0].Dist( HL.Start() );
	  
	  _shower_tree->Fill();
	  
	  _rce = _rce1;
	  _mce = _mce1;
	  _mcedep = _mcedep1;
	  _anglediff = _angle1;
	  _rcx = _rc1x;
	  _rcy = _rc1y;
	  _rcz = _rc1z;
	  _rcpx = rc1dir[0];
	  _rcpy = rc1dir[1];
	  _rcpz = rc1dir[2];
	  HL  = geoalgo::HalfLine( _rcx, _rcy, _rcz, _rcpx, _rcpy, _rcpz);
	  pts = _geoAlgo.Intersection(HL,_TPC);
	  _dwall = 1036.;
	  if (pts.size() == 1)
	    _dwall = pts[0].Dist( HL.Start() );
	  
	  _shower_tree->Fill();

	}// if we have MC info available for matching
	  
      }// if a pi0 was reconstructed

    }// if a shower-pair was found

    _tree->Fill();

    return true;
  }

  bool Pi0Selection::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();
    if (_pi0_tree) _pi0_tree->Write();
    if (_shower_tree) _shower_tree->Write();

    return true;
  }

  const std::vector<size_t> Pi0Selection::FilterShowers(const larlite::event_shower* shr_v) {

    // list of event_shower indices which pass selection
    std::vector<size_t> filtered_shower_idx_v;

    for (size_t i=0; i < shr_v->size(); i++) {

      auto const& shr = shr_v->at(i);

      if ( shr.Energy() < _emin) continue;

      ::geoalgo::Point_t strt(shr.ShowerStart().X(),shr.ShowerStart().Y(),shr.ShowerStart().Z());

      if ( strt.Dist( _rcvtx ) > _radlenmax ) continue;

      filtered_shower_idx_v.push_back( i );

    }// for all showers

    return filtered_shower_idx_v;
  }

  const std::vector< std::pair<size_t,size_t> > Pi0Selection::Combinatorics(const std::vector<size_t> idx_v) {

    std::vector< std::pair<size_t, size_t> > shr_pairs;

    for (size_t i=0; i < idx_v.size(); i++) {
      for (size_t j=i+1; j < idx_v.size(); j++) {
	shr_pairs.push_back( std::make_pair(i,j) );
      }
    }

    return shr_pairs;
  }

  // pair returned is [ rc, mc ]
  std::vector< std::pair<size_t,size_t> >  Pi0Selection::MCRCMatch(const std::vector<larlite::shower>& shr_v,
								   const std::vector<larlite::mcshower>& mcs_v) {
    
    
    std::vector< std::pair<size_t,size_t> > matched_indices;

    // match the two reco showers to the two true showers according to their direction
    // best match is determined by which pairing leads to the smallest sum of 3D opening angles.

    double angle00 = ( 180. / 3.14 ) * shr_v[0].Direction().Angle( mcs_v[0].Start().Momentum().Vect() );
    double angle01 = ( 180. / 3.14 ) * shr_v[0].Direction().Angle( mcs_v[1].Start().Momentum().Vect() );
    double angle11 = ( 180. / 3.14 ) * shr_v[1].Direction().Angle( mcs_v[1].Start().Momentum().Vect() );
    double angle10 = ( 180. / 3.14 ) * shr_v[1].Direction().Angle( mcs_v[0].Start().Momentum().Vect() );

    if ( (angle00 + angle11) < (angle01 + angle10) ) {
      matched_indices.push_back( std::make_pair(0,0) );
      matched_indices.push_back( std::make_pair(1,1) );
    }
    else {
      matched_indices.push_back( std::make_pair(0,1) );
      matched_indices.push_back( std::make_pair(1,0) );
    }
    
    return matched_indices;
  }

  double Pi0Selection::ContainmentCorr(double d) {

    double frac = -1.039E-08 * pow(d,4) + 4.34E-06 * pow(d,3) - 0.000647 * pow(d,2) + 0.04224 * d - 0.04999;

    return 1./frac;
  }

}
#endif
