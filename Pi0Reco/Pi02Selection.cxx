#ifndef LARLITE_PI02SELECTION_CXX
#define LARLITE_PI02SELECTION_CXX

#include "Pi02Selection.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  Pi02Selection::Pi02Selection()
    : _tree(nullptr)
    , _pi0_tree(nullptr)
    , _shower_tree(nullptr)
  { 
    _name="Pi02Selection";
    _fout=0;
    _edepmin = 0;
    _containmentcorrection = false;
  }

  bool Pi02Selection::initialize() {

    _ctr = -1;

    _TPC = ::geoalgo::AABox_t(0,-116,0,256,116,1036);

    if (_pi0_tree) delete _pi0_tree;

    _pi0_tree = new TTree("_pi0_tree","pi0 tree");

    // energies
    _pi0_tree->Branch("el",&_el,"el/D");
    _pi0_tree->Branch("eh",&_eh,"eh/D");
    _pi0_tree->Branch("elc",&_elc,"elc/D");
    _pi0_tree->Branch("ehc",&_ehc,"ehc/D");
    
    // conversion disntaces
    _pi0_tree->Branch("rl",&_rl,"rl/D");
    _pi0_tree->Branch("rh",&_rh,"rh/D");

    // start point
    _pi0_tree->Branch("hxs",&_hxs,"hxs/D");
    _pi0_tree->Branch("hys",&_hys,"hys/D");
    _pi0_tree->Branch("hzs",&_hzs,"hzs/D");
    _pi0_tree->Branch("lxs",&_lxs,"lxs/D");
    _pi0_tree->Branch("lys",&_lys,"lys/D");
    _pi0_tree->Branch("lzs",&_lzs,"lzs/D");

    // end point
    _pi0_tree->Branch("hxe",&_hxe,"hxe/D");
    _pi0_tree->Branch("hye",&_hye,"hye/D");
    _pi0_tree->Branch("hze",&_hze,"hze/D");
    _pi0_tree->Branch("lxe",&_lxe,"lxe/D");
    _pi0_tree->Branch("lye",&_lye,"lye/D");
    _pi0_tree->Branch("lze",&_lze,"lze/D");

    // dEdx
    _pi0_tree->Branch("dedxh",&_dedxh,"dedxh/D");
    _pi0_tree->Branch("dedxl",&_dedxl,"dedxl/D");

    // d Wall
    _pi0_tree->Branch("dwallh",&_dwallh,"dwallh/D");
    _pi0_tree->Branch("dwalll",&_dwalll,"dwalll/D");

    // pitch
    _pi0_tree->Branch("pitchh",&_pitchh,"pitchh/D");
    _pi0_tree->Branch("pitchl",&_pitchl,"pitchl/D");

    // pi0 reconstructed properties
    _pi0_tree->Branch("angle",&_angle,"angle/D");
    _pi0_tree->Branch("mass",&_mass,"mass/D");
    _pi0_tree->Branch("massc",&_massc,"massc/D");
    _pi0_tree->Branch("ip",&_ip,"ip/D");
    
    // event info
    _pi0_tree->Branch("run",&_run,"run/I");
    _pi0_tree->Branch("sub",&_sub,"sub/I");
    _pi0_tree->Branch("evt",&_evt,"evt/I");
    _pi0_tree->Branch("ctr",&_ctr,"ctr/I");
    
    return true;
  }
  
  bool Pi02Selection::analyze(storage_manager* storage) {

    _npi0     = 0;
    _ngamma   = 0;
    _mass     = -1;

    _ctr += 1;

    std::vector<larlite::mcshower> pi0_gamma_v;

    auto ev_shr = storage->get_data<event_shower>(_shr_producer);
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
      std::vector<double> mass_v;
      
      for (size_t pidx = 0; pidx < shr_pairs_idx_v.size(); pidx++ ) {
	
	auto const& shr_pair = shr_pairs_idx_v.at(pidx);
	
	auto const& shr1 = ev_shr->at(shr_pair.first);
	auto const& shr2 = ev_shr->at(shr_pair.second);
	
	auto const& srt1 = shr1.ShowerStart();
	auto const& dir1 = shr1.Direction();
	auto const& HL1  = geoalgo::HalfLine( srt1.X(), srt1.Y(), srt1.Z(),
					  -dir1.X(),-dir1.Y(),-dir1.Z());
	auto const& srt2 = shr2.ShowerStart();
	auto const& dir2 = shr2.Direction();
	auto const& HL2  = geoalgo::HalfLine( srt2.X(), srt2.Y(), srt2.Z(),
					  -dir2.X(),-dir2.Y(),-dir2.Z());
	

	::geoalgo::Point_t pt1, pt2;
	_ip = sqrt( _geoAlgo.SqDist(HL1,HL2,pt1,pt2) );
	_angle = HL1.Dir().Angle( HL2.Dir() );
	_ipvtx = _rcvtx.Dist( ((pt1+pt2)/2.) );

	// save mass value regardless
	_mass  = sqrt( 2 * shr1.Energy() * shr2.Energy() * ( 1 - cos(_angle) ) );

	mass_v.push_back(_mass);
	
	//_pi0_tree->Fill();
	
      }// for all shower-pairs

      std::cout << "pairs : " << mass_v.size() << std::endl;

      // if one pairs:
      if (mass_v.size() == 1) {
	auto const& shr_pair = shr_pairs_idx_v.at(0);
	FillTree(ev_shr->at(shr_pair.first),ev_shr->at(shr_pair.second));
      }// if 1 pair

      
      // if three pairs:
      if (mass_v.size() == 3) {
	
	// find one closest to pi0 mass:
	double dmin = 1000.;
	int   idx  = 0;
	for (size_t i=0; i < mass_v.size(); i++) {
	  auto m = mass_v[i];
	  if ( fabs(m - 134.*0.9) < dmin ) { dmin = fabs(m - 134.*0.9); idx = i; }
	}// for all pairs
	
	auto const& shr_pair = shr_pairs_idx_v.at(idx);
	FillTree(ev_shr->at(shr_pair.first),ev_shr->at(shr_pair.second));
	
      }// if 3 pairs

      // if 6 pairs:
      if (mass_v.size() == 6) {
	
	// 0 : (0,5) -> (0,1) & (2,3)
	// 1 : (1,4) -> (0,2) & (1,3)
	// 2 : (2,3) -> (0,3) & (1,2)

    	double dmin = 1000.;
	int   idx  =  0;
	double mm0 =  fabs(mass_v[0] - 134.*0.9) + fabs(mass_v[5] - 134.*0.9);
	double mm1 =  fabs(mass_v[1] - 134.*0.9) + fabs(mass_v[4] - 134.*0.9);
	double mm2 =  fabs(mass_v[2] - 134.*0.9) + fabs(mass_v[3] - 134.*0.9);

	std::vector<int> pair0 = {0,0};
	std::vector<int> pair1 = {0,0};
	
	if ( (mm0 < mm1) && (mm0 < mm2) ) {
	  pair0 = {0,1};
	  pair1 = {2,3};
	}

	if ( (mm1 < mm0) && (mm1 < mm2) ) {
	  pair0 = {0,2};
	  pair1 = {1,3};
	}

	if ( (mm2 < mm1) && (mm2 < mm0) ) {
	  pair0 = {0,3};
	  pair1 = {1,2};
	}

	FillTree(ev_shr->at(pair0[0]),ev_shr->at(pair0[1]));
	FillTree(ev_shr->at(pair1[0]),ev_shr->at(pair1[1]));
	
      }// if 6 pairs

    }// if at least 2 showers
      
      return true;
  }
  
  bool Pi02Selection::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();
    if (_pi0_tree) _pi0_tree->Write();
    if (_shower_tree) _shower_tree->Write();

    return true;
  }

  const std::vector<size_t> Pi02Selection::FilterShowers(const larlite::event_shower* shr_v) {

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


  void Pi02Selection::FillTree(const larlite::shower& shr1, const larlite::shower& shr2) {

    // transform into half-line objects
    auto const& srt1 = shr1.ShowerStart();
    auto const& dir1 = shr1.Direction();
    auto const& HL1  = geoalgo::HalfLine( srt1.X(), srt1.Y(), srt1.Z(),
					  -dir1.X(),-dir1.Y(),-dir1.Z());
    auto const& HL1C = geoalgo::HalfLine( srt1.X(), srt1.Y(), srt1.Z(),
					  dir1.X(),dir1.Y(),dir1.Z());
    auto const& srt2 = shr2.ShowerStart();
    auto const& dir2 = shr2.Direction();
    auto const& HL2  = geoalgo::HalfLine( srt2.X(), srt2.Y(), srt2.Z(),
					  -dir2.X(),-dir2.Y(),-dir2.Z());
    auto const& HL2C = geoalgo::HalfLine( srt2.X(), srt2.Y(), srt2.Z(),
					  dir2.X(),dir2.Y(),dir2.Z());

    geoalgo::Point_t end1 = srt1 + dir1 * shr1.Length();
    geoalgo::Point_t end2 = srt2 + dir2 * shr2.Length();
    
    auto pts1 = _geoAlgo.Intersection(HL1C,_TPC);
    auto pts2 = _geoAlgo.Intersection(HL2C,_TPC);
    
    auto geomH = larutil::GeometryHelper::GetME();
    
    _dwallh = 1036.;
    _dwalll = 1036.;

    if ( shr1.Energy() < shr2.Energy() ) {
      _el = shr1.Energy();
      _eh = shr2.Energy();
      _rl = HL1.Start().Dist( _rcvtx );
      _rh = HL2.Start().Dist( _rcvtx );
      _dedxl = shr1.dEdx_v().at(2);
      _dedxh = shr2.dEdx_v().at(2);
      if (pts1.size() == 1)
	_dwalll = pts1[0].Dist( HL1C.Start() );
      if (pts2.size() == 1)
	_dwallh = pts2[0].Dist( HL2C.Start() );
      _lxs = srt1[0];
      _lys = srt1[1];
      _lzs = srt1[2];
      _hxs = srt2[0];
      _hys = srt2[1];
      _hzs = srt2[2];
      _lxe = end1[0];
      _lye = end1[1];
      _lze = end1[2];
      _hxe = end2[0];
      _hye = end2[1];
      _hze = end2[2];
      _pitchh = geomH->GetPitch( shr2.Direction(), 2);
      _pitchl = geomH->GetPitch( shr1.Direction(), 2);
    }
    else {
      _el = shr2.Energy();
      _eh = shr1.Energy();
      _rl = HL2.Start().Dist( _rcvtx );
      _rh = HL1.Start().Dist( _rcvtx );
      _dedxl = shr2.dEdx_v().at(2);
      _dedxh = shr1.dEdx_v().at(2);
      if (pts2.size() == 1)
	_dwalll = pts2[0].Dist( HL2C.Start() );
      if (pts1.size() == 1)
	_dwallh = pts1[0].Dist( HL1C.Start() );
      _lxs = srt2[0];
      _lys = srt2[1];
      _lzs = srt2[2];
      _hxs = srt1[0];
      _hys = srt1[1];
      _hzs = srt1[2];
      _lxe = end2[0];
      _lye = end2[1];
      _lze = end2[2];
      _hxe = end1[0];
      _hye = end1[1];
      _hze = end1[2];
      _pitchh = geomH->GetPitch( shr1.Direction(), 2);
      _pitchl = geomH->GetPitch( shr2.Direction(), 2);
    }
    
    ::geoalgo::Point_t pt1, pt2;
    _ip = sqrt( _geoAlgo.SqDist(HL1,HL2,pt1,pt2) );
    _angle = HL1.Dir().Angle( HL2.Dir() );
    _ipvtx = _rcvtx.Dist( ((pt1+pt2)/2.) );
    _mass  = sqrt( 2 * shr1.Energy() * shr2.Energy() * ( 1 - cos(_angle) ) );
    _ehc = _eh * ContainmentCorr(_dwallh);
    _elc = _el * ContainmentCorr(_dwalll);
    _massc = _mass * sqrt( ContainmentCorr(_dwallh) * ContainmentCorr(_dwalll) );
    
    _pi0_tree->Fill();
    
    return;
  }

  const std::vector< std::pair<size_t,size_t> > Pi02Selection::Combinatorics(const std::vector<size_t> idx_v) {

    std::vector< std::pair<size_t, size_t> > shr_pairs;

    for (size_t i=0; i < idx_v.size(); i++) {
      for (size_t j=i+1; j < idx_v.size(); j++) {
	shr_pairs.push_back( std::make_pair(i,j) );
      }
    }


    return shr_pairs;
  }

  double Pi02Selection::ContainmentCorr(double d) {

    if (d > 100) return 1;

    double frac = -1.039E-08 * pow(d,4) + 4.34E-06 * pow(d,3) - 0.000647 * pow(d,2) + 0.04224 * d - 0.04999;

    return 1./frac;
  }

}
#endif
