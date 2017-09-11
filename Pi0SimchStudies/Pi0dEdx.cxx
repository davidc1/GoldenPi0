#ifndef LARLITE_PI0DEDX_CXX
#define LARLITE_PI0DEDX_CXX

#include "Pi0dEdx.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/DetectorProperties.h"

#include "DataFormat/mcpart.h"
#include "DataFormat/simch.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool Pi0dEdx::initialize() {

    auto geomH   = ::larutil::GeometryHelper::GetME();
    _t2cm = geomH->TimeToCm();
    _w2cm = geomH->WireToCm();

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");

    _tree->Branch("_event",&_event,"event/I");

    _tree->Branch("_etrue0",&_etrue0,"etrue0/D");
    _tree->Branch("_edep0", &_edep0, "edep0/D" );
    _tree->Branch("_qcol0", &_qcol0, "qcol0/D" );
    _tree->Branch("_ehit0", &_ehit0, "ehit0/D" );
    _tree->Branch("_qhit0", &_qhit0, "qhit0/D" );
    _tree->Branch("_mcsdedx0", &_mcsdedx0, "mcsdedx0/D" );
    _tree->Branch("_truededx0", &_truededx0, "truededx0/D" );
    _tree->Branch("_dedx0", &_dedx0, "dedx0/D" );
    _tree->Branch("_e00", &_e00, "e00/D" );
    _tree->Branch("_e01", &_e01, "e01/D" );
    _tree->Branch("_a0", &_a0, "a0/D" );
    _tree->Branch("_truededx0_v","std::vector<double>" , &_truededx0_v );
    _tree->Branch("_recodedx0_v","std::vector<double>" , &_recodedx0_v );
    _tree->Branch("_adcqdedx0_v","std::vector<double>" , &_adcqdedx0_v );
    _tree->Branch("_pitch0", &_pitch0, "pitch0/D" );
    _tree->Branch("_process0", &_process0);

    _tree->Branch("_etrue1",&_etrue1,"etrue1/D");
    _tree->Branch("_edep1", &_edep1, "edep1/D" );
    _tree->Branch("_qcol1", &_qcol1, "qcol1/D" );
    _tree->Branch("_ehit1", &_ehit1, "ehit1/D" );
    _tree->Branch("_qhit1", &_qhit1, "qhit1/D" );
    _tree->Branch("_mcsdedx1", &_mcsdedx1, "mcsdedx1/D" );
    _tree->Branch("_truededx1", &_truededx1, "truededx1/D" );
    _tree->Branch("_dedx1", &_dedx1, "dedx1/D" );
    _tree->Branch("_e10", &_e10, "e10/D" );
    _tree->Branch("_e11", &_e11, "e11/D" );
    _tree->Branch("_a1", &_a1, "a1/D" );
    _tree->Branch("_truededx1_v","std::vector<double>" , &_truededx1_v );
    _tree->Branch("_recodedx1_v","std::vector<double>" , &_recodedx1_v );
    _tree->Branch("_adcqdedx1_v","std::vector<double>" , &_adcqdedx1_v );
    _tree->Branch("_pitch1", &_pitch1, "pitch1/D" );
    _tree->Branch("_process1", &_process1);

    _tree->Branch("_wmin0",&_wmin0,"wmin0/I");
    _tree->Branch("_wmax0",&_wmax0,"wmax0/I");
    _tree->Branch("_wmin1",&_wmin1,"wmin1/I");
    _tree->Branch("_wmax1",&_wmax1,"wmax1/I");

    _tree->Branch("_angle", &_angle, "angle/D" );

    return true;
  }
  
  bool Pi0dEdx::analyze(storage_manager* storage) {

    auto geom    = ::larutil::Geometry::GetME();
    auto geomH   = ::larutil::GeometryHelper::GetME();
    auto detProp = ::larutil::DetectorProperties::GetME();

    _event = storage->event_id();

    _chTickMap.clear();
    _chTickMap = std::vector< std::vector< std::pair<int,int> > >(8500,std::vector<std::pair<int,int> >());

    _chIDEmap.clear();
    _chIDEmap = std::vector< std::map<int,double> >(8500,std::map<int,double>());
    
    // Retrieve mcshower data product
    auto ev_mcs   = storage->get_data<event_mcshower>("mcreco");
    // Retrieve hit data product
    auto ev_hit   = storage->get_data<event_hit>("gaushit");
    // Retrieve simch data product
    auto ev_simch = storage->get_data<event_simch>("largeant");
    // Retrieve mcparticle data product
    auto ev_mcp = storage->get_data<event_mcpart>("largeant");

    if (!ev_simch){
      std::cout << "No simch data-product -> exit" << std::endl;
      return false;
    }

    if (!ev_hit){
      std::cout << "No hit data-product -> exit" << std::endl;
      return false;
    }

    // used track id vector
    std::vector<unsigned int> used_trk_id;

    // map connecting 1st vs. 2nd shower filled and mcindex used
    std::map<size_t,size_t> shrmapidx;

    // Create G4 track ID vector for MCShowers & MCTracks in the event
    std::vector<std::vector<unsigned int> > g4_trackid_v;
    std::vector<unsigned int> mc_index_v;
    // keep track of all the shower and track trackIDs
    g4_trackid_v.reserve(ev_mcs->size());
    // keep track of track ID of the two photons
    int trkid0, trkid1;

    // for each mcshower, find the geant trackIDs associated with that shower
    for (size_t mc_index = 0; mc_index < ev_mcs->size(); ++mc_index) {
      auto const& mcs = (*ev_mcs)[mc_index];
      
      // only keep showers form pi0 decay
      if (mcs.MotherPdgCode() != 111) continue;


      std::vector<unsigned int> id_v;
      id_v.reserve(mcs.DaughterTrackID().size());

      for (auto const& id : mcs.DaughterTrackID()) {
	if (id == mcs.TrackID()) continue;
	used_trk_id.push_back(id);
	if (id == 1)
	  std::cout << "ENTERING TrackID == 1!" << std::endl;
	id_v.push_back(id);
      }

      id_v.push_back(mcs.TrackID());

      used_trk_id.push_back(mcs.TrackID());
      g4_trackid_v.push_back(id_v);
      mc_index_v.push_back(mc_index);

      shrmapidx[mc_index] = mc_index_v.size() - 1;
      
      if (mc_index_v.size() == 1) {
	_etrue0 = mcs.Start().E();
	_edep0  = mcs.DetProfile().E();
	_qcol0  = mcs.Charge(2);
	trkid0  = mcs.TrackID();
	_mcsdedx0 = mcs.dEdx();
      }
      if (mc_index_v.size() == 2) {
	_etrue1 = mcs.Start().E();
	_edep1  = mcs.DetProfile().E();
	_qcol1  = mcs.Charge(2);
	trkid1  = mcs.TrackID();
	_mcsdedx1 = mcs.dEdx();
      }	
    }// for all mcshowers

    if (mc_index_v.size() != 2) {
      std::cout << "did not find 2 MCShowers... quit " << std::endl;
      return false;
    }

    _process0 = _process1 = "";

    // loop through MCParticles and find end process of the photons
    std::vector<int> daughters0, daughters1;
    for (size_t p=0; p < ev_mcp->size(); p++) {
      auto const& mcp = ev_mcp->at(p);
      if (mcp.TrackId() == trkid0){
	//_process0 = mcp.Process();
	//std::cout << "process0 @ trackid " << trkid0 << " is " << _process0 << std::endl;
	//std::cout << "This gamma has " << mcp.Daughters().size() << " children " << std::endl;
	for (auto const& d : mcp.Daughters())
	  daughters0.push_back(d);
      }
      if (mcp.TrackId() == trkid1){
	//_process1 = mcp.Process();
	//std::cout << "process1 @ trackid " << trkid1 << " is " << _process1 << std::endl;
	//std::cout << "This gamma has " << mcp.Daughters().size() << " children " << std::endl;
	for (auto const& d : mcp.Daughters())
	  daughters1.push_back(d);
      }
    }// for all MCParticles

    for (size_t p=0; p < ev_mcp->size(); p++) {
      auto const& mcp = ev_mcp->at(p);
      if (daughters0.size()){
	if (daughters0[0] == mcp.TrackId() )
	  _process0 = mcp.Process();
      }
      if (daughters1.size()){
	if (daughters1[0] == mcp.TrackId() ){
	  _process1 = mcp.Process();
	}
      }
    }

    // save energy and opening angle of pair-production output
    if (daughters0.size() == 2) {
      larlite::mcstep step0, step1;
      for (size_t p=0; p < ev_mcp->size(); p++) {
	auto const& mcp = ev_mcp->at(p);
	if (daughters0[0] == mcp.TrackId() ) 
	  step0 = mcp.Trajectory().at(0);
    	if (daughters0[1] == mcp.TrackId() ) 
	  step1 = mcp.Trajectory().at(0);	
      }// for all particles
      _e00 = step0.E();
      _e01 = step1.E();
      _a0  = step0.Momentum().Vect().Angle(step1.Momentum().Vect());
    }// if part0 has 2 daughters
    if (daughters1.size() == 2) {
      larlite::mcstep step0, step1;
      for (size_t p=0; p < ev_mcp->size(); p++) {
	auto const& mcp = ev_mcp->at(p);
    	if (daughters1[0] == mcp.TrackId() ) 
	  step0 = mcp.Trajectory().at(0);
    	if (daughters1[1] == mcp.TrackId() ) 
	  step1 = mcp.Trajectory().at(0);	
      }// for all particles
      _e10 = step0.E();
      _e11 = step1.E();
      _a1  = step0.Momentum().Vect().Angle(step1.Momentum().Vect());
    }// if part1 has 2 daughters

    // get start point and direction for both showers
    auto srt0 = ev_mcs->at(mc_index_v[0]).DetProfile().Position().Vect();
    auto dir0 = ev_mcs->at(mc_index_v[0]).Start().Momentum().Vect().Unit();

    auto srt1 = ev_mcs->at(mc_index_v[1]).DetProfile().Position().Vect();
    auto dir1 = ev_mcs->at(mc_index_v[1]).Start().Momentum().Vect().Unit();

    _truededx0_v.clear();
    _truededx1_v.clear();

    _truededx0_v = std::vector<double>(_dmax*3*2,0.);
    _truededx1_v = std::vector<double>(_dmax*3*2,0.);

    _recodedx0_v.clear();
    _recodedx1_v.clear();

    _adcqdedx0_v.clear();
    _adcqdedx1_v.clear();

    //_recodedx0_v = std::vector<double>(_dmax*3*2,0.);
    //_recodedx1_v = std::vector<double>(_dmax*3*2,0.);

    // loop through simchannels and calculate "true" dedx
    _truededx0 = _truededx1 = 0.;
    for (size_t e=0; e < ev_simch->size(); e++) {
      
      if (ev_simch->at(e).Channel() < 4800) continue;
      
      auto const& ide_v = ev_simch->at(e).TrackIDsAndEnergies(0,30000);

      for (auto const& ide : ide_v) {

	TVector3 ideR(ide.x,ide.y,ide.z);
	
	double d0 = (ideR-srt0).Mag();
	double d1 = (ideR-srt1).Mag();
	
	double dz0 = fabs(ide.z - srt0.Z());
	double dz1 = fabs(ide.z - srt1.Z());

	if (d0 < _dmax) { 
	  _truededx0 += ide.energy; 
	  if (dz0*3 > _truededx0_v.size() ) continue;
	  _truededx0_v[dz0*3] += ide.energy;
	}
	if (d1 < _dmax) { _truededx1 += ide.energy; }
	if (dz1*3 > _truededx1_v.size() ) continue;
	  _truededx1_v[dz1*3] += ide.energy;
      }
      
    }// for all IDEs

    _truededx0 /= _dmax;
    _truededx1 /= _dmax;

    //std::cout << "start0 : " << srt0.X() << ", " << srt0.Y() << ", " << srt0.Z() << std::endl;
    //std::cout << "start1 : " << srt1.X() << ", " << srt1.Y() << ", " << srt1.Z() << std::endl;

    // calculate trigger-time offset
    auto tvtx = ev_mcs->at(mc_index_v[1]).DetProfile().T();
    // time in tick :
    auto vtxtick = (tvtx / 1000.) * 2.;
    // time in cm :
    auto vtxtimecm = vtxtick * _t2cm; 

    // pitch calculation
    _pitch0 = geomH->GetPitch(dir0, 2);
    _pitch1 = geomH->GetPitch(dir1, 2);

    _angle = ev_mcs->at(mc_index_v[0]).Start().Momentum().Vect().Unit().Dot( ev_mcs->at(mc_index_v[1]).Start().Momentum().Vect().Unit() );

    // reset charge integrators
    _ehit0 = _qhit0 = _ehit1 = _qhit1 = 0;

    _wmin0 = _wmin1 = 9000;
    _wmax0 = _wmax1 = 0;

    // reset MCBTAlg
    _bt_algo.Reset(g4_trackid_v,*ev_simch);

    // create a vector of vectors in which to store the hit associations
    // for each cluster produced
    std::vector<std::vector<unsigned int> > cluster_hit_v;
    cluster_hit_v.resize( 3 * mc_index_v.size() );
    // vector to save the plane information for each cluster
    std::vector<larlite::geo::View_t> cluster_plane_v(cluster_hit_v.size(),
						      larlite::geo::View_t::kUnknown);

    // whether to use a hit for the dE/dx or not
    bool use4dEdx0 = false;
    bool use4dEdx1 = false;

    // number of hits used in dedx calculation for each shower
    int ndedx0 = 0;
    int ndedx1 = 0;
    _dedx0 = _dedx1 = 0.;

    // min hit distance to the two shower start points
    double dmin0 = 1e9;
    double dmin1 = 1e9;

    double xmin0 = 1e9;
    double xmin1 = 1e9;
    double zmin0 = 1e9;
    double zmin1 = 1e9;

    // loop through hits, use the back-tracker to find which MCX object
    // it should belong to, and add that hit to the cluster that is
    // indended for that MCX object
    // only use hits from association to rawclusters
    for (size_t i=0; i < ev_hit->size(); i++) {

      auto const& hit_idx = i;
      auto const& hit     = ev_hit->at(hit_idx);
      
      auto const& ch     = hit.Channel();
      auto const& tstart = hit.PeakTime() - 2*hit.RMS();// + 2255;//3050;
      auto const& tend   = hit.PeakTime() + 2*hit.RMS();// + 2255;//3050;
      auto const& pl     = hit.View();

      double* origin;
      origin = new double[3];
      geom->PlaneOriginVtx( pl, origin);
      float planeOffset = origin[0];

      double hitw = hit.WireID().Wire * _w2cm;
      double hitt = ( hit.PeakTime() - detProp->TriggerOffset() ) * _t2cm + planeOffset - vtxtimecm;

      double d2D0 = sqrt( pow((hitw - srt0.Z()),2) + pow((hitt - srt0.X()),2) );
      double d2D1 = sqrt( pow((hitw - srt1.Z()),2) + pow((hitt - srt1.X()),2) );
      double d3D0 = d2D0 / (1 - dir0.Y() * dir0.Y() );
      double d3D1 = d2D1 / (1 - dir1.Y() * dir1.Y() );

      //std::cout << "\t hit wire/time : " << hitw << ", " << hitt << " has d3D0 of " << d3D0 << std::endl;
      //std::cout << "\t hit wire/time : " << hitw << ", " << hitt << " has d3D1 of " << d3D1 << std::endl;

      if (d3D0 < dmin0) { dmin0 = d2D0; xmin0 = hitt; zmin0 = hitw; }
      if (d3D1 < dmin1) { dmin1 = d2D1; xmin1 = hitt; zmin1 = hitw; }
      
      // is hit within dE/dx region of either shower
      if ( d3D0 > _dmax ) use4dEdx0 = false;
      else use4dEdx0 = true;
      if ( d3D1 > _dmax ) use4dEdx1 = false;
      else use4dEdx1 = true;

      // avoid duplicating time-intervals
      std::pair<double,double> timeinterval = std::make_pair(tstart,tend);

      if (_avoid_duplicate_ticks)
	auto timeinterval = getTimeSubset(ch,(int)tstart,(int)tend);

      if (timeinterval.first >= timeinterval.second) {
	//std::cout << "skipping hit " << std::endl;
	continue;
      }

      _chTickMap[ch].push_back( timeinterval );

      // further avoid duplication
      std::vector< std::pair<int,int> > timeinterval_v;

      timeinterval_v.push_back( timeinterval );

      for (auto const& T : timeinterval_v) {
      
	// create a wire-range object with channel + (start,end) time info for the hit
	::btutil::WireRange_t wr(ch,T.first,T.second);
	// use the back-tracker to return a vector of track IDs associated to this hit
	// contents of return vector are proportional to the fraction of the hit's
	// charge that is associated with the MCX id at that element's position
	// vector ordered such that
	// mcq_v [ index ] corresponds to index wihtin mc_index_v
	auto mcq_v = _bt_algo.MCQ(wr);
	// grab the same information but for the deposied energy
	auto mce_v = _bt_algo.MCE(wr);
	// find entry with largest value -> this lets us know which MCX object
	// to assign this hit to and thus which cluster this hit should
	// be associated with
	size_t idx  = 0;
	double max_edep = 0;
	for (size_t j=0; j < mce_v.size(); j++){
	  if (mce_v[j] > max_edep){
	    max_edep = mce_v[j];
	    idx  = j;
	  }
	}
	double max_qdep = mcq_v[idx];
	//std::cout << "Edep for this hit : " << max_edep << " associated with idx " << idx << std::endl;
	// if the maximum amount of charge is 0
	// ignore this hit
	if (max_edep == 0)
	  continue;
	// if the idx found is == to mcq_v.size() - 1
	// this means that most of the charge belongs to 
	// none of the MCX objects we are interested in
	// -> ignore this hit
	if (idx == mcq_v.size() -1 )
	  continue;
	
	if (pl == 2) {
	  if (shrmapidx[idx] == 0) {
	    _ehit0 += max_edep;
	    _qhit0 += max_qdep;
	    _ahit0 += hit.Integral();
	    if (use4dEdx0) {
	      //std::cout << "adding dedx @ shr 0" << std::endl;
	      _dedx0 += max_qdep / _pitch0;
	      ndedx0 += 1;
	      _recodedx0_v.push_back( max_qdep );
	      _adcqdedx0_v.push_back( hit.Integral() );
	    }
	    if (ch < _wmin0) _wmin0 = ch;
	    if (ch > _wmax0) _wmax0 = ch;
	  }
	  if (shrmapidx[idx] == 1) {
	    _ehit1 += max_edep;
	    _qhit1 += max_qdep;
	    _ahit1 += hit.Integral();
	    if (use4dEdx1) {
	      //std::cout << "adding dedx @ shr 1" << std::endl;
	      _dedx1 += max_qdep / _pitch1;
	      ndedx1 += 1;
	      _recodedx1_v.push_back( max_qdep );
	      _adcqdedx1_v.push_back( hit.Integral() );
	    }
	    if (ch < _wmin1) _wmin1 = ch;
	    if (ch > _wmax1) _wmax1 = ch;
	  }

	}// if on collection plane

	// if not, associate this hit with the appropriate cluster
	// ( "i" is the hit index )
	cluster_hit_v[ pl * mc_index_v.size() + idx ].push_back( hit_idx );
	cluster_plane_v[ pl * mc_index_v.size() + idx ] = pl;
	
      }// for all found time-intervals

    }// for all hits in cluster

    _dedx0 /= ndedx0;
    _dedx1 /= ndedx1;

    //std::cout << "dmin0 = " << dmin0 << "\t xmin0 = " << xmin0 << "\t zmin0 = " << zmin0 << std::endl;
    //std::cout << "dmin1 = " << dmin1 << "\t xmin1 = " << xmin1 << "\t zmin1 = " << zmin1 << std::endl;
    
    _tree->Fill();

    if (!_save_clusters) return true;

    // Retrieve cluster data product (output)
    auto ev_mccluster = storage->get_data<event_cluster>("mccluster");
    // association information to go here
    auto cluster_ass_v = storage->get_data<event_ass>(ev_mccluster->name());
    
    // since we are creating a new data product,
    // reload the event information
    storage->set_id(ev_hit->run(),ev_hit->subrun(),ev_hit->event_id()); 
    
    // now create the clusters and save them as data-products
    
    // create a cluster for each entry in g4_trackid_v -> the vector
    // of MCShower / MCTrack trackIDs that we want to keep
    // also create a new vector for cluster associations
    // to be filled only if the cluster for the corresponding
    // hit list will be created after meeting the minimum
    // num of hit requirement
    std::vector<std::vector<unsigned int> > cluster_hit_ass_v;
    for (size_t idx=0; idx < cluster_hit_v.size(); idx++){

      //std::cout << "saving cluster w/ "  << cluster_hit_v[idx].size() << std::endl;
      
      // create a new cluster
      cluster clus;
      clus.set_n_hits(cluster_hit_v[idx].size());
      clus.set_view(cluster_plane_v[idx]);
      ev_mccluster->push_back(clus);
      cluster_hit_ass_v.push_back(cluster_hit_v[idx]);

    }// for all clusters created

    // now save the associations for the cluster
    cluster_ass_v->set_association(ev_mccluster->id(),product_id(data::kHit,ev_hit->name()),
				   cluster_hit_ass_v);

  
    return true;
  }

  bool Pi0dEdx::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();

    return true;
  }


  std::pair<int,int> Pi0dEdx::getTimeSubset(const int& ch,const int& tstart,const int& tend) {

    if (_chTickMap[ ch ].size() == 0) {
      std::pair<int,int> tpair = std::make_pair(tstart,tend);
      return tpair;
    }

    int startunique = tstart;
    int endunique   = tend;

    //std::cout << "new interval @ channel "  << ch << " [ " << tstart << ", " << tend << " ]" << std::endl;

    // found, let's eliminate time-intervals
    for (auto const& tpair : _chTickMap[ch] ) {

      if (tpair.first > tend) continue;

      if (tpair.second < tstart) continue;
      
      // if interval fully contained in already scanned interval
      // return nothing
      if ( (tstart > tpair.first) && (tend < tpair.second) )
	return std::make_pair(0,0);

      if ( (tstart > tpair.first) && (tstart < tpair.second) ) {
	// only make interval smaller
	if (tpair.second > startunique)
	  startunique = tpair.second;
      }
      
      if ( (tend > tpair.first) && (tend < tpair.second) ) {
	// only make interval smaller
	if (tpair.first < endunique)
	  endunique = tpair.first;
      }

      // does the new it completely surround the old interval?
      //if ( (tend >= tpair.second) && (tstart <= tpair.first) )
      //std::cout << "Fully surrounded!" << std::endl;
    }// for all time intervals already saved
    
    //if ( (startunique != tstart) || (endunique != tend) )
    //std::cout << "time-interval @ channel " << ch << " changed from [ " << tstart << ", " << tend << " ]" 
    //<< " to [ " << startunique << ", " << endunique << " ]" << std::endl;
    
    return std::make_pair(startunique,endunique);
    
  }
  
}
#endif
