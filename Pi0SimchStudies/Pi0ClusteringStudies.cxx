#ifndef LARLITE_PI0CLUSTERINGSTUDIES_CXX
#define LARLITE_PI0CLUSTERINGSTUDIES_CXX

#include "Pi0ClusteringStudies.h"

#include "DataFormat/simch.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool Pi0ClusteringStudies::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");

    _tree->Branch("_etrue0",&_etrue0,"etrue0/D");
    _tree->Branch("_edep0", &_edep0, "edep0/D" );
    _tree->Branch("_qcol0", &_qcol0, "qcol0/D" );
    _tree->Branch("_eclus0", &_eclus0, "eclus0/D" );
    _tree->Branch("_qclus0", &_qclus0, "qclus0/D" );
    _tree->Branch("_ehit0", &_ehit0, "ehit0/D" );
    _tree->Branch("_qhit0", &_qhit0, "qhit0/D" );

    _tree->Branch("_etrue1",&_etrue1,"etrue1/D");
    _tree->Branch("_edep1", &_edep1, "edep1/D" );
    _tree->Branch("_qcol1", &_qcol1, "qcol1/D" );
    _tree->Branch("_eclus1", &_eclus1, "eclus1/D" );
    _tree->Branch("_qclus1", &_qclus1, "qclus1/D" );
    _tree->Branch("_ehit1", &_ehit1, "ehit1/D" );
    _tree->Branch("_qhit1", &_qhit1, "qhit1/D" );

    _tree->Branch("_wmin0",&_wmin0,"wmin0/I");
    _tree->Branch("_wmax0",&_wmax0,"wmax0/I");
    _tree->Branch("_wmin1",&_wmin1,"wmin1/I");
    _tree->Branch("_wmax1",&_wmax1,"wmax1/I");

    _tree->Branch("_angle", &_angle, "angle/D" );

    if (_hit_tree) delete _hit_tree;
    _hit_tree = new TTree("hit_tree","hit tree");
    _hit_tree->Branch("_ch",&_ch,"ch/I");
    _hit_tree->Branch("_q", &_q, "q/D" );
    _hit_tree->Branch("_eshr", &_eshr, "eshr/D" );
    _hit_tree->Branch("_qshr", &_qshr, "qshr/D" );
    _hit_tree->Branch("_eall", &_eall, "eall/D" );
    _hit_tree->Branch("_qall", &_qall, "qall/D" );
    _hit_tree->Branch("_adc", &_adc, "adc/D" );

    return true;
  }
  
  bool Pi0ClusteringStudies::analyze(storage_manager* storage) {


    
    // Retrieve mcshower data product
    auto ev_mcs   = storage->get_data<event_mcshower>("mcreco");
    // Retrieve simch data product
    auto ev_simch = storage->get_data<event_simch>("largeant");

    if (!ev_simch){
      std::cout << "No simch data-product -> exit" << std::endl;
      return false;
    }

    // grab showers and clusters associated to them
    auto ev_shr                     = storage->get_data<event_shower>("showerreco");
    larlite::event_pfpart*  ev_pfp  = nullptr;
    larlite::event_cluster* ev_clus = nullptr;
    larlite::event_hit*     ev_hit  = nullptr;

    if (!ev_shr) {
      std::cout << "No shower data-product -> exit" << std::endl;
      return false;
    }    

    auto const& ass_shr_pfp_v = storage->find_one_ass(ev_shr->id(), ev_pfp, ev_shr->name());
    
    if (!ev_pfp) {
      std::cout << "No pfp data-product -> exit" << std::endl;
      return false;
    }

    auto const& ass_pfp_clus_v = storage->find_one_ass(ev_pfp->id(), ev_clus, ev_pfp->name());
    
    if (!ev_clus) {
      std::cout << "No clus data-product -> exit" << std::endl;
      return false;
    } 
    
    auto const& ass_clus_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit) {
      std::cout << "No hit data-product -> exit" << std::endl;
      return false;
    } 

    // find cluster indices associated with the various showers on the collection plane
    std::vector< std::vector<size_t> > shower_hit_idx_v_v;
    for (auto const& pfp_idx_v : ass_shr_pfp_v) {
      std::vector<size_t> shr_hit_idx_v;
      auto const& clus_idx_v = ass_pfp_clus_v[ pfp_idx_v[0] ];
      for (auto const& clus_idx : clus_idx_v) {
	auto const& hit_idx_v = ass_clus_hit_v[clus_idx];
	for (auto const& hit_idx : hit_idx_v) {
	  auto const& hit = ev_hit->at(hit_idx);
	  if (hit.WireID().Plane == 2)
	    shr_hit_idx_v.push_back( hit_idx );
	}// for all hits in cluster
      }// for all clusters associated to the shower
      shower_hit_idx_v_v.push_back( shr_hit_idx_v );
    }// for all showers

    // used track id vector
    std::vector<unsigned int> used_trk_id;

    // map connecting 1st vs. 2nd shower filled and mcindex used
    std::map<size_t,size_t> shrmapidx;

    // Create G4 track ID vector for MCShowers & MCTracks in the event
    std::vector<std::vector<unsigned int> > g4_trackid_v;
    std::vector<unsigned int> mc_index_v;
    // keep track of all the shower and track trackIDs
    g4_trackid_v.reserve(ev_mcs->size());
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
      }
      if (mc_index_v.size() == 2) {
	_etrue1 = mcs.Start().E();
	_edep1  = mcs.DetProfile().E();
	_qcol1  = mcs.Charge(2);
      }	

    }// for all mcshowers

    if (mc_index_v.size() != 2) {
      std::cout << "did not find 2 MCShowers... quit " << std::endl;
      return false;
    }

    _angle = ev_mcs->at(mc_index_v[0]).Start().Momentum().Vect().Unit().Dot( ev_mcs->at(mc_index_v[1]).Start().Momentum().Vect().Unit() );

    // reset MCBTAlg
    _bt_algo.Reset(g4_trackid_v,*ev_simch);

    // create a vector of vectors in which to store the hit associations
    // for each cluster produced
    std::vector<std::vector<unsigned int> > cluster_hit_v;
    cluster_hit_v.resize( 3 * mc_index_v.size() );
    // vector to save the plane information for each cluster
    std::vector<larlite::geo::View_t> cluster_plane_v(cluster_hit_v.size(),
						      larlite::geo::View_t::kUnknown);

    std::cout << "There are " << shower_hit_idx_v_v.size() << " showers" << std::endl;


    _qhit0 = _qhit1 = _ehit0 = _ehit0 = 0;

    _chTickMap.clear();
    _chTickMap = std::vector< std::vector< std::pair<int,int> > >(8500,std::vector<std::pair<int,int> >());
    
    _chIDEmap.clear();
    _chIDEmap = std::vector< std::map<int,double> >(8500,std::map<int,double>());
    
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

      // avoid duplicating time-intervals
      std::pair<double,double> timeinterval = std::make_pair(tstart,tend);

      if (_avoid_duplicate_ticks)
	auto timeinterval = getTimeSubset(ch,(int)tstart,(int)tend);

      if (timeinterval.first >= timeinterval.second) {
	std::cout << "skipping hit " << std::endl;
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
	for (size_t j=0; j < mcq_v.size(); j++){
	  if (mcq_v[j] > max_edep){
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
	    if (ch < _wmin0) _wmin0 = ch;
	    if (ch > _wmax0) _wmax0 = ch;
	  }
	  if (shrmapidx[idx] == 1) {
	    _ehit1 += max_edep;
	    _qhit1 += max_qdep;
	    _ahit1 += hit.Integral();
	    if (ch < _wmin1) _wmin1 = ch;
	    if (ch > _wmax1) _wmax1 = ch;
	  }

	}// if on collection plane
	
	
      }// for all found time-intervals
      
    }// for all hits in cluster

    // loop through clusters associated to showers on the collection plane
    for (auto const& shr_hit_idx_v : shower_hit_idx_v_v) {

      _chTickMap.clear();
      _chTickMap = std::vector< std::vector< std::pair<int,int> > >(8500,std::vector<std::pair<int,int> >());
      
      _chIDEmap.clear();
      _chIDEmap = std::vector< std::map<int,double> >(8500,std::map<int,double>());
      
      // reset charge integrators
      _eclus0 = _qclus0 = _eclus1 = _qclus1 = 0;
      
      _wmin0 = _wmin1 = 9000;
      _wmax0 = _wmax1 = 0;
      
      // loop through hits in cluster, use the back-tracker to find which MCX object
      // it should belong to, and add that hit to the cluster that is
      // indended for that MCX object
      // only use hits from association to rawclusters
      for (auto const& i : shr_hit_idx_v) {

	auto const& hit_idx = i;
	auto const& hit     = ev_hit->at(hit_idx);
	auto const& ch      = hit.Channel();
	auto const& tstart  = hit.PeakTime() - 2*hit.RMS();// + 2255;//3050;
	auto const& tend    = hit.PeakTime() + 2*hit.RMS();// + 2255;//3050;
	auto const& pl      = hit.View();
	
	// avoid duplicating time-intervals
	std::pair<double,double> timeinterval = std::make_pair(tstart,tend);
	
	if (_avoid_duplicate_ticks)
	  auto timeinterval = getTimeSubset(ch,(int)tstart,(int)tend);
	
	if (timeinterval.first >= timeinterval.second) {
	  std::cout << "skipping hit " << std::endl;
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
	  for (size_t j=0; j < mcq_v.size(); j++){
	    if (mcq_v[j] > max_edep){
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
	  if ( idx == mcq_v.size() -1 )
	    continue;
	  
	  if (pl == 2) {

	    if (shrmapidx[idx] == 0) {
	      _eclus0 += max_edep;
	      _qclus0 += max_qdep;
	      if (ch < _wmin0) _wmin0 = ch;
	      if (ch > _wmax0) _wmax0 = ch;
	    }

	    if (shrmapidx[idx] == 1) {
	      _eclus1 += max_edep;
	      _qclus1 += max_qdep;
	      if (ch < _wmin1) _wmin1 = ch;
	      if (ch > _wmax1) _wmax1 = ch;
	    }

	  }// if on collection plane
	  
	}// for all found time-intervals
	
      }// for all hits in cluster

      _tree->Fill();

    }// for all showers
    
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

  bool Pi0ClusteringStudies::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();
    if (_hit_tree) _hit_tree->Write();

    return true;
  }


  std::pair<int,int> Pi0ClusteringStudies::getTimeSubset(const int& ch,const int& tstart,const int& tend) {

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
