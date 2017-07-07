#ifndef LARLITE_PI0HITTHRESHOLDSTUDIES_CXX
#define LARLITE_PI0HITTHRESHOLDSTUDIES_CXX

#include "Pi0HitThresholdStudies.h"

#include "DataFormat/simch.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool Pi0HitThresholdStudies::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_etrue0",&_etrue0,"etrue0/D");
    _tree->Branch("_edep0", &_edep0, "edep0/D" );
    _tree->Branch("_ehit0", &_ehit0, "ehit0/D" );
    _tree->Branch("_qhit0", &_qhit0, "qhit0/D" );

    _tree->Branch("_etrue1",&_etrue1,"etrue1/D");
    _tree->Branch("_edep1", &_edep1, "edep1/D" );
    _tree->Branch("_ehit1", &_ehit1, "ehit1/D" );
    _tree->Branch("_qhit1", &_qhit1, "qhit1/D" );

    _tree->Branch("_angle", &_angle, "angle/D" );

    return true;
  }
  
  bool Pi0HitThresholdStudies::analyze(storage_manager* storage) {

    // Retrieve mcshower data product
    auto ev_mcs   = storage->get_data<event_mcshower>("mcreco");
    // Retrieve hit data product
    auto ev_hit   = storage->get_data<event_hit>("gaushit");
    // Retrieve simch data product
    auto ev_simch = storage->get_data<event_simch>("largeant");

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
      }
      if (mc_index_v.size() == 2) {
	_etrue1 = mcs.Start().E();
	_edep1  = mcs.DetProfile().E();
      }	

    }// for all mcshowers

    if (mc_index_v.size() != 2) {
      std::cout << "did not find 2 MCShowers... quit " << std::endl;
      return false;
    }

    _angle = ev_mcs->at(mc_index_v[0]).Start().Momentum().Vect().Unit().Dot( ev_mcs->at(mc_index_v[1]).Start().Momentum().Vect().Unit() );

    // reset charge integrators
    _ehit0 = _qhit0 = _ehit1 = _qhit1 = 0;

    // reset MCBTAlg
    _bt_algo.Reset(g4_trackid_v,*ev_simch);

    // create a vector of vectors in which to store the hit associations
    // for each cluster produced
    std::vector<std::vector<unsigned int> > cluster_hit_v;
    cluster_hit_v.resize( 3 * mc_index_v.size() );
    // vector to save the plane information for each cluster
    std::vector<larlite::geo::View_t> cluster_plane_v(cluster_hit_v.size(),
						      larlite::geo::View_t::kUnknown);

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
      
      // create a wire-range object with channel + (start,end) time info for the hit
      ::btutil::WireRange_t wr(ch,tstart,tend);
      // use the back-tracker to return a vector of track IDs associated to this hit
      // contents of return vector are proportional to the fraction of the hit's
      // charge that is associated with the MCX id at that element's position
      // vector ordered such that
      // mcq_v [ index ] corresponds to index wihtin mc_index_v
      auto mcq_v = _bt_algo.MCQ(wr);
      // find entry with largest value -> this lets us know which MCX object
      // to assign this hit to and thus which cluster this hit should
      // be associated with
      size_t idx  = 0;
      double max_edep = 0;
      for (size_t j=0; j < mcq_v.size(); j++){
	if (mcq_v[j] > max_edep){
	  max_edep = mcq_v[j];
	  idx  = j;
	}
      }
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
	  _qhit0 += hit.Integral();
	}
	if (shrmapidx[idx] == 1) {
	  _ehit1 += max_edep;
	  _qhit1 += hit.Integral();
	}
      }// if on collection plane
      
      // if not, associate this hit with the appropriate cluster
      // ( "i" is the hit index )
      cluster_hit_v[ pl * mc_index_v.size() + idx ].push_back( hit_idx );
      cluster_plane_v[ pl * mc_index_v.size() + idx ] = pl;
    }// for all hits in cluster

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

  bool Pi0HitThresholdStudies::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();

    return true;
  }

}
#endif
