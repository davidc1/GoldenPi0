#ifndef LARLITE_CCTRUE_CXX
#define LARLITE_CCTRUE_CXX

#include "CCtrue.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool CCtrue::initialize() {

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
    _tree->Branch("_rc_mass"  ,&_rc_mass  ,"_rc_mass/D" );

    _tree->Branch("_dwallmin",&_dwallmin,"dwallmin/D");
    
    
    return true;
  }
  
  bool CCtrue::analyze(storage_manager* storage) {

    _event = storage->event_id();

    // start with mc info

    auto *ev_mctruth  = storage->get_data<event_mctruth>("generator");

      
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

    std::cout << "particle vector size : " << part_v.size() << std::endl;
    
    for (size_t i=0; i < part_v.size(); i++){
      auto const& part = part_v[i];
      std::cout << "PDG code = " << part.PdgCode() << "\t status : " << part.StatusCode() <<  std::endl;
      //std::cout << "found " << part.PdgCode() << std::endl;
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	//std::cout << "found " << part.PdgCode() << "\tw/ energy : " << part.Trajectory().at(0).E() <<  std::endl;
	n_pi0    += 1;
	pi0_idx   = i;
	pi0_trkid = part.TrackId();
	_mc_vtx_x = part.Trajectory().at( 0 ).X();
	_mc_vtx_y = part.Trajectory().at( 0 ).Y();
	_mc_vtx_z = part.Trajectory().at( 0 ).Z();
	_pi0_e    = part.Trajectory().at( 0 ).E() * 1000.;
      }
    }

    _tree->Fill();

    return true;
  }
  
  bool CCtrue::finalize() {

    if (_fout) _fout->cd();
    _tree->Write();
    
    return true;
  }

}
#endif
