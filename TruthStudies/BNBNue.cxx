#ifndef LARLITE_BNBNUE_CXX
#define LARLITE_BNBNUE_CXX

#include "BNBNue.h"

#include "DataFormat/mctruth.h"

namespace larlite {

  BNBNue::BNBNue()
    : _tree(nullptr)
  {
    _name="BNBNue";
    _fout=0;
  }

  bool BNBNue::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","shower study tree");
    _tree->Branch("_nu_e",&_nu_e,"nu_e/D");
    _tree->Branch("_mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/D");
    _tree->Branch("_mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/D");
    _tree->Branch("_mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/D");

    return true;
  }
  
  bool BNBNue::analyze(storage_manager* storage) {

    _mc_vtx_x = -999;
    _mc_vtx_y = -999;
    _mc_vtx_z = -999;

    auto *ev_mctruth  = storage->get_data<event_mctruth>("generator");

    if ( (!ev_mctruth) or (ev_mctruth->size() == 0) )
      return false; 
    
    // get all MCParticles
    auto mctruth = ev_mctruth->at(0);

    auto nu = mctruth.GetNeutrino().Nu();
    
    _nu_e = nu.Trajectory().at(0).E();

    auto part_v = mctruth.GetParticles();

    for (size_t i=0; i < part_v.size(); i++){
      auto const& part = part_v[i];
      if ( part.StatusCode() == 1 ){
	_mc_vtx_x = part.Trajectory().at( 0 ).X();
	_mc_vtx_y = part.Trajectory().at( 0 ).Y();
	_mc_vtx_z = part.Trajectory().at( 0 ).Z();
	break;
      }
    }
    
    _tree->Fill();
  
    return true;
  }

  bool BNBNue::finalize() {

    _tree->Write();
  
    return true;
  }

}
#endif
