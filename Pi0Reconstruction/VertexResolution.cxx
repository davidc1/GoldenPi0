#ifndef LARLITE_VERTEXRESOLUTION_CXX
#define LARLITE_VERTEXRESOLUTION_CXX

#include "VertexResolution.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  VertexResolution::VertexResolution() 
    : _tree(nullptr)
  {
    _fout = 0;
    _name = "VertexResolution";
  }


  bool VertexResolution::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("_mc_x",&_mc_x,"mc_x/D");
    _tree->Branch("_mc_y",&_mc_y,"mc_y/D");
    _tree->Branch("_mc_z",&_mc_z,"mc_z/D");
    _tree->Branch("_rc_x",&_rc_x,"rc_x/D");
    _tree->Branch("_rc_y",&_rc_y,"rc_y/D");
    _tree->Branch("_rc_z",&_rc_z,"rc_z/D");
    
    return true;
  }
  
  bool VertexResolution::analyze(storage_manager* storage) {

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    auto ev_mc  = storage->get_data<event_mctruth>("generator");

    if (!ev_mc or (ev_mc->size() == 0) ){
      std::cout << "No truth info..." << std::endl;
      return false;
    }

    if (!ev_vtx or (ev_vtx->size() == 0) ){
      std::cout << "No vertex info..." << std::endl;
      return false;
    }
    
    auto mctruth = ev_mc->front();
    auto const& mcvtx = mctruth.GetNeutrino().Lepton().Trajectory().front();

    _mc_x = mcvtx.X();
    _mc_y = mcvtx.Y();
    _mc_z = mcvtx.Z();

    auto const& vtx = ev_vtx->at(0);
    _rc_x = vtx.X();
    _rc_y = vtx.Y();
    _rc_z = vtx.Z();

    _tree->Fill();

    return true;
  }

  bool VertexResolution::finalize() {

    if (_fout) 
      _fout->cd();

    _tree->Write();

    return true;
  }

}
#endif
