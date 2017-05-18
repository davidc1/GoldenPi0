#ifndef LARLITE_MAKERANDOMVERTEX_CXX
#define LARLITE_MAKERANDOMVERTEX_CXX

#include "MakeRandomVertex.h"

#include "DataFormat/vertex.h"

#include <random>

namespace larlite {

  bool MakeRandomVertex::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("randomvtx","randomvtx");
    _tree->Branch("_x",&_x,"x/D");
    _tree->Branch("_y",&_y,"y/D");
    _tree->Branch("_z",&_z,"z/D");

    srand( time(NULL) );
    
    return true;
  }
  
  bool MakeRandomVertex::analyze(storage_manager* storage) {
  
    auto ev_vtx = storage->get_data<event_vertex>("randomvertex");

    storage->set_id(storage->run_id(), storage->subrun_id(), storage->event_id());

    // random x / y / z positions

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> xdis(20,236);
    std::uniform_real_distribution<> ydis(-96,96);
    std::uniform_real_distribution<> zdis(20,1016);


    _x = xdis(gen);
    _y = ydis(gen);
    _z = zdis(gen);
    
    _tree->Fill();

    double xyz[3] = {_x,_y,_z};
    vertex new_vtx(xyz);
    ev_vtx->push_back(new_vtx);
    
    return true;
  }

  bool MakeRandomVertex::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();

    return true;
  }

}
#endif
