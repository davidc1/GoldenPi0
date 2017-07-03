#ifndef LARLITE_SHOWERCONTAINMENT_CXX
#define LARLITE_SHOWERCONTAINMENT_CXX

#include "ShowerContainment.h"

#include "DataFormat/mcshower.h"

namespace larlite {

  bool ShowerContainment::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_edep",&_edep,"edep/D");
    _tree->Branch("_etot",&_etot,"etot/D");
    _tree->Branch("_r",&_r,"r/D");
    _tree->Branch("_x",&_x,"x/D");
    _tree->Branch("_y",&_y,"y/D");
    _tree->Branch("_z",&_z,"z/D");
    _tree->Branch("_px",&_px,"px/D");
    _tree->Branch("_py",&_py,"py/D");
    _tree->Branch("_pz",&_pz,"pz/D");

    _TPC = ::geoalgo::AABox(0,-116,0,256,116,1037);

    return true;
  }
  
  bool ShowerContainment::analyze(storage_manager* storage) {

    auto const& ev_mcs = storage->get_data<event_mcshower>("mcreco");

    for (auto const &mcs : *ev_mcs) {

      if (mcs.DetProfile().E() <= 0) continue;

      _edep = mcs.DetProfile().E();
      _etot = mcs.Start().E();
      _x    = mcs.DetProfile().X();
      _y    = mcs.DetProfile().Y();
      _z    = mcs.DetProfile().Z();
      _px   = mcs.DetProfile().Px();
      _py   = mcs.DetProfile().Py();
      _pz   = mcs.DetProfile().Pz();
      double mag = sqrt( _px*_px + _py*_py + _pz*_pz );
      _px /= mag;
      _py /= mag;
      _pz /= mag;
     

      ::geoalgo::Vector strt(_x,_y,_z);

      // create HalfLine for this shower
      ::geoalgo::HalfLine hline(_x,_y,_z,_px, _py, _pz);

      auto intersectionPts = _geoAlgo.Intersection(_TPC,hline);
      
      if (intersectionPts.size() != 1) continue;
      
      _r = intersectionPts[0].Dist( strt );

      _tree->Fill();

    }// for all MCShowers
  
    return true;
  }

  bool ShowerContainment::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();

    return true;
  }

}
#endif
