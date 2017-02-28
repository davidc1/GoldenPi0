#ifndef LARLITE_SIMCHCOLLECTOR_CXX
#define LARLITE_SIMCHCOLLECTOR_CXX

#include "SimchCollector.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include "DataFormat/simch.h"
#include "DataFormat/hit.h"

namespace larlite {
  
  bool SimchCollector::initialize() {
    
    _tree = new TTree("_ide_tree","IDE TREE");
    _tree->Branch("_mc_e",&_mc_e,"mc_e/D");
    _tree->Branch("_mc_edep",&_mc_edep,"mc_edep/D");
    _tree->Branch("_simch_e",&_simch_e,"simch_e/D");
    _tree->Branch("_adc_e",&_adc_e,"adc_e/D");
    _tree->Branch("_shr_z",&_shr_z,"shr_z/D");

    return true;
  }
  
  bool SimchCollector::analyze(storage_manager* storage) {

    std::cout << "Event : " << storage->event_id() << std::endl;
    
    auto *ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto *ev_simch    = storage->get_data<event_simch>("largeant");
    auto *ev_hit      = storage->get_data<event_hit>("gaushit");
    
    std::cout << "there are " << ev_mcshower->size() << " MCShowers" << std::endl;

    if (ev_mcshower->size() != 1) return true;

    auto const& mcs = ev_mcshower->at(0);

    _mc_e = mcs.Start().E();
    
    _mc_edep = mcs.DetProfile().E();
    
    _shr_z = mcs.DetProfile().Z();
    
    _simch_e = 0;

    for (size_t l=0; l < ev_simch->size(); l++){
      
      auto const& simch = ev_simch->at(l);
      
      if (simch.Channel() < 4800)
	continue;
      
      auto const& all_ide = simch.TrackIDsAndEnergies(0,19600);
      
      for (size_t j=0; j < all_ide.size(); j++){
	
	auto const& ide = all_ide[j];

	_simch_e += ide.energy;
	
      }// fr all IDEs

    }// for all simchannels

    _adc_e = 0;

    for (size_t h=0; h < ev_hit->size(); h++){

      auto const& hit = ev_hit->at(h);

      if (hit.WireID().Plane != 2) continue;
      
      _adc_e += hit.Integral();

    }

    _tree->Fill();
    
    return true;
  }

  bool SimchCollector::finalize() {

    _tree->Write();

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
