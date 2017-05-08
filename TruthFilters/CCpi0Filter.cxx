#ifndef LARLITE_CCPI0FILTER_CXX
#define LARLITE_CCPI0FILTER_CXX

#include "CCpi0Filter.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool CCpi0Filter::initialize() {

    return true;
  }
  
  bool CCpi0Filter::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator");

    if (!ev_mctruth) { std::cout << "no truth" << std::endl; }

    auto nu = ev_mctruth->at(0).GetNeutrino();
    auto parts = ev_mctruth->at(0).GetParticles();

    // CCNC() == 1 for NC and == 0 for CC
    if ( nu.CCNC() == 1 ) return false;
    
    int npi0 = 0.;
    
    for (size_t i=0; i < parts.size(); i++){
      auto const& part = parts[i];
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) )
	npi0 += 1;
    }// for all mcparticles

    if (npi0 != _npi0) return false;
    
    return true;
  }

  bool CCpi0Filter::finalize() {

    return true;
  }

}
#endif
