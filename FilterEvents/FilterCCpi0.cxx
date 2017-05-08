#ifndef LARLITE_FILTERCCPI0_CXX
#define LARLITE_FILTERCCPI0_CXX

#include "FilterCCpi0.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool FilterCCpi0::initialize() {

    return true;
  }
  
  bool FilterCCpi0::analyze(storage_manager* storage) {
    
    auto ev_mctruth= storage->get_data<event_mctruth>("generator");
    
    if (!ev_mctruth) { std::cout << "no truth" << std::endl; }
    
    auto nu = ev_mctruth->at(0).GetNeutrino();
    auto parts = ev_mctruth->at(0).GetParticles();
    
    if (nu.CCNC() == 1) return false;
    
    int npi0 = 0;
    
    for (size_t i=0; i < parts.size(); i++){
      auto const& part = parts[i];
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) )
	npi0 += 1;
    }// for all mcparticles                                                                                                                                                                               

    if (npi0 != 0)
	return false;

    return true;
  }

  bool FilterCCpi0::finalize() {

    return true;
  }

}
#endif
