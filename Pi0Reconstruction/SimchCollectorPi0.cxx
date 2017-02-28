#ifndef LARLITE_SIMCHCOLLECTORPI0_CXX
#define LARLITE_SIMCHCOLLECTORPI0_CXX

#include "SimchCollectorPi0.h"

#include "DataFormat/mcpart.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include "DataFormat/simch.h"
#include "DataFormat/hit.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool SimchCollectorPi0::initialize() {

    _entry = 0;

    _tree = new TTree("_tree","TREE");
    _tree->Branch("_mc_e_1",&_mc_e_1,"mc_e_1/D");
    _tree->Branch("_simch_e_1",&_simch_e_1,"simch_e_1/D");
    _tree->Branch("_simch_q_1",&_simch_q_1,"simch_q_1/D");
    _tree->Branch("_adc_e_1",&_adc_e_1,"adc_e_1/D");
    _tree->Branch("_mc_e_2",&_mc_e_2,"mc_e_2/D");
    _tree->Branch("_simch_e_2",&_simch_e_2,"simch_e_2/D");
    _tree->Branch("_simch_q_2",&_simch_q_2,"simch_q_2/D");
    _tree->Branch("_adc_e_2",&_adc_e_2,"adc_e_2/D");
    _tree->Branch("_gamma_angle",&_gamma_angle,"gamma_angle/D");
    _tree->Branch("_gamma_angle_Y",&_gamma_angle_Y,"gamma_angle_Y/D");
    _tree->Branch("_pi0_e",&_pi0_e,"pi0_e/D");
    _tree->Branch("_vx",&_vx,"vx/D"); 
    _tree->Branch("_vy",&_vy,"vy/D");
    _tree->Branch("_vz",&_vz,"vz/D");
    _tree->Branch("_mass_truth",&_mass_truth,"mass_truth/D");
    _tree->Branch("_mass_simch",&_mass_simch,"mass_simch/D");
    _tree->Branch("_mass_adc",&_mass_adc,"mass_adc/D");

    _hit_tree = new TTree("_hit_tree","HIT TREE");
    _hit_tree->Branch("_hit_x",&_hit_x,"hit_x/D");
    _hit_tree->Branch("_hit_z",&_hit_z,"hit_z/D");
    _hit_tree->Branch("_entry",&_entry,"entry/I");

    return true;
  }
  
  bool SimchCollectorPi0::analyze(storage_manager* storage) {

    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    double t2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    auto *ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto *ev_simch    = storage->get_data<event_simch>("largeant");
    auto *ev_hit      = storage->get_data<event_hit>("gaushit");

    if (!ev_mcshower) { std::cout << "NO MCSHOWER" << std::endl; return true; }
    if (!ev_simch)    { std::cout << "NO SIMCH"    << std::endl; return true; }
    if (!ev_hit)      { std::cout << "NO HIT"      << std::endl; return true; }

    // get the pi0 information and that of the two showers.
    double px1, py1, pz1; // 1st shower direction
    double px2, py2, pz2; // 2nd shower direction
    
    if (ev_mcshower->size() != 2){
      std::cout << "Number of EM showers != 2...skip" << std::endl;
      return true;
    }
    
    auto const& mcs1 = ev_mcshower->at(0);
    auto const& mcs2 = ev_mcshower->at(1);

    px1 = mcs1.Start().Px();
    py1 = mcs1.Start().Py();
    pz1 = mcs1.Start().Pz();
    double mag1_Y = sqrt ( (px1*px1) + (pz1*pz1) );
    double mag1   = sqrt ( (px1*px1) + (py1*py1) + (pz1*pz1) );
    double px1_Y = px1 / mag1_Y;
    double pz1_Y = pz1 / mag1_Y;
    px1 /= mag1;
    py1 /= mag1;
    pz1 /= mag1;

    px2 = mcs2.Start().Px();
    py2 = mcs2.Start().Py();
    pz2 = mcs2.Start().Pz();
    double mag2_Y = sqrt ( (px2*px2) + (pz2*pz2) );
    double mag2   = sqrt ( (px2*px2) + (py2*py2) + (pz2*pz2) );
    double px2_Y = px2 / mag2_Y;
    double pz2_Y = pz2 / mag2_Y;
    px2 /= mag2;
    py2 /= mag2;
    pz2 /= mag2;

    _vx = mcs1.MotherStart().X();
    _vy = mcs1.MotherStart().Y();
    _vz = mcs1.MotherStart().Z();

    // calculate angle between two momenta on collection plane
    _gamma_angle_Y = (px1_Y*px2_Y) + (pz1_Y*pz2_Y);

    // calculate direction of line separating the two momenta vectors in 2D
    double px_median_Y = px1_Y + px2_Y;
    double pz_median_Y = pz1_Y + pz2_Y;
    double mag_median_Y = sqrt( (px_median_Y*px_median_Y) + (pz_median_Y)*(pz_median_Y) );
    px_median_Y /= mag_median_Y;
    pz_median_Y /= mag_median_Y;

    // calculate equation of line in Z-X plane which separates the two photons.

    double slope = px_median_Y / pz_median_Y;
    double intercept = _vx - slope * _vz;

    // "1" refers to the shower above the line
    if (mcs1.Start().X() > (mcs1.Start().Z() * slope + intercept) ){
      _mc_e_1 = mcs1.Start().E();
      _mc_e_2 = mcs2.Start().E();
    }
    else{
      _mc_e_2 = mcs1.Start().E();
      _mc_e_1 = mcs2.Start().E();
    }
    double _gamma_angle = (px1*px2) + (py1*py2) + (pz1*pz2);

    _mass_truth = sqrt( 2 * _mc_e_1 * _mc_e_2 * (1 - _gamma_angle) );
    
    _simch_e_1 = 0;
    _simch_q_1 = 0;

    _simch_e_2 = 0;
    _simch_q_2 = 0;

    for (size_t l=0; l < ev_simch->size(); l++){
      
      auto const& simch = ev_simch->at(l);
      
      if (simch.Channel() < 4800)
	continue;
      
      auto const& all_ide = simch.TrackIDsAndEnergies(0,19600);
      
      for (size_t j=0; j < all_ide.size(); j++){
	
	auto const& ide = all_ide[j];

	// figure out which side of the line the charge lies on
	if ( ide.x > (ide.z * slope + intercept) ){
	  _simch_e_1 += ide.energy;
	  _simch_q_1 += ide.numElectrons;
	}
	else{
	  _simch_e_2 += ide.energy;
	  _simch_q_2 += ide.numElectrons;
	}
	
      }// fr all IDEs

    }// for all simchannels

    _mass_simch = sqrt( 2 * _simch_e_1 * _simch_e_2 * (1 - _gamma_angle) );

    _tree->Fill();

    _adc_e_1 = 0;
    _adc_e_2 = 0;

    for (size_t h=0; h < ev_hit->size(); h++){

      auto const& hit = ev_hit->at(h);

      if (hit.WireID().Plane != 2) continue;

      double w = hit.WireID().Wire * w2cm ;
      double t = hit.PeakTime() * t2cm - 45.; 

      _hit_x = t - _vx;
      _hit_z = w - _vz;
      //_hit_tree->Fill();
      
      double adc = hit.Integral();
      
      if ( t > (w * slope + intercept) )
	_adc_e_1 += adc;
      else
	_adc_e_2 += adc;

    }// for all hits

    _mass_adc = sqrt( 2 * _adc_e_1 * _adc_e_2 * (1 - _gamma_angle) );

    _tree->Fill();
    
    return true;
  }

  bool SimchCollectorPi0::finalize() {

    _entry += 1;

    _tree->Write();
    //_hit_tree->Write();
    
    return true;
  }

}
#endif
