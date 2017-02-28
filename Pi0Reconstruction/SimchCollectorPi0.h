/**
 * \file SimchCollectorPi0.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class SimchCollectorPi0
 *
 * @author david
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_SIMCHCOLLECTORPI0_H
#define LARLITE_SIMCHCOLLECTORPI0_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class SimchCollectorPi0
     User custom analysis class made by SHELL_USER_NAME
   */
  class SimchCollectorPi0 : public ana_base{
  
  public:

    /// Default constructor
    SimchCollectorPi0()
      : _tree(nullptr)
      { _name="SimchCollectorPi0"; _fout=0;}

    /// Default destructor
    virtual ~SimchCollectorPi0(){}

    /** IMPLEMENT in SimchCollectorPi0.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimchCollectorPi0.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimchCollectorPi0.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _tree;
    double _mc_e_1, _simch_e_1, _simch_q_1, _adc_e_1;
    double _mc_e_2, _simch_e_2, _simch_q_2, _adc_e_2;
    double _gamma_angle_Y, _gamma_angle;
    double _pi0_e;
    double _vx, _vy, _vz;
    double _mass_truth;
    double _mass_simch;
    double _mass_adc;

    TTree* _hit_tree;
    double _hit_x, _hit_z;
    int _entry;
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
