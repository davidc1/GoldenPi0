/**
 * \file BNBPi0Properties.h
 *
 * \ingroup TruthStudies
 * 
 * \brief Class def header for a class BNBPi0Properties
 *
 * @author davidc1
 */

/** \addtogroup TruthStudies

    @{*/

#ifndef LARLITE_BNBPI0PROPERTIES_H
#define LARLITE_BNBPI0PROPERTIES_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class BNBPi0Properties
     User custom analysis class made by SHELL_USER_NAME
   */
  class BNBPi0Properties : public ana_base{
  
  public:

    /// Default constructor
    BNBPi0Properties();

    /// Default destructor
    virtual ~BNBPi0Properties(){}

    /** IMPLEMENT in BNBPi0Properties.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BNBPi0Properties.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BNBPi0Properties.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _tree;

    int _n_reco_showers;

    int _event;

    double _nu_e, _pi0_e;
    
    double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;

    double _mc_shr1_x,  _mc_shr1_y,  _mc_shr1_z;
    double _mc_shr1_px, _mc_shr1_py, _mc_shr1_pz;
    double _mc_shr1_e;

    double _mc_shr2_x,  _mc_shr2_y,  _mc_shr2_z;
    double _mc_shr2_px, _mc_shr2_py, _mc_shr2_pz;
    double _mc_shr2_e;
    
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
