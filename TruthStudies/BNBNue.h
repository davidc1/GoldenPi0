/**
 * \file BNBNue.h
 *
 * \ingroup TruthStudies
 * 
 * \brief Class def header for a class BNBNue
 *
 * @author davidc1
 */

/** \addtogroup TruthStudies

    @{*/

#ifndef LARLITE_BNBNUE_H
#define LARLITE_BNBNUE_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class BNBNue
     User custom analysis class made by SHELL_USER_NAME
   */
  class BNBNue : public ana_base{
  
  public:

    /// Default constructor
    BNBNue();

    /// Default destructor
    virtual ~BNBNue(){}

    /** IMPLEMENT in BNBNue.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BNBNue.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BNBNue.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _tree;

    double _nu_e, _pi0_e;
    double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
    
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
