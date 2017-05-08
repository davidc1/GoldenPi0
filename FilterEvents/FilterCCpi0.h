/**
 * \file FilterCCpi0.h
 *
 * \ingroup FilterEvents
 * 
 * \brief Class def header for a class FilterCCpi0
 *
 * @author davidc1
 */

/** \addtogroup FilterEvents

    @{*/

#ifndef LARLITE_FILTERCCPI0_H
#define LARLITE_FILTERCCPI0_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FilterCCpi0
     User custom analysis class made by SHELL_USER_NAME
   */
  class FilterCCpi0 : public ana_base{
  
  public:

    /// Default constructor
    FilterCCpi0(){ _name="FilterCCpi0"; _fout=0;}

    /// Default destructor
    virtual ~FilterCCpi0(){}

    /** IMPLEMENT in FilterCCpi0.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FilterCCpi0.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FilterCCpi0.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
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
