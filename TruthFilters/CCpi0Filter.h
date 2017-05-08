/**
 * \file CCpi0Filter.h
 *
 * \ingroup TruthFilters
 * 
 * \brief Class def header for a class CCpi0Filter
 *
 * @author david
 */

/** \addtogroup TruthFilters

    @{*/

#ifndef LARLITE_CCPI0FILTER_H
#define LARLITE_CCPI0FILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CCpi0Filter
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCpi0Filter : public ana_base{
  
  public:

    /// Default constructor
    CCpi0Filter(){ _name="CCpi0Filter"; _fout=0;}

    /// Default destructor
    virtual ~CCpi0Filter(){}

    /** IMPLEMENT in CCpi0Filter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCpi0Filter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CCpi0Filter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setNPi0(int n) { _npi0 = n; }

  protected:

    int _npi0; 
    
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
