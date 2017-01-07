/**
 * \file SimchCollector.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class SimchCollector
 *
 * @author david
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_SIMCHCOLLECTOR_H
#define LARLITE_SIMCHCOLLECTOR_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class SimchCollector
     User custom analysis class made by SHELL_USER_NAME
   */
  class SimchCollector : public ana_base{
  
  public:

    /// Default constructor
    SimchCollector()
      : _tree(nullptr)
      { _name="SimchCollector"; _fout=0;}

    /// Default destructor
    virtual ~SimchCollector(){}

    /** IMPLEMENT in SimchCollector.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimchCollector.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimchCollector.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _tree;
    double _mc_e, _simch_e, _adc_e;
    
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
