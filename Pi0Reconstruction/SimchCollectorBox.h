/**
 * \file SimchCollectorBox.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class SimchCollectorBox
 *
 * @author david
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_SIMCHCOLLECTORBOX_H
#define LARLITE_SIMCHCOLLECTORBOX_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class SimchCollectorBox
     User custom analysis class made by SHELL_USER_NAME
   */
  class SimchCollectorBox : public ana_base{
  
  public:

    /// Default constructor
    SimchCollectorBox()
      : _ide_tree(nullptr)
      , _hit_tree(nullptr)
      { _name="SimchCollectorBox"; _fout=0;}

    /// Default destructor
    virtual ~SimchCollectorBox(){}

    /** IMPLEMENT in SimchCollectorBox.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimchCollectorBox.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimchCollectorBox.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree* _ide_tree;
    double _x, _y, _z;
    double _e;
		       
    TTree* _hit_tree;
    double _w;
    double _t;
    double _adc;
    
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
