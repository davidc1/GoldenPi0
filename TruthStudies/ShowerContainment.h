/**
 * \file ShowerContainment.h
 *
 * \ingroup TruthStudies
 * 
 * \brief Class def header for a class ShowerContainment
 *
 * @author dcaratelli
 */

/** \addtogroup TruthStudies

    @{*/

#ifndef LARLITE_SHOWERCONTAINMENT_H
#define LARLITE_SHOWERCONTAINMENT_H

#include "Analysis/ana_base.h"

#include "GeoAlgo/GeoAlgo.h"

#include "TTree.h"

namespace larlite {
  /**
     \class ShowerContainment
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerContainment : public ana_base{
  
  public:

    /// Default constructor
    ShowerContainment()
      : _tree(nullptr)
      { _name="ShowerContainment"; _fout=0;}

    /// Default destructor
    virtual ~ShowerContainment(){}

    /** IMPLEMENT in ShowerContainment.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerContainment.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerContainment.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    ::geoalgo::AABox _TPC;
    ::geoalgo::GeoAlgo _geoAlgo;

    TTree* _tree;
    double _x,_y,_z;
    double _px,_py,_pz;
    double _edep, _etot;
    double _r;
    
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
