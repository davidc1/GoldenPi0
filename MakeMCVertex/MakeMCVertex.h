/**
 * \file MakeMCVertex.h
 *
 * \ingroup MakeMCVertex
 * 
 * \brief Class def header for a class MakeMCVertex
 *
 * @author david
 */

/** \addtogroup MakeMCVertex

    @{*/

#ifndef LARLITE_MAKEMCVERTEX_H
#define LARLITE_MAKEMCVERTEX_H

#include "Analysis/ana_base.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

namespace larlite {
  /**
     \class MakeMCVertex
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeMCVertex : public ana_base{
  
  public:

    /// Default constructor
    MakeMCVertex(){ _name="MakeMCVertex"; _fout=0;}

    /// Default destructor
    virtual ~MakeMCVertex(){}

    /** IMPLEMENT in MakeMCVertex.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeMCVertex.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeMCVertex.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetXOffset(double d) { _offset = d; }

  protected:

    double _offset;

    double _time2cm;

    larutil::SpaceChargeMicroBooNE *_SCE;

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
