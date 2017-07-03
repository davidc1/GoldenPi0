/**
 * \file MakeSingleParticleVertex.h
 *
 * \ingroup MakeMCVertex
 * 
 * \brief Class def header for a class MakeSingleParticleVertex
 *
 * @author dcaratelli
 */

/** \addtogroup MakeMCVertex

    @{*/

#ifndef LARLITE_MAKESINGLEPARTICLEVERTEX_H
#define LARLITE_MAKESINGLEPARTICLEVERTEX_H

#include "Analysis/ana_base.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

namespace larlite {
  /**
     \class MakeSingleParticleVertex
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeSingleParticleVertex : public ana_base{
  
  public:

    /// Default constructor
    MakeSingleParticleVertex(){ _name="MakeSingleParticleVertex"; _fout=0;}

    /// Default destructor
    virtual ~MakeSingleParticleVertex(){}

    /** IMPLEMENT in MakeSingleParticleVertex.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeSingleParticleVertex.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeSingleParticleVertex.cc! 
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
