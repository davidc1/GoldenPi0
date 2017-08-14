/**
 * \file SaveSel2VtxTrkMCC82.h
 *
 * \ingroup NuMuCCFilterStudy
 * 
 * \brief Class def header for a class SaveSel2VtxTrkMCC82
 *
 * @author david
 */

/** \addtogroup NuMuCCFilterStudy

    @{*/

#ifndef LARLITE_SAVESEL2VTXTRKMCC82_H
#define LARLITE_SAVESEL2VTXTRKMCC82_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SaveSel2VtxTrkMCC82
     User custom analysis class made by SHELL_USER_NAME
   */
  class SaveSel2VtxTrkMCC82 : public ana_base{
  
  public:

    /// Default constructor
    SaveSel2VtxTrkMCC82();

    /// Default destructor
    virtual ~SaveSel2VtxTrkMCC82(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetVerbose(bool on) { _verbose = on; }

  protected:

    bool _verbose;
    
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
