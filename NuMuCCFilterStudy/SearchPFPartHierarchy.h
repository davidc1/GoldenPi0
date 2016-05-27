/**
 * \file SearchPFPartHierarchy.h
 *
 * \ingroup NuMuCCFilterStudy
 * 
 * \brief Class def header for a class SearchPFPartHierarchy
 *
 * @author david
 */

/** \addtogroup NuMuCCFilterStudy

    @{*/

#ifndef LARLITE_SEARCHPFPARTHIERARCHY_H
#define LARLITE_SEARCHPFPARTHIERARCHY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SearchPFPartHierarchy
     User custom analysis class made by SHELL_USER_NAME
   */
  class SearchPFPartHierarchy : public ana_base{
  
  public:

    /// Default constructor
    SearchPFPartHierarchy();

    /// Default destructor
    virtual ~SearchPFPartHierarchy(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetFilterShowers(bool on, int n=0) { _filter_showers = on; _n_showers = n; }
    void SetFilterTracks (bool on, int n=0) { _filter_tracks  = on; _n_tracks  = n; }
    void SetMinTracks    (bool on, int n=0) { _filter_tracks  = on; _n_tracks  = n; _min_tracks  = true; }
    void SetMinShowers   (bool on, int n=0) { _filter_showers = on; _n_showers = n; _min_showers = true; }
    void SetMinDaughters (bool on, int n=0) { _min_daughters  = on; _n_daughters = n; }
    void SetVerbose(bool on) { _verbose = on; }

  protected:

    bool _verbose;
    bool _filter_showers;
    bool _filter_tracks;
    bool _min_tracks;
    bool _min_showers;
    bool _min_daughters;
    
    int _n_tracks;
    int _n_showers;
    int _n_daughters;
    
    
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
