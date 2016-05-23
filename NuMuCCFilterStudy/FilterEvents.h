/**
 * \file FilterEvents.h
 *
 * \ingroup NuMuCCFilterStudy
 * 
 * \brief Class def header for a class FilterEvents
 *
 * @author david
 */

/** \addtogroup NuMuCCFilterStudy

    @{*/

#ifndef LARLITE_FILTEREVENTS_H
#define LARLITE_FILTEREVENTS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FilterEvents
     User custom analysis class made by SHELL_USER_NAME
   */
  class FilterEvents : public ana_base{
  
  public:

    /// Default constructor
    FilterEvents();

    /// Default destructor
    virtual ~FilterEvents(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetFilterShowers(bool on, int n=0) { _filter_showers = on; _n_showers = n; }
    void SetFilterTracks (bool on, int n=0) { _filter_tracks  = on; _n_tracks  = n; }
    void SetMinTracks    (bool on, int n=0) { _filter_tracks  = on; _n_tracks  = n; _min_tracks  = true; }
    void SetMinShowers   (bool on, int n=0) { _filter_showers = on; _n_showers = n; _min_showers = true; }
    void SetMaxDist(double d) { _max_dist = d; }
    void SetVerbose(bool on) { _verbose = on; }

  protected:

    bool _verbose;
    bool _filter_showers;
    bool _filter_tracks;
    bool _min_tracks;
    bool _min_showers;

    double _max_dist;
    
    int _n_tracks;
    int _n_showers;
    
    
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
