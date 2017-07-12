/**
 * \file Pi0Selection.h
 *
 * \ingroup Pi0Reco
 * 
 * \brief Class def header for a class Pi0Selection
 *
 * @author david caratelli
 */

/** \addtogroup Pi0Reco

    @{*/

#ifndef LARLITE_PI0SELECTION_H
#define LARLITE_PI0SELECTION_H

#include "Analysis/ana_base.h"

#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"

#include "GeoAlgo/GeoAlgo.h"

#include "TVector3.h"
#include "TTree.h"

namespace larlite {
  /**
     \class Pi0Selection
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi0Selection : public ana_base{
  
  public:

    /// Default constructor
    Pi0Selection()
      : _tree(nullptr)
      , _pi0_tree(nullptr)
      { _name="Pi0Selection"; _fout=0;}

    /// Default destructor
    virtual ~Pi0Selection(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVtxProducer(std::string s) { _vtx_producer = s; }
    void setEMin(double emin) { _emin = emin; }
    void setRadLenMax(double l) { _radlenmax = l; }
    void setIPMax(double ip) { _ipmax = ip; }
    void setAngleMin(double amin) { _anglemin = amin; }
    void doMCMatch(bool on) { _mcmatch = on; }

  protected:

    // perform MC <-> RC shower matching
    std::vector< std::pair<size_t,size_t> > MCRCMatch(const std::vector<larlite::shower>& shr_v,
						      const std::vector<larlite::mcshower>& mcs_v);

    // filter out showers based on various cuts
    const std::vector<size_t> FilterShowers(const larlite::event_shower* shr_v);

    // create pair-combinatorics
    const std::vector< std::pair<size_t,size_t> > Combinatorics(const std::vector<size_t> idx_v);

    std::string _vtx_producer;

    // reconstructed and true vertex coordinates
    ::geoalgo::Point_t _vtx, _mcvtx;

    // vector where to store true gamma MCshower objects
    std::vector<larlite::mcshower> pi0_gamma_v;

    bool _mcmatch;

    /// minimum shower energy
    double _emin;
    /// maximum radiation length
    double _radlenmax;
    /// maximum IP allowed
    double _ipmax;
    /// minimum opening angle allowed
    double _anglemin;

    ::geoalgo::GeoAlgo _geoAlgo;

    TTree *_tree;
    double _ip;
    double _angle;
    double _mass;
    double _el,_eh; // low and high shower energies
    double _rl,_rh; // low and hgh conversion distances
    double _rc0x, _rc0y, _rc0z, _rc1x, _rc1y, _rc1z;
    double _mc0x, _mc0y, _mc0z, _mc1x, _mc1y, _mc1z;
    double _angle0, _angle1;
    double _d0, _d1;
    double _rce0, _rce1, _mce0, _mce1, _mcedep0, _mcedep1;
    int    _npi0;
    double _mcvtxx, _mcvtxy, _mcvtxz;
    double _rcvtxx, _rcvtxy, _rcvtxz;
    int    _nrecoshr;
    int    _run, _sub, _evt, _ctr;

    TTree *_pi0_tree;
    
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
