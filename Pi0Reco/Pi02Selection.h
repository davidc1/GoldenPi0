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

#ifndef LARLITE_PI02SELECTION_H
#define LARLITE_PI02SELECTION_H

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
     \class Pi02Selection
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi02Selection : public ana_base{
  
  public:

    /// Default constructor
    Pi02Selection();

    /// Default destructor
    virtual ~Pi02Selection(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVtxProducer(std::string s) { _vtx_producer = s; }
    void setShrProducer(std::string s) { _shr_producer = s; }
    void setEMin(double emin) { _emin = emin; }
    void setRadLenMax(double l) { _radlenmax = l; }
    void setIPMax(double ip) { _ipmax = ip; }
    void setAngleMin(double amin) { _anglemin = amin; }
    void doMCMatch(bool on) { _mcmatch = on; }
    void EDepMin(double e) { _edepmin = e; }
    void ApplyContainmentCorrection(bool on) { _containmentcorrection = on; }

    /// shower containment correction
    double ContainmentCorr(double dwall);

  protected:

    // filter out showers based on various cuts
    const std::vector<size_t> FilterShowers(const larlite::event_shower* shr_v);

    // create pair-combinatorics
    const std::vector< std::pair<size_t,size_t> > Combinatorics(const std::vector<size_t> idx_v);

    // fill ttree
    void FillTree(const larlite::shower& shr1, const larlite::shower& shr2);

    std::string _vtx_producer, _shr_producer;

    // reconstructed and true vertex coordinates
    ::geoalgo::Point_t _rcvtx, _mcvtx;

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
    /// minimum truth-level E de
    double _edepmin;

    ::geoalgo::GeoAlgo _geoAlgo;

    TTree *_tree;
    double _ip;
    double _angle;
    double _mass, _massc;
    double _el,_eh, _elc, _ehc; // low and high shower energies
    double _rl,_rh; // low and hgh conversion distances
    double _rc0x, _rc0y, _rc0z, _rc1x, _rc1y, _rc1z;
    double _mc0x, _mc0y, _mc0z, _mc1x, _mc1y, _mc1z;
    double _dedxl, _dedxh;
    double _pitchl, _pitchh;
    double _angle0, _angle1;
    double _d0, _d1;
    double _rce0, _rce1, _mce0, _mce1, _mcedep0, _mcedep1;
    int    _npi0;   // number of pi0s coming out of the interaction
    int    _ngamma; // number of showers coming out of the interaction
    double _pi0px, _pi0py, _pi0pz;
    double _pi0e, _pi0a; // energy and angle of the pi0
    double _mcvtxx, _mcvtxy, _mcvtxz;
    double _rcvtxx, _rcvtxy, _rcvtxz;
    double _ipvtx;
    int    _nrecoshr, _nrecoshrcut;
    int    _run, _sub, _evt, _ctr;

    TTree *_pi0_tree;
    double _hxs, _hys, _hzs, _lxs, _lys, _lzs, _hxe, _hye, _hze, _lxe, _lye, _lze;

    TTree *_shower_tree;
    double _rce;
    double _mce;
    double _mcedep;
    double _anglediff;
    double _rcx, _rcy, _rcz, _rcpx, _rcpy, _rcpz;
    double _dwall, _dwallh, _dwalll;

    ::geoalgo::AABox_t _TPC;

    // apply containmnent correction
    bool _containmentcorrection;

    double _m0, _m1; // mass of pair of pi0s if 2 are reconstructed
    
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
