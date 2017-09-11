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
    Pi0Selection();

    /// Default destructor
    virtual ~Pi0Selection(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVtxProducer(std::string s) { _vtx_producer = s; }
    void setShrProducer(std::string s) { _shr_producer = s; }
    void setEMin(double emin) { _emin = emin; }
    void setEHighMin(double e) { _ehighmin = e; }
    void setRadLenMax(double l) { _radlenmax = l; }
    void setIPMax(double ip) { _ipmax = ip; }
    void setAngleMin(double amin) { _anglemin = amin; }
    void doMCMatch(bool on) { _mcmatch = on; }
    void EDepMin(double e) { _edepmin = e; }
    void ApplyContainmentCorrection(bool on) { _containmentcorrection = on; }

    /// shower containment correction
    double ContainmentCorr(double dwall);

  protected:

    // perform MC <-> RC shower matching
    std::vector< std::pair<size_t,size_t> > MCRCMatch(const std::vector<larlite::shower>& shr_v,
						      const std::vector<larlite::mcshower>& mcs_v);

    // fill tree
    void FillTree(const larlite::shower& shr1, const larlite::shower& shr2);

    // filter out showers based on various cuts
    const std::vector<size_t> FilterShowers(const larlite::event_shower* shr_v);

    // create pair-combinatorics
    const std::vector< std::pair<size_t,size_t> > Combinatorics(const std::vector<size_t> idx_v,
								const larlite::event_shower* shr_v);

    std::string _vtx_producer, _shr_producer;

    // reconstructed and true vertex coordinates
    ::geoalgo::Point_t _rcvtx, _mcvtx;

    // vector where to store true gamma MCshower objects
    std::vector<larlite::mcshower> pi0_gamma_v;

    bool _mcmatch;

    /// minimum shower energy
    double _emin;
    /// minimum energy for larger shower
    double _ehighmin;
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
    double _el, _eh; // low and high shower energies
    double _elc, _ehc; // after containment correction
    double _rl, _rh; // low and hgh conversion distances
    // start point of reconstructed showers
    double _rchx, _rchy, _rchz, _rclx, _rcly, _rclz;
    // start point of truth showers
    double _mchx, _mchy, _mchz, _mclx, _mcly, _mclz;
    double _dedxh, _dedxl;
    double _pitchh, _pitchl;
    double _anglediffh, _anglediffl;
    double _dwallh, _dwalll;
    double _mcel, _mceh, _mcedepl, _mcedeph;
    int    _npi0;   // number of pi0s coming out of the interaction
    int    _ngamma; // number of showers coming out of the interaction
    double _pi0px, _pi0py, _pi0pz;
    double _pi0e, _pi0a; // energy and angle of the pi0
    double _mcvtxx, _mcvtxy, _mcvtxz;
    double _rcvtxx, _rcvtxy, _rcvtxz;
    double _ipvtx;
    int    _nrecoshr, _nrecoshrcut;
    int    _run, _sub, _evt, _ctr;
    ::geoalgo::Point_t _rcldir, _rchdir;

    TTree *_pi0_tree;

    TTree *_shower_tree;
    double _rce;
    double _mce;
    double _mcedep;
    double _anglediff;
    double _rcx, _rcy, _rcz, _rcpx, _rcpy, _rcpz;
    double _dwall, _dwall0, _dwall1;

    ::geoalgo::AABox_t _TPC;

    // apply containmnent correction
    bool _containmentcorrection;
    
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
