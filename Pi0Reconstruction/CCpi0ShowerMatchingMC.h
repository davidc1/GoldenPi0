/**
 * \file CCpi0ShowerMatchingMC.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class CCpi0ShowerMatchingMC
 *
 * @author david caratelli
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_CCPI0SHOWERMATCHINGMC_H
#define LARLITE_CCPI0SHOWERMATCHINGMC_H

#include "Analysis/ana_base.h"

#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"

#include "TTree.h"

namespace larlite {
  /**
     \class CCpi0ShowerMatchingMC
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCpi0ShowerMatchingMC : public ana_base{
  
  public:

    /// Default constructor
    CCpi0ShowerMatchingMC()
      : _tree(nullptr)
      { _name="CCpi0ShowerMatchingMC"; _fout=0;}

    /// Default destructor
    virtual ~CCpi0ShowerMatchingMC(){}

    /** IMPLEMENT in CCpi0ShowerMatchingMC.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCpi0ShowerMatchingMC.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CCpi0ShowerMatchingMC.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    void MinDWall();

    double _dwallmin;

    std::vector<int> Match(const std::vector<larlite::mcshower>& mcs_v,
			   const std::vector<larlite::shower>&   shr_v);

    bool loadVertex(event_vertex* ev_vtx);

    void Reset();

    std::vector<double> _vtx_w_cm, _vtx_t_cm;

    double _w2cm, _t2cm;
    
    TTree* _tree;

    int _n_reco_showers;

    int _event;

    double _nu_e, _pi0_e;
    int _n_trk;
    
    double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
    double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

    double _mc_shr1_x,  _mc_shr1_y,  _mc_shr1_z;
    double _mc_shr1_px, _mc_shr1_py, _mc_shr1_pz;
    double _mc_shr1_e;
    double _mc_shr2_x,  _mc_shr2_y,  _mc_shr2_z;
    double _mc_shr2_px, _mc_shr2_py, _mc_shr2_pz;
    double _mc_shr2_e;
    double _mcradlen1, _mcradlen2;
    
    double _mc_shr_x,  _mc_shr_y,  _mc_shr_z;
    double _mc_shr_px, _mc_shr_py, _mc_shr_pz;
    double _mc_shr_e;
    double _mcradlen;

    double _mc_oangle;

    double _mc_mass;

    // MC -> RC shower comparisons
    double _dot;
    double _strt;
    double _erc;

    // cluster metrics
    double _ip;
    double _lin;
    double _ssv;
    double _slope;
    
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
