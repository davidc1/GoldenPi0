/**
 * \file CCpi0ShowerMatching.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class CCpi0ShowerMatching
 *
 * @author david caratelli
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_CCPI0SHOWERMATCHING_H
#define LARLITE_CCPI0SHOWERMATCHING_H

#include "Analysis/ana_base.h"

#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"

#include "TTree.h"

namespace larlite {
  /**
     \class CCpi0ShowerMatching
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCpi0ShowerMatching : public ana_base{
  
  public:

    /// Default constructor
    CCpi0ShowerMatching()
      : _tree(nullptr)
      { _name="CCpi0ShowerMatching"; _fout=0;}

    /// Default destructor
    virtual ~CCpi0ShowerMatching(){}

    /** IMPLEMENT in CCpi0ShowerMatching.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCpi0ShowerMatching.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CCpi0ShowerMatching.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    std::pair<int,int> Match(const std::vector<larlite::mcshower>& mcs_v,
			     const std::vector<larlite::shower>&   shr_v);

    void Reset();

    TTree* _tree;

    int _n_reco_showers;

    int _event;

    double _nu_e, _pi0_e;
    
    double _mc_vtx_x, _mc_vtx_y, _mc_vtx_z;
    double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

    double _mc_shr1_x,  _mc_shr1_y,  _mc_shr1_z;
    double _mc_shr1_px, _mc_shr1_py, _mc_shr1_pz;
    double _mc_shr1_e;
    double _rc_shr1_x,  _rc_shr1_y,  _rc_shr1_z;
    double _rc_shr1_px, _rc_shr1_py, _rc_shr1_pz;
    double _rc_shr1_e;

    double _mc_shr2_x,  _mc_shr2_y,  _mc_shr2_z;
    double _mc_shr2_px, _mc_shr2_py, _mc_shr2_pz;
    double _mc_shr2_e;
    double _rc_shr2_x,  _rc_shr2_y,  _rc_shr2_z;
    double _rc_shr2_px, _rc_shr2_py, _rc_shr2_pz;
    double _rc_shr2_e;

    double _mc_oangle, _rc_oangle;

    double _rc_mass;

    double _dot1, _dot2;
    
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
