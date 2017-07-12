/**
 * \file CCtrue.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class CCtrue
 *
 * @author david caratelli
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_CCTRUE_H
#define LARLITE_CCTRUE_H

#include "Analysis/ana_base.h"

#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"

#include "TTree.h"

namespace larlite {
  /**
     \class CCtrue
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCtrue : public ana_base{
  
  public:

    /// Default constructor
    CCtrue()
      : _tree(nullptr)
      { _name="CCtrue"; _fout=0;}

    /// Default destructor
    virtual ~CCtrue(){}

    /** IMPLEMENT in CCtrue.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCtrue.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);


    /** IMPLEMENT in CCtrue.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    double _dwallmin;

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

    double _rcradlen1, _rcradlen2, _mcradlen1, _mcradlen2;

    double _mc_mass, _rc_mass;

    double _dot1, _dot2;
    double _strt1, _strt2;

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
