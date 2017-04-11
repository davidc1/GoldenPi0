/**
 * \file ShowerStudy.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class ShowerStudy
 *
 * @author david
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_SHOWERSTUDY_H
#define LARLITE_SHOWERSTUDY_H

#include "Analysis/ana_base.h"

#include "TTree.h"

namespace larlite {
  /**
     \class ShowerStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerStudy : public ana_base{
  
  public:

    /// Default constructor
    ShowerStudy()
      : _tree(nullptr)
      { _name="ShowerStudy"; _fout=0;}

    /// Default destructor
    virtual ~ShowerStudy(){}

    /** IMPLEMENT in ShowerStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

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
