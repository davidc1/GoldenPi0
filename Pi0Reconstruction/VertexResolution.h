/**
 * \file VertexResolution.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class VertexResolution
 *
 * @author davidc1
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_VERTEXRESOLUTION_H
#define LARLITE_VERTEXRESOLUTION_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class VertexResolution
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexResolution : public ana_base{
  
  public:

    /// Default constructor
    VertexResolution();

    /// Default destructor
    virtual ~VertexResolution(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setResolutionFilter(double d) { _dmax = d; }
    

  protected:

    TTree* _tree;
    double _rc_x, _rc_y, _rc_z;
    double _mc_x, _mc_y, _mc_z;
    double _dist;

    double _dmax;
    
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
