/**
 * \file VertexExample.h
 *
 * \ingroup Pi0Reconstruction
 * 
 * \brief Class def header for a class VertexExample
 *
 * @author davidc1
 */

/** \addtogroup Pi0Reconstruction

    @{*/

#ifndef LARLITE_VERTEXEXAMPLE_H
#define LARLITE_VERTEXEXAMPLE_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class VertexExample
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexExample : public ana_base{
  
  public:

    /// Default constructor
    VertexExample();

    /// Default destructor
    virtual ~VertexExample(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    // function to set the vertex producer name
    void SetVertexProducer(std::string s) { _vertex_producer = s; }

  protected:

    std::string _vertex_producer;

    
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
