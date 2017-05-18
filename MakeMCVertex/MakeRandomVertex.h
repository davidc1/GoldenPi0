/**
 * \file MakeRandomVertex.h
 *
 * \ingroup MakeMCVertex
 * 
 * \brief Class def header for a class MakeRandomVertex
 *
 * @author david
 */

/** \addtogroup MakeMCVertex

    @{*/

#ifndef LARLITE_MAKERANDOMVERTEX_H
#define LARLITE_MAKERANDOMVERTEX_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MakeRandomVertex
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeRandomVertex : public ana_base{
  
  public:

    /// Default constructor
    MakeRandomVertex()
      : _tree(nullptr)
      { _name="MakeRandomVertex"; _fout=0;}

    /// Default destructor
    virtual ~MakeRandomVertex(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    TTree *_tree;
    double _x,_y,_z;
    
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
