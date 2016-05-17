/**
 * \file SearchPFPartHierarchy.h
 *
 * \ingroup NuMuCCFilterStudy
 * 
 * \brief Class def header for a class SearchPFPartHierarchy
 *
 * @author david
 */

/** \addtogroup NuMuCCFilterStudy

    @{*/

#ifndef LARLITE_SEARCHPFPARTHIERARCHY_H
#define LARLITE_SEARCHPFPARTHIERARCHY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SearchPFPartHierarchy
     User custom analysis class made by SHELL_USER_NAME
   */
  class SearchPFPartHierarchy : public ana_base{
  
  public:

    /// Default constructor
    SearchPFPartHierarchy(){ _name="SearchPFPartHierarchy"; _fout=0;}

    /// Default destructor
    virtual ~SearchPFPartHierarchy(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:
    
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
