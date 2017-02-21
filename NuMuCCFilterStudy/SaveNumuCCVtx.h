/**
 * \file SaveNumuCCVtx.h
 *
 * \ingroup Playground
 * 
 * \brief Class def header for a class SaveNumuCCVtx
 *
 * @author davidc1
 */

/** \addtogroup Playground

    @{*/

#ifndef LARLITE_SAVENUMUCCVTX_H
#define LARLITE_SAVENUMUCCVTX_H

#include "Analysis/ana_base.h"
#include <iostream>
#include <fstream>
#include <utility>

namespace larlite {
  /**
     \class SaveNumuCCVtx
     User custom analysis class made by SHELL_USER_NAME
   */
  class SaveNumuCCVtx : public ana_base{
  
  public:

    /// Default constructor
    SaveNumuCCVtx(){ _name="SaveNumuCCVtx"; _fout=0; _vtxmap.clear(); }

    /// Default destructor
    virtual ~SaveNumuCCVtx(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void addVertex(const int& run, const int& event, const double& x, const double& y, const double& z);

  protected:
    
    std::ofstream event_file;

    std::map< std::pair<int,int> , std::vector<double> > _vtxmap;
    
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
