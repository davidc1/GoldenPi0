/**
 * \file Pi0dEdx.h
 *
 * \ingroup Pi0SimchStudies
 * 
 * \brief Class def header for a class Pi0dEdx
 *
 * @author dcaratelli
 */

/** \addtogroup Pi0SimchStudies

    @{*/

#ifndef LARLITE_PI0DEDX_H
#define LARLITE_PI0DEDX_H

#include "Analysis/ana_base.h"

#include "MCComp/MCBTAlg.h"

#include "TTree.h"

#include <map>

namespace larlite {
  /**
     \class Pi0dEdx
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi0dEdx : public ana_base{
  
  public:

    /// Default constructor
    Pi0dEdx()
      : _tree(nullptr)
      { _name="Pi0dEdx"; _fout=0;}

    /// Default destructor
    virtual ~Pi0dEdx(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setClusterProducer(std::string s) { _cluster_producer = s; }
    void SaveClusters(bool on) { _save_clusters = on; }
    void AvoidDuplicateHits(bool on) { _avoid_duplicate_ticks = on; }
    void setDMAX(double d) { _dmax = d; }

  protected:

    std::pair<int,int> getTimeSubset(const int& ch, const int& tstart, const int& tend);

    bool _save_clusters;
    
    bool _avoid_duplicate_ticks;

    // hit brack-tracking tool
    ::btutil::MCBTAlg _bt_algo;

    // cluster producer to be used
    std::string _cluster_producer;

    double _dmax; // max integration distance for dE/dx

    double _w2cm, _t2cm;

    TTree* _tree;
    int    _event;
    double _etrue0, _edep0, _qcol0, _ehit0, _qhit0, _ahit0;
    double _etrue1, _edep1, _qcol1, _ehit1, _qhit1, _ahit1;
    double _e00, _e01, _a0, _e10, _e11, _a1;
    double _angle;
    double _dedx0, _dedx1;
    std::vector<double> _truededx0_v, _truededx1_v;
    std::vector<double> _recodedx0_v, _recodedx1_v;
    std::vector<double> _adcqdedx0_v, _adcqdedx1_v;
    double _mcsdedx0, _mcsdedx1; // dedx from shower
    double _truededx0, _truededx1;
    double _pitch0, _pitch1;
    int    _wmin0, _wmax0, _wmin1, _wmax1;
    std::string _process0, _process1;
    
    std::vector< std::vector<std::pair<int,int> > > _chTickMap;

    // map linking channel to ide.Energy for each channel
    std::vector< std::map< int, double> > _chIDEmap;
    
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
