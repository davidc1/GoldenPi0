/**
 * \file Pi0HitThresholdStudies.h
 *
 * \ingroup Pi0SimchStudies
 * 
 * \brief Class def header for a class Pi0HitThresholdStudies
 *
 * @author dcaratelli
 */

/** \addtogroup Pi0SimchStudies

    @{*/

#ifndef LARLITE_PI0HITTHRESHOLDSTUDIES_H
#define LARLITE_PI0HITTHRESHOLDSTUDIES_H

#include "Analysis/ana_base.h"

#include "MCComp/MCBTAlg.h"

#include "TTree.h"

#include <map>

namespace larlite {
  /**
     \class Pi0HitThresholdStudies
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi0HitThresholdStudies : public ana_base{
  
  public:

    /// Default constructor
    Pi0HitThresholdStudies()
      : _tree(nullptr)
      , _hit_tree(nullptr)
      { _name="Pi0HitThresholdStudies"; _fout=0;}

    /// Default destructor
    virtual ~Pi0HitThresholdStudies(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setClusterProducer(std::string s) { _cluster_producer = s; }
    void SaveClusters(bool on) { _save_clusters = on; }
    void AvoidDuplicateHits(bool on) { _avoid_duplicate_ticks = on; }

  protected:

    std::pair<int,int> getTimeSubset(const int& ch, const int& tstart, const int& tend);

    bool _save_clusters;
    
    bool _avoid_duplicate_ticks;

    // hit brack-tracking tool
    ::btutil::MCBTAlg _bt_algo;

    // cluster producer to be used
    std::string _cluster_producer;

    TTree* _tree;
    double _etrue0, _edep0, _ehit0, _qhit0;
    double _etrue1, _edep1, _ehit1, _qhit1;
    double _angle;
    
    TTree* _hit_tree;
    int    _ch;
    double _q, _eall, _eshr, _adc;

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
