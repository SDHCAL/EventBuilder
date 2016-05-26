#ifndef ANALYSIS_H
#define ANALYSIS_H
// marlin
#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "Geometry/Geometry.h"
#include "TestedPlan.h"
#include "TrackFitter.h"
#include "UTIL/CellIDDecoder.h"
#include <map>
#include <array>
#include <string>
#include <vector>
#include "Config/Config.h"
#include "THnSparse.h"
enum Threshold{Threshold_2=1,Threshold_1,Threshold_3};
//class hitsInPlan;


std::map<std::string,std::array<std::map<std::vector<int>,std::vector<int>>,6>>Efficiency_per_pad;
double _Chi2;
double _Chi2Rate;
int _ShortEfficiency;
int _NbrHitPerPlaneMax ;
int _NbrPlaneUseForTracking ;
int _NbrPlaneUseForTrackingRate ;
double _dlimforPad;
int _NbrRun;
double _dlimforStrip;
bool EstimateNoiseContamination;
std::map<int ,std::vector<double> >Delimiter;
std::string _Delimiters;
std::vector<std::string> _hcalCollections;
class AnalysisProcessor : public marlin::Processor
{
public:
    AnalysisProcessor();
    ~AnalysisProcessor();
    void PrintStat(std::string);
    void PrintStatShort();
    void PrintStatShort(std::string);
    void operator=(const AnalysisProcessor&);
    virtual marlin::Processor* newProcessor() {
      return new AnalysisProcessor;
    }
    void init();
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
protected:
 
    LCWriter* _EventWriter;
    LCWriter* _NoiseWriter;
    std::string _FileNameGeometry;
    unsigned int _eventNr;
    unsigned int _skip;
    unsigned int _maxRecord;
    unsigned int _GlobalEvents;
    unsigned int _GlobalEventsSc;
    
    std::string _Config_xml;   
    std::string _Elog_xml;
    ConfigInfos conf;
    Geometry geom;
    std::string _ReaderType;
    std::map<std::string,std::map<int,hitsInPlan>> Planss;
    std::map<int,geometryplan> geometryplans;
    std::vector<testedPlan> testedPlanList;

};
#endif
