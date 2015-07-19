#ifndef TRIVENT_PROCESSOR
#define TRIVENT_PROCESSOR
#include <map>
#include <vector>
#include "marlin/Processor.h"
#include "EVENT/RawCalorimeterHit.h"
#include "Geometry/Geometry.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/CellIDEncoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <EVENT/LCRunHeader.h>
#include "UTIL/LCTOOLS.h"
#include "Trivent/HistoPlane.h"
#include <set>
#include <string>


class TriventProcessor : public marlin::Processor
{
public:
    TriventProcessor();
    ~TriventProcessor();
    marlin::Processor *newProcessor()
    {
        return new TriventProcessor();
    }
    void init();
    void Writer(IO::LCWriter* file,const char* name,std::map<int,std::vector<EVENT::RawCalorimeterHit *>> & vec,EVENT::LCEvent * event,unsigned int& nbr,unsigned int IsNoise);
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
    inline void FillTimes();
    void FillIJK(std::vector<RawCalorimeterHit *>vec, LCCollectionVec* col,CellIDEncoder<CalorimeterHitImpl>& cd,int IsNoise);
    void FillIJK(std::vector<RawCalorimeterHit *>vec);
private:
    void processCollection(EVENT::LCEvent *evtP,LCCollection* col);    
protected:
    Geometry geom;
    std::map<int,std::vector<double>>SinCos;
    std::string _FileNameGeometry;
    std::string _ReaderType;
    std::string _noiseFileName;
    std::string _outFileName;
    std::vector<std::string> _hcalCollections;
    LCWriter* _EventWriter;
    LCWriter* _NoiseWriter;
    unsigned int TouchedEvents;
    int _timeWin;
    int _eventNr;
    int _trig_count;
    int _noiseCut;
    int _LayerCut;
    int _TriggerTime;
    int _GlobalEvents;
    std::string _Delimiters;
    unsigned int _maxRecord;
    unsigned int eventtotal;
    unsigned int _rolling;
    int _skip;
    double _efficiencyFrontScintillator;
    int _IgnorebeginningSpill;
    double _efficiencyBackScintillator;
    bool _WantDistribution;
    bool _WantCalibration;
    bool _Spill_Study;
    std::string _Database_name;
    std::map<int,int>Times;
    unsigned int _Front_scintillator;
    unsigned int _Back_scintillator;
    unsigned int _Both_scintillator;
    //std::vector<std::map<int,int> >Times_Plates_perRun;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawHits;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > BehondTrigger;
    std::vector<RawCalorimeterHit *>EventsGrouped;
    std::vector<RawCalorimeterHit *>EventsGroupedScin;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawTimeDifs;
    float pos[3];
    std::map<int,bool>Warningg;
    std::map<std::vector<unsigned int>,std::map< int, int>>Negative;
    std::map<int,HistoPlane*>HistoPlanes;
    void save_calibration(std::string filename);
    void read_calibration(std::string filename);
    std::map<int ,std::vector<double> >Delimiter;
};


#endif
