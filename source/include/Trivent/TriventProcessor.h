#ifndef TRIVENT_PROCESSOR
#define TRIVENT_PROCESSOR

#include "marlin/Processor.h"
#include "EVENT/RawCalorimeterHit.h"
#include "Geometry/Geometry.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/CellIDEncoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <EVENT/LCRunHeader.h>
//std::string _Delimiters;
//std::map<int ,std::vector<double> >Delimiter;
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
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
    void FillTimes();
    void FillIJK(std::vector<RawCalorimeterHit *>vec, LCCollectionVec* col,CellIDEncoder<CalorimeterHitImpl>& cd,bool IsNoise);
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
    int _timeWin;
    int _eventNr;
    int _trig_count;
    int _noiseCut;
    int _LayerCut;
    int _TriggerTime;
    std::map<int,int>Times;
    std::vector<std::map<int,int> >Times_Plates;
    std::vector<std::map<int,int> >Times_Plates_perRun;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawHits;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > BehondTrigger;
    std::vector<RawCalorimeterHit *>EventsGrouped;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawTimeDifs;
    float pos[3];
};



#endif
