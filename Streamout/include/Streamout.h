#ifndef Streamout_h
#define Streamout_h
#include "Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include "Geometry.h"
#include "Config.h"
#include "DIFSlowControl.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/LCTOOLS.h"
#include <cstdint>
#include "Global.h"
#include "BufferNavigator.h"


using namespace lcio ;
using namespace marlin ;



class Streamout : public Processor
{
  public:
  
    virtual Processor*  newProcessor()
    {
      return new Streamout ;
    }
  
    Streamout() ;
    ~Streamout(){};

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init() ;

    /** Called for every run.
    */
    virtual void processRunHeader( LCRunHeader* run ) ;

    /** Called for every event - the working horse.
    */
    virtual void processEvent( LCEvent * evt ) ;


    /** Called after data processing for clean up.
    */
    virtual void end() ;


  private:
    Geometry geom;
    IMPL::LCFlagImpl chFlag;
    std::string _FileNameGeometry;
    std::string _ReaderType;
    /** Flags to DEBUG : use it with care : this is plan to crash may crash the application (use of assert)
    */
    bool _debugMode;
    unsigned int _eventNr;
    unsigned int _NbrRun=0;
    /** Input collection name.
     */
    std::string  _XDAQCollectionNames ;

    /** Output collection name.
     */
    std::string _RawHitCollectionName;
    std::string _RawHitCollectionNameTime;
    std::string _TcherenkovSignal;
    std::string _ScintillatorSignal;
    int _TcherenkovSignalDuration;

    /** Parameters
     */
    int _BitsToSkip;
    unsigned int _GlobalEvents;
    unsigned int _maxRecord;
    //int _rolling;
    unsigned int _skip;
    //statistical counters
    unsigned int _nevt;
    std::map<int,int> _CollectionSizeCounter;
    int _nWrongObj,_nProcessedObject, _hasSlowControl, _hasBadSlowControl;
    std::map<int,int> _DIFStarter;
    std::map<int,int> _DIFPtrValueAtReturnedPos;
    std::map<int,int> _SizeAfterDIFPtr;
    std::map<int,int> _SizeAfterAllData;
    std::map<int,int> _NonZeroValusAtEndOfData;

    void printCounter(std::string description, std::map<int,int> &);
} ;
#endif
