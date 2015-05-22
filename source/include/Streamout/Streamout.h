#ifndef Streamout_h
#define Streamout_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "IMPL/LCGenericObjectImpl.h"
#include <string>
#include <map>
#include "Geometry/Geometry.h"
#include "DIFSlowControl.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/LCTOOLS.h"
#include <cstdint>
#include "marlin/Global.h"


using namespace lcio ;
using namespace marlin ;

class SDHCAL_buffer : public std::pair<uint8_t*, uint32_t>
{
public:
    SDHCAL_buffer(uint8_t* b, uint32_t i) : pair<uint8_t*, uint32_t>(b,i)
    {
        ;
    }
    uint8_t* buffer()
    {
        return first;
    }
    uint8_t* endOfBuffer()
    {
        return first+second;
    }
    uint32_t getsize()
    {
        return second;
    }
    void printBuffer(uint32_t start, uint32_t stop,std::ostream& flux=std::cout);
    void printBuffer(uint32_t start=0,std::ostream& flux=std::cout)
    {
        printBuffer(start,getsize());
    }
};


//From an original class/code by Laurent Mirabito
class LMGeneric: public IMPL::LCGenericObjectImpl
{
public:
    LMGeneric()
    {
        ;
    }
    std::vector<int>& getIntVector()
    {
        return _intVec;
    }
    int* getIntBuffer()
    {
        return _intVec.empty() ? NULL : &_intVec[0];
    }
    uint8_t* getCharBuffer()
    {
        return (uint8_t*) getIntBuffer();
    }
    unsigned int nBytes()
    {
        return getNInt()*sizeof(int32_t);   //4 bytes for each int
    }
    SDHCAL_buffer getSDHCALBuffer()
    {
        return SDHCAL_buffer(getCharBuffer(),nBytes());
    }
};

//class to navigate in the raw data buffer
class SDHCAL_RawBuffer_Navigator
{
public:
    SDHCAL_RawBuffer_Navigator(SDHCAL_buffer b,unsigned int BitsToSkip); //BitsToSkip=92 in 2012, 24 in 2014
    ~SDHCAL_RawBuffer_Navigator()
    {
        if (_theDIFPtr!=NULL) delete _theDIFPtr;
    }
    bool validBuffer()
    {
        return _DIFstartIndex != 0;
    }
    uint32_t getStartOfDIF()
    {
        return _DIFstartIndex;
    }
    uint8_t* getDIFBufferStart()
    {
        return &(_buffer.buffer()[_DIFstartIndex]);
    }
    uint32_t getDIFBufferSize()
    {
        return _buffer.getsize()-_DIFstartIndex;
    }
    SDHCAL_buffer getDIFBuffer()
    {
        return SDHCAL_buffer(getDIFBufferStart(),getDIFBufferSize());
    }
    DIFPtr* getDIFPtr();
    uint32_t getEndOfDIFData()
    {
        return getDIFPtr()->getGetFramePtrReturn()+3;
    }
    uint32_t getSizeAfterDIFPtr()
    {
        return getDIFBufferSize()-getDIFPtr()->getGetFramePtrReturn();
    }
    uint32_t getDIF_CRC();
    bool hasSlowControlData()
    {
        return getDIFBufferStart()[getEndOfDIFData()]==0xb1;
    }
    SDHCAL_buffer getSCBuffer()
    {
        setSCBuffer();
        return _SCbuffer;
    }
    bool badSCData()
    {
        setSCBuffer();
        return _badSCdata;
    }
    SDHCAL_buffer getEndOfAllData();

private:
    void setSCBuffer();
    SDHCAL_buffer _buffer,_SCbuffer;
    uint32_t _DIFstartIndex;
    DIFPtr* _theDIFPtr;
    bool _badSCdata;
};


class Streamout : public Processor
{

public:

    virtual Processor*  newProcessor()
    {
        return new Streamout ;
    }


    Streamout() ;

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
    
    /** Input collection name.
     */
    std::string  _XDAQCollectionNames ;

    /** Output collection name.
     */
    std::string _RawHitCollectionName;
    std::string _RawHitCollectionNameTime;
    std::string _TcherenkovSignal;
    int _TcherenkovSignalDuration;

    /** Parameters
     */
    int _BitsToSkip;
    int _GlobalEvents;
    int _maxRecord;
    int _rolling;
    int _skip;
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
