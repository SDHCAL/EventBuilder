#ifndef BUFFER_NAVIGATOR_H
#define BUFFER_NAVIGATOR_H 1
#include <utility>
#include <cstdint>
#include "DIFSlowControl.h"

using namespace std;
class SDHCAL_buffer : public std::pair<uint8_t*, uint32_t>
{
    public:
    SDHCAL_buffer(uint8_t* b, uint32_t i) : pair<uint8_t*, uint32_t>(b,i){}
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

#endif
