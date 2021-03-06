#include "Streamout/BufferNavigator.h"
#include <iostream>
#include <cstdint>
#include "Streamout/DIFSlowControl.h"
using namespace std;

void SDHCAL_buffer::printBuffer(unsigned int start, unsigned int stop,std::ostream& flux)
{
    flux << std::hex;
    for (unsigned int k=start; k<stop; k++) flux << (unsigned int)(first[k]) << " - ";
    flux << std::dec <<  std::endl;
}

SDHCAL_RawBuffer_Navigator::SDHCAL_RawBuffer_Navigator(SDHCAL_buffer b, unsigned int BitsToSkip) :_buffer(b),_SCbuffer(0,0)
{
    _DIFstartIndex=DIFUnpacker::getStartOfDIF(_buffer.buffer(),_buffer.getsize(),BitsToSkip); //92 was here
    _theDIFPtr=NULL;
    _badSCdata=false;
}

DIFPtr* SDHCAL_RawBuffer_Navigator::getDIFPtr()
{
    if (NULL==_theDIFPtr) _theDIFPtr=new DIFPtr(getDIFBufferStart(),getDIFBufferSize());
    return _theDIFPtr;
}

uint32_t SDHCAL_RawBuffer_Navigator::getDIF_CRC()
{
    uint32_t i=getEndOfDIFData();
    uint32_t ret=0;
    ret |= ( (_buffer.buffer()[i-2])<<8 );
    ret |= _buffer.buffer()[i-1];
    return ret;
}

void SDHCAL_RawBuffer_Navigator::setSCBuffer()
{
    if (! hasSlowControlData() ) return;
    if (_SCbuffer.getsize()!=0 ) return; //deja fait
    if (_badSCdata) return;
    _SCbuffer.first=&(getDIFBufferStart()[getEndOfDIFData()]);
    //compute Slow Control size
    uint32_t maxsize=_buffer.getsize()-_DIFstartIndex-getEndOfDIFData()+1; // should I +1 here ?
    uint32_t k=1; //SC Header
    uint32_t dif_ID=_SCbuffer.first[1];
    uint32_t chipSize=_SCbuffer.first[3];
    while ((dif_ID != 0xa1 && _SCbuffer.first[k] != 0xa1 && k <maxsize) ||(dif_ID == 0xa1 && _SCbuffer.first[k+2]==chipSize && k<maxsize)) 
    {
      k+=2; //DIF ID + ASIC Header
      uint32_t scsize=_SCbuffer.first[k];
      if (scsize != 74 && scsize != 109) 
      {
        std::cout << "PROBLEM WITH SC SIZE " << scsize << std::endl;
        k=0;
        _badSCdata=true;
        break;
      }
      k++; //skip size bit
      k+=scsize; // skip the data
    }
    if (_SCbuffer.first[k] == 0xa1 && !_badSCdata ) _SCbuffer.second=k+1; //add the trailer
    else 
    {
        _badSCdata=true;
        std::cout << "PROBLEM SC TRAILER NOT FOUND " << std::endl;
    }
}

SDHCAL_buffer SDHCAL_RawBuffer_Navigator::getEndOfAllData()
{
    setSCBuffer();
    if (hasSlowControlData() && !_badSCdata) {
        return SDHCAL_buffer( &(_SCbuffer.buffer()[_SCbuffer.getsize()]), getSizeAfterDIFPtr()-3-_SCbuffer.getsize() );
    } else return SDHCAL_buffer( &(getDIFBufferStart()[getEndOfDIFData()]), getSizeAfterDIFPtr()-3 ); //remove the 2 bytes for CRC and the DIF trailer
}
