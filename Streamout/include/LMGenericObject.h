#ifndef LM_GENERIC_OBJECT_H
#define LM_GENERIC_OBJECT_H 1 1
#include "IMPL/LCGenericObjectImpl.h"
#include "BufferNavigator.h"
//From an original class/code by Laurent Mirabito
class LMGeneric: public IMPL::LCGenericObjectImpl
{
public:
    LMGeneric(){};
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
#endif
