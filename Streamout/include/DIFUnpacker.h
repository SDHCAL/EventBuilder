#ifndef DIFUnpacker_h
#define DIFUnpacker_h
#include <string>
#include <vector>
#include <cstdint>

#define DU_DATA_FORMAT_VERSION 13

#define DU_START_OF_DIF      0xB0
#define DU_START_OF_DIF_TEMP 0xBB
#define DU_END_OF_DIF        0xA0
#define DU_START_OF_LINES    0xC4
#define DU_END_OF_LINES      0xD4

#define DU_START_OF_FRAME    0xB4
#define DU_END_OF_FRAME      0xA3


#define DU_ID_SHIFT      1
#define DU_DTC_SHIFT     2
#define DU_GTC_SHIFT    10
#define DU_ABCID_SHIFT  14
#define DU_BCID_SHIFT   20
#define DU_LINES_SHIFT  23
#define DU_TASU1_SHIFT  24
#define DU_TASU2_SHIFT  28
#define DU_TDIF_SHIFT   32

#define DU_FRAME_ASIC_HEADER_SHIFT 0
#define DU_FRAME_BCID_SHIFT        1
#define DU_FRAME_DATA_SHIFT        4
#define DU_FRAME_SIZE              20

class DIFUnpacker {
public:
  static unsigned long long GrayToBin(unsigned long long n);
  static uint32_t getStartOfDIF(uint8_t* cbuf,uint32_t size_buf,uint32_t start);
  static uint32_t getID(uint8_t* cb,uint32_t idx=0);
  static uint32_t getDTC(uint8_t* cb,uint32_t idx=0);
  static uint32_t getGTC(uint8_t* cb,uint32_t idx=0);
  static unsigned long long getAbsoluteBCID(uint8_t* cb,uint32_t idx=0);
  static uint32_t getBCID(uint8_t* cb,uint32_t idx=0);
  static uint32_t getLines(uint8_t* cb,uint32_t idx=0);
  static bool hasLine(uint32_t line,uint8_t* cb,uint32_t idx=0);
  static uint32_t getTASU1(uint8_t* cb,uint32_t idx=0);
  static uint32_t getTASU2(uint8_t* cb,uint32_t idx=0);
  static uint32_t getTDIF(uint8_t* cb,uint32_t idx=0);
  static bool hasTemperature(uint8_t* cb,uint32_t idx=0);
  static bool hasAnalogReadout(uint8_t* cb,uint32_t idx=0);
  static uint8_t getFrameData(uint8_t* framePtr,uint32_t ip);
  static uint32_t getFrameAsicHeader(uint8_t* framePtr);
  static uint32_t getFrameBCID(uint8_t* framePtr);
  static bool getFramePAD(uint8_t* framePtr,uint32_t ip);
  static bool getFrameLevel(uint8_t* framePtr,uint32_t ip,uint32_t level);
  static uint32_t getAnalogPtr(std::vector<uint8_t*> &vLines,uint8_t* cb,uint32_t idx=0) throw (std::string);
  static uint32_t getFramePtr(std::vector<uint8_t*> &vFrame,std::vector<uint8_t*> &vLines,uint32_t max_size,uint8_t* cb,uint32_t idx=0) throw (std::string);
  static void dumpFrameOld(uint8_t* buf);
  static unsigned long swap_bytes(unsigned int n,uint8_t* buf); //Stolen from DCBufferReader
};
#endif
