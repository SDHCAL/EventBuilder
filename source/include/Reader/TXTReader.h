#ifndef TXTREADER_H
#define TXTREADER_H
#include"Reader.h"
#include<iostream>
class TXTReader: public Reader
{
  public:
    void Read(std::string& FileName,Geometry& geom);
    void Read(std::string& FileName,ConfigInfos& config){};
    Reader *Clone() { return new TXTReader(); }
    ~TXTReader(){};
};
#endif
