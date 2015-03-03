#ifndef TXTREADER_H
#define TXTREADER_H
#include"Reader.h"
#include<iostream>
class TXTReader: public Reader
{
  public:
    void Read(std::string& FileName,Geometry& geom);
    Reader *Clone() { return new TXTReader(); }
    ~TXTReader(){};
};
#endif
