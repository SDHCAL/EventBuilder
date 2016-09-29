#ifndef READER_H
#define READER_H
#include<string>
#include"Geometry.h"
#include"Config.h"
class Reader
{
  public:
  virtual void Read(std::string& FileName,Geometry& geom) = 0;
  virtual void Read(std::string& FileName,ConfigInfos& config)=0;
  virtual void Read(std::string& FileName,ConfigInfos& config, unsigned int RunNumber)=0;
  virtual Reader *Clone() = 0;
  virtual ~Reader(){};
};
#endif
