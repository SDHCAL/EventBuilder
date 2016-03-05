#ifndef XMLREADER_H
#define XMLREADER_H
#include"Reader.h"
#include<iostream>
#include<string>
#include "tinyxml.h"


class XMLReader: public Reader
{
  public:
  void Read(std::string &FileName,Geometry& geom);
  void Read(std::string &FileName,ConfigInfos& conf){};
  void Read(std::string& FileName,ConfigInfos& config, unsigned int RunNumber){};
  void Write(TiXmlElement*,const char*,unsigned int&, unsigned int&,double& ,std::string& FileName);
  Reader *Clone() { return new XMLReader(); }
  ~XMLReader(){};
};
#endif
