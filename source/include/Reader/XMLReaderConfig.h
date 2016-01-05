#ifndef XMLREADERCONF_H
#define XMLREADERCONF_H
#include"Reader.h"
#include<iostream>
#include<string>
#include "marlin/tinyxml.h"

class XMLReaderConfig: public Reader
{
  public:
  void Read(std::string &FileName,ConfigInfos& config);
  void Read(std::string& FileName,Geometry& geom){};
  Reader *Clone() { return new XMLReaderConfig(); }
  ~XMLReaderConfig(){};
};
#endif
