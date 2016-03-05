#ifndef XMLREADERCONF_H
#define XMLREADERCONF_H
#include"Reader.h"
#include<iostream>
#include<string>
#include "tinyxml.h"

class XMLReaderConfig: public Reader
{
  public:
  void Read(std::string &FileName,ConfigInfos& config);
  void Read(std::string& FileName,Geometry& geom){};
  void Read(std::string& FileName,ConfigInfos& config, unsigned int RunNumber){};
  Reader *Clone() { return new XMLReaderConfig(); }
  ~XMLReaderConfig(){};
};
#endif
