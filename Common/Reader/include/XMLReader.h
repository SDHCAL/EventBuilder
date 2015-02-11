#ifndef XMLREADER_H
#define XMLREADER_H
#include"Reader.h"
#include<iostream>
#include<string>
#include "marlin/tinyxml.h"

class XMLReader: public Reader
{
  public:
  void Read(std::string FileName,Geometry& geom);
  Reader *Clone() { return new XMLReader(); }
  ~XMLReader(){};
};
#endif
