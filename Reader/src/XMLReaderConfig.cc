#include<iostream>
#include<vector>
#include "VerbosityLevels.h"
#include "XMLReaderConfig.h"
#include "Colors.h"
#include "Types.h"
#include <string>
#include <iostream>
#include <sstream>
#include <utility>

std::array<unsigned int,64> ParsePAGAIN(std::string& pagain)
{
      std::array<unsigned int,64>arr;
      std::stringstream iss(pagain);
      unsigned int gain=0;
      unsigned short index_pad=0;
      while(iss>>gain)
      {
        //std::cout<<(unsigned int )(gain)<<std::endl;
        arr[index_pad]=gain;
        //std::cout<<index_pad<<std::endl;
        index_pad++;
      }
      return arr;
}



void XMLReaderConfig::Read(std::string &FileName,ConfigInfos& Conf)
{
 TiXmlDocument doc(FileName.c_str());
  doc.LoadFile();
  if(!doc.LoadFile())
  {
    streamlog_out( WARNING )<<"Error loading file"<< std::endl;
    streamlog_out( WARNING )<<"Error #" << doc.ErrorId() << " : " << doc.ErrorDesc()<<std::endl;
  }
  else
  {
    DifDatabaseInfo dif =DifDatabaseInfo();
    streamlog_out( MESSAGE )<<"File : "<<FileName<<std::endl;
    TiXmlHandle hdl(&doc);
    TiXmlElement*  element = hdl.FirstChildElement().FirstChildElement().Element();
    while (element)
    {
      std::string name=element->Value();
      if(name=="DIF")
      {
        std::string DIF_TYPE=element->FirstChildElement( "DIF_TYPE" )->GetText();
        std::string NAME=element->FirstChildElement( "NAME" )->GetText();
        unsigned int ID=std::stoi(element->FirstChildElement( "ID" )->GetText());
        bool ENABLED=bool(element->FirstChildElement( "ENABLED" )->GetText()) ;
        Conf.AddDif(DifDatabaseInfo(NAME,DIF_TYPE,ID,ENABLED));
      }
      if(name=="ASIC")
      {
        std::string ASIC_TYPE=element->FirstChildElement( "ASIC_TYPE" )->GetText();
        unsigned int B0=std::stoi(element->FirstChildElement( "B0" )->GetText());
        unsigned int B1=std::stoi(element->FirstChildElement( "B1" )->GetText());
        unsigned int B2=std::stoi(element->FirstChildElement( "B2" )->GetText());
        unsigned int DIF_ID=std::stoi(element->FirstChildElement( "DIF_ID" )->GetText());
        unsigned int ID=std::stoi(element->FirstChildElement( "ID" )->GetText());
        unsigned int HEADER=std::stoi(element->FirstChildElement( "HEADER" )->GetText());
        unsigned int NUM=std::stoi(element->FirstChildElement( "NUM" )->GetText());
        std::string PAGAIN=element->FirstChildElement( "PAGAIN" )->GetText();  
        Conf.AddAsic(AsicInfo(ASIC_TYPE,ParsePAGAIN(PAGAIN),B0,B1,B2,DIF_ID,ID,HEADER,NUM));
      }
      element=element->NextSiblingElement(); 
    }
  }
}
