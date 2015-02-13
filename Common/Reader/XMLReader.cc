#include<iostream>
#include<vector>
#include "marlin/VerbosityLevels.h"
#include "Reader/XMLReader.h"
void XMLReader::Read(std::string FileName, Geometry& geom)
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
    std::vector<int>DifM;
    streamlog_out( MESSAGE )<<"File : "<<FileName<<std::endl;
    TiXmlHandle hdl(&doc);
    TiXmlElement* Platee = hdl.FirstChildElement().FirstChildElement().Element();
    TiXmlElement* Diff=NULL;
    unsigned int PlateNumber=0;
    while (Platee)
    {
      PlateNumber++;
      streamlog_out( MESSAGE ) << Platee->Value() << std::endl;
      double x= atof(Platee->Attribute("x"));
      double y=atof(Platee->Attribute("y"));
      double z=atof(Platee->Attribute("z"));
      double xy=atof(Platee->Attribute("alpha"));
      double xz=atof(Platee->Attribute("beta"));
      double yz=atof(Platee->Attribute("gamma"));
      
	    //std::cout<<x<<y<<z<<xy<<xz<<yz<<std::endl;
      Diff=Platee->FirstChildElement();
      while (Diff)
      {
        int DifType=0;
        int up_down=0;
        if(Diff->Attribute("up_down")!=NULL)
        {
         if(strcmp(Diff->Attribute("up_down"), "up") == 0) up_down=1;
         else if (strcmp(Diff->Attribute("up_down"), "down") == 0) up_down=0;
         else
         {
          
          up_down=-1;std::cout<<"Error defining the position of the Dif (up, down)"<<std::endl;
         }
        }if(Diff->Attribute("DifType")!=NULL)
        {
         if(strcmp(Diff->Attribute("DifType"), "temporal") == 0) DifType=2;
         else if (strcmp(Diff->Attribute("DifType"), "positional") == 0) DifType=1;
	else if (strcmp(Diff->Attribute("DifType"), "tcherenkov") == 0) DifType=3;
         else
         {
          
          DifType=-1;std::cout<<"Error defining the use of the Dif (temporal,posicional)"<<std::endl;
         }
        }
        
              streamlog_out( MESSAGE ) << Diff->Value() << std::endl;
	      std::cout<<Diff->Attribute("DifId")<<Diff->Attribute("x")<<"DifType :  "<<DifType<<" Position :  "<<up_down<<std::endl;
	      DifM.push_back(atof(Diff->Attribute("DifId")));
	      geom.AddDif(atof(Diff->Attribute("x")),atof(Diff->Attribute("y")),atof(Diff->Attribute("DifId")),xy,xz,yz,PlateNumber,up_down,DifType);
	      Diff= Diff->NextSiblingElement();
	    }
	    geom.AddPlate(x,y,z,xy,xz,yz,DifM);
	    DifM.clear();
	    Platee=Platee->NextSiblingElement(); 
    }
	  //std::cout<<"Fini"<<std::endl;
	}
}
