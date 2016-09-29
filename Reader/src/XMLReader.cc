#include<iostream>
#include<vector>
#include "VerbosityLevels.h"
#include "XMLReader.h"
#include "Colors.h"
#include"Types.h"

void XMLReader::Write(TiXmlElement* elem,const char* who,unsigned int &wherePlate,unsigned int &whereDif,double& var,std::string &FileName)
{
	if(elem->Attribute(who)!=NULL) var=atof(elem->Attribute(who));
  else
  {
		if(whereDif>0)std::cout<<red<<who<<" argument mising in Plate : "<<wherePlate<<" Dif number : "<<whereDif<<" !!!"<<normal<<std::endl;  
    else std::cout<<red<<who<<" argument mising in Plate : "<<wherePlate<<" !!!"<<normal<<std::endl;  
    std::cout<<red<<"Please check the file : "<<FileName<<normal<<std::endl;
		std::exit(1);	
	}
}

void XMLReader::Read(std::string& FileName, Geometry& geom)
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
    TiXmlElement* Diff=nullptr;
    unsigned int PlateNumber=0;
    unsigned int DifNumber=0;
    double x,y,z,xy,xz,yz,SizeX,SizeY,I,J,DifId;
    double glass_type=-1;
    double Gaz_number=-1;
    std::string HV_channel="";
    std::string Gaz_channel="";
    while (Platee)
    {
      DifNumber=0;
      PlateNumber++;
      Write(Platee,"x",PlateNumber,DifNumber,x,FileName);
      Write(Platee,"y",PlateNumber,DifNumber,y,FileName);
      Write(Platee,"z",PlateNumber,DifNumber,z,FileName);
      Write(Platee,"alpha",PlateNumber,DifNumber,xy,FileName);
      Write(Platee,"beta",PlateNumber,DifNumber,xz,FileName);
      Write(Platee,"gamma",PlateNumber,DifNumber,yz,FileName);
      Write(Platee,"SizeX",PlateNumber,DifNumber,SizeX,FileName);
      Write(Platee,"SizeY",PlateNumber,DifNumber,SizeY,FileName);
      if(Platee->Attribute("Glass_type")!=NULL)
	    {
	        if(strcmp(Platee->Attribute("Glass_type"), "standard") == 0) glass_type=0;
	        else if (strcmp(Platee->Attribute("Glass_type"), "chinese") == 0) glass_type=1;
	        else if (strcmp(Platee->Attribute("Glass_type"), "") == 0) glass_type=-1;
	        else
		      {
		        glass_type=-1;
		        std::cout<<"Error defining Glasse_type (standard,chinese)"<<std::endl;
		      }
	    }
	    else glass_type=-1;
	    if(Platee->Attribute("Gaz_number")!=NULL)
	    {
	      std::string e=Platee->Attribute("Gaz_number");
	      if(e.size()!=0) Gaz_number=atof(Platee->Attribute("Gaz_number"));
	    }
	    else Gaz_number=-1;
	    if(Platee->Attribute("HV_channel")!=NULL) HV_channel=Platee->Attribute("HV_channel");
	    else HV_channel="";
	    if(Platee->Attribute("Gaz_channel")!=NULL) Gaz_channel=Platee->Attribute("Gaz_channel");
	    else Gaz_channel="";
      Diff=Platee->FirstChildElement();
      while (Diff)
      {
        DifNumber++;
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
	      }
	      if(Diff->Attribute("DifType")!=NULL)
	      {
	        if(strcmp(Diff->Attribute("DifType"), "temporal") == 0) DifType=temporal;
	        else if (strcmp(Diff->Attribute("DifType"), "positional") == 0) DifType=positional;
	        else if (strcmp(Diff->Attribute("DifType"), "tcherenkov") == 0) DifType=tcherenkov;
	        else if (strcmp(Diff->Attribute("DifType"), "tricot") == 0) DifType=tricot;
          else if (strcmp(Diff->Attribute("DifType"), "scintillator") == 0) DifType=scintillator;
          else if (strcmp(Diff->Attribute("DifType"), "pad") == 0) DifType=pad;
	        else
		      {
		        DifType=-1;
		        std::cout<<"Error defining the use of the Dif (temporal,posicional,tcherenkov,tricot)"<<std::endl;
		      }
	      }
        Write(Diff,"DifId",PlateNumber,DifNumber,DifId,FileName);
        Write(Diff,"I",PlateNumber,DifNumber,I,FileName);
        Write(Diff,"J",PlateNumber,DifNumber,J,FileName);
	      DifM.push_back(atof(Diff->Attribute("DifId")));
	      geom.AddDif(I,J,DifId,xy,xz,yz,PlateNumber,up_down,DifType);
	      Diff= Diff->NextSiblingElement();
      }
      geom.AddPlate(x,y,z,xy,xz,yz,DifM,SizeX,SizeY,glass_type,Gaz_number,HV_channel,Gaz_channel);
      DifM.clear();
      Platee=Platee->NextSiblingElement(); 
    }
  }
}
