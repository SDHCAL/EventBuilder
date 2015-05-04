#include<iostream>
#include<vector>
#include "marlin/VerbosityLevels.h"
#include "Reader/XMLReader.h"
#include "Colors.h"

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

void XMLReader::Read(std::string &FileName, Geometry& geom)
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
    unsigned int DifNumber=0;
    double x,y,z,xy,xz,yz,SizeX,SizeY,I,J,DifId;
    while (Platee)
    {
      DifNumber=0;
      PlateNumber++;
      //streamlog_out( MESSAGE ) << Platee->Value() << std::endl;
      Write(Platee,"x",PlateNumber,DifNumber,x,FileName);
      Write(Platee,"y",PlateNumber,DifNumber,y,FileName);
      Write(Platee,"z",PlateNumber,DifNumber,z,FileName);
      Write(Platee,"alpha",PlateNumber,DifNumber,xy,FileName);
      Write(Platee,"beta",PlateNumber,DifNumber,xz,FileName);
      Write(Platee,"gamma",PlateNumber,DifNumber,yz,FileName);
      Write(Platee,"SizeX",PlateNumber,DifNumber,SizeX,FileName);
      Write(Platee,"SizeY",PlateNumber,DifNumber,SizeY,FileName);
      /*double y=atof(Platee->Attribute("y"));
      double z=atof(Platee->Attribute("z"));
      double xy=atof(Platee->Attribute("alpha"));
      double xz=atof(Platee->Attribute("beta"));
      double yz=atof(Platee->Attribute("gamma"));
      double SizeX=atof(Platee->Attribute("SizeX"));
      double SizeY=atof(Platee->Attribute("SizeY"));*/
      
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
	  }if(Diff->Attribute("DifType")!=NULL)
	     {
	       if(strcmp(Diff->Attribute("DifType"), "temporal") == 0) DifType=2;
	       else if (strcmp(Diff->Attribute("DifType"), "positional") == 0) DifType=1;
	       else if (strcmp(Diff->Attribute("DifType"), "tcherenkov") == 0) DifType=3;
	       else if (strcmp(Diff->Attribute("DifType"), "tricot") == 0) DifType=4;
	       else
		 {
		   DifType=-1;std::cout<<"Error defining the use of the Dif (temporal,posicional,tcherenkov,tricot)"<<std::endl;
		 }
	     }
        
	//streamlog_out( MESSAGE ) << Diff->Value() << std::endl;
        Write(Diff,"DifId",PlateNumber,DifNumber,DifId,FileName);
        Write(Diff,"I",PlateNumber,DifNumber,I,FileName);
        Write(Diff,"J",PlateNumber,DifNumber,J,FileName);
	DifM.push_back(atof(Diff->Attribute("DifId")));
        /*double I=atof(Diff->Attribute("I"));
        double J=atof(Diff->Attribute("J"));
        double DifId=atof(Diff->Attribute("DifId"));*/
	geom.AddDif(I,J,DifId,xy,xz,yz,PlateNumber,up_down,DifType);
	Diff= Diff->NextSiblingElement();
      }
      geom.AddPlate(x,y,z,xy,xz,yz,DifM,SizeX,SizeY);
      DifM.clear();
      Platee=Platee->NextSiblingElement(); 
    }
  }
}
