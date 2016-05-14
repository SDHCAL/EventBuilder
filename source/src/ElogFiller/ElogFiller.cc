#include <vector>
#include <string>
#include <iostream>
#include "tinyxml.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

/*int main()
{
    ifstream iFile("~/SDHCAL/python/Geom.xml");	// input.txt has integers, one per line

    while (!iFile.eof())
    {
    	int x;
    	iFile >> x;
    	cerr << x << endl;
    }

    return 0;
}*/
const std::string elogcommand="elog";
const std::string hostname="lyosvn.in2p3.fr";
//const std::string hostname="localhost";
const std::string logbook="GIF";
const std::string port="443";
//const std::string port="8080";
const std::string subdir="elog";
std::string username_password="";

enum Types {pad,positional,temporal,tcherenkov,tricot,scintillator};

void Write(TiXmlElement* elem,const char* who,double& var)
{
	if(elem->Attribute(who)!=NULL) var=atof(elem->Attribute(who));
  else
  {
		std::exit(1);	
	}
}
void Read(std::string& FileName,std::string& parameters)
{

  static std::vector<std::string>Numbers{"Zero","One","Two","Three","Four","Five","Six","Seven","Eight","Nine"};
  std::string Types="";
  std::string HVChannels="";
  std::string GasChannels="";
  std::string DifTypes="";
  std::string Dif_Ids="";
  std::string zs="";
  std::string HVs="";
  std::string TcheScinti="";
  std::string Scintillator=" -a Scintillator=\"Not Set\"";
  std::string Tcherenkov=" -a Tcherenkov=No";
  std::string PositionDet="";
  unsigned int PlateNumber=0;
  TiXmlDocument doc(FileName.c_str());
  doc.LoadFile();
  if(!doc.LoadFile())
  {
    std::cout<<"not"<<std::endl;
    std::exit(1);	
  }
  else
  {
    TiXmlHandle hdl(&doc);
    TiXmlElement* Detectorr = hdl.FirstChildElement().Element();
    if(Detectorr->Attribute("PositionX")!=NULL)
    {
	std::string e=Detectorr->Attribute("PositionX");
	PositionDet+=" -a \"Position Bati X\"="+e+" ";
    }
    if(Detectorr->Attribute("PositionY")!=NULL)
    {
	std::string e=Detectorr->Attribute("PositionY");
	PositionDet+=" -a \"Position Bati Y\"="+e+" ";
    }
    if(Detectorr->Attribute("PositionZ")!=NULL)
    {
	std::string e=Detectorr->Attribute("PositionZ");
	PositionDet+=" -a \"Position Bati Z\"="+e+" ";
    }
    TiXmlElement* Platee = hdl.FirstChildElement().FirstChildElement().Element();
    TiXmlElement* Diff=nullptr;
    unsigned int DifNumber=0;
    double x,y,z,xy,xz,yz,SizeX,SizeY,I,J,DifId;
    double Gaz_number=-1;
    while (Platee)
    {
      DifNumber=0;
      PlateNumber++;
      Write(Platee,"x",x);
      Write(Platee,"y",y);
      Write(Platee,"z",z);
      std::string zs_to_test="";
      zs_to_test=" -a \"Position"+std::to_string(PlateNumber)+" Z\"="+std::to_string(int(z))+" ";
      Write(Platee,"alpha",xy);
      Write(Platee,"beta",xz);
      Write(Platee,"gamma",yz);
      Write(Platee,"SizeX",SizeX);
      Write(Platee,"SizeY",SizeY);
      std::string type_to_test="";
      if(Platee->Attribute("Glass_type")!=NULL)
	    {
	        if(strcmp(Platee->Attribute("Glass_type"), "standard") == 0) type_to_test+=" -a Type"+std::to_string(PlateNumber)+"=\"Standard\"";
	        else if (strcmp(Platee->Attribute("Glass_type"), "chinese") == 0) type_to_test+=" -a Type"+std::to_string(PlateNumber)+"=\"Chineese\"";
	        else if (strcmp(Platee->Attribute("Glass_type"), "") == 0) type_to_test+=" -a Type"+std::to_string(PlateNumber)+"=\"Unknown\"";
	        else
		      {
		        type_to_test+=" -a Type"+std::to_string(PlateNumber)+"=\"Unknown\"";
		      }
	    }
	    else type_to_test+=" -a Type"+std::to_string(PlateNumber)+"=\"Unknown\"";
	    if(Platee->Attribute("Gaz_number")!=NULL)
	    {
	      std::string e=Platee->Attribute("Gaz_number");
	      if(e.size()!=0) Gaz_number=atof(Platee->Attribute("Gaz_number"));
	    }
	    else Gaz_number=-1;
	    std::string hvchannel_to_test="";
	    std::string gazchannel_to_test="";
	    if(Platee->Attribute("HV_channel")!=NULL) hvchannel_to_test+=" -a \"HV Channel"+std::to_string(PlateNumber)+"\"="+Platee->Attribute("HV_channel")+" ";
	    if(Platee->Attribute("Gaz_channel")!=NULL) gazchannel_to_test+=" -a \"Gas Channel"+std::to_string(PlateNumber)+"\"="+Platee->Attribute("Gaz_channel")+" ";
      Diff=Platee->FirstChildElement();
      while (Diff)
      {
        DifNumber++;
        int DifType=0;
        int up_down=0;
	      if(Diff->Attribute("DifType")!=NULL)
	      {
	        if (strcmp(Diff->Attribute("DifType"), "tcherenkov") == 0) 
	        {
	          DifTypes+=" -a Tcherenkov=\"Yes\"";
	          Write(Diff,"DifId",DifId);
	          TcheScinti+=" -a Dif_Tcherenkov="+std::to_string(int(DifId));
	          Tcherenkov=" -a Tcherenkov=Yes";
	          PlateNumber--;
	        }
          else if (strcmp(Diff->Attribute("DifType"), "scintillator") == 0) 
          {
            DifTypes+=" -a Scintillator=\"Ours\"";
            Write(Diff,"DifId",DifId);
            TcheScinti+=" -a Dif_Scintillator="+std::to_string(int(DifId));
            Scintillator=" -a Scintillator=Ours";
            PlateNumber--;
          }
	      }
	      else
          {
            HVChannels+=hvchannel_to_test;
            Types+=type_to_test;
            GasChannels+=gazchannel_to_test;
            Write(Diff,"DifId",DifId);
            Dif_Ids+=" -a Dif_Id"+std::to_string(PlateNumber)+"="+std::to_string(int(DifId));
            zs+=zs_to_test;
            Write(Diff,"I",I);
            Write(Diff,"J",J);
          }
	      Diff= Diff->NextSiblingElement();
      }
      Platee=Platee->NextSiblingElement(); 
    }
  }
  for(unsigned int i=1;i!=PlateNumber+1;++i) HVs+=" -a HV"+std::to_string(i)+"=\"Not Set\"";
  parameters+=Types+" "+PositionDet+" "+HVChannels+" "+GasChannels+" "+DifTypes+" -a \"Number Dif\"="+Numbers[PlateNumber]+" "+Dif_Ids+" "+zs+" "+HVs+" "+ Tcherenkov+" "+Scintillator+" "+TcheScinti+" " ;
}

void FillElogPlease(std::string Filename,std::string run,std::string daqname)
{
  std::string FileNam=Filename;
  std::string command_part=elogcommand+" -h "+hostname+" -d "+subdir+" -l "+logbook+" -p "+port+" -u "+username_password+" -x -s ";
  std::string defaults="-a Good=\"Not Set\"  -a Beam=\"Not Set\" -a Source=\"Not Set\" -a Gaz=CMS";
  defaults+=" -a Author=DAQ -a Run="+run+" -a \"DAQ Name\"="+daqname+" -a Type=Data ";
  std::string parameters="";
  Read(FileNam,parameters);
  std::string command=command_part+" "+defaults+" "+parameters+" -m ~/SDHCAL/python/Geom.xml \" \"";
  std::cout<<command_part+" "+defaults+" "+parameters+" -m ~/SDHCAL/python/Geom.xml \" \""<<std::endl;
  int good =std::system(command.c_str());
}


 

int main(int argc, char *argv[])
{
  
  char* cpath=getenv("EventBuilder");
  if (cpath == NULL) std::cout<<"Please add EventBuilder path"<<std::endl;
  std::string confdb(cpath);
  size_t ipass = confdb.find(":");
  username_password=confdb.substr(0,ipass)+" "+confdb.substr(ipass+1,confdb.size()); 
  std::string FileName="";
  std::string run="";
  std::string daqname="";
  if(argc==4) 
  {
    FileName=argv[1];
    run=argv[2];
    daqname=argv[3];
    std::cout<<FileName<<"  "<<run<<"  "<<daqname<<std::endl;
  }
  else
  {
    std::cout<<FileName<<"  "<<run<<"  "<<daqname<<std::endl;
    std::exit(1);
  }
  FillElogPlease(FileName,run,daqname);
  return 0;
}

