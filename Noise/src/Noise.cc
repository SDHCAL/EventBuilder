#include "Noise.h"
#include <iostream>
#include <string> 
#include "marlin/Processor.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDDecoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "../../Common/Reader/include/ReaderFactory.h"
#include "../../Common/Reader/include/Reader.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TApplication.h>
#include "TObject.h"
#include "../../Common/Colors.h"
#define degtorad 0.0174532925
#include <cstdlib>
#include <cmath>
#ifndef COLORS_H
#define normal " "
#define red "  "
#define vert " "
#define blanc " "
#define green " "
#define blue " "
#endif
#include <fstream>
#include <algorithm>
#include <numeric>
using namespace std;
std::vector<TH1F*>Histos;
std::vector<TH1F*>Histos3;
std::vector<TH2F*>Noise_Histos;
TH3F* All = new TH3F("aaa","aaa",100,0,100,100,0,100,100,0,100);
std::vector<std::map<int,int>>NumberOfHits;
std::vector<double>NumberOfHitsTotal;
double TotalTime=0;
using namespace marlin;
NoiseProcessor aNoiseProcessor;
NoiseProcessor::NoiseProcessor() : Processor("NoiseProcessorType")
{
  std::vector<std::string> hcalCollections(1,"SDHCAL_HIT");    
  registerInputCollections( LCIO::RAWCALORIMETERHIT,"HCALCollections","HCAL Collection Names",_hcalCollections,hcalCollections); 
  _FileNameGeometry="";
  registerProcessorParameter("FileNameGeometry","Name of the Geometry File",_FileNameGeometry,_FileNameGeometry);
  _ReaderType="";
  registerProcessorParameter("ReaderType","Type of the Reader needed to read InFileName",_ReaderType ,_ReaderType);
  
}
NoiseProcessor::~NoiseProcessor() {}
void NoiseProcessor::init()
{
  _EVENT=0;
  ReaderFactory readerFactory;
  Reader* myReader = readerFactory.CreateReader(_ReaderType); 
  if(myReader)
  {
    myReader->Read(_FileNameGeometry,geom);
    std::map<int, Dif > Difs=geom.GetDifs();
    std::map<int,int> PlansType;
    for(std::map<int, Dif >::iterator it=Difs.begin();it!=Difs.end();++it)
    {
      if(geom.GetDifType(it->first)!=temporal)
      {
        PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
      }
    }
    for(std::map<int, int >::iterator it=PlansType.begin();it!=PlansType.end();++it)
    {
      testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second));
      std::string a="Time_Distribution"+ std::to_string( (long long int) it->first +1 );
      std::string b="number_of_hits_distribution"+std::to_string( (long long int)it->first +1 );
      std::string c="number_of_hits_distribution2"+std::to_string( (long long int)it->first +1 );
      Histos.push_back(new TH1F(a.c_str(),a.c_str(),50,0,50));
      Histos3.push_back(new TH1F(c.c_str(),c.c_str(),50,0,50));
      if(it->second==positional) {std::cout<<"Creation positional"<<std::endl;Noise_Histos.push_back(new TH2F(b.c_str(),b.c_str(),128,0,128,1,0,50));}
      else {std::cout<<"Creation one"<<std::endl;Noise_Histos.push_back(new TH2F(b.c_str(),b.c_str(),100,0,100,100,0,100));}
      Noise.push_back(std::map<int ,int>());
      NumberOfHitsTotal.push_back(0);
      Noise_vector.push_back(std::map<int ,std::vector<CalorimeterHit*>>());
      NumberOfHits.push_back(std::map<int ,int>());
    }
  }
  else
  {
    std::cout << "Reader type n'existe pas !!" << std::endl;
    std::exit(1);
  }
  delete myReader;
}




void NoiseProcessor::processEvent( LCEvent* evtP ) 
{	
  
  if (evtP != NULL)
  {
    _NbrRun=evtP->getRunNumber();
    _eventNr=evtP->getEventNumber();
    for(unsigned int i=0; i< _hcalCollections.size(); ++i)
    {
			LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
			if(col == NULL)
			{
				std::cout<< red << "TRIGGER SKIPED ..."<< normal <<std::endl;
	      break;
	    } 
	    CellIDDecoder<CalorimeterHit> cd(col);
	    int numElements = col->getNumberOfElements();
	    _EVENT+=numElements;
			for(unsigned int i=0; i<Noise_vector.size();++i)
			{
				Noise_vector[i].clear();
			}
	    for (int ihit=0; ihit < numElements; ++ihit) 
	    {
	    	CalorimeterHit *raw_hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(ihit)) ;
	      if (raw_hit != NULL)
	      {
					Noise[geom.GetDifNbrPlate(cd(raw_hit)["Dif_id"])-1][raw_hit->getTime()]++;
					if(raw_hit->getTime()>0)Noise_vector[geom.GetDifNbrPlate(cd(raw_hit)["Dif_id"])-1][raw_hit->getTime()].push_back(raw_hit);
					//std::cout<<red<<evtP->getTimeStamp()<<normal<<std::endl;
					TotalTime+=evtP->getTimeStamp();
					//std::cout<<magenta<<TotalTime<<normal<<std::endl;
	      }
	    }
    
  		
		for(unsigned int i=0; i<Noise_vector.size();++i)
		{
			
 			std::cout<<red<<i<<"  "<<yellow<<"  "<<Noise_vector[i].size()<<"  "<<normal<<std::endl;
 			for(std::map<int,std::vector<CalorimeterHit*>>::iterator it=Noise_vector[i].begin(); it!=Noise_vector[i].end();++it)
 			{
   			
				//if(it->second.size()>1)
				//{
			NumberOfHits[i][it->second.size()]++;
			NumberOfHitsTotal[i]+=it->second.size();
                        
				//std::cout<<blue<<it->first<<"  "<<it->second.size()<<normal<<"  ";
   			for(std::vector<CalorimeterHit*>::iterator itt=(it->second).begin();itt!=(it->second).end();++itt)//unsigned int k=0;k<it->second.size();++k)
   			{ 			
                                //std::cout<<red<<cd(*itt)["I"]<<normal<<std::endl;
				Noise_Histos[i]->Fill(cd(*itt)["I"],cd(*itt)["J"]);
				All->Fill(cd(*itt)["I"],cd(*itt)["J"],cd(*itt)["K"]);
				
					//std::cout<<"*"<<red<<(*itt)->getPosition()[0]<<"   "<<(*itt)->getPosition()[1]<<"   "<<(*itt)->getPosition()[2]<<normal<<"*";
			}
  			//std::cout<<std::endl;
				//}
 			}
		}
		
		//for(unsigned int i=0; i<Noise_vector.size();++i)
		//{
			
		//	if((Noise_vector[i].end()--)->first>max) max=(Noise_vector[i].end()--)->first;
			
		//}
}
		
	
	}

}

 


void NoiseProcessor::end()
{
//delete h;
  std::string b="Results_Noise"+ std::to_string( (long long int) _NbrRun)+".root";
  TFile *hfile = new TFile(b.c_str(),"RECREATE","Results");

  for(unsigned int i=0; i<Noise_vector.size();++i)
    {
      std::cout<<"Plane : "<<i<<std::endl;
      for(std::map<int,int>::iterator it=NumberOfHits[i].begin(); it!=NumberOfHits[i].end();++it)
	{
				std::cout<<"*"<<blue<<it->first<<"  "<<it->second<<normal<<"*";
				Histos3[i]->Fill(it->first,it->second);
	}
      Histos3[i]->Write();
      Noise_Histos[i]->Scale(TotalTime*200e-9);
      Noise_Histos[i]->Write();
    }

  All->Write();


  for(unsigned int i=0; i<Noise_vector.size();++i)
    {
      std::cout<<green<<NumberOfHitsTotal[i]/(TotalTime*200e-9*33*30)<<" "<<normal;
    }


  std::vector<int> number;
  for(unsigned int i=0;i!= Noise.size();++i)
    {
	Histos[i]->SetBit(TH1::kCanRebin);
	number.push_back(0);
	for(std::map<int,int>::iterator it = Noise[i].begin();it!=Noise[i].end();++it)
	  {
	    Histos[i]->Fill(it->first,it->second*1.0/_EVENT);
	  }
	
	Histos[i]->Write();
    }
  for(unsigned int i=0;i!= Histos.size();++i)
    {
      delete Histos[i];
      delete Histos3[i];
    }
  Histos.clear();
  Histos3.clear();
  delete hfile;
}
