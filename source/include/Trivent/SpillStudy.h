#ifndef SPillStudy
#define SPillStudy
#include<cstdint>
#include"IO/LCReader.h"
#include<vector>
#include<string>
#include<iostream>
#include "Trivent/HistoPlane.h"
#include<map>
#include<iterator>
#include"TH1.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDEncoder.h"
#include <EVENT/LCGenericObject.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <IMPL/LCRunHeaderImpl.h>
#include "EVENT/RawCalorimeterHit.h"
#include "Geometry/Geometry.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/CellIDEncoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <EVENT/LCRunHeader.h>
#include "UTIL/LCTOOLS.h"
#include "marlin/Global.h"
TH1D* timestamp=nullptr;
TH1D* timestamps=nullptr;
TH1D* time2read=nullptr;
TH1D* time2readtime=nullptr;
TH1D* TimeSpill=new TH1D("Duration Spill","Duration Spill",50000,0,500);
TH1D* TimeRamFull=new TH1D("Duration Frame","Duration Frame",1000,0,1);
TH1D* UsefullTime=new TH1D("UselessTime","UselessTime",1000,0,1);
std::vector<TH1*>typeee;
unsigned long int total_time_file=0;
unsigned long int total_time_by_user=0;
void SpillStudy(unsigned int& _skip, unsigned int& _GlobalEvents,unsigned int& _maxRecord,Geometry& geom,std::map<int,HistoPlane*>& HistoPlanes,std::vector<std::string>& LCIOFiles,bool& _Spill_Study,std::map<std::string,std::vector<double>> & Types)
{
  if(_Spill_Study)
  {
    std::map<unsigned long long int,int>Time_stamps;
    std::map<unsigned long long int,int>Time_stampss;
    int32_t timetime=0;
    LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
    std::vector<std::string>colee{"DHCALRawHits"};
    lcReader->setReadCollectionNames( colee ) ;
    unsigned eventnumber=-1;
    unsigned long long int min_by_user=999999999999;
    unsigned long long int  max_by_user=1;
    for(unsigned int i=0;i!=LCIOFiles.size();++i)
    {
      std::cout<<"I'm Readind The DATA in "<<LCIOFiles[i]<<std::endl;
      LCEvent* evt(0) ;
      lcReader->open( LCIOFiles[i] ) ;
      std::cout<<lcReader->getNumberOfRuns()<<" "<<lcReader->getNumberOfEvents()<<std::endl;
      evt=lcReader->readNextEvent();
      unsigned int counter=0;
      std::cout<<"Spill Study progression"<<std::endl;
      do
      {
        counter++;
        Progress(_skip,_GlobalEvents,_maxRecord,counter);
        LCCollection* col=evt->getCollection("DHCALRawHits");
        if(col!=nullptr)
		    {
          eventnumber++;
          RawCalorimeterHit * myhit = dynamic_cast<RawCalorimeterHit*>(col->getElementAt(0)) ;
         
          unsigned int dif_id=myhit->getCellID0()&0xFF;
          if (dif_id==0) return;
          std::string name="DIF"+patch::to_string(dif_id)+"_Triggers";
          
          lcio::IntVec vTrigger;
          col->getParameters().getIntVals(name,vTrigger);
          unsigned long long _bcid=0;
          if (vTrigger.size()>=5)
          {
            unsigned long long Shift=16777216ULL;
  	        _bcid=vTrigger[4]*Shift+vTrigger[3];
          }
          Time_stamps[_bcid]++;
          for(unsigned int i=0;i<col->getNumberOfElements();++i)
          {
            myhit=dynamic_cast<RawCalorimeterHit*>(col->getElementAt(i));
            unsigned int dif_id=myhit->getCellID0()&0xFF;
        
       
            if(geom.GetDifType(int(dif_id))<=1)if(geom.GetDifNbrPlate(int(dif_id))!=-1) HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("timestamp")->Fill(_bcid*200e-9,1);
        }
        Time_stampss[_bcid]=col->getNumberOfElements();
        if(eventnumber>=_skip&&eventnumber<=_maxRecord)
				{
          
					if(_bcid<min_by_user)
          {
				    min_by_user=_bcid;
			 			
					}
        	if(_bcid>max_by_user)
					{
		       
						max_by_user=_bcid;
					}
				}
		  }
		  
		  if(counter<=_maxRecord) evt=lcReader->readNextEvent();
      else evt =nullptr;
	  }
	  while(evt!=nullptr);
    lcReader->close();
    delete lcReader ;
  }
   //////////////////////////////////////////////
  unsigned int total_Ram_Full=0;
  unsigned int total_hits=0;
  for(std::map<unsigned long long int,int>::iterator it=Time_stampss.begin();it!=Time_stampss.end();++it) total_hits+=it->second;
  for(std::map<unsigned long long int,int>::iterator it=Time_stamps.begin();it!=Time_stamps.end();++it) total_Ram_Full+=it->second;
  double moyenne=total_hits/(total_Ram_Full*1.0);
  std::cout<<total_Ram_Full<<"  "<<total_hits<<"  "<<moyenne<<std::endl;
  //////////////////////////////////////////////
  std::cout<<Time_stampss.size()<<std::endl;
  unsigned long long int min=9999999999;
  unsigned long long int  max=1;
  double min_between=99999;
  double max_between=0;
  double min_betweentime=99999;
  double max_betweentime=0;
  std::vector<double>Vec_timebetween;
  std::vector<double>timebetweentime_toplot;
  std::map<double,std::vector<double>>Vec_timebetweentime;
  Vec_timebetween.reserve(Time_stamps.size()-1);
  std::map<unsigned long long int,int>::iterator ultimate=--Time_stamps.end();
  std::map<unsigned long long int,int>::iterator ultimateless=ultimate; --ultimateless;
  bool skip=false;
  for(std::map<unsigned long long int,int>::iterator it=Time_stamps.begin();it!=Time_stamps.end();++it)
  {
    if(skip) 
    {
      skip=false; 
      continue;
    }
    if(it!=ultimate)
    {
      std::map<unsigned long long int,int>::iterator itt=it;
		  ++itt;
      double timebetween=(itt->first-it->first)*200e-9;
      if(timebetween<min_between)min_between=timebetween;
      if(timebetween>max_between)max_between=timebetween;
      Vec_timebetween.push_back(timebetween);
      if(it==Time_stamps.begin())
      {
        Vec_timebetweentime[(itt->first-it->first)*200e-9].push_back(it->first);
			  Vec_timebetweentime[(itt->first-it->first)*200e-9].push_back(itt->first);
	    }
      if(it!=ultimateless)
		  {
			  std::map<unsigned long long int,int>::iterator ittt=itt;
			  ++ittt;
        double a= (ittt->first-itt->first);
			  double b=(itt->first-it->first);
			  double timebetweentime=(a-b)*200e-9;
        
        timebetweentime_toplot.push_back((b-a)*200e-9);
			
        if(timebetweentime<min_betweentime)min_betweentime=timebetweentime;
        if(timebetweentime>max_betweentime)max_betweentime=timebetweentime;
        std::map<double,std::vector<double>>::iterator toadd=Vec_timebetweentime.end();
			  --toadd;
        for(std::map<double,std::vector<double>>::iterator j=Vec_timebetweentime.begin();j!=Vec_timebetweentime.end();++j)
			  {
          
          if(fabs(fabs(timebetweentime)-j->first) <=0.50)
				  {
					  (j->second).push_back(ittt->first);
            
            break;
				  }
          else if(j==toadd&&it!=Time_stamps.begin()) 
          {
					  Vec_timebetweentime[(ittt->first-itt->first)*200e-9].push_back(ittt->first);
					 
            skip=true;
   					break;
				  }
			  }
		  }
	  }
	  if(it->first<min&&it->first!=0)min=it->first;
    if(it->first>max)max=it->first;
    
  }
   for(std::map<double,std::vector<double>>::iterator j=Vec_timebetweentime.begin();j!=Vec_timebetweentime.end();++j)
  {
	  if(j->first<=0.5) 
    {
		  //Between_spill[j->first]=j->second;
		  Types["RamFull"].insert(Types["RamFull"].end(),j->second.begin(),j->second.end());
	  }
    else if (j->first>=1) 
	  {
		  //Spills[j->first]=j->second;
		  Types["Spill"].insert(Types["Spill"].end(),j->second.begin(),j->second.end());
	  }
    else 
	  {
		  //Others[j->first]=j->second;
		  Types["Other"].insert(Types["Other"].end(),j->second.begin(),j->second.end());
	  }
  }
  int diffbetweentime=(int((max_betweentime-min_betweentime))+1)*10;
  time2readtime = new TH1D("Difference in time bettwen two time","Difference in time bettwen two time",10000000,-10,10);
  for(unsigned int i =0;i<timebetweentime_toplot.size();++i)
  {
	  time2readtime->Fill(timebetweentime_toplot[i]);
    
  }
  int diffbetween=(int((max_between-min_between))+1)*100;
  time2read = new TH1D("Distribution time between 2 read out","Distribution time between 2 read out",5000000,min_between,100);
  for(unsigned int i =0;i<Vec_timebetween.size();++i)
  {
	  time2read->Fill(Vec_timebetween[i]);
   
  }
  double diff=(max-min);
  if(diff>9999999)diff=9999999;
  timestamp= new TH1D("timestamp","timestamp",diff,min*200e-9,max*200e-9);
  timestamps= new TH1D("timestampS","timestampS",diff,min*200e-9,max*200e-9);
  for(std::map<double,std::vector<double>>::iterator j=Vec_timebetweentime.begin();j!=Vec_timebetweentime.end();++j)
  {
    std::string name = "type"+patch::to_string(j->first);
	  typeee.push_back(new TH1D(name.c_str(),name.c_str(),1000000,0.,10000.));
    for(unsigned int i =0;i<j->second.size();++i)
	  {
		  typeee[typeee.size()-1]->Fill(j->second[i]*200e-9);
		  
	  }
  }

  for(std::map<unsigned long long int ,int>::iterator it=Time_stamps.begin();it!=Time_stamps.end();++it)
  {
	  timestamp->Fill(it->first*200e-9,it->second);
	  
    
  }
  for(std::map<unsigned long long int ,int>::iterator it=Time_stampss.begin();it!=Time_stampss.end();++it)
  {
	  timestamps->Fill(it->first*200e-9,it->second);
    
  }
  total_time_by_user=(max_by_user-min_by_user)*200e-9;
  total_time_file=(max-min)*200e-9;
  }
}
#endif
