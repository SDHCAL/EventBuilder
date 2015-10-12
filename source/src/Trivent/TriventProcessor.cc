#include "Trivent/TriventProcessor.h"
#include "Progress.h"
#include "Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include "marlin/Processor.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDEncoder.h"
#include <EVENT/LCGenericObject.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "TH2F.h"
#include <iomanip>
#include <iostream>
#include <climits>
#include <ctime>
#include "TFile.h"
#include "TTree.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <IMPL/LCRunHeaderImpl.h>
#include "Colors.h"
#include "TBranch.h"
#include "TObject.h"
#include "TPDF.h"
#include "TCanvas.h"
#include "Trivent/SpillStudy.h"
#ifndef COLORS_H
#define normal " "
#define red "  "
#define vert " "
#define blanc " "
#endif
#include <fstream>
#include "Reader/ReaderFactory.h"
#include "Reader/Reader.h"
#include "Trivent/Mapping.h"
#include "Trivent/HistoPlane.h"
#include "TStyle.h"
#include "TF1.h"
#include "TObject.h"
#include "TList.h"
#include "TH3.h"
#include "TColor.h"
#include "TMath.h"
#include "Patch.h"
#include "TMarker.h"
#include "TColor.h"
#include "TNamed.h"
#include "THnSparse.h"
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCRunHeader.h" 
#include <algorithm> 
#include "EVENT/SimCalorimeterHit.h" 
#include "EVENT/CalorimeterHit.h" 
#include "EVENT/RawCalorimeterHit.h" 
//#include "THClass.h"
bool pdf=false;
bool HasScintiSignal=false;



enum Threshold{Threshold2=1,Threshold1,Threshold3};

using namespace marlin;
std::vector<int>ScintillatorCoincidence;
#define degtorad 0.0174532925
unsigned int EventsNoise=0;

unsigned int EventsSelected=0;
unsigned int EventsSelectedt=0;

unsigned int _eventNr=0;
#define size_pad 10.4125
#define size_strip 2.5
//unsigned long int total_time_file=0;
//unsigned long int total_time_by_user=0;
unsigned long int min_by_user=0;
unsigned long int max_by_user=0;
unsigned long long int HistoPlane::global_total_time =0;
//std::map<int,TGraphTime*>time_graph;

TH1F diff("diff","diff",1000,-500,500);
std::map<std::string,std::vector<double>>Types;
//class ToTree
//{
//public:
//int pI,pJ,pK,pAsic,pDifId,pAsicChannel;
//unsigned long long int pTime;
//double pX,pY,pZ;
//bool pEvent;
//};
//ToTree totree;
//std::string name="Tree";
//TTree* t= new TTree(name.c_str(), name.c_str());
//TBranch* Branch1 =  t->Branch("X",&(totree.pX));
//TBranch* Branch2 =  t->Branch("Y",&(totree.pY));
//TBranch* Branch3 =  t->Branch("Z",&(totree.pZ));
//TBranch* Branch4 =  t->Branch("I",&(totree.pI));
//TBranch* Branch5 =  t->Branch("J",&(totree.pJ));
//TBranch* Branch6 =  t->Branch("K",&(totree.pK));
//TBranch* Branch7 =  t->Branch("Time",&(totree.pTime));
//TBranch* Branch8 =  t->Branch("Asic",&(totree.pAsic));
//TBranch* Branch9 =  t->Branch("DifId",&(totree.pDifId));
//TBranch* Branch10 =  t->Branch("AsicChannel",&(totree.pAsicChannel));
//TBranch* Branch11 = t->Branch("Event",&(totree.pEvent));
std::vector<std::string  > th1{
                               "Time_Distr","Hits_Distr",
                               "Time_Distr_Events","Time_Distr_Events_S1","Time_Distr_Events_S2","Time_Distr_Events_S3",
                               "Hits_Distr_Events","Hits_Distr_Events_S1","Hits_Distr_Events_S2","Hits_Distr_Events_S3",
                               "Time_Distr_Noise","Time_Distr_Noise_S1","Time_Distr_Noise_S2","Time_Distr_Noise_S3",
                               "Hits_Distr_Noise","Hits_Distr_Noise_S1","Hits_Distr_Noise_S2","Hits_Distr_Noise_S3",
                               "timestamp"
                               };
std::vector<std::string> th2 {
                              "Flux_Noise","Flux_Events","EffiScintiOnly",
                              "Flux_Noise_S1","Flux_Noise_S2","Flux_Noise_S3","Flux_Noise_Pon",
                              "Flux_Events_S1","Flux_Events_S2","Flux_Events_S3","Flux_Events_Pon",
                             };
std::vector<std::string>th2_Asic{"Flux_Noise_Asic","Flux_Events_Asic"};
int _NbrRun=0;


TriventProcessor aTriventProcessor;

TriventProcessor::TriventProcessor() : Processor("TriventProcessorType")
{
  
  std::vector<std::string> hcalCollections(1,"DHCALRawHits");
  registerInputCollections( LCIO::RAWCALORIMETERHIT,"HCALCollections","HCAL Collection Names",_hcalCollections,hcalCollections);
  _FileNameGeometry="";
  registerProcessorParameter("FileNameGeometry","Name of the Geometry File",_FileNameGeometry,_FileNameGeometry);
  _ReaderType="";
  registerProcessorParameter("ReaderType","Type of the Reader needed to read InFileName",_ReaderType ,_ReaderType);
  _outFileName="LCIO_clean_run.slcio";
  registerProcessorParameter("LCIOOutputFile","LCIO file",_outFileName,_outFileName);
  _noiseFileName="";
  registerProcessorParameter("NOISEOutputFile" ,"NOISE file" ,_noiseFileName ,_noiseFileName);
  _timeWin = 2;
  registerProcessorParameter("timeWin" ,"time window = 2 in default",_timeWin ,_timeWin);
  _noiseCut = 7;
  registerProcessorParameter("noiseCut" ,"noise cut in time spectrum 7 in default",_noiseCut ,_noiseCut);
  _LayerCut = 3;
  registerProcessorParameter("LayerCut" ,"cut in number of layer 3 in default",_LayerCut ,_LayerCut);
  _TriggerTime = 0;
  registerProcessorParameter("TriggerTime" ,"All Events with Time greater than this number will be ignored (0) in case of Triggerless",_TriggerTime ,_TriggerTime);
  _WantDistribution = false;
  registerProcessorParameter("Distribution" ,"Create Distribution of hits for Plates, Difs, Asics, and Pads",_WantDistribution ,_WantDistribution);
  _WantCalibration = false;
  registerProcessorParameter("Calibration" ,"Create Calibration file for the Detector",_WantCalibration ,_WantCalibration);
  _Database_name ="";
  registerProcessorParameter("Database_name" ,"Name of the Database for Calibration",_Database_name ,_Database_name);
  _Spill_Study =false;
  registerProcessorParameter("Spill_Study" ,"Study the Spill",_Spill_Study ,_Spill_Study);
  _efficiencyFrontScintillator =1;
  registerProcessorParameter("efficiencyFrontScintillator" ,"efficiency Front Scintillator",_efficiencyFrontScintillator ,_efficiencyFrontScintillator);
  _efficiencyBackScintillator =1;
  registerProcessorParameter("efficiencyBackScintillator" ,"efficiency Back Scintillator",_efficiencyBackScintillator ,_efficiencyBackScintillator);
  _IgnorebeginningSpill=0;
  registerProcessorParameter("IgnorebeginingSpill" ,"Ignore begining of the Spills ",_IgnorebeginningSpill,_IgnorebeginningSpill);
  _Delimiters="";
  registerProcessorParameter("Delimiters" ,"Delimiters",_Delimiters,_Delimiters);  
}




int counttt=0;
TriventProcessor::~TriventProcessor() {}

void TriventProcessor::processRunHeader( LCRunHeader* run)
{
    LCTOOLS::dumpRunHeader(run);
}

void TriventProcessor::init()
{   
  TOTALNUMBERHITSCINTI=0;
  ReaderFactory readerFactory;
  Reader* myReader = readerFactory.CreateReader(_ReaderType);
  if(myReader) 
	{
    myReader->Read(_FileNameGeometry,geom);
    geom.PrintGeom();
    std::map<int, Dif > Difs=geom.GetDifs();
    unsigned int NbrPlate =0;
    for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it) 
		{
      if(geom.GetDifType(it->first)!=temporal&&geom.GetDifType(it->first)!=tcherenkov&&geom.GetDifType(it->first)!=scintillator) 
	    {
        SinCos[it->first]=std::vector<double>{cos(geom.GetDifAlpha(it->first)*degtorad),sin(geom.GetDifAlpha(it->first)*degtorad),cos(geom.GetDifBeta(it->first)*degtorad),sin(geom.GetDifBeta(it->first)*degtorad),cos(geom.GetDifGamma(it->first)*degtorad),sin(geom.GetDifGamma(it->first)*degtorad)};
        NbrPlate=geom.GetDifNbrPlate(it->first)-1;
        if(HistoPlanes.find(NbrPlate)==HistoPlanes.end()) 
        {
          EffiwithDiscri.push_back(0);
          HistoPlanes.insert(std::make_pair(NbrPlate, new HistoPlane(_WantDistribution,NbrPlate,geom.GetDifsInPlane(NbrPlate),geom.GetSizeX(NbrPlate),geom.GetSizeY(NbrPlate),th1,th2,th2_Asic)));
				}
      }
    }
    FillDelimiter(_Delimiters,HistoPlanes.size(),Delimiter);
  } 
	else 
	{
    std::cout << "Reader type n'existe pas !!" << std::endl;
    std::exit(1);
  }
  delete myReader;
  TouchedEvents=0;
  eventtotal=0;
  _maxRecord= Global::parameters->getIntVal("MaxRecordNumber")-1;
  _skip= Global::parameters->getIntVal("SkipNEvents");
  std::vector<std::string>LCIOFiles;
  Global::parameters->getStringVals("LCIOInputFiles" ,LCIOFiles );
  
  for(unsigned int i=0;i!=LCIOFiles.size();++i)
  {
    LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
    std::vector<std::string>colee{"DHCALRawHits"};
    lcReader->open( LCIOFiles[i] ) ;
    lcReader->setReadCollectionNames( colee ) ;
    _GlobalEvents+=lcReader->getNumberOfEvents()-1;
    delete lcReader;
  }
  std::cout<<yellow<<_GlobalEvents<<normal<<std::endl;
  SpillStudy(_skip,_GlobalEvents,_maxRecord,geom,HistoPlanes,LCIOFiles,_Spill_Study);
  printParameters();
  if(_WantCalibration==true&&_Database_name ==""){std::cout<<red<<"Name's Database is unknown from the xml file"<<normal<<std::endl;std::exit(1);}
  if(_LayerCut==-1){std::cout<<red<<"LayerCut set to -1, assuming that you want to use trigger to see events"<<normal<<std::endl;}
  _EventWriter = LCFactory::getInstance()->createLCWriter() ;
  _EventWriter->setCompressionLevel( 2 ) ;
  _EventWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ;
  if(_noiseFileName!="") 
  {
    _NoiseWriter = LCFactory::getInstance()->createLCWriter() ;
    _NoiseWriter->setCompressionLevel( 2 ) ;
    _NoiseWriter->open(_noiseFileName.c_str(),LCIO::WRITE_NEW) ;
  }
}
 
unsigned long long bcid_spill=0;
void TriventProcessor::processEvent( LCEvent * evtP )
{
  if (evtP== nullptr ) return;
  _NbrRun=evtP->getRunNumber();
  _eventNr=evtP->getEventNumber()+1;
  time_t date(evtP->getTimeStamp());
  Progress(_skip,_GlobalEvents,_maxRecord,_eventNr);
  LCCollection* col2=nullptr;
  LCCollection* col3=nullptr;
  std::vector<std::string>names=*evtP->getCollectionNames();
  for(unsigned int i=0;i<names.size();++i)
  {
    if(names[i]=="DHCALRawTimes")
  	{
  		col2 = evtP ->getCollection("DHCALRawTimes");
  		RawTimeDifs.clear();
      for (int ihit=0; ihit < col2->getNumberOfElements(); ++ihit) 
		  {
	      EVENT::CalorimeterHit* raw_time = dynamic_cast<EVENT::CalorimeterHit* >( col2->getElementAt(ihit)) ;
	  		// std::cout<<raw_time->getTime()<<"  "<<raw_time->getEnergyError()<<std::endl;
	  		//RawTimeDifs[raw_time->getTimeStamp()].push_back(raw_time);
		  }
	  }
	  if(names[i]=="Scintillator")
  	{
  	  HasScintiSignal=true;
  	  ScintillatorCoincidence.clear();
  		col3 = evtP ->getCollection("Scintillator");
      for (int ihit=0; ihit < col3->getNumberOfElements(); ++ihit) 
		  {
	      EVENT::LCGenericObject* raw_scin = dynamic_cast<EVENT::LCGenericObject* >( col3->getElementAt(ihit)) ;
        _Front_scintillator+=raw_scin->getIntVal(0);
        _Back_scintillator+=raw_scin->getIntVal(1);
        _Both_scintillator+=raw_scin->getIntVal(2);
        ScintillatorCoincidence.push_back(raw_scin->getIntVal(3));
      }
    }
  }
  for(unsigned int i=0; i< _hcalCollections.size(); i++) 
  {
    LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
    if(col==nullptr) 
    {
	    std::cout << "TRIGGER SKIPED ..."<<std::endl;
	    _trig_count++;
	    break;
    }
    RawCalorimeterHit * hittt = dynamic_cast<RawCalorimeterHit*>(col->getElementAt(/*hit*/0)) ;
    //std::cout<<blue<<myhit->getCellID0()<<"  "<<std::endl;
    unsigned int dif_id=hittt->getCellID0()&0xFF;
    if (dif_id==0) return;
    std::string name="DIF"+patch::to_string(dif_id)+"_Triggers";
    //std::cout<<name<<std::endl;
    lcio::IntVec vTrigger;
    col->getParameters().getIntVals(name,vTrigger);
    unsigned long long _bcid=0;
    if (vTrigger.size()>=5)
    {
      unsigned long long Shift=16777216ULL;
  	  _bcid=vTrigger[4]*Shift+vTrigger[3];
    }
    bool to_skip=false;
    
    for(unsigned int i=0;i<Types["Spill"].size();++i)
    if(Types["Spill"][i]==_bcid)
    {
        bcid_spill=_bcid;
    }
    //if(_bcid-bcid_spill<=2500000)std::cout<<green<<_bcid<<"  "<<bcid_spill<<"  "<<_bcid-bcid_spill<<normal<<std::endl;
    //else std::cout<<red<<_bcid<<"  "<<bcid_spill<<"  "<<_bcid-bcid_spill<<normal<<std::endl;
    if(_IgnorebeginningSpill>0)
    {
      if(_bcid-bcid_spill<=_IgnorebeginningSpill)
      {
	      std::cout<<"ignoring : "<<_bcid*200e-9<<std::endl;
        to_skip=true;
      }
    }
    if(to_skip==false)  processCollection(evtP,col);
  }
} 
unsigned long int debut_RamFull;
unsigned long int fin_RamFull;
unsigned long int debut_spill;
unsigned long int fin_spill;
bool NeverStarted = true;
void TriventProcessor::processCollection(EVENT::LCEvent *evtP,LCCollection* col)
{
  if(NeverStarted==true)
  {
	  NeverStarted=false;
	  if(Types["Spill"].size()!=0)fin_spill=Types["Spill"][0];
	  debut_spill=0;
	  debut_RamFull=0;
	  if(Types["RamFull"].size()!=0)fin_RamFull=Types["RamFull"][0];
  }
  Times.clear();
  RawHits.clear();
  for(unsigned int i =0;i<HistoPlanes.size();++i)
  {
	  HistoPlanes[i]->Init_local_min_max();
    if(_WantDistribution==true)HistoPlanes[i]->Init_Hit_In_Pad_Per_RamFull();
  }
  BehondTrigger.clear();
  int numElements = col->getNumberOfElements();

  for (int ihit=0; ihit < numElements; ++ihit) 
  {
    RawCalorimeterHit *raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit)) ;
    if (raw_hit != nullptr) 
    {
	    unsigned int dif_id  = (raw_hit)->getCellID0() & 0xFF ;
	    if(geom.GetDifType(dif_id)==scintillator)continue;
	    if(geom.GetDifNbrPlate(dif_id)==-1) 
	    {
	      if(Warningg[dif_id]!=true) 
	      {
	        Warningg[dif_id]=true;
	        std::cout<<"Please add DIF "<<dif_id<<" to your geometry file; I'm Skipping its data."<<std::endl;
	      }
	      continue;
	    }
	    if(raw_hit->getTimeStamp()<0)
	    {
	      std::vector<unsigned int>b{dif_id,(unsigned int)((raw_hit->getCellID0() & 0xFF00)>>8),(unsigned int)((raw_hit->getCellID0() & 0xFF00)>>16)}; Negative[b][raw_hit->getTimeStamp()]++;
	    }
	    if(_TriggerTime==0 || (raw_hit->getTimeStamp()<=_TriggerTime&&raw_hit->getTimeStamp()>=0))
	    {
        ////ADDEDD TO TEST///////////////////////
        ///MAYBE DETECT A PROBLEM 
	      /*
		    983601575  936784058  1170124157  46.668  983601575  984314249  0.142535
        [ VERBOSE "Trivent"] 984314249  936784058  1170124157  46.668  984314249  984314249  0
        [ VERBOSE "Trivent"] 1170124157  1170124157  1170809308  0.13703  984314249  984314249  0
        [ VERBOSE "Trivent"] 1170809308  1170809308  1404126445  46.6634  984314249  984314249  0
        [ VERBOSE "Trivent"] 1171504528  1170809308  1404126445  46.6634  1171504528  1172167416  0.132578
	      */
        if(ihit==0)
	      {
	        std::string name="DIF"+patch::to_string(dif_id)+"_Triggers";
          lcio::IntVec vTrigger;
          col->getParameters().getIntVals(name,vTrigger);
          unsigned long long _bcid=0;
          if (vTrigger.size()>=5)
          {
            unsigned long long Shift=16777216ULL;
  	        _bcid=vTrigger[4]*Shift+vTrigger[3];
          }
          for(std::map<string,std::vector<double>>::iterator it=Types.begin();it!=Types.end();++it)
	        {  	
            for(unsigned int i=1;i<it->second.size();++i)
	    		  {	
              if(it->second[i]==_bcid)
				      {
    			      if(it->first=="Spill")
					      {
						      debut_spill=it->second[i];
						      if(i!=it->second.size()-2) fin_spill=it->second[i+1];
         			    break;
					      }
					      if(it->first=="RamFull")
					      {
						      debut_RamFull=it->second[i];
						      if(i!=it->second.size()-2)if(it->second[i+1]<=fin_spill) fin_RamFull=it->second[i+1];
  						    break;
					      }
				      }
			      }
			    }
          //std::cout<<red<<_bcid<<"  "<<debut_spill<<"  "<<fin_spill<<"  "<<(fin_spill-debut_spill)*200e-9<<"  "<<debut_RamFull<<"  "<<fin_RamFull<<"  "<<(fin_RamFull-debut_RamFull)*200e-9<<normal<<std::endl;
          TimeSpill->Fill((fin_spill-debut_spill)*200e-9);
	        TimeRamFull->Fill((fin_RamFull-debut_RamFull)*200e-9);
	      }
        /////////////////////////////////////////
	      /////supress this in case of emergency
	      //int a,b,c,d;
	      //if(Delimiter.find(dif_id)==Delimiter.end()){a=Delimiter[1][0];b=Delimiter[1][1];c=Delimiter[1][2];d=Delimiter[1][3];}
	      //else {a=Delimiter[dif_id][0];b=Delimiter[dif_id][1];c=Delimiter[dif_id][2];d=Delimiter[dif_id][3];}
	      //std::cout<<Delimiter.size()<<std::endl;
	      //std::cout<<Delimiter[dif_id][0]<<"  "<<std::endl;//<<Delimiter[dif_id][1]<<"  "<<Delimiter[dif_id][2]<<"  "<<Delimiter[dif_id][3]<<std::endl;
	      //if(a<=I&&b>=I&&c<=J&&d>=J)
	      //{
	      ///////////////////////////////////////
	      //if(geom.GetDifType(dif_id)!=scintillator)
	      //{
	      
	      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_hit_trigger();
	      Times[raw_hit->getTimeStamp()]++;
	      RawHits[raw_hit->getTimeStamp()].push_back(raw_hit);//}
	      ////////////////////////////////////////
	      //}
	      /////////////////////////////////////////
	      //if(raw_hit->getTimeStamp()<0)std::cout<<yellow<<raw_hit->getTimeStamp()<<"  "<<((raw_hit)->getCellID0() & 0xFF)<<normal<<std::endl;
	    }
	    else
	    {
	      //if(geom.GetDifType(dif_id)!=scintillator)
	      //{
	      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_hit_other();
	      BehondTrigger[raw_hit->getTimeStamp()].push_back(raw_hit);
	      //}
	      //std::cout<<blue<<raw_hit->getTimeStamp()<<"  "<<BehondTrigger.size()<<normal<<std::endl;
	    }
	    //if(geom.GetDifType(dif_id)!=scintillator)
	      //{
	    if(raw_hit->getTimeStamp()>HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Get_local_max())HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_local_max(raw_hit->getTimeStamp());
	    if(raw_hit->getTimeStamp()<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Get_local_min()&&raw_hit->getTimeStamp()>=0)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_local_min(raw_hit->getTimeStamp());
	    HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Time_Plates(raw_hit->getTimeStamp());
	    HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Time_Plates_perRun(raw_hit->getTimeStamp());
	    //}
    }
  }
  if(Times.size()==0) std::cout<<red<<" 0 hits within the the TriggerTime given... You should verify your TriggerTime or your run is triggerless "<<normal<<std::endl;
  unsigned int  long long global_min=HistoPlanes[0]->Get_local_min();
  unsigned int long long  global_max=HistoPlanes[0]->Get_local_max();
  for(unsigned int i=0;i<HistoPlanes.size();++i) 
  { 
	  HistoPlanes[i]->Set_Total_Time(); 
	  HistoPlanes[i]->Set_Nbrof0Hits();
    if(HistoPlanes[i]->Get_local_max()>global_max) global_max=HistoPlanes[i]->Get_local_max();
	  if(HistoPlanes[i]->Get_local_min()<global_min) global_min=HistoPlanes[i]->Get_local_min();
  }
  UsefullTime->Fill(((fin_RamFull-debut_RamFull)-global_max-global_min)*200e-9);
  HistoPlanes[0]->Set_Global_Total_Time(global_max-global_min);
  if(_LayerCut!=-1)
  {
    FillTimes();
    for(std::map< int,int>::iterator itt=Times.begin(); itt!=Times.end(); ++itt) 
	  {
	    static int iii=0;
	    bool Scintillatorseeittoo=false;
	    int iiii=0;
	    for(unsigned int i=0;i!=ScintillatorCoincidence.size();++i)
	    if(ScintillatorCoincidence[i]-itt->first>=-4&&ScintillatorCoincidence[i]-itt->first<=0)
	    {
	      iiii++;
	      diff.Fill(ScintillatorCoincidence[i]-itt->first);
	      //std::cout<<ScintillatorCoincidence[i]<<"  "<<itt->first<<"  "<<abs(ScintillatorCoincidence[i]-itt->first)<<std::endl;
	      //if(iiii==1) std::cout<<yellow<<++iii<<normal<<std::endl;
	      Scintillatorseeittoo=true;
	    }
	    EventsGrouped.clear();
	    std::map<int,std::vector<RawCalorimeterHit *> >::iterator middle=RawHits.find(itt->first);
	    std::map<int,std::vector<RawCalorimeterHit *> >::iterator after=middle;
	    std::map<int,std::vector<RawCalorimeterHit *> >::iterator before=middle;
	    while(fabs(middle->first-before->first)<=_timeWin && before!=RawHits.begin()) --before;
	    ++before;
	    while(fabs(after->first-middle->first)<=_timeWin && after!=RawHits.end()) ++after;
	    std::map<int,int> nbrPlanestouched;
	      for(middle=before; middle!=after; ++middle ) 
	      {  
	        for(int unsigned i=0; i<(middle->second).size(); ++i) 
		      {
		        int dif_id=((middle->second)[i])->getCellID0() & 0xFF;
		        nbrPlanestouched[geom.GetDifNbrPlate(dif_id)]++;
		      }
	      }
	      if(Scintillatorseeittoo)
	      {
	      //////////////////////////Just Scintillator
	      std::cout<<red<<counttt++<<normal<<std::endl;
	      for(middle=before; middle!=after;++middle ) 
		    {
		      EventsGroupedScin.insert(EventsGroupedScin.end(),middle->second.begin(),middle->second.end());
		    }
		    FillIJK(EventsGroupedScin);
		    
		    //////////////////////////////////////////////
		    ////////////////Scintillator & Timewin etc
		    if(nbrPlanestouched.size()>=(unsigned int)(_LayerCut)) 
	      {
	        EventsSelectedt++;
	        LCEventImpl*  evtt = new LCEventImpl() ;
	        LCCollectionVec* col_eventt = new LCCollectionVec(LCIO::CALORIMETERHIT);
		      CellIDEncoder<CalorimeterHitImpl> cdt( "I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" ,col_eventt) ;
	        col_eventt->setFlag(col_eventt->getFlag()|( 1 << LCIO::RCHBIT_LONG));
	        col_eventt->setFlag(col_eventt->getFlag()|( 1 << LCIO::RCHBIT_TIME));
	        FillIJK(EventsGroupedScin, col_eventt,cdt,0);
	        evtt->addCollection(col_eventt, "SDHCAL_HIT_SC");
	        evtt->setEventNumber(EventsSelectedt);
	        evtt->setTimeStamp(evtP->getTimeStamp());
	        evtt->setRunNumber(evtP->getRunNumber());
	        _EventWriter->writeEvent( evtt ) ;
	        delete evtt;
		    }
		    EventsGroupedScin.clear();
		    }
	      
	      
	      /////////////////////////////////////////////////////////////////////////
	      ////////////////////////////////////////////////////////////////////////
	      if(nbrPlanestouched.size()>=(unsigned int)(_LayerCut)) 
	      {
	        EventsSelected++;
	        for(middle=before; middle!=after; ) 
		      {
		        EventsGrouped.insert(EventsGrouped.end(),middle->second.begin(),middle->second.end());
		        RawHits.erase(middle++);
		      }
	        LCEventImpl*  evt = new LCEventImpl() ;
	        LCCollectionVec* col_event = new LCCollectionVec(LCIO::CALORIMETERHIT);
		      CellIDEncoder<CalorimeterHitImpl> cd( "I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" ,col_event) ;
	        col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
	        col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
	        
	       FillIJK(EventsGrouped, col_event,cd,0);
	        
	        evt->addCollection(col_event, "SDHCAL_HIT");
	        evt->setEventNumber(EventsSelected);
	        evt->setTimeStamp(evtP->getTimeStamp());
	        evt->setRunNumber(evtP->getRunNumber());
	        _EventWriter->writeEvent( evt ) ;
	        delete evt;
	      }
	  }
	}
  else
  {
    EventsSelected++;
    Writer(_EventWriter,"SDHCAL_HIT",RawHits, evtP,EventsSelected,0);
  }
  if(_noiseFileName!=""&&_LayerCut!=-1) 
  {
    EventsNoise++;
    Writer(_NoiseWriter,"SDHCAL_HIT_NOISE_IN_TRIGGER_TIME",RawHits, evtP,EventsNoise,1);
    Writer(_NoiseWriter,"SDHCAL_HIT_NOISE",BehondTrigger, evtP,EventsNoise,1);
  }
  if(_noiseFileName!=""&&_LayerCut==-1)
  {
    EventsNoise++;
    Writer(_NoiseWriter,"SDHCAL_HIT_NOISE",BehondTrigger, evtP,EventsNoise,1);
  }
  if(_WantDistribution==true) for(unsigned int i=0;i<HistoPlanes.size();++i)HistoPlanes[i]->Fill_TH1_Hit_In_Pad_Per_RamFull();
}

void TriventProcessor::end()
{  
  
  std::string name="Results_Trivent_"+ patch::to_string(_NbrRun)+".root";
  TFile *hfile = new TFile(name.c_str(),"RECREATE","Results");
 
  /*for(std::map<int,TGraphTime*>::iterator it=time_graph.begin();it!=time_graph.end();++it)
  {
      std::string name = "a"+to_string(it->first);
      it->second->SetSleepTime(200)
      it->second->Write(name.c_str());
  }*/
  //t->Write();
  TF1 * tf = new TF1("TransferFunction", transfer_function);
  
  for(unsigned int i=0;i<typeee.size();++i)
  {
	  typeee[i]->Write();
    delete typeee[i];
  }
  hs.Write();
  diff.Write();
  hs2.Write();
  int coord[3];
  double total=0;
  double max=0;
  TH3D* h1= hs.Projection(2,1,0);
    TH3D* h2= hss.Projection(2,1,0);
    Make_good_TH3(h1);
    Make_good_TH3(h2);
    h1->GetListOfFunctions()->Add(tf);
    h2->GetListOfFunctions()->Add(tf);
    h1->Write();
    if(pdf)h1->Draw("glcolz");
    std::string namepdf="plots"+patch::to_string(_NbrRun)+".pdf";
    if(pdf)canvas->Print((namepdf+"(").c_str());
    h2->Write();
    if(pdf)h2->Draw("glcolz");
    if(pdf)canvas->Print(namepdf.c_str());
    TH3D* h12= hs2.Projection(2,1,0);
    TH3D* h22= hss2.Projection(2,1,0);
    Make_good_TH3(h12);
    Make_good_TH3(h22);
    h12->GetListOfFunctions()->Add(tf);
    h22->GetListOfFunctions()->Add(tf);
    h12->Write();
    if(pdf)h12->Draw("glcolz");
    if(pdf)canvas->Print(namepdf.c_str());
    if(_Spill_Study)
  {
  timestamp->Write();
  if(pdf){timestamp->Draw();
  canvas->Print(namepdf.c_str());}
  delete timestamp;
  timestamps->Write();
  if(pdf){timestamps->Draw();
  canvas->Print(namepdf.c_str());}
  delete timestamps;
  time2read->Write();
  if(pdf){time2read->Draw();
  canvas->Print(namepdf.c_str());}
  delete time2read;
  time2readtime->Write();
  if(pdf){time2readtime->Draw();
  canvas->Print(namepdf.c_str());}
  delete time2readtime;
  TimeSpill->Write();
  if(pdf){TimeSpill->Draw();
  canvas->Print(namepdf.c_str());}
  delete TimeSpill;
  TimeRamFull->Write();
  if(pdf){TimeRamFull->Draw();
  canvas->Print(namepdf.c_str());}
  delete TimeRamFull;
  UsefullTime->Write();
  if(pdf){UsefullTime->Draw();
  canvas->Print(namepdf.c_str());}
  delete UsefullTime;
  }
    for(unsigned int i=0; i<HistoPlanes.size(); ++i) 
    {
    	HistoPlanes[i]->Save(hfile,namepdf);
    }
    h22->Write();
    h22->Draw("glcolz");
    if(pdf)canvas->Print((namepdf+")").c_str(),"Title:Results");
    delete tf;
    delete h1;
    delete h2;
    delete h12;
    delete h22;
    //delete Branch1;
    //delete Branch2;
    //delete Branch3;
    //delete Branch4;
    //delete Branch5;
    //delete Branch6;
    //delete Branch7;
    //delete Branch8;
    //delete Branch9;
    //delete Branch10;
    //hist->Write();
    //histt->Write();
    //histtt->Write();
    //tf->Write();
    //delete Branch11;
    //delete t;
    hfile->Close();
    delete hfile;
    _EventWriter->close();
    if(_noiseFileName!="") _NoiseWriter->close();
    std::cout << "TriventProcess::end() !! "<<_trig_count<<" Events Trigged"<< std::endl;
    std::cout <<TouchedEvents<<" Events were overlaping "<<"("<<(TouchedEvents*1.0/(TouchedEvents+eventtotal))*100<<"%)"<<std::endl;
    std::cout <<"Total nbr Events : "<<eventtotal<<" Events with nbr of plates >="<<_LayerCut<<" : "<<EventsSelected<<" ("<<EventsSelected*1.0/eventtotal*100<<"%)"<< std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Total Time "<<i+1<<" : "<<HistoPlanes[i]->Get_Total_Time()*200e-9<<"  "; std::cout<<std::endl;
    std::cout<<green<<"Total time of the File : "<<total_time_file<<"s Total time of selected events : "<<total_time_by_user<<"s ; Time global ( Usefull Time ) : "<<HistoPlanes[0]->Get_Global_Total_Time()*200e-9<<"s; In percent : "<<HistoPlanes[0]->Get_Global_Total_Time()*100.0*200e-9/total_time_by_user<<"% of usefull time and "<<total_time_by_user*100.0/total_time_file<<"% of time selected in file"<<normal<<std::endl;
    for(std::map<int,bool>::iterator it=Warningg.begin(); it!=Warningg.end(); it++) std::cout<<red<<"REMINDER::Data from Dif "<<it->first<<" are skipped !"<<normal<<std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Mean noise in plane "<<i+1<<" : "<<HistoPlanes[i]->Get_Means()<<" Hz.cm-2 "; std::cout<<std::endl;
    if(_LayerCut==-1) for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Efficiency "<<i<<" : "<<HistoPlanes[i]->Efficiency()<<"  "; std::cout<<std::endl;
    
    if(hitinit.size()!=0)
    {
      std::cout<<green<<"Counter Front scintillator : "<<_Front_scintillator<<" Counter Back scintillator : "<< _Back_scintillator<<" Counter Both scintillator : "<<_Both_scintillator<<"  "<<normal<<std::endl;
      std::cout<<green<<"Counter Front scintillator : "<<_Front_scintillator*1.0/_efficiencyFrontScintillator<<" Counter Back scintillator : "<< _Back_scintillator*1.0/_efficiencyBackScintillator<<" Counter Both scintillator : "<<_Both_scintillator*1.0/(_efficiencyFrontScintillator*_efficiencyBackScintillator)<<"  "<<normal<<std::endl;
      std::cout<<green<<"Efficiency Calculated with Scintillators : "<<std::endl;
     
      for(unsigned int i=0;i!=hitinit.size();++i)
      {
        std::cout<<"Plate : "<<i+1<<"  "<<EffiwithDiscri[i]*100.0/counttt<<"  "<<normal;
      }
      std::cout<<std::endl;
   
   
    ofstream fichier;
fichier.open("ResultsScinti.txt", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
   if(fichier) { // si l'ouverture a réussi
        fichier<<_NbrRun<<"   ";
        for(unsigned int i=0;i!=hitinit.size();++i)
     {
      fichier<<"Plate : "<<i+1<<"  "<<EffiwithDiscri[i]*100.0/counttt<<"  ";
     }
        fichier<<std::endl;
        fichier.close();  // on referme le fichier
    }
    }
   
   
   
    
    if(_WantCalibration==true)
    {
        std::string calib="Calibration"+patch::to_string(_NbrRun)+".py";
        std::ofstream file(calib,std::ios_base::out); 
        std::cout << "first pass ? enter y"<< std::endl;
    	std::string c;
    	std::cin>>c;
    	bool firstPass=(c=="y" || c=="Y" || c=="Yes");
    	//Name of last DataBase
    	std::cout<<"Number of last DataBase ?"<<std::endl;
    	std::string reponse;
    	std::cin>>reponse;
        std::string OracleDB = "s=oa.OracleAccess(\""+_Database_name+"_"+reponse+"\")";
        if (firstPass)
      	for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i]->Get_Calibration(1,254,10,false);
    	else
      	{
		read_calibration("Calib.txt");
		for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i]->Get_Calibration(0,1,10,true);
      	}
        //for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i]->Get_Calibration();
	//std::ofstream file( "Calibration.py", std::ios_base::out ); 
    	file<<"import OracleAccess as oa"<<std::endl;
    	//file<<"s=oa.OracleAccess(\"T9_AOUT2014_76\")"<<std::endl;
        file<<OracleDB<<std::endl;
    	for(unsigned int i=0;i<HistoPlanes.size();++i)
    	{
      	HistoPlanes[i]->Print_Calibration(file);
    	} 
        file<<"s.uploadChanges()"<<std::endl;
    	file.close();
        save_calibration("Calib.txt");
    }
    if(Negative.size()!=0)
    {
	std::cout<<red<<"WARNING !!! : Negative Value(s) of timeStamp found. They are written in Negative_Values.txt"<<normal<<std::endl;
        std::ofstream fileNeg( "Negative_Values"+patch::to_string(_NbrRun)+".txt", std::ios_base::out ); 
	for(std::map<std::vector<unsigned int>,std::map< int, int>>::iterator it=Negative.begin();it!=Negative.end();++it)
    	{
		fileNeg<<"Dif_Id : "<<it->first[0]<<" Asic_Id : "<<it->first[1]<<" Channel_Id : "<<it->first[2];
                for(std::map< int, int>::iterator itt=it->second.begin();itt!=it->second.end();++itt)fileNeg<<" Value : "<<itt->first<<" - "<<itt->second<<" Times; ";
                fileNeg<<std::endl;
    	}
        fileNeg.close();
    }
    for(std::map<int,HistoPlane*>::iterator it=HistoPlanes.begin();it!=HistoPlanes.end();++it)
    {
	delete(it->second);
    }
    delete canvas;
}


