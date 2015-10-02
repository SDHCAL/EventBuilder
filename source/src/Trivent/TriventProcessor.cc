#include "Trivent/TriventProcessor.h"
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
#include "marlin/Global.h"
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
std::vector<int>EffiwithDiscri;


enum Threshold{Threshold2=1,Threshold1,Threshold3};

using namespace marlin;
std::vector<int>ScintillatorCoincidence;
#define degtorad 0.0174532925
unsigned int EventsNoise=0;

unsigned int EventsSelected=0;
unsigned int EventsSelectedt=0;
unsigned int TOTALNUMBERHITSCINTI=0;
unsigned int _eventNr=0;
#define size_pad 10.4125
#define size_strip 2.5
unsigned long int total_time_file=0;
unsigned long int total_time_by_user=0;
unsigned long int min_by_user=0;
unsigned long int max_by_user=0;
unsigned long long int HistoPlane::global_total_time =0;
//std::map<int,TGraphTime*>time_graph;
int bin[3]={110,110,700};
int bin2[3]={100,100,1000};
double xmin[3]={0,0,0};
double xmax[3]={110,110,700};
double xmin2[3]={0,0,0};
double xmax2[3]={1000,1000,1000};
extern TCanvas* canvas;
THnSparseI hs("Noise", "Noise", 3, bin, xmin, xmax);
THnSparseI hs2("Events", "Events", 3, bin, xmin, xmax);
THnSparseD hss("Noise_2", "Noise_2", 3, bin2, xmin2, xmax2);
THnSparseD hss2("Events_2", "Events_2", 3, bin2, xmin2, xmax2);
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
TH1D* timestamp=nullptr;
TH1D* timestamps=nullptr;
TH1D* time2read=nullptr;
TH1D* time2readtime=nullptr;
TH1D* TimeSpill=new TH1D("Duration Spill","Duration Spill",50000,0,500);
TH1D* TimeRamFull=new TH1D("Duration Frame","Duration Frame",1000,0,1);
TH1D* UsefullTime=new TH1D("UselessTime","UselessTime",1000,0,1);
std::vector<TH1*>typeee;
void TriventProcessor::FillTimes()
{
  bool eraseFirst=false;
  bool nextHasBeenErased=false;
  for (std::map<int,int>::iterator it=Times.begin(); it!= Times.end(); ++it) 
  {
    if (nextHasBeenErased) --it;
    nextHasBeenErased=false;
    bool eraseIt=(it->second<_noiseCut);
    if (!eraseIt) 
    {
      std::map<int,int>::iterator itnext=it;
      ++itnext;
      if (fabs(itnext->first-it->first)<=_timeWin) 
      {
        if (itnext->second >= it->second)	eraseIt=true;
        else 
        {
          Times.erase(itnext);
          nextHasBeenErased=true;
        }
      }
    }
    if (eraseIt) 
    {
      std::map<int,int>::iterator itprev=it;
      --itprev;
      if (it == Times.begin()) eraseFirst=true;
      else 
      {
        Times.erase(it);
        it=itprev;
      }
    }
  }
  if (eraseFirst) Times.erase(Times.begin());
  std::set<int>touched;
  for(std::map< int,int>::iterator firstit=Times.begin(); firstit!=Times.end(); ++firstit) 
  {
    eventtotal++;
    if(firstit!=--(Times.end())) 
    {
      std::map< int,int>::iterator secondit=firstit;
      secondit++;
      if(secondit->first-firstit->first<2*_timeWin)
      {
        streamlog_message(DEBUG,std::cout<<magenta<<secondit->first<<"  "<<firstit->first<<normal<<std::endl;touched.insert(firstit->first);touched.insert(secondit->first);,"";);
      } 
      else streamlog_message(DEBUG,std::cout<<green<<secondit->first<<"  "<<firstit->first<<normal<<std::endl; ,"";);
    }
  }
  TouchedEvents+=touched.size();
  for(std::set< int>::iterator it=touched.begin(); it!=touched.end(); ++it) 
  {
    Times.erase(*it);
  }
  streamlog_message(MESSAGE0,if(touched.size()!=0)std::cout<<touched.size()<<" Events are touched !"<<std::endl; ,"";);
}

void TriventProcessor::FillIJK(std::vector<RawCalorimeterHit *>vec, LCCollectionVec* col,CellIDEncoder<CalorimeterHitImpl>& cd, int IsNoise)
{
  std::vector<std::map<int,int> >Times_Plates_S1;
  std::vector<std::map<int,int> >Times_Plates_S2;
  std::vector<std::map<int,int> >Times_Plates_S3;
  std::vector<std::vector<std::map<int,int> >>Times_Plates{Times_Plates_S1,Times_Plates_S2,Times_Plates_S3};
  for(unsigned int j=0; j<HistoPlanes.size(); ++j) 
  {
    for(unsigned int i=0; i<Times_Plates.size(); ++i)
    {
      Times_Plates[i].emplace_back(std::map<int,int>());
    }
  }
  for(std::vector<RawCalorimeterHit *>::iterator it=vec.begin(); it!=vec.end(); ++it) 
  {
    CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
    int dif_id  = (*it)->getCellID0() & 0xFF ;
    int asic_id = ((*it)->getCellID0() & 0xFF00)>>8;
    int chan_id = ((*it)->getCellID0() & 0x3F0000)>>16;
    int Seuil = ((*it)->getAmplitude());
    double ca=SinCos[dif_id][0];
	  double sa=SinCos[dif_id][1];
    double cb=SinCos[dif_id][2];
	  double sb=SinCos[dif_id][3];
    double cg=SinCos[dif_id][4];
	  double sg=SinCos[dif_id][5];
    /*float ca=cos(geom.GetDifAlpha(dif_id)*degtorad);
    float sa=sin(geom.GetDifAlpha(dif_id)*degtorad);
    float cb=cos(geom.GetDifBeta(dif_id)*degtorad);
    float sb=sin(geom.GetDifBeta(dif_id)*degtorad);
    float cg=cos(geom.GetDifGamma(dif_id)*degtorad);
    float sg=sin(geom.GetDifGamma(dif_id)*degtorad);*/
    unsigned int NbrPlate =geom.GetDifNbrPlate(dif_id)-1;
    float Z= geom.GetPlatePositionZ(NbrPlate);
    cd["Dif_id"]=dif_id;
    cd["Asic_id"]=asic_id;
    cd["Chan_id"]=chan_id;
    caloHit->setTime(float((*it)->getTimeStamp()));
    caloHit->setEnergy(float((*it)->getAmplitude()&3));
    unsigned int K =geom.GetDifNbrPlate(dif_id);
    unsigned int I=0;
    unsigned int J=0;
    if(geom.GetDifType(dif_id)==pad) 
    {
      I =(1+MapILargeHR2[chan_id]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id);
      J =(32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id);
      pos[0] = cg*cb*I*size_pad+(-sg*ca+cg*sb*sa)*J*size_pad+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
      pos[1] = sg*cb*I*size_pad+(cg*ca+sg*sb*sa)*J*size_pad+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
      pos[2] = -sb*I*size_pad+cb*sa*J*size_pad+cb*ca*Z;
    }
    if(geom.GetDifType(dif_id)==positional) 
    {
      if(asic_id%2==0) Z= geom.GetPlatePositionZ(NbrPlate)+2;
      //if((asic_id%2==0&&geom.GetDifUpDown(dif_id)==1)||(asic_id%2==1&&geom.GetDifUpDown(dif_id)==0))
      if(geom.GetDifUpDown(dif_id)==1) 
      {
        I =(2*chan_id)+geom.GetDifPositionX(dif_id);
      } 
      else I =2*(64-chan_id)-1+geom.GetDifPositionX(dif_id);
      J =0;
      pos[0] = cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
      if(asic_id%2==1) 
      {
        pos[0]=cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate)+1;
      }
        pos[1] = sg*cb*I*size_strip+(cg*ca+sg*sb*sa)*J*size_strip+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
        pos[2] = -sb*I*size_strip+cb*sa*J*size_strip+cb*ca*Z;
    }
    Times_Plates[Seuil-1][geom.GetDifNbrPlate(dif_id)-1][(*it)->getTimeStamp()]++;
    if(IsNoise==1) 
    {
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise")->Fill((*it)->getTimeStamp(),1);
      if(Seuil==1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S1")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S2")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S3")->Fill((*it)->getTimeStamp(),1);
    } 
    else if (IsNoise==0) 
    {
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events")->Fill((*it)->getTimeStamp(),1);
      if(Seuil==1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S1")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S2")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S3")->Fill((*it)->getTimeStamp(),1);
    }
    if(IsNoise==1) 
    {
      //TMarker *m = new TMarker(I*10.4125,J*10.4125,21);
      //m->SetMarkerColor(kRed);
      //m->SetMarkerSize(1.4125);
      //time_graph[geom.GetDifNbrPlate(dif_id)-1]->Add(m,_eventNr);
      //time_graph[geom.GetDifNbrPlate(dif_id)-1]->Add(new TPaveLabel(.90,.92,.98,.97,Form("%d",_eventNr),"brNDC"),_eventNr);
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise")->Fill(I,J);
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_Pon")->Fill(I,J,Seuil);
      if(Seuil==1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S1")->Fill(I,J);
      else if (Seuil==2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S2")->Fill(I,J);
      else if (Seuil==3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S3")->Fill(I,J);
      if(geom.GetDifType(dif_id)==positional)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_Asic")->Fill(asic_id,asic_id);
      else HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_Asic")->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
      ///////////////
      if(_WantDistribution==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Hit_In_Pad_Per_RamFull(dif_id,asic_id,chan_id);
	    //////////
      //if(_WantCalibration==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Calibration(dif_id,asic_id,chan_id);
      if(_WantCalibration==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_NumberHitsDistribution(dif_id,asic_id,chan_id);
      //std::cout<<green<<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Get_Calibration(dif_id,asic_id,chan_id)<<normal<<std::endl;
    } 
    else if (IsNoise==0) 
    {
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events")->Fill(I,J);
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_Pon")->Fill(I,J,Seuil);
      if(Seuil==1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S1")->Fill(I,J);
      else if (Seuil==2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S2")->Fill(I,J);
      else if (Seuil==3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S3")->Fill(I,J);
      if(geom.GetDifType(dif_id)==positional)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_Asic")->Fill(asic_id,asic_id);
      else HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_Asic")->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
    }
    cd["I"] = I ;
    cd["J"] = J ;
    cd["K"] = K ;
    //totree.pI=I;
    //totree.pJ=J;
    //totree.pK=K;
    //totree.pX=pos[0];
    //totree.pY=pos[1];
    //totree.pZ=pos[2];
    //totree.pAsic=asic_id;
    //totree.pDifId=dif_id;
    //totree.pAsicChannel=chan_id;
    //totree.pTime=(*it)->getTimeStamp();
    //if(IsNoise==1)totree.pEvent=0;
    //else totree.pEvent=1;
    caloHit->setPosition(pos);
    cd.setCellID( caloHit ) ;
    double fill[3]={double(I),double(J),double(10*K)};
    double fill2[3]={pos[0],pos[1],pos[2]};
    if(IsNoise==1) 
    {            
      hs.Fill(fill,1);
      hss.Fill(fill2,1);
      //hist->Fill(I,J,10*K,1);
    }
    //else histt->Fill(I,J,10*K,1);
    else if (IsNoise==0) 
    {
      hs2.Fill(fill,1);
		  hss2.Fill(fill2,1);
	  }
    //int a,b,c,d;
    //if(Delimiter.find(dif_id)==Delimiter.end()){a=Delimiter[1][0];b=Delimiter[1][1];c=Delimiter[1][2];d=Delimiter[1][3];}
    //else {a=Delimiter[dif_id][0];b=Delimiter[dif_id][1];c=Delimiter[dif_id][2];d=Delimiter[dif_id][3];}
    //std::cout<<Delimiter.size()<<std::endl;
    //std::cout<<Delimiter[dif_id][0]<<"  "<<std::endl;//<<Delimiter[dif_id][1]<<"  "<<Delimiter[dif_id][2]<<"  "<<Delimiter[dif_id][3]<<std::endl;
    //if(a<=I&&b>=I&&c<=J&&d>=J)
    //{
    col->addElement(caloHit);
    //}
    //std::cout<<magenta<<totree.pI<<"  "<<totree.pJ<<"  "<<red<<(*it)->getTimeStamp()<<"  "<<totree.pTime<<normal<<std::endl;
    //t->Fill();
  }
  if(IsNoise==1) 
  {
    for(unsigned int k=0; k<Times_Plates.size(); ++k)for(unsigned int i=0; i<Times_Plates[k].size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[k][i].begin(); it!=Times_Plates[k][i].end(); ++it) 
    {
      HistoPlanes[i]->Return_TH1F("Hits_Distr_Noise")->Fill(it->second,1);
      if(k==0)HistoPlanes[i]->Return_TH1F("Hits_Distr_Noise_S1")->Fill(it->second,1);
      else if (k==1)HistoPlanes[i]->Return_TH1F("Hits_Distr_Noise_S2")->Fill(it->second,1);
      else if (k==2)HistoPlanes[i]->Return_TH1F("Hits_Distr_Noise_S3")->Fill(it->second,1);  
    }
  } 
  else if (IsNoise==0) 
  {
    for(unsigned int k=0; k<Times_Plates.size(); ++k)for(unsigned int i=0; i<Times_Plates[k].size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[k][i].begin(); it!=Times_Plates[k][i].end(); ++it)
    { 
      HistoPlanes[i]->Return_TH1F("Hits_Distr_Events")->Fill(it->second,1);
      if(k==0)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S1")->Fill(it->second,1);
      else if (k==1)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S2")->Fill(it->second,1);
      else if (k==2)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S3")->Fill(it->second,1);  
    }
  }
}

std::vector<bool>hitinit;

void TriventProcessor::FillIJK(std::vector<RawCalorimeterHit *>vec)
{
  std::vector<std::map<int,int> >Times_Plates;
  hitinit.clear();
  for(unsigned int j=0; j<HistoPlanes.size(); ++j) 
  {
    Times_Plates.emplace_back(std::map<int,int>());
    hitinit.push_back(false);
  }
  for(std::vector<RawCalorimeterHit *>::iterator it=vec.begin(); it!=vec.end(); ++it) 
  {
    TOTALNUMBERHITSCINTI++;
    int dif_id  = (*it)->getCellID0() & 0xFF ;
    int asic_id = ((*it)->getCellID0() & 0xFF00)>>8;
    int chan_id = ((*it)->getCellID0() & 0x3F0000)>>16;
    double ca=SinCos[dif_id][0];
	  double sa=SinCos[dif_id][1];
    double cb=SinCos[dif_id][2];
	  double sb=SinCos[dif_id][3];
    double cg=SinCos[dif_id][4];
	  double sg=SinCos[dif_id][5];
    unsigned int NbrPlate =geom.GetDifNbrPlate(dif_id)-1;
    float Z= geom.GetPlatePositionZ(NbrPlate);
    unsigned int K =geom.GetDifNbrPlate(dif_id);
    unsigned int I=0;
    unsigned int J=0;
    if(geom.GetDifType(dif_id)==pad) 
    {
      I =(1+MapILargeHR2[chan_id]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id);
      J =(32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id);
      pos[0] = cg*cb*I*size_pad+(-sg*ca+cg*sb*sa)*J*size_pad+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
      pos[1] = sg*cb*I*size_pad+(cg*ca+sg*sb*sa)*J*size_pad+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
      pos[2] = -sb*I*size_pad+cb*sa*J*size_pad+cb*ca*Z;
    }
    if(geom.GetDifType(dif_id)==positional) 
    {
      if(asic_id%2==0) Z= geom.GetPlatePositionZ(NbrPlate)+2;
      //if((asic_id%2==0&&geom.GetDifUpDown(dif_id)==1)||(asic_id%2==1&&geom.GetDifUpDown(dif_id)==0))
      if(geom.GetDifUpDown(dif_id)==1) 
      {
        I =(2*chan_id)+geom.GetDifPositionX(dif_id);
      } 
      else I =2*(64-chan_id)-1+geom.GetDifPositionX(dif_id);
      J =0;
      pos[0] = cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
      if(asic_id%2==1) 
      {
        pos[0]=cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate)+1;
      }
        pos[1] = sg*cb*I*size_strip+(cg*ca+sg*sb*sa)*J*size_strip+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
        pos[2] = -sb*I*size_strip+cb*sa*J*size_strip+cb*ca*Z;
    }
    Times_Plates[geom.GetDifNbrPlate(dif_id)-1][(*it)->getTimeStamp()]++; 
      int a,b,c,d;
      
	    if(Delimiter.find(dif_id)==Delimiter.end()){a=Delimiter[1][0];b=Delimiter[1][1];c=Delimiter[1][2];d=Delimiter[1][3];}
	    else {a=Delimiter[dif_id][0];b=Delimiter[dif_id][1];c=Delimiter[dif_id][2];d=Delimiter[dif_id][3];}
	    
      if(a<=I&&b>=I&&c<=J&&d>=J)
	    {
	    
	    hitinit[geom.GetDifNbrPlate(dif_id)-1]=true; 
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("EffiScintiOnly")->Fill(I,J);
      
      
      }
    }
  for(unsigned int i=0;i!=hitinit.size();++i)
  {
    if(hitinit[i]==true)
    {
      EffiwithDiscri[i]+=1;
    }
      std::cout<<green<<i<<"  "<<EffiwithDiscri[i]<<"  ";
  }
  std::cout<<normal<<std::endl;
}

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

void TriventProcessor::Writer(IO::LCWriter* file,const char * name,std::map<int,std::vector<EVENT::RawCalorimeterHit *> >& vec, EVENT::LCEvent* event,unsigned int & nbr,unsigned int IsNoise)
{
  LCEventImpl*  evt = new LCEventImpl() ;
  LCCollectionVec* col_event = new LCCollectionVec(LCIO::CALORIMETERHIT);
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
  CellIDEncoder<CalorimeterHitImpl> cd( "I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" ,col_event) ;
  for(std::map<int,std::vector<RawCalorimeterHit *> >::iterator itt=vec.begin(); itt!=vec.end(); ++itt) 
	{
    FillIJK((itt->second),col_event,cd,IsNoise);
  }
  evt->addCollection(col_event,name);
  evt->setEventNumber(nbr);
  evt->setTimeStamp(event->getTimeStamp());
  evt->setRunNumber(event->getRunNumber());
  file->writeEvent( evt ) ;
  delete evt;
}
int counttt=0;
TriventProcessor::~TriventProcessor() {}

void TriventProcessor::processRunHeader( LCRunHeader* run)
{
    LCTOOLS::dumpRunHeader(run);
}
std::map<unsigned long long int,int>Time_stamps;
std::map<unsigned long long int,int>Time_stampss;
void TriventProcessor::init()
{   
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
  if(_Spill_Study)
  {
  
  int32_t timetime=0;
  LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
  std::vector<std::string>colee{"DHCALRawHits"};
  std::string namea="DHCALRawHits";
  lcReader->setReadCollectionNames( colee ) ;
  unsigned eventnumber=-1;
  unsigned long long int min_by_user=999999999999;
  unsigned long long int  max_by_user=1;
  for(unsigned int i=0;i!=LCIOFiles.size();++i)
  {
    std::cout<<"I'm Readind The DATA in "<<LCIOFiles[i]<<std::endl;
    LCEvent* evt(0) ;
    lcReader->open( LCIOFiles[i] ) ;
    _GlobalEvents=lcReader->getNumberOfEvents()-1;
    std::cout<<lcReader->getNumberOfRuns()<<" "<<lcReader->getNumberOfEvents()<<std::endl;
    evt=lcReader->readNextEvent();
    int counter=0;
    do
    {
      counter++;
		  if(counter%1000==0)std::cout<<counter<<std::endl;
      LCCollection* col=evt->getCollection("DHCALRawHits");
      if(col!=nullptr)
      //for(unsigned int hit=0;hit<col->getNumberOfElements();++hit)
	    //{
		  {
        eventnumber++;
        RawCalorimeterHit * myhit = dynamic_cast<RawCalorimeterHit*>(col->getElementAt(/*hit*/0)) ;
        //std::cout<<blue<<myhit->getCellID0()<<"  "<<std::endl;
        unsigned int dif_id=myhit->getCellID0()&0xFF;
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
        Time_stamps[_bcid]++;
        for(unsigned int i=0;i<col->getNumberOfElements();++i)
        {
        myhit=dynamic_cast<RawCalorimeterHit*>(col->getElementAt(i));
        unsigned int dif_id=myhit->getCellID0()&0xFF;
        
        //std::cout<<yellow<<dif_id<<"  "<<geom.GetDifNbrPlate(int(dif_id))<<"  "<<geom.GetDifType(int(dif_id))<<normal<<std::endl;
        if(geom.GetDifType(int(dif_id))<=1)if(geom.GetDifNbrPlate(int(dif_id))!=-1) HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("timestamp")->Fill(_bcid*200e-9,1);
        }
        Time_stampss[_bcid]=col->getNumberOfElements();
        //std::cout<<red<<col->getNumberOfElements()<<normal<<std::endl;
        if(eventnumber>=_skip&&eventnumber<=_maxRecord)
				{
          //std::cout<<eventnumber<<" "<<_skip<<"  "<<_maxRecord<<"  "<<_bcid<<std::endl;
					if(_bcid<min_by_user)
          {
				    min_by_user=_bcid;
			 			//std::cout<<red<<_bcid<<normal<<std::endl;
					}
        	if(_bcid>max_by_user)
					{
		        //std::cout<<green<<_bcid<<normal<<std::endl;
						max_by_user=_bcid;
					}
				}
		  }
		  //}
		  if(counter<=_maxRecord) evt=lcReader->readNextEvent();
      else evt =nullptr;
	  }
	  while(evt!=nullptr);
    lcReader->close();
    delete lcReader ;
  }
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
        //std::cout<<red<<fabs((b-a))*200e-9;
        timebetweentime_toplot.push_back((b-a)*200e-9);
			  //std::cout<<yellow<<ittt->first<<"  "<<itt->first<<"  "<<it->first<<"  "<<(ittt->first-itt->first)<<"  "<<(itt->first-it->first)<<"  "<<timebetweentime<<normal<<std::endl;
        if(timebetweentime<min_betweentime)min_betweentime=timebetweentime;
        if(timebetweentime>max_betweentime)max_betweentime=timebetweentime;
        std::map<double,std::vector<double>>::iterator toadd=Vec_timebetweentime.end();
			  --toadd;
        for(std::map<double,std::vector<double>>::iterator j=Vec_timebetweentime.begin();j!=Vec_timebetweentime.end();++j)
			  {
          //std::cout<<"Vector size : "<<Vec_timebetweentime.size()<<std::endl;
          //std::cout<<"fabs(timebetweentime-j->first) : "<<fabs(timebetweentime-j->first) << "  "<< timebetweentime<< " "<<j->first<<"  "<<std::endl;
          if(fabs(fabs(timebetweentime)-j->first) <=0.50)
				  {
					  (j->second).push_back(ittt->first);
            //std::cout<<"ittt->first : "<<ittt->first<<"  "<<std::endl;
            break;
				  }
          else if(j==toadd&&it!=Time_stamps.begin()) 
          {
					  Vec_timebetweentime[(ittt->first-itt->first)*200e-9].push_back(ittt->first);
					  //std::cout<<red<<"Added to vec "<<(ittt->first-itt->first)*200e-9<<"  "<<ittt->first<<normal<<std::endl;
            skip=true;
   					break;
				  }
			  }
		  }
	  }
	  if(it->first<min&&it->first!=0)min=it->first;
    if(it->first>max)max=it->first;
    //std::cout<<it->first<<std::endl;
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
  //Vec_timebetweentime.clear();
  int diffbetweentime=(int((max_betweentime-min_betweentime))+1)*10;
  time2readtime = new TH1D("Difference in time bettwen two time","Difference in time bettwen two time",10000000,-10,10);
  for(unsigned int i =0;i<timebetweentime_toplot.size();++i)
  {
	  time2readtime->Fill(timebetweentime_toplot[i]);
    //std::cout<<Vec_timebetween[i]<<"  "<<diffbetween<<"  "<<min_between<<"  "<<max_between<<std::endl;
  }
  int diffbetween=(int((max_between-min_between))+1)*100;
  time2read = new TH1D("Distribution time between 2 read out","Distribution time between 2 read out",5000000,min_between,100);
  for(unsigned int i =0;i<Vec_timebetween.size();++i)
  {
	  time2read->Fill(Vec_timebetween[i]);
    //std::cout<<Vec_timebetween[i]<<"  "<<diffbetween<<"  "<<min_between<<"  "<<max_between<<std::endl;
  }
  double diff=(max-min);
  if(diff>9999999)diff=9999999;
  std::cout<<red<<"rrrrrrrrrrrrrrrrrrrrrrrrrr"<<diff<<"  "<<INT_MAX<<normal<<std::endl;
  timestamp= new TH1D("timestamp","timestamp",diff,min*200e-9,max*200e-9);
  timestamps= new TH1D("timestampS","timestampS",diff,min*200e-9,max*200e-9);
  for(std::map<double,std::vector<double>>::iterator j=Vec_timebetweentime.begin();j!=Vec_timebetweentime.end();++j)
  {
    std::string name = "type"+patch::to_string(j->first);
	  typeee.push_back(new TH1D(name.c_str(),name.c_str(),1000000,0.,10000.));
    for(unsigned int i =0;i<j->second.size();++i)
	  {
		  typeee[typeee.size()-1]->Fill(j->second[i]*200e-9);
		  //timestamp->Fill(j->second[i]*200e-9);
	  }
  }

  for(std::map<unsigned long long int ,int>::iterator it=Time_stamps.begin();it!=Time_stamps.end();++it)
  {
	  timestamp->Fill(it->first*200e-9,it->second);
	  
    //std::cout<<yellow<<it->first<<"  "<<it->second<<normal<<std::endl;
  }
  for(std::map<unsigned long long int ,int>::iterator it=Time_stampss.begin();it!=Time_stampss.end();++it)
  {
	  timestamps->Fill(it->first*200e-9,it->second);
    //std::cout<<red<<it->first<<"  "<<it->second<<normal<<std::endl;
  }
  total_time_by_user=(max_by_user-min_by_user)*200e-9;
  total_time_file=(max-min)*200e-9;
  //std::cout<<min<<"   "<<max<<"  "<<(max-min)*200e-9<<std::endl;
  
  }
  _rolling=Every(_maxRecord);
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
  //std::cout<<asctime(localtime(&date))<<std::endl;
  int skip=0;
  if(_skip!=0)skip=_skip+1;
  int maxRecordplusskip=0;
  if(_maxRecord+skip>=_GlobalEvents) 
  {
    maxRecordplusskip=_GlobalEvents;
  }else maxRecordplusskip=_maxRecord+skip;
  if(_maxRecord>=_GlobalEvents)_maxRecord=_GlobalEvents ;
  if(_eventNr %_rolling==0 || _eventNr==_GlobalEvents || _eventNr==maxRecordplusskip || _eventNr==1)
  {
    if(_maxRecord==-1)
	  {
		  std::cout<<red<<"[";
      int percent=int((_eventNr-skip)*100.0/(_GlobalEvents-skip));
		  std::cout<<Shift(percent)<<percent<<"%]"<<normal<<" Event Number : "<<_eventNr<<"/"<<_GlobalEvents<<std::endl;
	  }
    else 
	  {
		  std::cout<<red<<"[";
		  int percent=int((_eventNr-skip)*100.0/(_maxRecord));
		  std::cout<<Shift(percent)<<percent<<"%]"<<normal<<" Event Number : "<<_eventNr<<"/"<<maxRecordplusskip<<" Total : "<<_GlobalEvents<<std::endl;
	  }
  }
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
    if(geom.GetDifType(int(dif_id))<=1)if(to_skip==false)  processCollection(evtP,col);
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
  for(unsigned int i=0; i<HistoPlanes.size(); ++i)HistoPlanes[i]->Clear_Time_Plates_perRun();
  for (int ihit=0; ihit < numElements; ++ihit) 
  {
    RawCalorimeterHit *raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit)) ;
    if (raw_hit != nullptr) 
    {
	    unsigned int dif_id  = (raw_hit)->getCellID0() & 0xFF ;
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
	      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_hit_trigger();
	      Times[raw_hit->getTimeStamp()]++;
	      RawHits[raw_hit->getTimeStamp()].push_back(raw_hit);
	      ////////////////////////////////////////
	      //}
	      /////////////////////////////////////////
	      //if(raw_hit->getTimeStamp()<0)std::cout<<yellow<<raw_hit->getTimeStamp()<<"  "<<((raw_hit)->getCellID0() & 0xFF)<<normal<<std::endl;
	    }
	    else
	    {
	      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_hit_other();
	      BehondTrigger[raw_hit->getTimeStamp()].push_back(raw_hit);
	      //std::cout<<blue<<raw_hit->getTimeStamp()<<"  "<<BehondTrigger.size()<<normal<<std::endl;
	    }
	    if(raw_hit->getTimeStamp()>HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Get_local_max())HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_local_max(raw_hit->getTimeStamp());
	    if(raw_hit->getTimeStamp()<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Get_local_min()&&raw_hit->getTimeStamp()>=0)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Set_local_min(raw_hit->getTimeStamp());
	    HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Time_Plates(raw_hit->getTimeStamp());
	    HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Time_Plates_perRun(raw_hit->getTimeStamp());
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
    canvas->Print((namepdf+")").c_str(),"Title:Results");
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

void TriventProcessor::save_calibration(std::string filename)
{
  std::ofstream f(filename);
  for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i]->SaveCalibration(f);
}


void TriventProcessor::read_calibration(std::string filename)
{
  std::ifstream f(filename);
  int difN,asicN,channelN; double v;
  while (f.good())
    {
      f >> difN >> asicN >> channelN >> v;
      if (f.good() and !f.eof()) {HistoPlanes[geom.GetDifNbrPlate(difN)-1]->Set_Calibration(difN,asicN,channelN,v);}
    }
}
