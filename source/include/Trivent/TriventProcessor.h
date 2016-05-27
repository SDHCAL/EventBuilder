#ifndef TRIVENT_PROCESSOR
#define TRIVENT_PROCESSOR
#include <map>
#include <vector>
#include "marlin/Processor.h"
#include "EVENT/RawCalorimeterHit.h"
#include "Geometry/Geometry.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/CellIDEncoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <EVENT/LCRunHeader.h>
#include "UTIL/LCTOOLS.h"
#include "Trivent/HistoPlane.h"
#include "IMPL/LCEventImpl.h"
#include <set>
#include <string>
#include "Trivent/Mapping.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#define size_pad 10.4125
#define size_strip 2.5

int bin[3]={31,31,70};
int bin2[3]={310,310,100};
double xmin[3]={0,0,0};
double xmax[3]={31,31,100};
double xmin2[3]={0,0,0};
double xmax2[3]={310,310,100};
extern TCanvas* canvas;
THnSparseI hs("Noise", "Noise", 3, bin, xmin, xmax);
THnSparseI hs2("Events", "Events", 3, bin, xmin, xmax);
THnSparseD hss("Noise_2", "Noise_2", 3, bin2, xmin2, xmax2);
THnSparseD hss2("Events_2", "Events_2", 3, bin2, xmin2, xmax2);
enum Threshold{Threshold_2=1,Threshold_1,Threshold_3};

class TriventProcessor : public marlin::Processor
{
public:
    TriventProcessor();
    ~TriventProcessor();
    marlin::Processor *newProcessor()
    {
        return new TriventProcessor();
    }
    void init();
    void Writer(IO::LCWriter* file,const char* name,std::map<int,std::vector<EVENT::RawCalorimeterHit *>> & vec,EVENT::LCEvent * event,unsigned int& nbr,unsigned int IsNoise);
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
    inline void FillTimes();
    void FillIJK(std::vector<RawCalorimeterHit *>vec, LCCollectionVec* col,CellIDEncoder<CalorimeterHitImpl>& cd,int IsNoise);
    void CalculateEfficiencyScinti(std::vector<RawCalorimeterHit *>vec);
private:
    void processCollection(EVENT::LCEvent *evtP,LCCollection* col);    
protected:
    Geometry geom;
    std::map<int,std::vector<double>>SinCos;
    std::string _FileNameGeometry;
    std::string _ReaderType;
    std::string _noiseFileName;
    std::string _outFileName;
    std::vector<std::string> _hcalCollections;
    LCWriter* _EventWriter;
    LCWriter* _NoiseWriter;
    unsigned int TouchedEvents;
    int _timeWin;
    int _TimeWin_Noise;
    int _Time_from_Track;
    unsigned int _eventNr;
    int _trig_count;
    int _noiseCut;
    int _LayerCut;
    int _TriggerTime;
    unsigned int _GlobalEvents;
    std::string _Delimiters;
    unsigned int _maxRecord;
    unsigned int eventtotal;
    unsigned int _skip;
    double _efficiencyFrontScintillator;
    int _IgnorebeginningSpill;
    double _efficiencyBackScintillator;
    bool _WantDistribution;
    bool _WantCalibration;
    bool _Spill_Study;
    unsigned int TOTALNUMBERHITSCINTI;
    std::vector<int>EffiwithDiscri;
    std::string _Database_name;
    std::map<int,int>Times;
    unsigned int _Front_scintillator;
    unsigned int _Back_scintillator;
    unsigned int _Both_scintillator;
    //std::vector<std::map<int,int> >Times_Plates_perRun;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawHits;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > BehondTrigger;
    std::vector<RawCalorimeterHit *>EventsGrouped;
    std::vector<RawCalorimeterHit *>EstimationNoiseAfter;
    std::vector<RawCalorimeterHit *>EstimationNoiseBefore;
    std::vector<RawCalorimeterHit *>EventsGroupedScin;
    std::map< int,std::vector<EVENT::RawCalorimeterHit*> > RawTimeDifs;
    float pos[3];
    std::map<int,bool>Warningg;
    std::map<std::vector<unsigned int>,std::map< int, int>>Negative;
    std::map<int,HistoPlane*>HistoPlanes;
    void save_calibration(std::string filename);
    void read_calibration(std::string filename);
    std::map<int ,std::vector<double> >Delimiter;
 private:
    //helper methods
    void processEvent_ProcessCollectionNamed_DHCALRawTimes(LCEvent *evtP);
    void processEvent_ProcessCollectionNamed_Scintillator(LCEvent *evtP);
    void Loop_On_Scintillator_Collection_unused(LCCollection *col3);
    void Loop_On_Scintillator_Collection(LCCollection *col3);
    unsigned int getDifId_of_first_hit_in_collection(LCCollection* col);
    bool data_is_at_beginning_of_spill(LCCollection* col, unsigned int dif_id);
    unsigned long long bcid_of_current_spill=0;
    
    static unsigned long long Shift_by_one_long;

};





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
      if(Seuil==Threshold_1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S1")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==Threshold_2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S2")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==Threshold_3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise_S3")->Fill((*it)->getTimeStamp(),1);
    } 
    else if (IsNoise==0) 
    {
      HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events")->Fill((*it)->getTimeStamp(),1);
      if(Seuil==Threshold_1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S1")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==Threshold_2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S2")->Fill((*it)->getTimeStamp(),1);
      else if (Seuil==Threshold_3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events_S3")->Fill((*it)->getTimeStamp(),1);
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
      if(Seuil==Threshold_1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S1")->Fill(I,J);
      else if (Seuil==Threshold_2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S2")->Fill(I,J);
      else if (Seuil==Threshold_3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_S3")->Fill(I,J);
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
      if(Seuil==Threshold_1)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S1")->Fill(I,J);
      else if (Seuil==Threshold_2)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S2")->Fill(I,J);
      else if (Seuil==Threshold_3)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events_S3")->Fill(I,J);
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
      HistoPlanes[i]->Return_TH1F("Hits_Distr")->Fill(it->second,1);
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
      HistoPlanes[i]->Return_TH1F("Hits_Distr")->Fill(it->second,1);
      if(k==0)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S1")->Fill(it->second,1);
      else if (k==1)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S2")->Fill(it->second,1);
      else if (k==2)HistoPlanes[i]->Return_TH1F("Hits_Distr_Events_S3")->Fill(it->second,1);  
    }
  }
}


void TriventProcessor::CalculateEfficiencyScinti(std::vector<RawCalorimeterHit *>vec)
{
  unsigned roll=0;
  for(std::vector<RawCalorimeterHit *>::iterator it=vec.begin(); it!=vec.end(); ++it) 
  {
    
    TOTALNUMBERHITSCINTI++;
    int dif_id  = (*it)->getCellID0() & 0xFF ;
    int asic_id = ((*it)->getCellID0() & 0xFF00)>>8;
    int chan_id = ((*it)->getCellID0() & 0x3F0000)>>16;
    unsigned int NbrPlate =geom.GetDifNbrPlate(dif_id)-1;
    float Z= geom.GetPlatePositionZ(NbrPlate);
    unsigned int K =geom.GetDifNbrPlate(dif_id);
    unsigned int I=0;
    unsigned int J=0;
    if(geom.GetDifType(dif_id)==pad) 
    {
      I =(1+MapILargeHR2[chan_id]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id);
      J =(32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id);
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
    }
      unsigned int a,b,c,d;
      
	    if(Delimiter.find(dif_id)==Delimiter.end()){a=Delimiter[1][0];b=Delimiter[1][1];c=Delimiter[1][2];d=Delimiter[1][3];}
	    else {a=Delimiter[dif_id][0];b=Delimiter[dif_id][1];c=Delimiter[dif_id][2];d=Delimiter[dif_id][3];}
	    
      if(a<=I&&b>=I&&c<=J&&d>=J)
	    {
	   
	    roll++;
	    if(roll==1)EffiwithDiscri[geom.GetDifNbrPlate(dif_id)-1]+=1;
     // if(roll==1)std::cout<<geom.GetDifNbrPlate(dif_id)-1<<"  "<<EffiwithDiscri[geom.GetDifNbrPlate(dif_id)-1]<<std::endl;
      
      }
    }
    }

void TriventProcessor::Writer(IO::LCWriter* file,const char * name,std::map<int,std::vector<EVENT::RawCalorimeterHit *> >& vec, EVENT::LCEvent* event,unsigned int & nbr,unsigned int IsNoise)
{
  LCEventImpl*  evt = new LCEventImpl() ;
  LCCollectionVec* col_event = new LCCollectionVec(LCIO::CALORIMETERHIT);
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
  CellIDEncoder<CalorimeterHitImpl> cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7",col_event) ;
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
#endif
