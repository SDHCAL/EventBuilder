#include "Trivent/TriventProcessor.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "marlin/Processor.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDEncoder.h"
#include <EVENT/LCGenericObject.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "IMPL/CalorimeterHitImpl.h"
#include <IMPL/LCRunHeaderImpl.h>
#include "Colors.h"
#include "TBranch.h"
#include "TObject.h"
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

/*const UInt_t Number = 4;
Double_t Red[Number]   = { 0.0, 1.0,0.0, 1.0 };
Double_t Green[Number] = { 0.0, 0.0,1.0, 1.0 };
Double_t Blue[Number]  = { 1.0, 0.0,0.0, 1.0 };
Double_t Stops[Number] = { 0.0, 0.50,0.99, 1.0 };
int nb=1000000;
Int_t a=TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);*/
using namespace marlin;
#define degtorad 0.0174532925
unsigned int EventsNoise=0;
unsigned int eventtotal=0;
unsigned int EventsSelected=0;
unsigned int TouchedEvents=0;
unsigned int _eventNr=0;
#define size_pad 10.4125
#define size_strip 2.5

unsigned long long int HistoPlane::global_total_time =0;
//std::map<int,TGraphTime*>time_graph;
int bin[3]={110,110,700};
int bin2[3]={100,100,1000};
double xmin[3]={0,0,0};
double xmax[3]={110,110,700};
double xmin2[3]={0,0,0};
double xmax2[3]={1000,1000,1000};
Double_t transfer_function(const Double_t *x, const Double_t * /*param*/)
{
   if (*x <=0)return 0.00;
   else if(*x<=0.1) return 0.1;
   else if(*x<=0.2) return 0.2;
   else if(*x<=0.3) return 0.3;
   else if(*x<=0.4) return 0.4;
   else if(*x<=0.5) return 0.5;
   else if(*x<=0.6) return 0.6;
   else if(*x<=0.7) return 0.7;
   else if(*x<=0.8) return 0.8;
   else if(*x<=0.9) return 0.9;
   else return 1.0;
}

void Make_good_TH3(TH3* t)
{
 for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,TMath::Log(t->GetBinContent(i, j, k)));
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    long int val=0;
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		val += t->GetBinContent(i, j, k);
            		}
        	}
    }
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,t->GetBinContent(i, j, k)*1.0/val);
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    double min=99999;
    double max=-1;
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) 
				{
					if(t->GetBinContent(i, j, k)<min)min=t->GetBinContent(i, j, k);
					if(t->GetBinContent(i, j, k)>max)max=t->GetBinContent(i, j, k);
				}
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
    for(int i = 1; i <= t->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= t->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= t->GetNbinsZ(); ++k)
            		{
                                
                		if(t->GetBinContent(i, j, k)!=0) t->SetBinContent(i, j, k,t->GetBinContent(i, j, k)*1.0/max);
                                //std::cout<<hist->GetBinContent(i, j, k)<<std::endl;
            		}
        	}
    }
}

THnSparseI hs("Noise", "Noise", 3, bin, xmin, xmax);
THnSparseI hs2("Events", "Events", 3, bin, xmin, xmax);
THnSparseD hss("Noise_2", "Noise_2", 3, bin2, xmin2, xmax2);
THnSparseD hss2("Events_2", "Events_2", 3, bin2, xmin2, xmax2);
//TH3F* hist=NULL;
//TH3F* histt=NULL; 
//TH3F* histtt=NULL;
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
std::vector<std::string  > th1 {"Time_Distr","Hits_Distr","Time_Distr_Events","Hits_Distr_Events","Time_Distr_Noise","Hits_Distr_Noise"};
std::vector<std::string> th2 {"Flux_Noise","Flux_Events"};
std::vector<std::string>th2_Asic{"Flux_Noise_Asic","Flux_Events_Asic"};
int _NbrRun=0;
void TriventProcessor::FillTimes()
{
    bool eraseFirst=false;
    bool nextHasBeenErased=false;
    for (std::map<int,int>::iterator it=Times.begin(); it!= Times.end(); ++it) {
        if (nextHasBeenErased) --it;
        nextHasBeenErased=false;
        bool eraseIt=(it->second<_noiseCut);
        if (!eraseIt) {
            std::map<int,int>::iterator itnext=it;
            ++itnext;
            if (fabs(itnext->first-it->first)<=_timeWin) {
                if (itnext->second >= it->second)	eraseIt=true;
                else {
                    Times.erase(itnext);
                    nextHasBeenErased=true;
                }
            }
        }
        if (eraseIt) {
            std::map<int,int>::iterator itprev=it;
            --itprev;
            if (it == Times.begin()) eraseFirst=true;
            else {
                Times.erase(it);
                it=itprev;
            }
        }
    }
    if (eraseFirst) Times.erase(Times.begin());
    std::set<int>touched;
    for(std::map< int,int>::iterator firstit=Times.begin(); firstit!=Times.end(); ++firstit) {
        eventtotal++;
        if(firstit!=--(Times.end())) {
            std::map< int,int>::iterator secondit=firstit;
            secondit++;
            if(secondit->first-firstit->first<2*_timeWin) {
                streamlog_message(DEBUG,std::cout<<magenta<<secondit->first<<"  "<<firstit->first<<normal<<std::endl;
                                  touched.insert(firstit->first);
                                  touched.insert(secondit->first); ,"";);
            } else streamlog_message(DEBUG,std::cout<<green<<secondit->first<<"  "<<firstit->first<<normal<<std::endl; ,"";);
        }
    }
    TouchedEvents+=touched.size();
    for(std::set< int>::iterator it=touched.begin(); it!=touched.end(); ++it) {
        Times.erase(*it);
    }
    streamlog_message(MESSAGE0,if(touched.size()!=0)std::cout<<touched.size()<<" Events are touched !"<<std::endl; ,"";);
}

void TriventProcessor::FillIJK(std::vector<RawCalorimeterHit *>vec, LCCollectionVec* col,CellIDEncoder<CalorimeterHitImpl>& cd, bool IsNoise)
{
    std::vector<std::map<int,int> >Times_Plates;
    for(unsigned int j=0; j<HistoPlanes.size(); ++j) Times_Plates.emplace_back(std::map<int,int>());
    for(std::vector<RawCalorimeterHit *>::iterator it=vec.begin(); it!=vec.end(); ++it) {

        CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
        int dif_id  = (*it)->getCellID0() & 0xFF ;
        
        int asic_id = ((*it)->getCellID0() & 0xFF00)>>8;
        int chan_id = ((*it)->getCellID0() & 0x3F0000)>>16;
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

        if(geom.GetDifType(dif_id)==pad) {
            I =(1+MapILargeHR2[chan_id]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id);
            J =(32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id);

            pos[0] = cg*cb*I*size_pad+(-sg*ca+cg*sb*sa)*J*size_pad+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
            pos[1] = sg*cb*I*size_pad+(cg*ca+sg*sb*sa)*J*size_pad+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
            pos[2] = -sb*I*size_pad+cb*sa*J*size_pad+cb*ca*Z;

        }

        if(geom.GetDifType(dif_id)==positional) {
            if(asic_id%2==0) Z= geom.GetPlatePositionZ(NbrPlate)+2;
            //if((asic_id%2==0&&geom.GetDifUpDown(dif_id)==1)||(asic_id%2==1&&geom.GetDifUpDown(dif_id)==0))
            if(geom.GetDifUpDown(dif_id)==1) {
                I =(2*chan_id)+geom.GetDifPositionX(dif_id);
            } else I =2*(64-chan_id)-1+geom.GetDifPositionX(dif_id);
            J =0;
            pos[0] = cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
            if(asic_id%2==1) {
                pos[0]=cg*cb*I*size_strip+(-sg*ca+cg*sb*sa)*J*size_strip+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate)+1;
            }
            pos[1] = sg*cb*I*size_strip+(cg*ca+sg*sb*sa)*J*size_strip+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
            pos[2] = -sb*I*size_strip+cb*sa*J*size_strip+cb*ca*Z;
        }
         Times_Plates[geom.GetDifNbrPlate(dif_id)-1][(*it)->getTimeStamp()]++;
        if(IsNoise==1) {
              HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Noise")->Fill((*it)->getTimeStamp(),1);
        } else {
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH1F("Time_Distr_Events")->Fill((*it)->getTimeStamp(),1);
        }
        if(IsNoise==1) {
            //TMarker *m = new TMarker(I*10.4125,J*10.4125,21);
            //m->SetMarkerColor(kRed);
            //m->SetMarkerSize(1.4125);
            //time_graph[geom.GetDifNbrPlate(dif_id)-1]->Add(m,_eventNr);
            //time_graph[geom.GetDifNbrPlate(dif_id)-1]->Add(new TPaveLabel(.90,.92,.98,.97,Form("%d",_eventNr),"brNDC"),_eventNr);
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise")->Fill(I,J);
            if(geom.GetDifType(dif_id)==positional)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_Asic")->Fill(asic_id,asic_id);
            else HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Noise_Asic")->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
            ///////////////
            if(_WantDistribution==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Hit_In_Pad_Per_RamFull(dif_id,asic_id,chan_id);
	    //////////
            //if(_WantCalibration==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_Calibration(dif_id,asic_id,chan_id);
              if(_WantCalibration==true)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Fill_NumberHitsDistribution(dif_id,asic_id,chan_id);
            //std::cout<<green<<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Get_Calibration(dif_id,asic_id,chan_id)<<normal<<std::endl;
        } else {
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1]->Return_TH2F("Flux_Events")->Fill(I,J);
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
        else 
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
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) HistoPlanes[i]->Return_TH1F("Hits_Distr_Noise")->Fill(it->second,1);
    } else 
    {
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) HistoPlanes[i]->Return_TH1F("Hits_Distr_Events")->Fill(it->second,1);
    }
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
    //_Delimiters="";
    //registerProcessorParameter("Delimiters" ,"Delimiters",_Delimiters,_Delimiters);
    
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

TriventProcessor::~TriventProcessor() {}

void TriventProcessor::processRunHeader( LCRunHeader* run)
{
    LCTOOLS::dumpRunHeader(run);
}

void TriventProcessor::init()
{
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
    ReaderFactory readerFactory;
    Reader* myReader = readerFactory.CreateReader(_ReaderType);
    if(myReader) {
        myReader->Read(_FileNameGeometry,geom);
        geom.PrintGeom();
        //hist=new TH3F("g","g",110,0,110,110,0,110,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        //histt=new TH3F("gg","gg",110,0,110,110,0,110,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        //histt=new TH3F("fg","fg",100,0,100,100,0,100,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        //histtt=new TH3F("ffg","ffg",100,0,100,100,0,100,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        //histt->GetListOfFunctions()->Add(tf);
        //histtt->GetListOfFunctions()->Add(tf);
        //TF1 * tf = new TF1("TransferFunction", transfer_function);
        //hist->GetListOfFunctions()->Add(tf);
        //hs.GetListOfFunctions()->Add(tf);
        //histt->GetListOfFunctions()->Add(tf);
        //hs2.GetListOfFunctions()->Add(tf);
        std::map<int, Dif > Difs=geom.GetDifs();
        //std::map<int,int> PlansType;
        unsigned int NbrPlate =0;
        for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it) 
	{
            if(geom.GetDifType(it->first)!=temporal) 
	    {
                SinCos[it->first]=std::vector<double>{cos(geom.GetDifAlpha(it->first)*degtorad),sin(geom.GetDifAlpha(it->first)*degtorad),cos(geom.GetDifBeta(it->first)*degtorad),sin(geom.GetDifBeta(it->first)*degtorad),cos(geom.GetDifGamma(it->first)*degtorad),sin(geom.GetDifGamma(it->first)*degtorad)};
                //PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
                NbrPlate=geom.GetDifNbrPlate(it->first)-1;
                if(HistoPlanes.find(NbrPlate)==HistoPlanes.end()) 
                {
		  //HistoPlane a(_WantDistribution,NbrPlate,geom.GetDifsInPlane(NbrPlate),geom.GetSizeX(NbrPlate),geom.GetSizeY(NbrPlate),th1,th2,th2_Asic);
                  HistoPlanes.insert(std::make_pair(NbrPlate, new HistoPlane(_WantDistribution,NbrPlate,geom.GetDifsInPlane(NbrPlate),geom.GetSizeX(NbrPlate),geom.GetSizeY(NbrPlate),th1,th2,th2_Asic)));
                  //time_graph.insert(std::make_pair(NbrPlate,new TGraphTime(10000,0,0,600,600)));
		}
               
            }
        }
        //FillDelimiter(_Delimiters,PlansType.size());
    } else {
        std::cout << "Reader type n'existe pas !!" << std::endl;
        std::exit(1);
    }
    delete myReader;
}

void TriventProcessor::processEvent( LCEvent * evtP )
{
  if (nullptr == evtP) return;

  _NbrRun=evtP->getRunNumber();
  _eventNr=evtP->getEventNumber();
  if(_eventNr %1000 ==0)std::cout<<"Event Number : "<<_eventNr<<std::endl;

  LCCollection* col2 = evtP ->getCollection("DHCALRawTimes");
  RawTimeDifs.clear();
  if(col2!=nullptr) 
    {
      for (int ihit=0; ihit < col2->getNumberOfElements(); ++ihit) 
	{
	  EVENT::CalorimeterHit* raw_time = dynamic_cast<EVENT::CalorimeterHit* >( col2->getElementAt(ihit)) ;
	  // std::cout<<raw_time->getTime()<<"  "<<raw_time->getEnergyError()<<std::endl;
	  //RawTimeDifs[raw_time->getTimeStamp()].push_back(raw_time);
	}
    }


  for(unsigned int i=0; i< _hcalCollections.size(); i++) {
    LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
    if(col2 == nullptr || col==nullptr) 
      {
	std::cout<< "TRIGGER SKIPED ..."<<std::endl;
	_trig_count++;
	break;
      }
    processCollection(evtP,col);
  }
} 
    
void TriventProcessor::processCollection(EVENT::LCEvent *evtP,LCCollection* col)
{
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
      if (raw_hit != nullptr) {
	unsigned int dif_id  = (raw_hit)->getCellID0() & 0xFF ;
	if(geom.GetDifNbrPlate(dif_id)==-1) {
	  if(Warningg[dif_id]!=true) {
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
  HistoPlanes[0]->Set_Global_Total_Time(global_max-global_min);
  if(_LayerCut!=-1)
    {
      FillTimes();
      
      for(std::map< int,int>::iterator itt=Times.begin(); itt!=Times.end(); ++itt) 
	{
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
	      col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
	      col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
	      CellIDEncoder<CalorimeterHitImpl> cd( "I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" ,col_event) ;
	      FillIJK(EventsGrouped, col_event,cd,0);
	      evt->addCollection(col_event, "SDHCAL_HIT");
	      evt->setEventNumber(EventsSelected);
	      evt->setTimeStamp(evtP->getTimeStamp());
	      evt->setRunNumber(evtP->getRunNumber());
	      _EventWriter->writeEvent( evt ) ;
	      delete evt;
	    }
	}}
  else
    {
      EventsSelected++;
      Writer(_EventWriter,"SDHCAL_HIT",RawHits, evtP,EventsSelected,0);
    }
  
  
  if(_noiseFileName!=""&&_LayerCut!=-1) {
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
        it->second->SetSleepTime(200);
        it->second->Write(name.c_str());
    }*/
    //t->Write();
   
       TF1 * tf = new TF1("TransferFunction", transfer_function);
    hs.Write();
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
    h2->Write();
    /*    for(unsigned int i=0;i!=hs2.GetNbins();++i)
    {
         double value=TMath::Log(hs2.GetBinContent(i,coord));
	 hs2.SetBinContent(coord,value);
    }
    for(unsigned int i=0;i!=hss2.GetNbins();++i)
    {
         double value=TMath::Log(hss2.GetBinContent(i,coord));
	 hss2.SetBinContent(coord,value);
    }*/
    TH3D* h12= hs2.Projection(2,1,0);
    TH3D* h22= hss2.Projection(2,1,0);
    Make_good_TH3(h12);
    Make_good_TH3(h22);
    h12->GetListOfFunctions()->Add(tf);
   h22->GetListOfFunctions()->Add(tf);
    h12->Write();
   h22->Write();
    delete tf;
    delete h1;
    delete h2;
    delete h12;
    delete h22;

    for(unsigned int i=0; i<HistoPlanes.size(); ++i) 
    {
    	HistoPlanes[i]->Save(hfile);
    }
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
    std::cout<<green<<"Time global : "<<HistoPlanes[0]->Get_Global_Total_Time()*200e-9<<normal<<std::endl;
    for(std::map<int,bool>::iterator it=Warningg.begin(); it!=Warningg.end(); it++) std::cout<<red<<"REMINDER::Data from Dif "<<it->first<<" are skipped !"<<normal<<std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Mean noise in plane "<<i+1<<" : "<<HistoPlanes[i]->Get_Means()<<" Hz.cm-2 "; std::cout<<std::endl;
    if(_LayerCut==-1) for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Efficiency "<<i<<" : "<<HistoPlanes[i]->Efficiency()<<"  "; std::cout<<std::endl;
    
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
        std::ofstream fileNeg( "Negative_Values.txt", std::ios_base::out ); 
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
