#include "Trivent/TriventProcessor.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "marlin/Processor.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDEncoder.h"
#include "TGraph2D.h"
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
#include "TGraph.h"
#include "Trivent/HistoPlane.h"
#include "TStyle.h"
#include "TF1.h"
#include "TObject.h"
#include "TList.h"
#include "TH3.h"
#include "TColor.h"
#include "TMath.h"
using namespace marlin;
#define degtorad 0.0174532925
unsigned int EventsNoise=0;
unsigned int eventtotal=0;
unsigned int EventsSelected=0;
unsigned int TouchedEvents=0;
unsigned int _eventNr=0;
#define size_pad 10.4125
#define size_strip 2.5
std::map<int,bool>Warningg;
std::map<std::vector<int>,std::map<int,int>>Negative;
Double_t my_transfer_function(const Double_t *x, const Double_t * /*param*/)
{
   if (*x <=0)return 0.00;
   else return 0.99;
}

TF1* tf =new TF1("TransferFunction", my_transfer_function);

TH3F* hist=NULL;
TH3F* histt=NULL; 
TH3F* histtt=NULL;
class ToTree
{
public:
int pI,pJ,pK,pAsic,pDifId,pAsicChannel;
unsigned long long int pTime;
double pX,pY,pZ;
bool pEvent;
};
std::map<std::vector<int>,double>Calibration;
ToTree totree;
std::string name="Tree";
TTree* t= new TTree(name.c_str(), name.c_str());
TBranch* Branch1 =  t->Branch("X",&(totree.pX));
TBranch* Branch2 =  t->Branch("Y",&(totree.pY));
TBranch* Branch3 =  t->Branch("Z",&(totree.pZ));
TBranch* Branch4 =  t->Branch("I",&(totree.pI));
TBranch* Branch5 =  t->Branch("J",&(totree.pJ));
TBranch* Branch6 =  t->Branch("K",&(totree.pK));
TBranch* Branch7 =  t->Branch("Time",&(totree.pTime));
TBranch* Branch8 =  t->Branch("Asic",&(totree.pAsic));
TBranch* Branch9 =  t->Branch("DifId",&(totree.pDifId));
TBranch* Branch10 =  t->Branch("AsicChannel",&(totree.pAsicChannel));
TBranch* Branch11 = t->Branch("Event",&(totree.pEvent));
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
              HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH1F("Time_Distr_Noise")->Fill((*it)->getTimeStamp(),1);
        } else {
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH1F("Time_Distr_Events")->Fill((*it)->getTimeStamp(),1);
        }
        if(IsNoise==1) {
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Noise")->Fill(I,J);
            if(geom.GetDifType(dif_id)==positional)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Noise_Asic")->Fill(asic_id,asic_id);
            else HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Noise_Asic")->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Fill_Calibration(dif_id,asic_id,chan_id);
            //std::cout<<green<<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Get_Calibration(dif_id,asic_id,chan_id)<<normal<<std::endl;
        } else {
            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Events")->Fill(I,J);
            if(geom.GetDifType(dif_id)==positional)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Events_Asic")->Fill(asic_id,asic_id);
            else HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Return_TH2F("Flux_Events_Asic")->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
        }
         
        cd["I"] = I ;
        cd["J"] = J ;
        cd["K"] = K ;
        totree.pI=I;
        totree.pJ=J;
        totree.pK=K;
        totree.pX=pos[0];
        totree.pY=pos[1];
        totree.pZ=pos[2];
        totree.pAsic=asic_id;
        totree.pDifId=dif_id;
        totree.pAsicChannel=chan_id;
        totree.pTime=(*it)->getTimeStamp();
        if(IsNoise==1)totree.pEvent=0;
        else totree.pEvent=1;
        caloHit->setPosition(pos);
        cd.setCellID( caloHit ) ;
        if(IsNoise==1) 
        {
        	hist->Fill(I,J,10*K,1);
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
        t->Fill();

    }
    if(IsNoise==1) 
    {
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) HistoPlanes[i].Return_TH1F("Hits_Distr_Noise")->Fill(it->second,1);
    } else 
    {
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) HistoPlanes[i].Return_TH1F("Hits_Distr_Events")->Fill(it->second,1);
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
    if(_LayerCut==-1){std::cout<<red<<"LayerCut set to -1, assuming that you want use trigger to see events"<<normal<<std::endl;}
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
        hist=new TH3F("g","g",100,0,100,100,0,100,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        histt=new TH3F("fg","fg",100,0,100,100,0,100,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        histtt=new TH3F("ffg","ffg",100,0,100,100,0,100,10*geom.GetNumberPlates()+1,10,10*geom.GetNumberPlates()+1);
        histt->GetListOfFunctions()->Add(tf);
        histtt->GetListOfFunctions()->Add(tf);
        
        std::map<int, Dif > Difs=geom.GetDifs();
        std::map<int,int> PlansType;

        for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it) {
            if(geom.GetDifType(it->first)!=temporal) {
                SinCos[it->first]=std::vector<double>{cos(geom.GetDifAlpha(it->first)*degtorad),sin(geom.GetDifAlpha(it->first)*degtorad),cos(geom.GetDifBeta(it->first)*degtorad),sin(geom.GetDifBeta(it->first)*degtorad),cos(geom.GetDifGamma(it->first)*degtorad),sin(geom.GetDifGamma(it->first)*degtorad)};
                PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
                unsigned int NbrPlate=geom.GetDifNbrPlate(it->first)-1;
                HistoPlanes.insert(std::pair<int,HistoPlane>(NbrPlate,HistoPlane(NbrPlate,geom.GetSizeX(NbrPlate),geom.GetSizeY(NbrPlate),th1,th2,th2_Asic)));
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
    _NbrRun=evtP->getRunNumber();
    if (evtP != NULL) {
        _eventNr=evtP->getEventNumber();
        if(_eventNr %1000 ==0)std::cout<<"Event Number : "<<_eventNr<<std::endl;
        for(unsigned int i=0; i< _hcalCollections.size(); i++) {
            Times.clear();
            RawHits.clear();
            for(unsigned int i =0;i<HistoPlanes.size();++i)HistoPlanes[i].Init_local_min_max();
            BehondTrigger.clear();
            RawTimeDifs.clear();
            LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
            LCCollection* col2 = evtP ->getCollection("DHCALRawTimes");
            if(col2 == NULL || col==NULL) 
            {
                std::cout<< "TRIGGER SKIPED ..."<<std::endl;
                _trig_count++;
                break;
            }
	    if(col2!=NULL) 
	    {
                for (int ihit=0; ihit < col2->getNumberOfElements(); ++ihit) 
		{
                    EVENT::CalorimeterHit* raw_time = dynamic_cast<EVENT::CalorimeterHit* >( col2->getElementAt(ihit)) ;
                   // std::cout<<raw_time->getTime()<<"  "<<raw_time->getEnergyError()<<std::endl;
                    //RawTimeDifs[raw_time->getTimeStamp()].push_back(raw_time);
                }
            }
            int numElements = col->getNumberOfElements();
            for(unsigned int i=0; i<HistoPlanes.size(); ++i)HistoPlanes[i].Clear_Time_Plates_perRun();
            for (int ihit=0; ihit < numElements; ++ihit) 
            {
                RawCalorimeterHit *raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit)) ;
                if (raw_hit != NULL) {
                    int dif_id  = (raw_hit)->getCellID0() & 0xFF ;
                    if(geom.GetDifNbrPlate(dif_id)==-1) {
                        if(Warningg[dif_id]!=true) {
                            Warningg[dif_id]=true;
                            std::cout<<"Please add DIF "<<dif_id<<" to your geometry file; I'm Skipping its data."<<std::endl;
                        }
                        continue;
                    }
                    if(raw_hit->getTimeStamp()<0)
                    {
			std::vector<int>b{dif_id,(raw_hit->getCellID0() & 0xFF00)>>8,(raw_hit->getCellID0() & 0xFF00)>>16}; Negative[b][raw_hit->getTimeStamp()]++;
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
                      	HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Set_hit_trigger();
		     	Times[raw_hit->getTimeStamp()]++;
                         RawHits[raw_hit->getTimeStamp()].push_back(raw_hit);
                     ////////////////////////////////////////
              	     //}
                    /////////////////////////////////////////
                    //if(raw_hit->getTimeStamp()<0)std::cout<<yellow<<raw_hit->getTimeStamp()<<"  "<<((raw_hit)->getCellID0() & 0xFF)<<normal<<std::endl;
                    }
                    else
                    {
                       HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Set_hit_other();
                       BehondTrigger[raw_hit->getTimeStamp()].push_back(raw_hit);
                       //std::cout<<blue<<raw_hit->getTimeStamp()<<"  "<<BehondTrigger.size()<<normal<<std::endl;
                    }
                    if(raw_hit->getTimeStamp()>HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Get_local_max())HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Set_local_max(raw_hit->getTimeStamp());
                    if(raw_hit->getTimeStamp()<HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Get_local_min()&&raw_hit->getTimeStamp()>=0)HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Set_local_min(raw_hit->getTimeStamp());
                    HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Fill_Time_Plates(raw_hit->getTimeStamp());
	            HistoPlanes[geom.GetDifNbrPlate(dif_id)-1].Fill_Time_Plates_perRun(raw_hit->getTimeStamp());
                    
                }
            }
            if(Times.size()==0) std::cout<<red<<" 0 hits within the the TriggerTime given... You should verify your TriggerTime or your run is triggerless "<<normal<<std::endl;
            for(unsigned int i=0;i<HistoPlanes.size();++i) { HistoPlanes[i].Set_Total_Time(); HistoPlanes[i].Set_Nbrof0Hits();}
           
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
        }
    }
}

void TriventProcessor::end()
{  
    std::ofstream file( "Calibration.py", std::ios_base::out ); 
    std::string name="Results_Trivent_"+ std::to_string( (long long int) _NbrRun)+".root";
    TFile *hfile = new TFile(name.c_str(),"RECREATE","Results");
    t->Write();
    for(unsigned int i=0; i<HistoPlanes.size(); ++i) 
    {
    	HistoPlanes[i].Save(hfile);
    }
    delete Branch1;
    delete Branch2;
    delete Branch3;
    delete Branch4;
    delete Branch5;
    delete Branch6;
    delete Branch7;
    delete Branch8;
    delete Branch9;
    delete Branch10;
    hist->Write();
   
    for(int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        	for(int j = 1; j <= hist->GetNbinsY(); ++j)
        	{
            		for(int k = 1; k <= hist->GetNbinsZ(); ++k)
            		{
                                
                		float val = hist->GetBinContent(i, j, k);
                		 if(val>=1000) histt->SetBinContent(i, j, k, val);
            		}
        	}
    }   
    histt->Write();
    //histtt->Write();
    //tf->Write();
    delete Branch11;
    //delete t;
    hfile->Close();
    delete hfile;
    _EventWriter->close();
    if(_noiseFileName!="") _NoiseWriter->close();
    std::cout << "TriventProcess::end() !! "<<_trig_count<<" Events Trigged"<< std::endl;
    std::cout <<TouchedEvents<<" Events were overlaping "<<"("<<(TouchedEvents*1.0/(TouchedEvents+eventtotal))*100<<"%)"<<std::endl;
    std::cout <<"Total nbr Events : "<<eventtotal<<" Events with nbr of plates >="<<_LayerCut<<" : "<<EventsSelected<<" ("<<EventsSelected*1.0/eventtotal*100<<"%)"<< std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Total Time "<<i+1<<" : "<<HistoPlanes[i].Get_Total_Time()*200e-9<<"  "; std::cout<<std::endl;
    for(std::map<int,bool>::iterator it=Warningg.begin(); it!=Warningg.end(); it++) std::cout<<red<<"REMINDER::Data from Dif "<<it->first<<" are skipped !"<<normal<<std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Mean noise in plane "<<i+1<<" : "<<HistoPlanes[i].Get_Means()<<" Hz.cm-2 "; std::cout<<std::endl;
    if(_LayerCut==-1) for(unsigned int i=0;i<HistoPlanes.size();++i)std::cout <<"Efficiency "<<i<<" : "<<HistoPlanes[i].Efficiency()<<"  "; std::cout<<std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i].Get_Flux();
    file<<"import OracleAccess as oa"<<std::endl;
    file<<"s=oa.OracleAccess(\"T9_AOUT2014_76\")"<<std::endl;
    for(unsigned int i=0;i<HistoPlanes.size();++i) HistoPlanes[i].Print_Calibration(file);
    if(Negative.size()!=0)
    {
	std::cout<<red<<"WARNING !!! : Negative Value(s) of timeStamp found"<<normal<<std::endl;
	for(std::map<std::vector<int>,std::map<int,int>>::iterator it=Negative.begin();it!=Negative.end();++it)
    	{
		std::cout<<red<<"Dif_Id : "<<it->first[0]<<" Asic_Id : "<<it->first[1]<<" Channel_Id : "<<it->first[2]<<normal;
                for(std::map<int,int>::iterator itt=it->second.begin();itt!=it->second.end();++itt)std::cout<<" Value : "<<itt->first<<","<<itt->second<<" Times; ";
                std::cout<<std::endl;
    	}
    }
    file.close();
}
