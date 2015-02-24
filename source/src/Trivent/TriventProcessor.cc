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
#include "TH3F.h"
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
#include "Reader/ReaderFactory.h"
#include "Reader/Reader.h"
#include "Trivent/Mapping.h"
#include "TGraph.h"
using namespace marlin;
#define degtorad 0.0174532925
unsigned int EventsNoise=0;
unsigned int eventtotal=0;
unsigned int EventsSelected=0;
unsigned int TouchedEvents=0;
unsigned int _eventNr=0;
#define size_pad 10.4125
#define size_strip 2.5
std::map<int,bool>Warning;
/*class ToTree
{
public:
int pI,pJ,pK,pAsic,pDifId,pAsicChannel;
unsigned int pTime;
double pX,pY,pZ;
};*/

/*ToTree totree;
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
TBranch* Branch10 =  t->Branch("AsicChannel",&(totree.pAsicChannel));*/
std::vector<TH1F*>Time_Distr;
std::vector<TH1F*>Hits_Distr;
std::vector<TH1F*>Time_Distr_Events;
std::vector<TH1F*>Hits_Distr_Events;
std::vector<TH1F*>Time_Distr_Noise;
std::vector<TH1F*>Hits_Distr_Noise;
std::vector<TH2F*>Flux_Noise;
std::vector<TH2F*>Flux_Events;
std::vector<TH2F*>Flux_Noise_Asic;
std::vector<TH2F*>Flux_Events_Asic;
std::vector<long>Nbrof0Hits;
int _NbrRun=0;
long long int total_time=0;
double timemax=0;
double timemin=2147483647;


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
    for(unsigned int j=0; j<Time_Distr.size(); ++j) Times_Plates.emplace_back(std::map<int,int>());
    for(std::vector<RawCalorimeterHit *>::iterator it=vec.begin(); it!=vec.end(); ++it) {


        CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
        int dif_id  = (*it)->getCellID0() & 0xFF ;
        Times_Plates[geom.GetDifNbrPlate(dif_id)-1][(*it)->getTimeStamp()]++;

        if(IsNoise==1) {

            Time_Distr_Noise[geom.GetDifNbrPlate(dif_id)-1]->Fill((*it)->getTimeStamp(),1);
        } else {
            Time_Distr_Events[geom.GetDifNbrPlate(dif_id)-1]->Fill((*it)->getTimeStamp(),1);
        }
        int asic_id = ((*it)->getCellID0() & 0xFF00)>>8;
        int chan_id = ((*it)->getCellID0() & 0x3F0000)>>16;
        float ca=cos(geom.GetDifAlpha(dif_id)*degtorad);
        float sa=sin(geom.GetDifAlpha(dif_id)*degtorad);
        float cb=cos(geom.GetDifBeta(dif_id)*degtorad);
        float sb=sin(geom.GetDifBeta(dif_id)*degtorad);
        float cg=cos(geom.GetDifGamma(dif_id)*degtorad);
        float sg=sin(geom.GetDifGamma(dif_id)*degtorad);

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
        if(IsNoise==1) {
            Flux_Noise[geom.GetDifNbrPlate(dif_id)-1]->Fill(I,J);
            if(geom.GetDifType(dif_id)==positional)Flux_Noise_Asic[geom.GetDifNbrPlate(dif_id)-1]->Fill(asic_id,asic_id);
            else Flux_Noise_Asic[geom.GetDifNbrPlate(dif_id)-1]->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
        } else {
            Flux_Events[geom.GetDifNbrPlate(dif_id)-1]->Fill(I,J);
            if(geom.GetDifType(dif_id)==positional)Flux_Events_Asic[geom.GetDifNbrPlate(dif_id)-1]->Fill(asic_id,asic_id);
            else Flux_Events_Asic[geom.GetDifNbrPlate(dif_id)-1]->Fill((AsicShiftI[asic_id]+geom.GetDifPositionX(dif_id))/8,(32-AsicShiftJ[asic_id]+geom.GetDifPositionY(dif_id))/8);
        }
        cd["I"] = I ;
        cd["J"] = J ;
        cd["K"] = K ;
        caloHit->setPosition(pos);
        cd.setCellID( caloHit ) ;
        col->addElement(caloHit);
        /*totree.pI=I;
        totree.pJ=J;
        totree.pK=K;
        totree.pX=pos[0];
        totree.pY=pos[1];
        totree.pZ=pos[2];
        totree.pAsic=asic_id;
        totree.pDifId=dif_id;
        totree.pAsicChannel=chan_id;
        totree.pTime=(*it)->getTimeStamp();*/

        //std::cout<<magenta<<totree.pI<<"  "<<totree.pJ<<"  "<<red<<(*it)->getTimeStamp()<<"  "<<totree.pTime<<normal<<std::endl;
        /*t->Fill();*/

    }
    if(IsNoise==1) {
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) Hits_Distr_Noise[i]->Fill(it->second,1);
    } else {
        for(unsigned int i=0; i<Times_Plates.size(); ++i)for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) Hits_Distr_Events[i]->Fill(it->second,1);
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
    _noiseFileName="0";
    registerProcessorParameter("NOISEOutputFile" ,"NOISE file" ,_noiseFileName ,_noiseFileName);
    _timeWin = 2;
    registerProcessorParameter("timeWin" ,"time window = 2 in default",_timeWin ,_timeWin);
    _noiseCut = 7;
    registerProcessorParameter("noiseCut" ,"noise cut in time spectrum 7 in default",_noiseCut ,_noiseCut);
    _LayerCut = 3;
    registerProcessorParameter("LayerCut" ,"cut in number of layer 3 in default",_LayerCut ,_LayerCut);
}

TriventProcessor::~TriventProcessor() {}


void TriventProcessor::processRunHeader( LCRunHeader* run)
{
    LCTOOLS::dumpRunHeader(run);
}



void TriventProcessor::init()
{
    printParameters();
    _EventWriter = LCFactory::getInstance()->createLCWriter() ;
    _EventWriter->setCompressionLevel( 0 ) ;
    _EventWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ;
    if(_noiseFileName!="0") {
        _NoiseWriter = LCFactory::getInstance()->createLCWriter() ;
        _NoiseWriter->setCompressionLevel( 0 ) ;
        _NoiseWriter->open(_noiseFileName.c_str(),LCIO::WRITE_NEW) ;
    }
    ReaderFactory readerFactory;
    Reader* myReader = readerFactory.CreateReader(_ReaderType);
    if(myReader) {
        myReader->Read(_FileNameGeometry,geom);
        geom.PrintGeom();
        std::map<int, Dif > Difs=geom.GetDifs();
        std::map<int,int> PlansType;
        for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it) {
            if(geom.GetDifType(it->first)!=temporal) {
                PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
            }
        }
        for(std::map<int, int >::iterator it=PlansType.begin(); it!=PlansType.end(); ++it) {
            std::string b="Number_hits_Events"+ std::to_string( (long long int) it->first +1 );
            Times_Plates.emplace_back(std::map<int ,int>());
            Times_Plates_perRun.emplace_back(std::map<int ,int>());

            if(it->second==positional) {
                Flux_Events.emplace_back(new TH2F(b.c_str(),b.c_str(),128,0,128,1,0,50));
            } else {
                Flux_Events.emplace_back(new TH2F(b.c_str(),b.c_str(),100,0,100,100,0,100));
            }
            std::string c="Number_hits_Noise"+ std::to_string( (long long int) it->first +1 );
            if(it->second==positional) {
                Flux_Noise.emplace_back(new TH2F(c.c_str(),c.c_str(),128,0,128,1,0,50));
            } else {
                Flux_Noise.emplace_back(new TH2F(c.c_str(),c.c_str(),100,0,100,100,0,100));
            }
            std::string d="Number_hits_Events_Asic"+ std::to_string( (long long int) it->first +1 );
            if(it->second==positional) {
                Flux_Events_Asic.emplace_back(new TH2F(d.c_str(),d.c_str(),2,0,3,1,0,50));
            } else {
                Flux_Events_Asic.emplace_back(new TH2F(d.c_str(),d.c_str(),100,0,100,100,0,100));
            }
            std::string e="Number_hits_Noise_Asic"+ std::to_string( (long long int) it->first +1 );
            if(it->second==positional) {
                Flux_Noise_Asic.emplace_back(new TH2F(e.c_str(),e.c_str(),2,1,3,1,0,50));
            } else {
                Flux_Noise_Asic.emplace_back(new TH2F(e.c_str(),e.c_str(),100,0,100,100,0,100));
            }
            std::string h="Time_Distribution"+ std::to_string( (long long int) it->first +1 );
            Time_Distr.emplace_back(new TH1F(h.c_str(),h.c_str(),25000,0,25000));
            std::string i="Number_hit_per_clock"+ std::to_string( (long long int) it->first +1 );
            Hits_Distr.emplace_back(new TH1F(i.c_str(),i.c_str(),25000,0,25000));
            std::string j="Time_Distribution_Events"+ std::to_string( (long long int) it->first +1 );
            Time_Distr_Events.emplace_back(new TH1F(j.c_str(),j.c_str(),25000,0,25000));
            std::string k="Number_hit_per_clock_Events"+ std::to_string( (long long int) it->first +1 );
            Hits_Distr_Events.emplace_back(new TH1F(k.c_str(),k.c_str(),25000,0,25000));
            std::string l="Time_Distribution_Noise"+ std::to_string( (long long int) it->first +1 );
            Time_Distr_Noise.emplace_back(new TH1F(l.c_str(),l.c_str(),25000,0,25000));
            std::string m="Number_hit_per_clock_Noise"+ std::to_string( (long long int) it->first +1 );
            Hits_Distr_Noise.emplace_back(new TH1F(m.c_str(),m.c_str(),25000,0,25000));
            Nbrof0Hits.push_back(0);
        }

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

            RawTimeDifs.clear();
            LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
            LCCollection* col2 = evtP ->getCollection("DHCALRawTimes");
            if(col2!=NULL) {
                for (int ihit=0; ihit < col2->getNumberOfElements(); ++ihit) {
                    EVENT::CalorimeterHit* raw_time = dynamic_cast<EVENT::CalorimeterHit* >( col2->getElementAt(ihit)) ;
                    std::cout<<raw_time->getTime()<<"  "<<raw_time->getEnergyError()<<std::endl;
                    //RawTimeDifs[raw_time->getTimeStamp()].push_back(raw_time);
                }
            }

            if(col == NULL) {
                std::cout<< "TRIGGER SKIPED ..."<<std::endl;
                _trig_count++;
                break;
            }
            int numElements = col->getNumberOfElements();
            long timemax_local=0;
            long timemin_local=17976931348620;
            for(unsigned int i=0; i<Times_Plates_perRun.size(); ++i)Times_Plates_perRun[i].clear();
            for (int ihit=0; ihit < numElements; ++ihit) {

                RawCalorimeterHit *raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit)) ;
                if (raw_hit != NULL) {
                    int dif_id  = (raw_hit)->getCellID0() & 0xFF ;
                    if(geom.GetDifNbrPlate(dif_id)==-1) {
                        if(Warning[dif_id]!=true) {
                            Warning[dif_id]=true;
                            std::cout<<"Please add DIF "<<dif_id<<" to your geometry file; I'm Skipping its data."<<std::endl;
                        }
                        continue;
                    }
                    //std::cout<<red<<dif_id<<blue<<geom.GetDifNbrPlate(dif_id)-1<<"  "<<dif_id<<normal<<std::endl;
                    Times[raw_hit->getTimeStamp()]++;
                    Times_Plates[geom.GetDifNbrPlate(dif_id)-1][raw_hit->getTimeStamp()]++;
                    Times_Plates_perRun[geom.GetDifNbrPlate(dif_id)-1][raw_hit->getTimeStamp()]++;
                    if(raw_hit->getTimeStamp()>timemax)timemax=raw_hit->getTimeStamp();
                    if(raw_hit->getTimeStamp()>timemax_local)timemax_local=raw_hit->getTimeStamp();
                    if(raw_hit->getTimeStamp()<timemin_local)timemin_local=raw_hit->getTimeStamp();
                    RawHits[raw_hit->getTimeStamp()].push_back(raw_hit);
                }
            }

            double delta=timemax_local-timemin_local;
            total_time+=delta;
            for(unsigned int j =0; j<Nbrof0Hits.size(); ++j) {
                Nbrof0Hits[j]+=(delta-Times_Plates_perRun[j].size());
            }
            //std::cout<<red<<timemax_local<<"  "<<timemin_local<<"  "<<delta<<"  "<<green<<Times_Plates_perRun[0].size()<<yellow<<"  "<<delta-Times_Plates_perRun[0].size()<<red<<"  "<<Nbrof0Hits[0]<<normal<<std::endl;

            FillTimes();

            for(std::map< int,int>::iterator itt=Times.begin(); itt!=Times.end(); ++itt) {
                EventsGrouped.clear();
                std::map<int,std::vector<RawCalorimeterHit *> >::iterator middle=RawHits.find(itt->first);
                std::map<int,std::vector<RawCalorimeterHit *> >::iterator after=middle;
                std::map<int,std::vector<RawCalorimeterHit *> >::iterator before=middle;
                while(fabs(middle->first-before->first)<=_timeWin && before!=RawHits.begin()) --before;
                ++before;
                while(fabs(after->first-middle->first)<=_timeWin && after!=RawHits.end()) ++after;
                std::map<int,int> nbrPlanestouched;
                for(middle=before; middle!=after; ++middle ) {

                    for(int unsigned i=0; i<(middle->second).size(); ++i) {
                        int dif_id=((middle->second)[i])->getCellID0() & 0xFF;
                        nbrPlanestouched[geom.GetDifNbrPlate(dif_id)]++;
                    }

                }

                if(nbrPlanestouched.size()>=(unsigned int)(_LayerCut)) {
                    EventsSelected++;
                    for(middle=before; middle!=after; ) {
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
            }

            if(_noiseFileName!="0") {
                EventsNoise++;
                LCEventImpl*  evt = new LCEventImpl() ;
                LCCollectionVec* col_event = new LCCollectionVec(LCIO::CALORIMETERHIT);
                col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
                col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));

                CellIDEncoder<CalorimeterHitImpl> cd( "I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" ,col_event) ;
                for(std::map<int,std::vector<RawCalorimeterHit *> >::iterator itt=RawHits.begin(); itt!=RawHits.end(); ++itt) {
                    FillIJK((itt->second),col_event,cd,1);
                }
                evt->addCollection(col_event, "SDHCAL_HIT");
                evt->setEventNumber(EventsNoise);
                evt->setRunNumber(evtP->getRunNumber());
                //std::cout<<(timemax-timemin)*200e-9<<std::endl;
                evt->setTimeStamp(timemax-timemin);
                _NoiseWriter->writeEvent( evt ) ;
                delete evt;
            }
        }
    }
}

void TriventProcessor::end()
{
    std::string name="Results_Trivent_"+ std::to_string( (long long int) _NbrRun)+".root";
    TFile *hfile = new TFile(name.c_str(),"RECREATE","Results");
    for(unsigned int i=0; i<Flux_Noise.size(); ++i) {
        //Flux_Hits[i]->Scale(total_time*200e-9);
        //Flux_Hits[i]->Write();
        //Flux_Noise[i]->Scale(total_time*200e-9);
        Flux_Noise[i]->Write();
        Flux_Noise_Asic[i]->Write();
        //Flux_Events[i]->Scale(total_time*200e-9);
        Flux_Events[i]->Write();
        Flux_Events_Asic[i]->Write();
        //Noise3D->Scale(total_time*200e-9);
        //Events3D->Scale(total_time*200e-9);
        //Noise3D->Write();
        //Events3D->Write();
        for(std::map<int,int>::iterator it = Times_Plates[i].begin(); it!=Times_Plates[i].end(); ++it) {
            Time_Distr[i]->Fill(it->first,it->second);
            Hits_Distr[i]->Fill(it->second,1);
        }
        Hits_Distr[i]->Fill(0.0,Nbrof0Hits[i]);
        Hits_Distr_Noise[i]->Fill(0.0,Nbrof0Hits[i]);
        Time_Distr[i]->GetXaxis()->SetRange(0,Times_Plates[i].size()+10);
        Hits_Distr[i]->GetXaxis()->SetRange(0,150);
        Time_Distr[i]->Write();
        Hits_Distr[i]->Write();
        Time_Distr_Events[i]->Write();
        Time_Distr_Noise[i]->Write();
        Hits_Distr_Noise[i]->Write();
        Hits_Distr_Events[i]->Write();
    }
    for(unsigned int i=0; i<Flux_Noise.size(); ++i) {
        delete Flux_Noise[i];
        delete Flux_Events[i];
        delete Flux_Noise_Asic[i];
        delete Flux_Events_Asic[i];
        delete Time_Distr[i];
        delete Hits_Distr[i];
        delete Time_Distr_Events[i];
        delete Time_Distr_Noise[i];
        delete Hits_Distr_Events[i];
        delete Hits_Distr_Noise[i];
    }
    /*t->Write();
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
    delete t;*/
    hfile->Close();
    delete hfile;
    _EventWriter->close();
    if(_noiseFileName!="0") _NoiseWriter->close();
    std::cout << "TriventProcess::end() !! "<<_trig_count<<" Events Trigged"<< std::endl;
    std::cout <<TouchedEvents<<" Events were overlaping "<<"("<<(TouchedEvents*1.0/(TouchedEvents+eventtotal))*100<<"%)"<<std::endl;
    std::cout <<"Total nbr Events : "<<eventtotal<<" Events with nbr of plates >="<<_LayerCut<<" : "<<EventsSelected<<" ("<<EventsSelected*1.0/eventtotal*100<<"%)"<< std::endl;
    std::cout <<"Total Time : "<<total_time <<std::endl;
    for(std::map<int,bool>::iterator it=Warning.begin(); it!=Warning.end(); it++) std::cout<<red<<"REMINDER::Data from Dif "<<it->first<<" are skipped !"<<normal<<std::endl;
}
