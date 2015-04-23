#include "Analysis/Analysis.h"
#include <iostream>
#include <string>
#include <iomanip>
#include "marlin/Processor.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/CellIDDecoder.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "Reader/ReaderFactory.h"
#include "Reader/Reader.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "UTIL/CellIDDecoder.h"
#include "Colors.h"
#define degtorad 0.0174532925
#include <cstdlib>
#include <cmath>
#include "TH2F.h"
#include "TTree.h"
#ifndef COLORS_H
#define normal " "
#define red "  "
#define vert " "
#define blanc " "
#define green " "
#define blue " "
#endif
#include <fstream>
#define size_pad 10.4125
#define size_strip 2.5
int totalTrace=0;
int eventnbrr=0;
std::ofstream Verif( "Verif.txt", std::ios_base::out ); 

class ToTreee
{
public:
double ChiXZ, ChiYZ, CDXZ, CDYZ, OrdXZ, OrdYZ;
};
ToTreee totreee;
std::string name="Treee";
TTree* tt= new TTree(name.c_str(), name.c_str());
TBranch* Branch1 =  tt->Branch("ChiXZ",&(totreee.ChiXZ));
TBranch* Branch2 =  tt->Branch("ChiYZ",&(totreee.ChiYZ));
TBranch* Branch3 =  tt->Branch("CDXZ",&(totreee.CDXZ));
TBranch* Branch4 =  tt->Branch("CDYZ",&(totreee.CDYZ));
TBranch* Branch5 =  tt->Branch("OrdXZ",&(totreee.OrdXZ));
TBranch* Branch6 =  tt->Branch("OrdYZ",&(totreee.OrdYZ));


using namespace std;
std::vector<TH2F*>Distribution_hits;
std::vector<TH2F*>Distribution_exp;
//std::vector<TH2F*>Correlations;
unsigned NbrReadOut=0;
void AnalysisProcessor::processRunHeader( LCRunHeader* run)
{
}
unsigned NbrRunn=0;
void AnalysisProcessor::PrintStatShort()
{
    ofstream fichier;
    fichier.open("ResultsShorts"+std::to_string((long long int)_NbrRun)+".txt", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
    if(fichier) { // si l'ouverture a réussi
        fichier<<_NbrRun<<"   ";
        for(unsigned int i=0; i!=testedPlanList.size(); ++i) {
            fichier<<testedPlanList[i].efficiencyShort()<<" "<<sqrt(testedPlanList[i].GetNumberOKShort()*testedPlanList[i].efficiencyShort()*(1-testedPlanList[i].efficiencyShort()))*1.0/testedPlanList[i].GetNumberOKShort()<<" "<<testedPlanList[i].multiplicityShort()<<" 0 "<<"  ";
        }
        fichier<<std::endl;
        fichier.close();  // on referme le fichier
    }

    for(unsigned int i=0; i!=testedPlanList.size(); ++i) {
        // testedPlanList[i].print();
        std::cout<<red<<_NbrRun<<normal<<std::endl;
        std::cout<<red<<"Efficacite "<<testedPlanList[i].efficiencyShort()<<" Multiplicite "<<testedPlanList[i].multiplicityShort()<<normal<<std::endl;
    }
}



void FillDelimiter(std::string ToParse,int size)
{
    std::string delimiter_Dif = "|";
    std::string delimiter_Difs = "*";
    std::string delimiter_others=":";
    size_t pos = 0;
    std::string token;
    std::vector<std::string>a;
    bool findstar =ToParse.find(delimiter_Difs);
    if(findstar) {
        pos=ToParse.find(delimiter_Dif);
        a.push_back(ToParse.substr(0, pos));
        std::vector<double> tab(4);
        int Dif=0;
        int j =0;
        size_t posi =0;
        while ((posi = a[0].find(delimiter_others)) != std::string::npos) {
            std::string token;
            token = a[0].substr(0, posi);
            tab[j]=atof(token.c_str());/*std::cout <<green<< token <<normal<< std::endl;*/



            a[0].erase(0, posi + delimiter_others.length());
            ++j;
        }
        tab[3]=atof(a[0].c_str());
        //std::cout<<green<<tab[3]<<normal<<std::endl;
        for(unsigned int i=0; i<size; ++i) Delimiter[i+1]=tab;
    } else {
        while ((pos = ToParse.find(delimiter_Dif)) != std::string::npos) {
            token = ToParse.substr(0, pos);
            //std::cout << token << std::endl;
            a.push_back(token);
            ToParse.erase(0, pos + delimiter_Dif.length());
        }
        for(unsigned int i=0; i<a.size(); ++i) {
            std::vector<double> tab(4);
            int Dif=0;
            int j =0;
            size_t posi =0;
            while ((posi = a[i].find(delimiter_others)) != std::string::npos) {
                token = a[i].substr(0, posi);

                if (j==0) {
                    Dif=atof(token.c_str());/*std::cout <<blue<< token <<normal<< std::endl;*/
                } else {
                    tab[j-1]=atof(token.c_str());/*std::cout <<green<< token <<normal<< std::endl;*/
                }
                a[i].erase(0, posi + delimiter_others.length());
                ++j;
            }
            tab[3]=atof(a[i].c_str());
            //std::cout<<green<<tab[3]<<normal<<std::endl;
            Delimiter[Dif]=tab;
        }
    }
    if(a.size()==0)std::cout<<red<<"Warning:No Delimiters given "<<normal<<std::endl;
    if(a.size()!=size&&a.size()!=0&&!findstar) {
        std::cout<<red<<"Error:Delimiters no well set ! "<<normal<<std::endl;
        std::exit(2);
    }

    //std::cout<<red<<a.size()<<normal<<std::endl;

    std::vector<string>word {"Imin : "," Imax : "," Jmin : "," Jmax : "};
    std::cout<<green<<"Delimiters"<<normal<<std::endl ;
    for(std::map<int,std::vector<double>>::iterator it=Delimiter.begin(); it!=Delimiter.end(); ++it) {

        std::cout<<green<<"Plane "<<it->first<<" : "<<normal;
        for(int i=0; i<it->second.size(); ++i) {
            std::cout<<green<<word[i]<<(it->second)[i]<<normal;
        }
        std::cout<<normal<<std::endl;
    }

}

int plan::countHitAt(double& x, double& y, double dlim)
{
    CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
    int n=0;
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) {
        if(fabs(x-(*it)->getPosition()[0])<dlim&&fabs(y-(*it)->getPosition()[1])<dlim) {
            n++;
            Distribution_hits[cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
        }
    }
    return n;
}

int plan::countHitAtStrip(double& x, double dlim)
{
    int n=0;
    CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) {
        if(fabs(x-(*it)->getPosition()[0])<dlim) {
            n++;
            //std::cout<<fabs(x-(*it)->getPosition()[0])<<"  "<<dlim<<std::endl;
            Distribution_hits[cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
        }
    }
    return n;
}

void plan::computeBarycentre( )
{
    for (int i=0; i<3; i++) barycentre[i]=0;
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it) {
        for (int i=0; i<3; i++) {
            barycentre[i]+=(*it)->getPosition()[i];
            // std::cout<<green<<(*it)->getPosition()[i]<<normal<<"  ";
        }
    }
    //std::cout<<std::endl;
    if (nHits() != 0)
        for (int i=0; i<3; i++) barycentre[i]/=nHits();
    //std::cout<<red<<"Barycentre:"<<    barycentre[0]<<"  "<<barycentre[1]<<"  "<<barycentre[2]<<normal<<std::endl;
}

void plan::computeMaxima()
{
    for (int i=0; i<3; i++) {
        min[i]=10000000;
        max[i]=-10000000;
    }
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it) {
        for (int i=0; i<3; i++) {
            if((*it)->getPosition()[i]<min[i])min[i]=(*it)->getPosition()[i];
            if((*it)->getPosition()[i]>max[i])max[i]=(*it)->getPosition()[i];
        }
    }
}
void plan::GivePoint()
{
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it) 
    {
        for (int i=0; i<2; i++) 
        {
           Verif<<(*it)->getPosition()[i]<<"  ";
        }
    }
}

void testedPlan::testYou(std::map<int,plan>& mapDIFplan)
{
    counts[TESTYOUCALLED]++;
    std::vector<plan*> plansUsedForTrackMaking;
    plan* thisPlan=nullptr;
    for (std::map<int,plan>::iterator it=mapDIFplan.begin(); it!=mapDIFplan.end(); ++it) {
        if (Nbr!=it->first) plansUsedForTrackMaking.push_back(&(it->second));
        else thisPlan=&(it->second);
    }

    for (std::vector<plan*>::iterator it=plansUsedForTrackMaking.begin(); it != plansUsedForTrackMaking.end(); ++it) if ((*it)->nHits()>=_NbrHitPerPlaneMax ) return;
    //std::cout<<blue<<plansUsedForTrackMaking.size()<<std::endl;
    if(plansUsedForTrackMaking.size()<_NbrPlaneUseForTracking) return;
    counts[NOTOOMUCHHITSINPLAN]++;
    ////////////////////////////////////////////////////////////////////////////////////
    TGraphErrors grxz(plansUsedForTrackMaking.size());
    TGraphErrors gryz(plansUsedForTrackMaking.size());
    for (unsigned int i=0; i < plansUsedForTrackMaking.size(); ++i) {
        plan &p=*(plansUsedForTrackMaking[i]);
        p.computeBarycentre();
        p.computeMaxima();
        grxz.SetPoint(i,p.barycentreZ(),p.barycentreX());
        if(p.GetType()==pad) {
            gryz.SetPoint(i,p.barycentreZ(),p.barycentreY());
            gryz.SetPointError(i,p.ErrorZ(),p.ErrorY());
        }
        grxz.SetPointError(i,p.ErrorZ(),p.ErrorX());
    }
    grxz.Fit("pol1","QRO","",-50000.0,50000.0);
    TF1 *myfitxz = (TF1*) grxz.GetFunction("pol1");
    double  kxz = myfitxz->GetChisquare();
    if (kxz>= _Chi2) return;
    double pxz0 = myfitxz->GetParameter(0);
    double  pxz1 = myfitxz->GetParameter(1);
    counts[XZTRACKFITPASSED]++;
    gryz.Fit("pol1","QRO","",-50000.0,50000.0);
    TF1 *myfityz = (TF1*) gryz.GetFunction("pol1");
    double  kyz = myfityz->GetChisquare();
    if(this->GetType()==positional) kyz=0;
    if (kyz>= _Chi2) return;
    double pyz0 = myfityz->GetParameter(0);
    double  pyz1 = myfityz->GetParameter(1);
    counts[YZTRACKFITPASSED]++;
    double Zexp=this->GetZexp(pxz0,pyz0,pxz1,pyz1);
    //std::cout<<blue<<pxz0<<"  "<<myfit.GetParErrors()[0]<<"  "<<pxz1<<"  "<<myfit.GetParErrors()[1]<<"  "<<pyz0<<"  "<<myfit2.GetParErrors()[0]<<"  "<<pyz1<<"  "<<myfit2.GetParErrors()[1]<<normal<<std::endl;
    
    ///////////////////////////////
    //double Xexp = pxz0+pxz1*Zexp;
    //double Yexp = pyz0+pyz1*Zexp;
    ///////////////////////////////
    double Projectioni=GetProjectioni(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
    double Projectionj=GetProjectionj(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
    
    //std::cout<<red<<this->GetIp()<<"  "<<this->GetIm()<<"  "<<this->GetJp()<<"  "<<this->GetJm()<<normal<<std::endl;
    bool Pass;
    if(Delimiter.size()==0)Pass=1;
    else Pass=Projectioni<=this->GetIp()&&Projectioni>=this->GetIm()&&Projectionj<=this->GetJp()&&Projectionj>=this->GetJm();
    if(Pass) {
        nombreTests++;
        
        if (nullptr==thisPlan) return;
        Distribution_exp[this->NbrPlate()]->Fill(Projectioni,Projectionj);
        counts[NOHITINPLAN]++;
        int nhit;
        if(this->GetType()==pad) {
            nhit=thisPlan->countHitAt(Projectioni,Projectionj,/*6*10.4125*/_dlimforPad);
        } else {
            nhit=thisPlan->countHitAtStrip(Projectioni,_dlimforStrip);
        }
        if (nhit>0) nombreTestsOK++;
        sommeNombreHits+=nhit;
        totalTrace++;
	Verif<<NbrRunn<<"  "<<eventnbrr<<"  "<<totalTrace<<"  "<<kxz<<"  "<<kyz<<"  "<<pxz0<<"  "<<pyz0<<"  "<<pxz1<<"  "<<pyz1<<"  "<<plansUsedForTrackMaking.size()<<"  ";
        totreee.ChiXZ=kxz;
        totreee.ChiYZ=kyz;
        totreee.CDXZ=pxz1;
        totreee.CDYZ=pyz1;
        totreee.OrdXZ=pxz0;
        totreee.OrdYZ=pyz0;
        tt->Fill();
	for (unsigned int i=0; i < plansUsedForTrackMaking.size(); ++i){plan &p=*(plansUsedForTrackMaking[i]);p.computeBarycentre();Verif<<p.barycentreX()/10<<"  ";}
        for (unsigned int i=0; i < plansUsedForTrackMaking.size(); ++i){plan &p=*(plansUsedForTrackMaking[i]);p.computeBarycentre();Verif<<p.barycentreY()/10<<"  ";}
	Verif<<std::endl;
        
    }
    delete myfityz;
    delete myfitxz;
    ///////////////////////////////////////////////////////////////////////////////////////////
}

void testedPlan::print()
{
    std::cout<<blue<<"Plane Number (in geometry file): "<<Nbr+1<<" Z = "<<Z0<<" : "<<normal<<std::endl;
    std::cout<<blue<<"Number of Test : "<<counts[0]<<"; with >="<<_NbrPlaneUseForTracking<<" planes for tracking : "<<counts[1]<<"; with ChiXZ <"<<_Chi2<<" : "<<counts[2]<<"; with ChiYZ <"<<_Chi2<<" : "<<counts[3]<<" ; with track in the Delimiters "<<nombreTests<<"; with hits in it : "<<counts[4]<<" ; with hits in dlim : "<<nombreTestsOK<<normal<<std::endl;
    std::cout<<blue<<"Sum of hits in dlim : "<<sommeNombreHits<<normal<<std::endl;
//<<" NombreTests = "<<nombreTests<<" nombreTestsOk = "<<nombreTestsOK<<"  sommeNHits = "<<sommeNombreHits<< "  ( type = " <<GetType()<<" )"<<normal<<std::endl;
  //  std::vector<std::string>Text{"TESTYOUCALLED","NOTOOMUCHHITSINPLAN","XZTRACKFITPASSED","YZTRACKFITPASSED","PRESENCEOFHITINPLAN"};
    //for (int i=0; i<NCOUNTERS; i++) std::cout <<blue<< Text[i]<<" : "<<counts[i]<<"  "<<normal;
    std::cout<<std::endl;
}

void AnalysisProcessor::PrintStat()
{
    ofstream fichier;
    fichier.open("Results.txt", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
    if(fichier) { // si l'ouverture a réussi
        fichier<<_NbrRun<<"   ";
        for(unsigned int i=0; i!=testedPlanList.size(); ++i) {
            fichier<<testedPlanList[i].efficiency()<<" "<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" "<<testedPlanList[i].multiplicity()<<" 0 "<<"  ";
        }
        fichier<<std::endl;
        fichier.close();  // on referme le fichier
    }
    std::cout<<"Run Number : "<<_NbrRun<<std::endl;
    for(unsigned int i=0; i!=testedPlanList.size(); ++i) {
        // testedPlanList[i].print();
        
        std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanList[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanList[i].efficiency()<<" Error : "<<setw(6)<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" Multiplicity : "<<setw(6)<<testedPlanList[i].multiplicity()<<normal<<std::endl;
    }
}

using namespace marlin;
AnalysisProcessor aAnalysisProcessor;
AnalysisProcessor::AnalysisProcessor() : Processor("AnalysisProcessorType")
{
    std::vector<std::string> hcalCollections(1,"SDHCAL_HIT");
    registerInputCollections( LCIO::RAWCALORIMETERHIT,"HCALCollections","HCAL Collection Names",_hcalCollections,hcalCollections);
    _FileNameGeometry="";
    registerProcessorParameter("FileNameGeometry","Name of the Geometry File",_FileNameGeometry,_FileNameGeometry);
    _ReaderType="";
    registerProcessorParameter("ReaderType","Type of the Reader needed to read InFileName",_ReaderType ,_ReaderType);
    _Chi2 = 1.0;
    registerProcessorParameter("Chi2" ,"Value of the Chi2  ",_Chi2 ,_Chi2);
    _NbrHitPerPlaneMax = 6;
    registerProcessorParameter("NbrHitPerPlaneMax" ,"Maximal number of Hit in each Plane (<=6 by default)  ",_NbrHitPerPlaneMax ,_NbrHitPerPlaneMax);
    _NbrPlaneUseForTracking = 3;
    registerProcessorParameter("NbrPlaneUseForTracking" ,"Number minimal of PLanes used for tracking (>=3 by default)",_NbrPlaneUseForTracking,_NbrPlaneUseForTracking);
    _dlimforPad=40;
    registerProcessorParameter("dlimforPad" ,"dlim for Pad ",_dlimforPad,_dlimforPad);
    _dlimforStrip=40;
    registerProcessorParameter("dlimforStrip" ,"dlim for Strip ",_dlimforStrip,_dlimforStrip);
    _Delimiters="";
    registerProcessorParameter("Delimiters" ,"Delimiters",_Delimiters,_Delimiters);
    _ShortEfficiency=0;
    registerProcessorParameter("ShortEfficiency" ,"ShortEfficiency",_ShortEfficiency,_ShortEfficiency);
}

AnalysisProcessor::~AnalysisProcessor() {}

void AnalysisProcessor::init()
{
    printParameters();
    Verif<<"Run Event Num trace ChiXZ ChiYZ CDXZ CDYZ OrdXZ OrdYZ LayTouch PosX1 PosX2 PosX3 PosX4 PosX5 PosX6 PosX7 PosX8 PosY1 PosY2 PosY3 PosY4PosY5 PosY6 PosY7 PosY8"<<std::endl;
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
                //PlansType[geom.GetDifNbrPlate(it->first)-1]=geom.GetDifType(it->first);
            }
        }
        FillDelimiter(_Delimiters,PlansType.size());
        for(std::map<int, int >::iterator it=PlansType.begin(); it!=PlansType.end(); ++it) {
            if(Delimiter.size()!=0)
                testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second,Delimiter[it->first+1][1],Delimiter[it->first+1][0],Delimiter[it->first+1][3],Delimiter[it->first+1][2]));
            else testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second,0,0,0,0));
            //std::string b="Correlations"+ std::to_string( (long long int) it->first +1 );
            std::string a="Distribution hit selectionner par analysis"+ std::to_string( (long long int) it->first +1 );
            if(it->second==positional) {
                Distribution_hits.push_back(new TH2F(a.c_str(),a.c_str(),700,0,700,700,0,700));
            } else {
                Distribution_hits.push_back(new TH2F(a.c_str(),a.c_str(),700,0,700,700,0,700));
            }
             std::string b="Hit expected"+ std::to_string( (long long int) it->first +1 );
            if(it->second==positional) {
                Distribution_exp.push_back(new TH2F(b.c_str(),b.c_str(),700,0,700,700,0,700));
            } else {
                Distribution_exp.push_back(new TH2F(b.c_str(),b.c_str(),700,0,700,700,0,700));
            }
            //Correlations.push_back(new TH2F(b.c_str(),b.c_str(),200,0,200,200,0,200));

        }



    } else {
        std::cout << "Reader type n'existe pas !!" << std::endl;
        std::exit(1);
    }
    delete myReader;
}

void AnalysisProcessor::processEvent( LCEvent * evtP )
{ 
    _eventNr=evtP->getEventNumber();
    eventnbrr=_eventNr;
    if(_eventNr %1000 ==0)std::cout<<"Event Number : "<<_eventNr<<std::endl;
    _NbrRun=evtP->getRunNumber();
    NbrRunn=_NbrRun;
    Plans.clear();
    if (evtP != nullptr) {
        _eventNr=evtP->getEventNumber();
        for(unsigned int i=0; i< _hcalCollections.size(); i++) {
            LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
            if(col == nullptr) {
                std::cout<< red << "TRIGGER SKIPED ..."<< normal <<std::endl;
                break;
            }
            CellIDDecoder<CalorimeterHit> cd(col);
            int numElements = col->getNumberOfElements();
            for (int ihit=0; ihit < numElements; ++ihit) {
                CalorimeterHit *raw_hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(ihit)) ;
                if (raw_hit != nullptr) {
                    int dif_id=cd(raw_hit)["Dif_id"];
                    int I=cd(raw_hit)["I"];
                    int J=cd(raw_hit)["J"];
                    Plans[geom.GetDifNbrPlate(dif_id)-1].addHit(raw_hit);
                    Plans[geom.GetDifNbrPlate(dif_id)-1].SetType(geom.GetDifType(dif_id));
                    /*for(int jhit=ihit+1; jhit<numElements; ++jhit) {
                        CalorimeterHit *raw_hit2 = dynamic_cast<CalorimeterHit*>( col->getElementAt(jhit)) ;
                        int dif_id2=cd(raw_hit2)["Dif_id"];
                        int I2=cd(raw_hit2)["I"];
                        int J2=cd(raw_hit2)["J"];
                        //if((geom.GetDifNbrPlate(dif_id)==1 && geom.GetDifNbrPlate(dif_id2)==4)&&(raw_hit->getTime()==raw_hit2->getTime())) Correlations[4]->Fill((I-1),I2);
                    }*/
                }

	     
            }
        }
        for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter) iter->testYou(Plans);
        if(_ShortEfficiency!=0 && NbrReadOut%_ShortEfficiency==0) 
        {
        	PrintStatShort();
            	for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
		{
                	testedPlanList[i].ClearShort();
            	}
    	}
}
}

void AnalysisProcessor::end()
{
    
    std::cout<<"List of counters : "<<std::endl;
    for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter) iter->print();
    std::string b="Results_Analysis_"+std::to_string( (long long int) _NbrRun)+".root";
    TFile *hfile = new TFile(b.c_str(),"RECREATE");
    tt->Write();
    for(unsigned int i=0; i<Distribution_hits.size(); ++i) {
        Distribution_hits[i]->Write();
        Distribution_exp[i]->Write();
        //Correlations[i]->Write();
    }
    for(unsigned int i=0; i<Distribution_hits.size(); ++i) {
        delete Distribution_hits[i];
        delete Distribution_exp[i];
        //delete Correlations[i];
    }
    hfile->Close();
    delete hfile;
    PrintStat();
}
