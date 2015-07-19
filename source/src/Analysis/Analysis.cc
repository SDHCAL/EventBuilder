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
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TGLTH3Composition.h"
#include "TGraph2D.h"
#include "Utilities.h"
#include <map>
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
std::map<std::string,std::vector<TH1F*>>Short_Efficiency;
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
std::vector<std::string> names={"NORMAL"};
///////////////////////////////////
//std::vector<TH2F*>Distribution_hits;
//std::vector<TH2F*>Efficiency_pads;
//std::vector<TH2F*>Multiplicity_pads;
//std::vector<TH2F*>Efficiency_asics;
//std::vector<TH1F*>HowLongFromExpectedX;
//std::vector<TH1F*>HowLongFromExpectedY;
//std::vector<TGraph2D*>Distribution_hits_tgraph;
//TGraph2D* Distribution_hits_tgraph=new TGraph2D(1);
//TGraph2D* Distribution_exp_tgraph= new TGraph2D(1);
//std::vector<TH2F*>Distribution_exp;
///////////////////////////////////////////////
std::map<std::string,std::vector<TH2F*>>Distribution_hits;
std::map<std::string,std::vector<TH2F*>>Efficiency_pads;
std::map<std::string,std::vector<TH2F*>>Multiplicity_pads;
std::map<std::string,std::vector<TH2F*>>Efficiency_asics;
std::map<std::string,std::vector<TH1F*>>HowLongFromExpectedX;
std::map<std::string,std::vector<TH1F*>>HowLongFromExpectedY;
//std::vector<TGraph2D*>Distribution_hits_tgraph;
std::map<std::string,TGraph2D*> Distribution_hits_tgraph={{"NORMAL",new TGraph2D(1)},{"SCINTILLATOR",new TGraph2D(1)}};
//Distribution_hits_tgraph.insert(std::pair<std::string,TGraph2D*>("NORMAL",new TGraph2D(1)));
//Distribution_hits_tgraph.insert("SCINTILLATOR",new TGraph2D(1));
std::map<std::string,TGraph2D*> Distribution_exp_tgraph={{"NORMAL",new TGraph2D(1)},{"SCINTILLATOR",new TGraph2D(1)}};
//Distribution_exp_tgraph.insert("NORMAL",new TGraph2D(1));
//Distribution_exp_tgraph.insert("SCINTILLATOR",new TGraph2D(1));
std::map<std::string,std::vector<TH2F*>>Distribution_exp={{"NORMAL",{}},{"SCINTILLATOR",{}}};
//std::vector<TGraph2D*>Distribution_exp_tgraph;
//std::vector<TH2F*>Correlations;
unsigned NbrRunn=0;
void AnalysisProcessor::PrintStatShort(bool IsScinti)
{
    ofstream fichier;
    std::string namee="ResultsShorts"+std::to_string((long long int)_NbrRun);
    if(IsScinti==true)namee+="Scinti";
    namee+=".txt";
    fichier.open(namee.c_str(), ios::out | ios::app);  //déclaration du flux et ouverture du fichier
    if(fichier) { // si l'ouverture a réussi
        fichier<<_NbrRun<<"   ";
        for(unsigned int i=0; i!=testedPlanList.size(); ++i) {
            if(IsScinti==true)Short_Efficiency["SCINTILLATOR"][i]->Fill(_eventNr,testedPlanList[i].efficiencyShort());
            else Short_Efficiency["NORMAL"][i]->Fill(_eventNr,testedPlanList[i].efficiencyShort());
            fichier<<testedPlanList[i].efficiencyShort()<<" "<<sqrt(testedPlanList[i].GetNumberOKShort()*testedPlanList[i].efficiencyShort()*(1-testedPlanList[i].efficiencyShort()))*1.0/testedPlanList[i].GetNumberOKShort()<<" "<<testedPlanList[i].multiplicityShort()<<" 0 "<<"  ";
        }
        fichier<<std::endl;
        fichier.close();  // on referme le fichier
    }
    if(IsScinti==true)
    {
      for(unsigned int i=0; i!=testedPlanListScinti.size(); ++i) 
      {
      // testedPlanList[i].print();
      std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanListScinti[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanListScinti[i].efficiency()<<" Error : "<<setw(6)<<sqrt(testedPlanListScinti[i].GetNumberOK()*testedPlanListScinti[i].efficiency()*(1-testedPlanListScinti[i].efficiency()))*1.0/testedPlanListScinti[i].GetNumberOK()<<" Multiplicity : "<<setw(6)<<testedPlanListScinti[i].multiplicity()<<normal<<std::endl;
      }
    }
    else
    {
      for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
      {
      // testedPlanList[i].print();
      std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanList[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanList[i].efficiency()<<" Error : "<<setw(6)<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" Multiplicity : "<<setw(6)<<testedPlanList[i].multiplicity()<<normal<<std::endl;
      }
    }
}

unsigned int Number_hits=-1;
int plan::countHitAt(double& x, double& y, double dlim,int Iexpected,int Jexpected,int Kexpected,double Imax,double Imin,double Jmax,double Jmin,bool IsScinti)
{     
    CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
    int n=0;
    std::vector<int>IJKexpected={Iexpected,Jexpected,Kexpected};
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) {
        //if((*it)->getPosition()[0]>=Imin&&(*it)->getPosition()[0]<=Imax&&(*it)->getPosition()[1]>=Jmin&&(*it)->getPosition()[1]<=Jmax&&fabs(x-(*it)->getPosition()[0])<dlim&&fabs(y-(*it)->getPosition()[1])<dlim) {
     if(fabs(x-(*it)->getPosition()[0])<dlim&&fabs(y-(*it)->getPosition()[1])<dlim) {
            n++;
            Number_hits++;
            if(IsScinti==true)
            {
              Distribution_hits["SCINTILLATOR"][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
              Distribution_hits_tgraph["SCINTILLATOR"]->SetPoint(Number_hits,(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
              HowLongFromExpectedX["SCINTILLATOR"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
              HowLongFromExpectedY["SCINTILLATOR"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[1]-y);
            }
            else
            {
              Distribution_hits["NORMAL"][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
              Distribution_hits_tgraph["NORMAL"]->SetPoint(Number_hits,(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
              HowLongFromExpectedX["NORMAL"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
              HowLongFromExpectedY["NORMAL"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[1]-y);
            }
        }
    }
    if(IsScinti==true) Efficiency_per_padScinti[IJKexpected].push_back(n);
    else Efficiency_per_pad[IJKexpected].push_back(n);
    
    
    return n;
}

int plan::countHitAtStrip(double& x, double dlim, bool IsScinti)
{
    int n=0;
    CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
    for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) {
        if(fabs(x-(*it)->getPosition()[0])<dlim) {
            n++;
            Number_hits++;
            if(IsScinti==true)
            {
              //std::cout<<fabs(x-(*it)->getPosition()[0])<<"  "<<dlim<<std::endl;
              Distribution_hits["SCINTILLATOR"][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
              //Distribution_hits[cd(*it)["K"]-1]->Fill((*it)->getPosition()[0],(*it)->getPosition()[1]);
              Distribution_hits_tgraph["SCINTILLATOR"]->SetPoint(Number_hits,(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
              HowLongFromExpectedX["SCINTILLATOR"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
              HowLongFromExpectedY["SCINTILLATOR"][cd(*it)["K"]-1]->Fill(0);
            }
            else
            {
              //std::cout<<fabs(x-(*it)->getPosition()[0])<<"  "<<dlim<<std::endl;
              Distribution_hits["NORMAL"][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
              //Distribution_hits[cd(*it)["K"]-1]->Fill((*it)->getPosition()[0],(*it)->getPosition()[1]);
              Distribution_hits_tgraph["NORMAL"]->SetPoint(Number_hits,(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
              HowLongFromExpectedX["NORMAL"][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
              HowLongFromExpectedY["NORMAL"][cd(*it)["K"]-1]->Fill(0);
            }
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

void testedPlan::testYou(std::map<int,plan>& mapDIFplan,bool IsScinti)
{
    counts[TESTYOUCALLED]++;
    std::vector<plan*> plansUsedForTrackMaking;
    plan* thisPlan=nullptr;
    for (std::map<int,plan>::iterator it=mapDIFplan.begin(); it!=mapDIFplan.end(); ++it) 
    {
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
    for (unsigned int i=0; i < plansUsedForTrackMaking.size(); ++i) 
    {
        plan &p=*(plansUsedForTrackMaking[i]);
        p.computeBarycentre();
        p.computeMaxima();
        grxz.SetPoint(i,p.barycentreZ(),p.barycentreX());
        if(p.GetType()==pad) 
        {
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
    if(Pass) 
    {
        nombreTests++;
        nombreTestsShort++;
        if (nullptr==thisPlan) return;
        int I,J,K;
        ca=this->get_ca();
        sa=this->get_sa();
        cb=this->get_cb();
        sb=this->get_sb();
        cg=this->get_cg();
        sg=this->get_sg();
        
        if(IsScinti==true) Distribution_exp_tgraph["SCINTILLATOR"]->SetPoint(nombreTests,Projectioni,Projectionj,Zexp);
        else Distribution_exp_tgraph["NORMAL"]->SetPoint(nombreTests,Projectioni,Projectionj,Zexp);
        
        
        counts[NOHITINPLAN]++;
        int nhit;
        if(this->GetType()==pad) 
        {
          I=cg*cb*1.0/size_pad*(Projectioni-this->get_X0())+sg*cb*1.0/size_pad*(Projectionj-this->get_Y0())+-sb*this->get_Z0();
        	J=(-sg*ca+cg*sb*sa)*1.0/size_pad*(Projectioni-this->get_X0())+(cg*ca+sg*sb*sa)*1.0/size_pad*(Projectionj-this->get_Y0())+cb*sa*this->get_Z0();
        	//K=(sg*sa+cg*sb*ca)*1.0/size_pad*(Projectioni-this->get_X0())+(-cg*sa+sg*sb*ca)*1.0/size_pad*(Projectionj-this->get_Y0())+cb*ca*this->get_Z0();
        	K=this->NbrPlate()+1;
          if(IsScinti==true)Distribution_exp["SCINTILLATOR"][this->NbrPlate()]->Fill(ceil(I),ceil(J));
          else Distribution_exp["NORMAL"][this->NbrPlate()]->Fill(ceil(I),ceil(J));
          nhit=thisPlan->countHitAt(Projectioni,Projectionj,/*6*10.4125*/_dlimforPad,ceil(I),ceil(J),K,this->GetIp(),this->GetIm(),this->GetJp(),this->GetJm(),IsScinti);
        } 
        else 
        {  
            nhit=thisPlan->countHitAtStrip(Projectioni,_dlimforStrip,IsScinti);
        }
        if (nhit>0) 
	{
		nombreTestsOK++;	
		nombreTestsOKShort++;
	}
        sommeNombreHits+=nhit;
        sommeNombreHitsShort+=nhit;
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
    std::cout<<std::endl;
}

void AnalysisProcessor::PrintStat(bool IsScinti)
{
    ofstream fichier;
    std::string namme ="Results";
    if(IsScinti==true) namme+="Scinti";
    namme+=".txt";
    fichier.open(namme, ios::out | ios::app);  //déclaration du flux et ouverture du fichier
    if(fichier) 
    { // si l'ouverture a réussi
        fichier<<_NbrRun<<"   ";
        if(IsScinti==false)
        {
          for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
          {
            fichier<<testedPlanList[i].efficiency()<<" "<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" "<<testedPlanList[i].multiplicity()<<" 0 "<<"  ";
          }
        }
        else
        {
          for(unsigned int i=0; i!=testedPlanListScinti.size(); ++i) 
          {
            fichier<<testedPlanListScinti[i].efficiency()<<" "<<sqrt(testedPlanListScinti[i].GetNumberOK()*testedPlanListScinti[i].efficiency()*(1-testedPlanListScinti[i].efficiency()))*1.0/testedPlanListScinti[i].GetNumberOK()<<" "<<testedPlanListScinti[i].multiplicity()<<" 0 "<<"  ";
          }
        }
        fichier<<std::endl;
        fichier.close();  // on referme le fichier
    }
    std::cout<<"Run Number : "<<_NbrRun<<std::endl;
    if(IsScinti==false)
    {
      for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
      {
        // testedPlanList[i].print();
        std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanList[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanList[i].efficiency()<<" Error : "<<setw(6)<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" Multiplicity : "<<setw(6)<<testedPlanList[i].multiplicity()<<normal<<std::endl;
      }
    }
    else
    {
      for(unsigned int i=0; i!=testedPlanListScinti.size(); ++i) 
      {
        // testedPlanList[i].print();
        std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanListScinti[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanListScinti[i].efficiency()<<" Error : "<<setw(6)<<sqrt(testedPlanListScinti[i].GetNumberOK()*testedPlanListScinti[i].efficiency()*(1-testedPlanListScinti[i].efficiency()))*1.0/testedPlanListScinti[i].GetNumberOK()<<" Multiplicity : "<<setw(6)<<testedPlanListScinti[i].multiplicity()<<normal<<std::endl;
      }
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

    if(myReader) 
    {
      myReader->Read(_FileNameGeometry,geom);
      geom.PrintGeom();
      std::map<int, Dif > Difs=geom.GetDifs();
      std::map<int,int> PlansType;
      for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it) 
      {
        if(geom.GetDifType(it->first)==scintillator) names.push_back("SCINTILLATOR");
        if(geom.GetDifType(it->first)!=temporal&&geom.GetDifType(it->first)!=scintillator&&geom.GetDifType(it->first)!=tcherenkov) 
        {
          PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
        //PlansType[geom.GetDifNbrPlate(it->first)-1]=geom.GetDifType(it->first);
        }
      }
      FillDelimiter(_Delimiters,PlansType.size(),Delimiter);
      for(std::map<int, int >::iterator it=PlansType.begin(); it!=PlansType.end(); ++it) 
	    {
	      if(Delimiter.size()!=0) testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second,Delimiter[it->first+1][1],Delimiter[it->first+1][0],Delimiter[it->first+1][3],Delimiter[it->first+1][2]));
        else testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second,0,0,0,0));
        std::string a="Distribution hit selectionner par analysis"+ std::to_string( (long long int) it->first +1 );
        std::string b="Hit expected"+ std::to_string( (long long int) it->first +1 );
        std::string c="Efficiency of the voisinage of the pad"+ std::to_string( (long long int) it->first +1 );
        std::string d="Efficiency of the Asic"+ std::to_string( (long long int) it->first +1 );
        std::string e="Distance to the expected hit in X "+ std::to_string( (long long int) it->first +1 );
        std::string f="Distance to the expected hit in Y "+ std::to_string( (long long int) it->first +1 );
        std::string g="Multiplicity of the voisinage of the pad "+ std::to_string( (long long int) it->first +1 );
        for(unsigned int i=0;i<names.size();++i)
        {
          if(it->second==positional) 
          {
            Distribution_hits[names[i]].push_back(new TH2F((a+names[i]).c_str(),(a+names[i]).c_str(),100,0,100,100,0,100));
            Distribution_exp[names[i]].push_back(new TH2F((b+names[i]).c_str(),(b+names[i]).c_str(),100,0,100,100,0,100));
            Efficiency_pads[names[i]].push_back(new TH2F((c+names[i]).c_str(),(c+names[i]).c_str(),100,0,100,100,0,100));
            Efficiency_asics[names[i]].push_back(new TH2F((d+names[i]).c_str(),(d+names[i]).c_str(),100,0,100,100,0,100));
                
            Multiplicity_pads[names[i]].push_back(new TH2F((g+names[i]).c_str(),(g+names[i]).c_str(),100,0,100,100,0,100));
		        HowLongFromExpectedX[names[i]].push_back(new TH1F((e+names[i]).c_str(),(e+names[i]).c_str(),2*(_dlimforStrip),-_dlimforStrip,_dlimforStrip));
            HowLongFromExpectedY[names[i]].push_back(new TH1F((f+names[i]).c_str(),(f+names[i]).c_str(),2*(_dlimforStrip),-_dlimforStrip,_dlimforStrip));
          } 
          else 
          {
            Distribution_hits[names[i]].push_back(new TH2F((a+names[i]).c_str(),(a+names[i]).c_str(),100,0,100,100,0,100));
            Distribution_exp[names[i]].push_back(new TH2F((b+names[i]).c_str(),(b+names[i]).c_str(),100,0,100,100,0,100));
            Efficiency_pads[names[i]].push_back(new TH2F((c+names[i]).c_str(),(c+names[i]).c_str(),100,0,100,100,0,100));
            Efficiency_asics[names[i]].push_back(new TH2F((d+names[i]).c_str(),(d+names[i]).c_str(),100,0,100,100,0,100));
            HowLongFromExpectedX[names[i]].push_back(new TH1F((e+names[i]).c_str(),(e+names[i]).c_str(),2*(_dlimforPad),-_dlimforPad,_dlimforPad));
            HowLongFromExpectedY[names[i]].push_back(new TH1F((f+names[i]).c_str(),(f+names[i]).c_str(),2*(_dlimforPad),-_dlimforPad,_dlimforPad));
            Multiplicity_pads[names[i]].push_back(new TH2F((g+names[i]).c_str(),(g+names[i]).c_str(),100,0,100,100,0,100));
          }
            //Correlations.push_back(new TH2F(b.c_str(),b.c_str(),200,0,200,200,0,200));
        }
      }
      if(_ShortEfficiency!=0)
	    {
		    for(unsigned int i=0;i<testedPlanList.size();++i)
		    {
			    std::string name="Short_Efficiency"+patch::to_string(i+1);
			    for(unsigned int i=0;i<names.size();++i)
          {
			      Short_Efficiency[names[i]].push_back(new TH1F((name+names[i]).c_str(),name.c_str(),1000,0,_ShortEfficiency*1000));
			    }
		    }
	    }
	  } 
	  else 
	  {
        std::cout << "Reader type n'existe pas !!" << std::endl;
        std::exit(1);
    }
    delete myReader;
    testedPlanListScinti=testedPlanList;
}

void AnalysisProcessor::processEvent( LCEvent * evtP )
{ 
    _NbrRun=evtP->getRunNumber();
    NbrRunn=_NbrRun;
    Plans.clear();
    PlansScintillator.clear();
    if (evtP != nullptr) {
        
        std::vector<std::string>names=*evtP->getCollectionNames();
        for(unsigned int i=0; i< _hcalCollections.size(); i++) 
        {
            IsScinti=false;
            LCCollection* col=nullptr;
            for(unsigned int i=0;i<names.size();++i) 
            {
              if(names[i]=="SDHCAL_HIT") 
              {
                col = evtP ->getCollection(names[i].c_str());
                _eventNr=evtP->getEventNumber();
                eventnbrr=_eventNr;
                if(_eventNr %1000 ==0)std::cout<<"Event Number : "<<_eventNr<<std::endl;
              }
              else if(names[i]=="SDHCAL_HIT_SC") 
              {
                //std::cout<<"SC one"<<std::endl;
                col = evtP ->getCollection(names[i].c_str());
                _eventNrSC=evtP->getEventNumber();
                if(_eventNrSC %1000 ==0)std::cout<<"Event Number Scintillator : "<<_eventNrSC<<std::endl;
                IsScinti=true;
              }
            }
            if(col == nullptr) 
            {
                if(IsScinti==true) std::cout<< red << "TRIGGER WITH SCINTILLATOR SKIPED ..."<< normal <<std::endl;
                else  std::cout<< red << "TRIGGER WITH SKIPED ..."<< normal <<std::endl;
                break;
            }
            CellIDDecoder<CalorimeterHit> cd(col);
            int numElements = col->getNumberOfElements();
            for (int ihit=0; ihit < numElements; ++ihit) {
                CalorimeterHit *raw_hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(ihit)) ;
                if (raw_hit != nullptr) 
                {
                    int dif_id=cd(raw_hit)["Dif_id"];
                    int I=cd(raw_hit)["I"];
                    int J=cd(raw_hit)["J"];
                    if(IsScinti==true)
                    {
                      PlansScintillator[geom.GetDifNbrPlate(dif_id)-1].addHit(raw_hit);
                      PlansScintillator[geom.GetDifNbrPlate(dif_id)-1].SetType(geom.GetDifType(dif_id));
                    }
                    else
                    {
                      Plans[geom.GetDifNbrPlate(dif_id)-1].addHit(raw_hit);
                      Plans[geom.GetDifNbrPlate(dif_id)-1].SetType(geom.GetDifType(dif_id));
                    }
                }
            }
        }
        if(IsScinti==true) for (std::vector<testedPlan>::iterator iter=testedPlanListScinti.begin(); iter != testedPlanListScinti.end(); ++iter) iter->testYou(PlansScintillator,true);
        else for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter) iter->testYou(Plans,false);
        if(_ShortEfficiency!=0) 
        {
        	if(_eventNr%_ShortEfficiency==0)
        	{
        	  PrintStatShort(0);
        	  for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
		        {
              testedPlanList[i].ClearShort();
            }
          }
        	if(_eventNrSC%_ShortEfficiency==0)
        	{
        	  PrintStatShort(1);
            for(unsigned int i=0; i!=testedPlanList.size(); ++i) 
		        {
              testedPlanListScinti[i].ClearShort();
            }
          }
    	  }
      }
}

void AnalysisProcessor::end()
{
    for(std::map<std::vector<int>,std::vector<int>>::iterator it=Efficiency_per_pad.begin();it!=Efficiency_per_pad.end();++it)
    {
        unsigned int was_at_least_a_hit=0;
        unsigned int number_of_hits=0;
        for(unsigned int i =0;i!=(it->second).size();++i)
	      {
                
		      if(it->second[i]>0)
          {
			      was_at_least_a_hit++;
			      number_of_hits+=it->second[i];
                         //std::cout<<red<<was_at_least_a_hit<<"  "<<number_of_hits<<normal<<std::endl;
		      }
	      }
        //std::cout<<it->first[0]<<"  "<<it->first[1]<<"  "<<it->first[2]-1<<"  "<<was_at_least_a_hit*1.0/(it->second).size()<<"  "<<number_of_hits*1.0/was_at_least_a_hit<<std::endl;
        Efficiency_pads["NORMAL"][it->first[2]-1]->Fill(it->first[0],it->first[1],was_at_least_a_hit*1.0/(it->second).size());
        if(was_at_least_a_hit!=0)Multiplicity_pads["NORMAL"][it->first[2]-1]->Fill(it->first[0],it->first[1],number_of_hits*1.0/was_at_least_a_hit);
    }
    
    for(std::map<std::vector<int>,std::vector<int>>::iterator it=Efficiency_per_padScinti.begin();it!=Efficiency_per_padScinti.end();++it)
    {
        unsigned int was_at_least_a_hit=0;
        unsigned int number_of_hits=0;
        for(unsigned int i =0;i!=(it->second).size();++i)
	      {
                
		      if(it->second[i]>0)
          {
			      was_at_least_a_hit++;
			      number_of_hits+=it->second[i];
                         //std::cout<<red<<was_at_least_a_hit<<"  "<<number_of_hits<<normal<<std::endl;
		      }
	      }
        //std::cout<<it->first[0]<<"  "<<it->first[1]<<"  "<<it->first[2]-1<<"  "<<was_at_least_a_hit*1.0/(it->second).size()<<"  "<<number_of_hits*1.0/was_at_least_a_hit<<std::endl;
        Efficiency_pads["SCINTILLATOR"][it->first[2]-1]->Fill(it->first[0],it->first[1],was_at_least_a_hit*1.0/(it->second).size());
        if(was_at_least_a_hit!=0)Multiplicity_pads["SCINTILLATOR"][it->first[2]-1]->Fill(it->first[0],it->first[1],number_of_hits*1.0/was_at_least_a_hit);
    }
    
    
    
    std::cout<<"List of counters : "<<std::endl;
    for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter) iter->print();
    std::cout<<"List of counters with Scintillator: "<<std::endl;
    for (std::vector<testedPlan>::iterator iter=testedPlanListScinti.begin(); iter != testedPlanListScinti.end(); ++iter) iter->print();
    std::string b="Results_Analysis_"+std::to_string( (long long int) _NbrRun)+".root";
    TFile *hfile = new TFile(b.c_str(),"RECREATE");
    tt->Write();
    for(unsigned int naa=0;naa<names.size();++naa)
    {
      Distribution_exp_tgraph[names[naa]]->Write();
      Distribution_hits_tgraph[names[naa]]->Write();
      for(unsigned int i =0 ;i!=HowLongFromExpectedX.size();++i)
      {
        std::string plate="Plate "+ patch::to_string(i+1);
   	    if(naa==0)hfile->mkdir(plate.c_str(),plate.c_str());
    	  hfile->cd(plate.c_str());
	      HowLongFromExpectedX[names[naa]][i]->Write();
    	  HowLongFromExpectedY[names[naa]][i]->Write();
        Distribution_hits[names[naa]][i]->Write();
        if(_ShortEfficiency!=0)Short_Efficiency[names[naa]][i]->Write();
        Distribution_exp[names[naa]][i]->Write();
        Efficiency_pads[names[naa]][i]->Write();
        Multiplicity_pads[names[naa]][i]->Write();
      }
    
        for(unsigned int i=0; i<Distribution_hits.size(); ++i) 
        {
          delete Distribution_hits[names[naa]][i];
          delete Distribution_exp[names[naa]][i];
          delete HowLongFromExpectedX[names[naa]][i];
          delete HowLongFromExpectedY[names[naa]][i];
          delete Efficiency_pads[names[naa]][i];
          delete Multiplicity_pads[names[naa]][i];
          if(_ShortEfficiency!=0) delete Short_Efficiency[names[naa]][i];
          //delete Correlations[i];
        }
        delete Distribution_hits_tgraph[names[naa]];
      delete Distribution_exp_tgraph[names[naa]];
    }
    hfile->Close();
    delete hfile;
    delete Branch1;
    delete Branch2;
    delete Branch3;
    delete Branch4;
    delete Branch5;
    delete Branch6;
    std::cout<<"Results"<<std::endl;
    PrintStat(0);
    std::cout<<"Results with scintillators"<<std::endl;
    PrintStat(1);
    /*for(std::map<std::vector<int>,std::vector<int>>::iterator it=Efficiency_per_pad.begin();it!=Efficiency_per_pad.end();++it)
	{
		std::cout<<it->first[0]<<"  "<<it->first[1]<<"  "<<it->first[2]<<"  "<<it->second.size()<<std::endl;
	}*/
}
