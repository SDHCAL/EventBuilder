#include "Analysis/Analysis.h"
#include <iostream>
#include <string>
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
#ifndef COLORS_H
#define normal " "
#define red "  "
#define vert " "
#define blanc " "
#define green " "
#define blue " "
#endif
#include <fstream>

using namespace std;
std::vector<TH2F*>Distribution_hits;
std::vector<TH2F*>Correlations;
void AnalysisProcessor::processRunHeader( LCRunHeader* run) 
{ 

} 

int plan::countHitAt(double& x, double& y, double dlim)
{
	CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
  int n=0;
  for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) 
	{        
	  if(fabs(x-(*it)->getPosition()[0])<dlim&&fabs(y-(*it)->getPosition()[1])<dlim)
		{ 
			n++;
      Distribution_hits[cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);}
		}
  	return n;
}

int plan::countHitAtStrip(double& x, double dlim)
{
  int n=0;
  CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
  for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!= hits.end(); ++it) 
	{
		if(fabs(x-(*it)->getPosition()[0])<dlim) 
    {
    	n++;
	  Distribution_hits[cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);}
	}
  return n;
}

void plan::computeBarycentre( )
{
  for (int i=0; i<3; i++) barycentre[i]=0;
  for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it)
  {
	for (int i=0; i<3; i++)barycentre[i]+=(*it)->getPosition()[i];
  }
  if (nHits() != 0)
  for (int i=0; i<3; i++) barycentre[i]/=nHits();      
}

void plan::computeMaxima()
{
  for (int i=0; i<3; i++)
  { 
    min[i]=10000000;
    max[i]=-10000000;
  }
  for (std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it)
  {
    for (int i=0; i<3; i++)
    {
	if((*it)->getPosition()[i]<min[i])min[i]=(*it)->getPosition()[i];
	if((*it)->getPosition()[i]>max[i])max[i]=(*it)->getPosition()[i];
    }
  }     
}

TF1 myfit = TF1("myfit","[1]*x + [0]", 0, 50);
TF1 myfit2 = TF1("myfit2","[1]*x + [0]", 0, 50);
void testedPlan::testYou(std::map<int,plan>& mapDIFplan)
{
  counts[TESTYOUCALLED]++;
  std::vector<plan*> plansUsedForTrackMaking;
  plan* thisPlan=nullptr;
  for (std::map<int,plan>::iterator it=mapDIFplan.begin(); it!=mapDIFplan.end(); ++it)
  {
    if (Nbr!=it->first) plansUsedForTrackMaking.push_back(&(it->second));
    else thisPlan=&(it->second);
  }
  
  for (std::vector<plan*>::iterator it=plansUsedForTrackMaking.begin();it != plansUsedForTrackMaking.end(); ++it) if ((*it)->nHits()>=_NbrHitPerPlaneMax ) return;
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
    if(p.GetType()==pad){gryz.SetPoint(i,p.barycentreZ(),p.barycentreY());gryz.SetPointError(i,p.ErrorZ(),p.ErrorY());}
    grxz.SetPointError(i,p.ErrorZ(),p.ErrorX());
  }
 
  myfit.SetParameter(1, 0.1);
  myfit.SetParameter(0, 0.1);
  grxz.Fit("myfit","Q");	
  TF1 *myfitxz = (TF1*) grxz.GetFunction("myfit");
  double  kxz = myfitxz->GetChisquare();
  if (kxz>= _Chi2) return;
  double pxz0 = myfitxz->GetParameter(0);           
  double  pxz1 = myfitxz->GetParameter(1);
  counts[XZTRACKFITPASSED]++;
  myfit2.SetParameter(1, 0.1);
  myfit2.SetParameter(0, 0.1);
  gryz.Fit("myfit2","Q");	
  TF1 *myfityz = (TF1*) gryz.GetFunction("myfit2");
  double  kyz = myfityz->GetChisquare();
  if (kyz>= _Chi2) return;
  double pyz0 = myfityz->GetParameter(0);           
  double  pyz1 = myfityz->GetParameter(1);
  counts[YZTRACKFITPASSED]++;
  double Zexp=this->GetZexp(pxz0,pyz0,pxz1,pyz1);
  ///////////////////////////////
  //double Xexp = pxz0+pxz1*Zexp;
  //double Yexp = pyz0+pyz1*Zexp;
  ///////////////////////////////
  double Projectioni=GetProjectioni(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
  double Projectionj=GetProjectionj(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
  if(Projectioni<=187.425&&Projectioni>=135.3625&&Projectionj<=218.6625&&Projectionj>=166.6)
  {
    nombreTests++;
    if (nullptr==thisPlan) return;
    counts[NOHITINPLAN]++;
    int nhit;
    if(this->GetType()==pad)
    {
      nhit=thisPlan->countHitAt(Projectioni,Projectionj,6*10.4125);
    }
    else
    {
      nhit=thisPlan->countHitAtStrip(Projectioni,40);
    }
    if (nhit>0) nombreTestsOK++;
    sommeNombreHits+=nhit;
  }
  delete myfityz;
  delete myfitxz;
  ///////////////////////////////////////////////////////////////////////////////////////////
}

void testedPlan::print() 
{ 
  std::cout<<"Plane Number (in geometry file): "<<Nbr<<" Z="<<Z0<<" NombreTests="<<nombreTests<<" pass="<<nombreTestsOK<<"  sommeNHits="<<sommeNombreHits<< "  (type=" <<GetType()<<" )"<<std::endl;
  for (int i=0; i<NCOUNTERS; i++) std::cout << i<<":"<<counts[i]<<"  ";
  std::cout<<std::endl;
}

void AnalysisProcessor::PrintStat()
{
 ofstream fichier;
 fichier.open("Results.txt", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
 if(fichier)  // si l'ouverture a réussi
 { 
		fichier<<_NbrRun<<"   ";
    for(unsigned int i=0; i!=testedPlanList.size();++i) 
    {
    	fichier<<testedPlanList[i].efficiency()<<" "<<sqrt(testedPlanList[i].GetNumberOK()*testedPlanList[i].efficiency()*(1-testedPlanList[i].efficiency()))*1.0/testedPlanList[i].GetNumberOK()<<" "<<testedPlanList[i].multiplicity()<<" 0 "<<"  ";
    }
    fichier<<std::endl;
    fichier.close();  // on referme le fichier
  }
 
  for(unsigned int i=0; i!=testedPlanList.size();++i) 
  {
 		// testedPlanList[i].print();
 	 	std::cout<<_NbrRun<<std::endl;
		std::cout<<"Efficacite "<<testedPlanList[i].efficiency()<<" Multiplicite "<<testedPlanList[i].multiplicity()<<std::endl;
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
}

AnalysisProcessor::~AnalysisProcessor() {}

void AnalysisProcessor::init()
{
  printParameters();
  ReaderFactory readerFactory;
  Reader* myReader = readerFactory.CreateReader(_ReaderType);
  
  if(myReader)
  {
    myReader->Read(_FileNameGeometry,geom);
    geom.PrintGeom();
    std::map<int, Dif > Difs=geom.GetDifs();
    std::map<int,int> PlansType;
    for(std::map<int, Dif >::iterator it=Difs.begin();it!=Difs.end();++it)
    {
      if(geom.GetDifType(it->first)!=temporal)
      {
        PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
        //PlansType[geom.GetDifNbrPlate(it->first)-1]=geom.GetDifType(it->first);
      }
    }
    for(std::map<int, int >::iterator it=PlansType.begin();it!=PlansType.end();++it)
    {
      testedPlanList.push_back(testedPlan(it->first,geom.GetPlatePositionX(it->first),geom.GetPlatePositionY(it->first),geom.GetPlatePositionZ(it->first),geom.GetDifPlateAlpha(it->first),geom.GetDifPlateBeta(it->first),geom.GetDifPlateGamma(it->first),it->second));
      std::string b="Correlations"+ std::to_string( (long long int) it->first +1 );
      std::string a="Distribution hit selectionner par analysis"+ std::to_string( (long long int) it->first +1 );
      if(it->second==positional) {Distribution_hits.push_back(new TH2F(a.c_str(),a.c_str(),128,0,128,1,0,50));}
      else {Distribution_hits.push_back(new TH2F(a.c_str(),a.c_str(),100,0,100,100,0,100));}
     Correlations.push_back(new TH2F(b.c_str(),b.c_str(),200,0,200,200,0,200));
			
    }
  }
  else
  {
    std::cout << "Reader type n'existe pas !!" << std::endl;
    std::exit(1);
  }
  delete myReader;
}
void AnalysisProcessor::processEvent( LCEvent * evtP ) 
{	
	_NbrRun=evtP->getRunNumber();  
	Plans.clear();
  if (evtP != nullptr)
  {
    _eventNr=evtP->getEventNumber();
    for(unsigned int i=0; i< _hcalCollections.size(); i++)
    {
	    LCCollection* col = evtP ->getCollection(_hcalCollections[i].c_str());
	    if(col == nullptr)
	    {
	      std::cout<< red << "TRIGGER SKIPED ..."<< normal <<std::endl;
	      break;
	    } 
	    CellIDDecoder<CalorimeterHit> cd(col);
	    int numElements = col->getNumberOfElements();
	    for (int ihit=0; ihit < numElements; ++ihit) 
	    {
	      CalorimeterHit *raw_hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(ihit)) ;
	      if (raw_hit != nullptr)
	      {
	        int dif_id=cd(raw_hit)["Dif_id"];
		int I=cd(raw_hit)["I"];
		int J=cd(raw_hit)["J"];
	        Plans[geom.GetDifNbrPlate(dif_id)-1].addHit(raw_hit);
	        Plans[geom.GetDifNbrPlate(dif_id)-1].SetType(geom.GetDifType(dif_id));
                for(int jhit=ihit+1;jhit<numElements;++jhit)
		{
			CalorimeterHit *raw_hit2 = dynamic_cast<CalorimeterHit*>( col->getElementAt(jhit)) ;
			int dif_id2=cd(raw_hit2)["Dif_id"];
			int I2=cd(raw_hit2)["I"];
			int J2=cd(raw_hit2)["J"];
			if((geom.GetDifNbrPlate(dif_id)==1 && geom.GetDifNbrPlate(dif_id2)==4)&&(raw_hit->getTime()==raw_hit2->getTime())) Correlations[4]->Fill((I-1),I2);}
	      }
	    }
    }
    for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter)
    {
      iter->testYou(Plans);
    }
  }
}

void AnalysisProcessor::end()
{
  std::string b="Results_Analysis_"+std::to_string( (long long int) _NbrRun)+".root";
	TFile *hfile = new TFile(b.c_str(),"RECREATE");
  for(unsigned int i=0; i<Distribution_hits.size();++i)
	{
		Distribution_hits[i]->Write();
		Correlations[i]->Write();
	}
	for(unsigned int i=0; i<Distribution_hits.size();++i)
	{
		delete Distribution_hits[i];
		delete Correlations[i];
	}
  hfile->Close();
	delete hfile;
  PrintStat();
}