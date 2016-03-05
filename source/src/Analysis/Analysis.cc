#include "Analysis/Analysis.h"
#include "Intro.h"
#include <iostream>
#include <string>
#include <iomanip>
#include "Progress.h"
#include "marlin/Processor.h"
#include "Utilities.h"
#include "THnSparse.h"
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
#include <cstdlib>
#include "Version.h"
#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "marlin/Global.h"
#include "TH3F.h"
#include "TTree.h"
#include "Trivent/Mapping.h"
#include "TGraph2D.h"
#include "Utilities.h"
#include <map>
#include <cmath>
#include <IMPL/RawCalorimeterHitImpl.h>
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
#define degtorad 0.0174532925
using namespace std;
ReaderFactory readerFactory;
std::vector<std::vector<std::array<double,7>>>useforrealrate;
std::map<std::string,std::vector<std::vector<TH1F*>>>Short_Efficiency;
std::map<std::string,std::vector<std::vector<TH1F*>>>Short_Multiplicity;
std::map<std::string,unsigned int>Numbers;
unsigned int SumCombinaison(unsigned int n,unsigned int kmin)
{
  unsigned int sumcomb=0;
  for(unsigned int i=kmin;i<=n;++i) sumcomb+=TMath::Binomial(n,i);
  return sumcomb;
}

class ToTreee
{
  public:
  double ChiXZ, ChiYZ, CDXZ, CDYZ, OrdXZ, OrdYZ;
};

ToTreee totreee;
TTree* tt= new TTree("Tree", "Tree");
TBranch* Branch1 =  tt->Branch("ChiXZ",&(totreee.ChiXZ));
TBranch* Branch2 =  tt->Branch("ChiYZ",&(totreee.ChiYZ));
TBranch* Branch3 =  tt->Branch("CDXZ",&(totreee.CDXZ));
TBranch* Branch4 =  tt->Branch("CDYZ",&(totreee.CDYZ));
TBranch* Branch5 =  tt->Branch("OrdXZ",&(totreee.OrdXZ));
TBranch* Branch6 =  tt->Branch("OrdYZ",&(totreee.OrdYZ));

THnSparseD* hss2=nullptr;
THnSparseD* hss=nullptr;
std::vector<TGraphErrors>xzaxis;
std::map<int,int long>RealNumberPlane;
std::vector<TGraphErrors>yzaxis;
std::vector<std::string> names={"SDHCAL_HIT"};
enum Tresholds {Threshold1,Threshold2,Threshold3,Threshold12,Threshold23,Thresholdall};
std::array<std::string,6>Thresholds_name{"Threshold1","Threshold2","Threshold3","Threshold12","Threshold23","Thresholdall"};
std::map<std::string,std::vector<TH2F*>>Distribution_hits;
std::array<std::map<std::string,std::vector<TH2F*>>,6>All_The_Hits_From_Trivent;
std::array<std::map<std::string,std::vector<TH2F*>>,6>Efficiency_pads;
std::array<std::map<std::string,std::vector<TH2F*>>,6>Multiplicity_pads;
std::array<std::map<std::string,std::vector<TH2F*>>,6>Efficiency_asics;
std::map<std::string,std::vector<TH1F*>>HowLongFromExpectedX;
std::map<std::string,std::vector<TH1F*>>HowLongFromExpectedY;
std::map<std::string,std::vector<TH2D*>>difxy;
std::map<int,TH1D*> estimateur;
std::map<int,TH1D*> estimateur2;
TH1D* Rate_Estimation = new TH1D("Nbr_Tracks_divided_by_time","Nbr_Tracks_divided_by_time",1,0,1);
std::map<int,std::map<int,TH2F*>> Nbrviewedtracks;
std::map<int,std::map<int,TH2F*>> eff;
std::map<int,std::map<int,TH2F*>> effsquared;
std::map<int,std::map<int,TH2F*>> estimatedrate;
std::map<int,TH1D*> estimateurmax;
std::map<int,TH1D*> rate;
std::map<std::string,std::vector<TH1D*>>difr;
std::array<std::map<std::string,std::vector<TH2F*>>,3>ThresholdMap;
std::map<std::string,std::vector<TH2F*>>Gain;
std::map<int,std::map<int,TH2F*>>RealRateWithSelectedZone;
std::map<int,std::map<int,TH2F*>>RealRateWithSelectedZonemin;
std::map<int,TH2F*>RateSeen;
std::map<int,TH2F*>RateSeenPon;
std::map<int,TH2F*>EfficacityVsRate;
std::map<std::string,TGraph2D*> Distribution_exp_tgraph;
std::map<std::string,std::vector<TH2F*>>Distribution_exp;
std::map<std::string,TGraph2D*> Distribution_hits_tgraph;
std::map<std::string,TGraph2D*> RateSeen3D;
std::vector<std::string>Maketracks{"SDHCAL_HIT","SDHCAL_HIT_SC"};
std::vector<std::string>Makeeffi{"NOISE_ESTIMATION_BEFORE","NOISE_ESTIMATION_AFTER"};
unsigned int numbertracks;

void AnalysisProcessor::PrintStatShort()
{
  for(unsigned f=0;f!=_hcalCollections.size();++f)
    {
      if(_ShortEfficiency!=0)
	{
	  if(Numbers[_hcalCollections[f]]%_ShortEfficiency==0&&Numbers[_hcalCollections[f]]!=0)
	    {
	      PrintStatShort(_hcalCollections[f]);
	      for(unsigned int i=0; i!=testedPlanList.size(); ++i)
		{
		  testedPlanList[i].ClearShort(_hcalCollections[f]);
		}
            }
	}
    }
}


void testedPlan::print(std::string name)
{
    std::cout<<red<<"Plane Number (in geometry file): "<<geomplan.NbrPlate()+1<<" Z = "<<geomplan.GetZ0()<<" : "<<normal<<std::endl;
    std::cout<<blue<<"Number of Test : "<<Counts[name][0]<<"; with >="<<_NbrPlaneUseForTracking<<" planes for tracking : "<<Counts[name][1]<<"; with ChiXZ <"<<_Chi2<<" : "<<Counts[name][2]<<"; with ChiYZ <"<<_Chi2<<" : "<<Counts[name][3]<<" ; with track in the Delimiters "<<nombreTests[name]<<"; with hits in it : "<<Counts[name][4]<<" ; "<<std::endl;
    std::cout<<red<<"with hits in dlim : "<<normal<<std::endl;
    for(unsigned int i=0;i!=Thresholds_name.size();++i)
    {
    std::cout<<Thresholds_name[i]<<" : "<<GetNumberOK(i,name)<<" "<<normal;
    }
    std::cout<<std::endl;
    std::cout<<red<<"Sum of hits in dlim : "<<normal<<std::endl;
    for(unsigned int i=0;i!=Thresholds_name.size();++i)
    {
    std::cout<<Thresholds_name[i]<<" : "<<GetNombreHits(i,name)<<" "<<normal;
    }
    std::cout<<std::endl;
}

void testedPlan::testYou(std::map<std::string,std::map<int,hitsInPlan>>&mapDIFplan,std::vector<testedPlan>& tested)
{
  for(std::map<std::string,std::map<int,hitsInPlan>>::iterator itt=mapDIFplan.begin();itt!=mapDIFplan.end();++itt)
    {
      bool Doit=false;
      for(unsigned int mm=0;mm!=Maketracks.size();++mm)
	{
	  if(itt->first!=Maketracks[mm]) continue;
	  else Doit=true;
	}
      if(Doit==true)
	{
	  std::vector<std::string>ToComputeEffi=Makeeffi;
	  ToComputeEffi.push_back(itt->first);
	  //std::cout<<itt->first<<std::endl;
	  for(unsigned int i=0;i!=ToComputeEffi.size();++i) Counts[ToComputeEffi[i]][TESTYOUCALLED]++;
	  std::vector<hitsInPlan*> plansUsedForTrackMaking;
	  hitsInPlan* thisPlan=nullptr;
	  std::vector<int>PlaneNbr;
	  for (std::map<int,hitsInPlan>::iterator it=mapDIFplan[itt->first].begin(); it!=mapDIFplan[itt->first].end(); ++it)
	    {
	      if (geomplan.NbrPlate()!=it->first) 
		{
		  if(it->second.nHits()>0)//Verify is hits are present
		    {
		      plansUsedForTrackMaking.push_back(&(it->second));
		      PlaneNbr.push_back(it->first);
		    }
		}
	      else thisPlan=&(it->second);
	    }

	  for (std::vector<hitsInPlan*>::iterator it=plansUsedForTrackMaking.begin(); it != plansUsedForTrackMaking.end(); ++it) if ((int)(*it)->nHits()>=_NbrHitPerPlaneMax ) return;
	  if((int)plansUsedForTrackMaking.size()<_NbrPlaneUseForTracking) return;
	  for(unsigned int i=0;i!=ToComputeEffi.size();++i) Counts[ToComputeEffi[i]][NOTOOMUCHHITSINPLAN]++;
	  ////////////////////////////////////////////////////////////////////////////////////
	  trackFitter tfit;
	  bool trackfound=tfit.Find(plansUsedForTrackMaking,_Chi2,this->GetType(),itt->first);
	  if (! tfit.foundTracksXZ() ) return;
	  for(unsigned int i=0;i!=ToComputeEffi.size();++i) Counts[ToComputeEffi[i]][XZTRACKFITPASSED]++;
	  if (! trackfound) return;
	  for(unsigned int i=0;i!=ToComputeEffi.size();++i) Counts[ToComputeEffi[i]][YZTRACKFITPASSED]++;

	  double  kxz = tfit.getXZChisquare();
	  double pxz0 = tfit.getXZParameter(0);
	  double  pxz1 = tfit.getXZParameter(1);
	  double  kyz = tfit.getYZChisquare();
	  if(this->GetType()==positional) kyz=0;
	  double pyz0 = tfit.getYZParameter(0);
	  double  pyz1 = tfit.getYZParameter(1);

	  double Zexp=this->GetZexp(pxz0,pyz0,pxz1,pyz1);
	  //xzaxis.push_back(grxz);
	  //yzaxis.push_back(gryz);
	  ///////////////////////////////
	  //double Xexp = pxz0+pxz1*Zexp;
	  //double Yexp = pyz0+pyz1*Zexp;
	  ///////////////////////////////
	  double Projectioni=GetProjectioni(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
	  double Projectionj=GetProjectionj(pxz0+pxz1*Zexp,pyz0+pyz1*Zexp,Zexp);
	  bool Pass;
	  if(Delimiter.size()==0)Pass=1;
	  else Pass=Projectioni<=this->GetIp()&&Projectioni>=this->GetIm()&&Projectionj<=this->GetJp()&&Projectionj>=this->GetJm();
	  if(Pass)
	    {
	      //std::cout<<std::endl;
	      for(unsigned int i=0;i!=ToComputeEffi.size();++i) nombreTests[ToComputeEffi[i]]++;
	      //std::cout<<nombreTests[itt->first]<<"  "<<itt->first<<std::endl;
	      for(unsigned int i=0;i!=ToComputeEffi.size();++i) nombreTestsShort[ToComputeEffi[i]]++;

	      if (nullptr==thisPlan) return;
	      int I,J,K;
	      double ca=this->get_ca();
	      double sa=this->get_sa();
	      double cb=this->get_cb();
	      double sb=this->get_sb();
	      double cg=this->get_cg();
	      double sg=this->get_sg();
	      Distribution_exp_tgraph[itt->first]->SetPoint(nombreTests[itt->first],Projectioni,Projectionj,Zexp);
        
	      for(unsigned int i=0;i!=ToComputeEffi.size();++i) Counts[ToComputeEffi[i]][NOHITINPLAN]++;
	      //for(unsigned int i=0;i!=ToComputeEffi.size();++i) std::cout<<yellow<<ToComputeEffi[i]<<normal<<std::endl;
	      for(unsigned int i=0;i!=ToComputeEffi.size();++i)
		{
		  std::array<double,6>Thresholds;
		  std::map<std::string,int> nhit;
		  hitsInPlan *planLookedFor;
		  if (ToComputeEffi[i]==itt->first)
		    planLookedFor=&(mapDIFplan[ToComputeEffi[i]][geomplan.NbrPlate()]);
		  else
		    planLookedFor=NULL; 
		  //countHitAt fills some counters and need to be called anyway for the moment
		  hitsInPlan dummy;
		  if (NULL==planLookedFor) planLookedFor=&dummy; 

		  if(this->GetType()==pad)
		    {
		      I=cg*cb*1.0/size_pad*(Projectioni-this->GetX0())+sg*cb*1.0/size_pad*(Projectionj-this->GetY0())+-sb*this->GetZ0();
		      J=(-sg*ca+cg*sb*sa)*1.0/size_pad*(Projectioni-this->GetX0())+(cg*ca+sg*sb*sa)*1.0/size_pad*(Projectionj-this->GetY0())+cb*sa*this->GetZ0();
		      K=this->NbrPlate()+1;
		      Distribution_exp[itt->first][this->NbrPlate()]->Fill(ceil(I),ceil(J));

		      Thresholds=planLookedFor->countHitAt(Projectioni,Projectionj,_dlimforPad,ceil(I),ceil(J),K,this->GetIp(),this->GetIm(),this->GetJp(),this->GetJm(),ToComputeEffi[i]);
		    }
		  else
		    {
		      nhit=planLookedFor->countHitAtStrip(Projectioni,_dlimforStrip,ToComputeEffi[i]);
		      Thresholds[5]=nhit[itt->first];
		    }
          
		  for(int kk=0;kk!=Thresholds.size();++kk)
		    {
		      if (Thresholds[kk]==0)
			{
			  nombreTestsOK[ToComputeEffi[i]][kk]=nombreTestsOK[ToComputeEffi[i]][kk];
			  sommeNombreHits[ToComputeEffi[i]][kk]+=Thresholds[kk];
			}
		      if (Thresholds[kk]>0)
			{
			  nombreTestsOK[ToComputeEffi[i]][kk]++;
			  nombreTestsOKShort[ToComputeEffi[i]][kk]++;
			  sommeNombreHits[ToComputeEffi[i]][kk]+=Thresholds[kk];
			  sommeNombreHitsShort[ToComputeEffi[i]][kk]+=Thresholds[kk];
  
			}
		    }
		}
	      totreee.ChiXZ=kxz;
	      totreee.ChiYZ=kyz;
	      totreee.CDXZ=pxz1;
	      totreee.CDYZ=pyz1;
	      totreee.OrdXZ=pxz0;
	      totreee.OrdYZ=pyz0;
	      tt->Fill();
	    } //end if (Pass)
	} //end if Doit=true
    }
  ///////////////////////////////////////////////////////////////////////////////////////////
}

void AnalysisProcessor::PrintStatShort(std::string name)
{
  static std::map<std::string,std::string>Mess;
  for(unsigned int i=0;i!=_hcalCollections.size();++i)Mess[_hcalCollections[i]]="List of counters "+ _hcalCollections[i]+" : ";
  ofstream fichier;
  for(unsigned int hhh=0;hhh!=Thresholds_name.size();++hhh)
  {
    std::string namee="ResultsShorts"+std::to_string((long long int)_NbrRun);
    namee+="_"+name+"_";
    namee+=Thresholds_name[hhh];
    namee+=".txt";
    fichier.open(namee.c_str(), ios::out | ios::app);  //déclaration du flux et ouverture du fichier
    if(fichier) 
    { // si l'ouverture a réussi
      fichier<<_NbrRun<<"   ";
      for(unsigned int i=0; i!=testedPlanList.size(); ++i)
      {
        if(isfinite(testedPlanList[i].efficiencyShort(hhh,name)))Short_Efficiency[name][i][hhh]->Fill(Numbers[name],testedPlanList[i].efficiencyShort(hhh,name));
        if(isfinite(testedPlanList[i].multiplicityShort(hhh,name)))Short_Multiplicity[name][i][hhh]->Fill(Numbers[name],testedPlanList[i].multiplicityShort(hhh,name));
        fichier<<testedPlanList[i].efficiencyShort(hhh,name)<<" "<<sqrt(testedPlanList[i].GetNumberOKShort(hhh,name)*testedPlanList[i].efficiencyShort(hhh,name)*(1-testedPlanList[i].efficiencyShort(hhh,name)))*1.0/testedPlanList[i].GetNumberOKShort(hhh,name)<<" "<<testedPlanList[i].multiplicityShort(hhh,name)<<" 0 "<<"  ";
      }
      fichier<<std::endl;
    }
    fichier.close();
  }
  for(unsigned int hhh=0;hhh!=Thresholds_name.size();++hhh)
  {
    if(hhh==0)std::cout<<Mess[name]<<std::endl;
    std::cout<<"Run Number : "<<_NbrRun<<" Thresholds "<<Thresholds_name[hhh]<<std::endl;
    for(unsigned int i=0; i!=testedPlanList.size(); ++i)
    {
      std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanList[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanList[i].efficiency(hhh,name)<<" Error : "<<setw(6)<<sqrt(testedPlanList[i].GetNumberOK(hhh,name)*testedPlanList[i].efficiency(hhh,name)*(1-testedPlanList[i].efficiency(hhh,name)))*1.0/testedPlanList[i].GetNumberOK(hhh,name)<<" Multiplicity : "<<setw(6)<<testedPlanList[i].multiplicity(hhh,name)<<normal<<std::endl;
    }
  }
}

std::array<double,6> hitsInPlan::countHitAt(double& x, double& y, double dlim,int Iexpected,int Jexpected,int Kexpected,double Imax,double Imin,double Jmax,double Jmin,std::string collectionName)
{
  static std::map<std::string,unsigned int>Number_hits;
  std::array<double,6>Threshold_Counters;
  for(unsigned int i=0;i!=Threshold_Counters.size();++i)Threshold_Counters[i]=0;
  std::vector<int>IJKexpected={Iexpected,Jexpected,Kexpected};
  std::vector<CalorimeterHit*>Hits=GetHits();
  for (std::vector<CalorimeterHit*>::iterator it=Hits.begin(); it!=Hits.end(); ++it)
  {
    CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
    difxy[collectionName][cd(*it)["K"]-1]->Fill(x-(*it)->getPosition()[0],y-(*it)->getPosition()[1]);
    difr[collectionName][cd(*it)["K"]-1]->Fill(sqrt((x-(*it)->getPosition()[0])*(x-(*it)->getPosition()[0])+(y-(*it)->getPosition()[1])*(y-(*it)->getPosition()[1])));
    if(fabs(x-(*it)->getPosition()[0])<dlim&&fabs(y-(*it)->getPosition()[1])<dlim)
    {
      Number_hits[collectionName]++;
      int Threshold_Hit=(*it)->getEnergy();
      if(Threshold_Hit==Threshold_3)
      {
        Threshold_Counters[Threshold3]++;
        Threshold_Counters[Threshold23]++;
        Threshold_Counters[Thresholdall]++;
      }
      else if(Threshold_Hit==Threshold_2)
      {
        Threshold_Counters[Threshold2]++;
        Threshold_Counters[Threshold12]++;
        Threshold_Counters[Threshold23]++;
        Threshold_Counters[Thresholdall]++;
      }
      else if(Threshold_Hit==Threshold_1)
      {
        Threshold_Counters[Threshold1]++;
        Threshold_Counters[Threshold12]++;
        Threshold_Counters[Thresholdall]++;
      }
      Distribution_hits[collectionName][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
      Distribution_hits_tgraph[collectionName]->SetPoint(Number_hits[collectionName],(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
      HowLongFromExpectedX[collectionName][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
      HowLongFromExpectedY[collectionName][cd(*it)["K"]-1]->Fill((*it)->getPosition()[1]-y);
    }
  }
  for(int i=0;i!=Threshold_Counters.size();++i)
  {
    Efficiency_per_pad[collectionName][i][IJKexpected].push_back(Threshold_Counters[i]);
  }
  return Threshold_Counters;
}

std::map<std::string,int> hitsInPlan::countHitAtStrip(double& x, double dlim,std::string collectionName)
{
  std::array<double,6>Threshold_Counters;
  for(unsigned int i=0;i!=Threshold_Counters.size();++i)Threshold_Counters[i]=0;
  //std::vector<int>IJKexpected={Iexpected,Jexpected,Kexpected};
  std::vector<CalorimeterHit*>Hits=GetHits();
  std::map<std::string,int>N;
  N.clear();
  static std::map<std::string,unsigned int>Number_hits;
  CellIDDecoder<CalorimeterHit>cd("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7" );
  for (std::vector<CalorimeterHit*>::iterator it=Hits.begin(); it!= Hits.end(); ++it) 
  {
    if(fabs(x-(*it)->getPosition()[0])<dlim) 
    {
      N[collectionName]++;
      Number_hits[collectionName]++;
      Distribution_hits[collectionName][cd(*it)["K"]-1]->Fill(cd(*it)["I"],cd(*it)["J"]);
      Distribution_hits_tgraph[collectionName]->SetPoint(Number_hits[collectionName],(*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2]);
      HowLongFromExpectedX[collectionName][cd(*it)["K"]-1]->Fill((*it)->getPosition()[0]-x);
      HowLongFromExpectedY[collectionName][cd(*it)["K"]-1]->Fill(0);
    }
  }
  return N;  
}


void Tracks(std::map<std::string,std::map<int,hitsInPlan>>& mapDIFplan,Geometry geom,std::map<int,geometryplan> geometryplans,std::vector<std::vector<std::array<double,7>>>& useforrealrate)
{
  
  double NbrPlateTotal=geom.GetNumberPlates();
  
  for(std::map<std::string,std::map<int,hitsInPlan>>::iterator itt=mapDIFplan.begin();itt!=mapDIFplan.end();++itt)
  {
    std::vector<std::array<double,7>>XYZExpected;
    bool Doit=false;
    for(unsigned int mm=0;mm!=Maketracks.size();++mm)
    {
      if(itt->first!=Maketracks[mm]) continue;
      else Doit=true;
      
    }
    if(Doit==true)
    {
      std::vector<hitsInPlan*> plansUsedForTrackMaking;
      std::vector<int>PlaneNbr;
      std::vector<int>othersNbr;
      for (std::map<int,hitsInPlan>::iterator it=mapDIFplan[itt->first].begin(); it!=mapDIFplan[itt->first].end(); ++it)
      {
        plansUsedForTrackMaking.push_back(&(it->second));
        PlaneNbr.push_back(it->first);
      }
      for (std::vector<hitsInPlan*>::iterator it=plansUsedForTrackMaking.begin(); it != plansUsedForTrackMaking.end(); ++it) if ((int)(*it)->nHits()>=_NbrHitPerPlaneMax ) return;
      if(plansUsedForTrackMaking.size()<(unsigned int )(_NbrPlaneUseForTrackingRate)) return;
      for(unsigned int i=0;i!=NbrPlateTotal-1;++i)
      {
        int notinPlanNbr=1;
        for(unsigned int j=0;j!=PlaneNbr.size();++j)
        {
          if(PlaneNbr[j]==(int)i) 
          {
            notinPlanNbr=0;
          }
        }
        if(notinPlanNbr==1)
        {
          othersNbr.push_back(i);
          //std::cout<<yellow<<i<<normal<<std::endl;
        }
      }
      TGraphErrors grxz(plansUsedForTrackMaking.size());
      TGraphErrors gryz(plansUsedForTrackMaking.size());
      for (unsigned int i=0; i < plansUsedForTrackMaking.size(); ++i)
      {
        hitsInPlan &hp=*(plansUsedForTrackMaking[i]);
        hp.computeBarycentre();
        double ca=geometryplans[PlaneNbr[i]].get_ca();
        double sa=geometryplans[PlaneNbr[i]].get_sa();
        double cb=geometryplans[PlaneNbr[i]].get_cb();
        double sb=geometryplans[PlaneNbr[i]].get_sb();
        double cg=geometryplans[PlaneNbr[i]].get_cg();
        double sg=geometryplans[PlaneNbr[i]].get_sg();
        int istouched=1;
        double I=cg*cb*1.0/size_pad*(hp.barycentreX()-geometryplans[PlaneNbr[i]].GetX0())+sg*cb*1.0/size_pad*(hp.barycentreY()-geometryplans[PlaneNbr[i]].GetY0())+-sb*geometryplans[PlaneNbr[i]].GetZ0();
        double J=(-sg*ca+cg*sb*sa)*1.0/size_pad*(hp.barycentreX()-geometryplans[PlaneNbr[i]].GetX0())+(cg*ca+sg*sb*sa)*1.0/size_pad*(hp.barycentreY()-geometryplans[PlaneNbr[i]].GetY0())+cb*sa*geometryplans[PlaneNbr[i]].GetZ0();
        int K=geometryplans[PlaneNbr[i]].NbrPlate();
	std::array<double,7>Vec={ceil(I),ceil(J),double(K),double(istouched),hp.barycentreX(),hp.barycentreY(),geometryplans[PlaneNbr[i]].GetZ0()};
        XYZExpected.push_back(Vec);
        hp.computeMaxima(); //NB : computation results are not used
        grxz.SetPoint(i,hp.barycentreZ(),hp.barycentreX());
        if(hp.GetType()==pad)
        {
            gryz.SetPoint(i,hp.barycentreZ(),hp.barycentreY());
            gryz.SetPointError(i,hp.ErrorZ(),hp.ErrorY());
        }
        grxz.SetPointError(i,hp.ErrorZ(),hp.ErrorX());
      }
      grxz.Fit("pol1","QRO","",-50000.0,50000.0);
      TF1 *myfitxz = (TF1*) grxz.GetFunction("pol1");
      double  kxz = myfitxz->GetChisquare();
      if (kxz>= _Chi2Rate) return;
      double pxz0 = myfitxz->GetParameter(0);
      double  pxz1 = myfitxz->GetParameter(1);
      gryz.Fit("pol1","QRO","",-50000.0,50000.0);
      TF1 *myfityz = (TF1*) gryz.GetFunction("pol1");
      double  kyz = myfityz->GetChisquare();
      if (kyz>= _Chi2Rate) return;
      double pyz0 = myfityz->GetParameter(0);
      double  pyz1 = myfityz->GetParameter(1);
      if(itt->first=="SDHCAL_HIT")numbertracks++;
      for (unsigned int i=0; i < othersNbr.size(); ++i)
      {
        int istouched=0;
        double ca=geometryplans[othersNbr[i]].get_ca();
        double sa=geometryplans[othersNbr[i]].get_sa();
        double cb=geometryplans[othersNbr[i]].get_cb();
        double sb=geometryplans[othersNbr[i]].get_sb();
        double cg=geometryplans[othersNbr[i]].get_cg();
        double sg=geometryplans[othersNbr[i]].get_sg();
        double Zexpp=geometryplans[othersNbr[i]].GetZexp(pxz0,pyz0,pxz1,pyz1);
        double iexp=geometryplans[othersNbr[i]].GetProjectioni(pxz0+pxz1*Zexpp,pyz0+pyz1*Zexpp,Zexpp);
        double jexp=geometryplans[othersNbr[i]].GetProjectionj(pxz0+pxz1*Zexpp,pyz0+pyz1*Zexpp,Zexpp);
        double I=cg*cb*1.0/size_pad*(iexp-geometryplans[othersNbr[i]].GetX0())+sg*cb*1.0/size_pad*(jexp-geometryplans[othersNbr[i]].GetY0())+-sb*geometryplans[othersNbr[i]].GetZ0();
        double J=(-sg*ca+cg*sb*sa)*1.0/size_pad*(iexp-geometryplans[othersNbr[i]].GetX0())+(cg*ca+sg*sb*sa)*1.0/size_pad*(jexp-geometryplans[othersNbr[i]].GetY0())+cb*sa*geometryplans[othersNbr[i]].GetZ0();
        int K=geometryplans[othersNbr[i]].NbrPlate();
       //std::cout<<red<<Zexpp<<"  "<<"  "<<"  "<<I<<"  "<<J<<"  "<<K<<"  "<<int(ceil(I))<<"  "<<int(ceil(J))<<"  "<<istouched<<normal<<std::endl;
        std::array<double,7>Vec={ceil(I),ceil(J),double(K),double(istouched),iexp,jexp,Zexpp};
        XYZExpected.push_back(Vec);
        
        
      }
    delete myfityz;
    delete myfitxz;
    }
   if(itt->first=="SDHCAL_HIT")
    {
      useforrealrate.push_back(XYZExpected);
    }
  }
   
}

void AnalysisProcessor::PrintStat(std::string name)
{
    static std::map<std::string,std::string>Mess;
    for(unsigned int i=0;i!=_hcalCollections.size();++i)Mess[_hcalCollections[i]]="List of counters "+ _hcalCollections[i]+" : ";
    ofstream fichier;
    for(unsigned int hhh=0;hhh!=Thresholds_name.size();++hhh)
    {
      std::string namme ="Results";
      namme+="_"+name+"_";
      namme+=Thresholds_name[hhh];
      namme+=".txt";
      fichier.open(namme, ios::out | ios::app);  //déclaration du flux et ouverture du fichier
      if(fichier)
      { // si l'ouverture a réussi
        fichier<<_NbrRun<<";#;";
          for(unsigned int i=0; i!=testedPlanList.size(); ++i)
          {
            fichier<<testedPlanList[i].efficiency(hhh,name)<<";"<<sqrt(testedPlanList[i].GetNumberOK(hhh,name)*testedPlanList[i].efficiency(hhh,name)*(1-testedPlanList[i].efficiency(hhh,name)))*1.0/testedPlanList[i].GetNumberOK(hhh,name)<<";"<<testedPlanList[i].multiplicity(hhh,name)<<";0;";
          }
        fichier<<std::endl;
          // on referme le fichier
      }
      fichier.close();
    }
    
    for(unsigned int hhh=0;hhh!=Thresholds_name.size();++hhh)
    {
    std::cout<<"Run Number : "<<_NbrRun<<" Thresholds "<<Thresholds_name[hhh]<<std::endl;
      for(unsigned int i=0; i!=testedPlanList.size(); ++i)
      {
        std::cout<<green<<setprecision(3)<<"Plane Number (in geometry file) : "<<testedPlanList[i].NbrPlate()+1<< " Efficiency : "<<setw(6)<<testedPlanList[i].efficiency(hhh,name)<<" Error : "<<setw(6)<<sqrt(testedPlanList[i].GetNumberOK(hhh,name)*testedPlanList[i].efficiency(hhh,name)*(1-testedPlanList[i].efficiency(hhh,name)))*1.0/testedPlanList[i].GetNumberOK(hhh,name)<<" Multiplicity : "<<setw(6)<<testedPlanList[i].multiplicity(hhh,name)<<normal<<std::endl;
      }
    }

}

using namespace marlin;
AnalysisProcessor aAnalysisProcessor;
AnalysisProcessor::AnalysisProcessor() : Processor("AnalysisProcessorType")
{
    std::vector<std::string> hcalCollections{"SDHCAL_HIT","SDHCAL_HIT_SC","NOISE_ESTIMATION_BEFORE","NOISE_ESTIMATION_AFTER"};
    registerInputCollections( LCIO::RAWCALORIMETERHIT,"HCALCollections","HCAL Collection Names",_hcalCollections,hcalCollections);
    _FileNameGeometry="";
    registerProcessorParameter("FileNameGeometry","Name of the Geometry File",_FileNameGeometry,_FileNameGeometry);
    _ReaderType="";
    registerProcessorParameter("ReaderType","Type of the Reader needed to read InFileName",_ReaderType ,_ReaderType);
    _Chi2 = 1.0;
    registerProcessorParameter("Chi2" ,"Value of the Chi2  ",_Chi2 ,_Chi2);
    _Chi2Rate = 100.0;
    registerProcessorParameter("Chi2" ,"Value of the Chi2 for Rate estimation  ",_Chi2Rate ,_Chi2Rate);
    _NbrHitPerPlaneMax = 6;
    registerProcessorParameter("NbrHitPerPlaneMax" ,"Maximal number of Hit in each Plane (<=6 by default)  ",_NbrHitPerPlaneMax ,_NbrHitPerPlaneMax);
    _NbrPlaneUseForTracking = 3;
    registerProcessorParameter("NbrPlaneUseForTracking" ,"Number minimal of PLanes used for tracking (>=3 by default)",_NbrPlaneUseForTracking,_NbrPlaneUseForTracking);
    _NbrPlaneUseForTrackingRate = 3;
    registerProcessorParameter("NbrPlaneUseForTrackingRate" ,"Number minimal of PLanes used for rate estimation (>=3 by default)",_NbrPlaneUseForTrackingRate,_NbrPlaneUseForTrackingRate);
    _dlimforPad=40;
    registerProcessorParameter("dlimforPad" ,"dlim for Pad ",_dlimforPad,_dlimforPad);
    _dlimforStrip=40;
    registerProcessorParameter("dlimforStrip" ,"dlim for Strip ",_dlimforStrip,_dlimforStrip);
    _Delimiters="";
    registerProcessorParameter("Delimiters" ,"Delimiters",_Delimiters,_Delimiters);
    _ShortEfficiency=0;
    registerProcessorParameter("ShortEfficiency" ,"ShortEfficiency",_ShortEfficiency,_ShortEfficiency);
    _Config_xml="";
    registerProcessorParameter("Config_xml" ,"Config_xml",_Config_xml,_Config_xml);
    EstimateNoiseContamination=true;
    registerProcessorParameter("NoiseContaminationEstimation" ,"Estimation Noise Contamination",EstimateNoiseContamination,EstimateNoiseContamination);
    _Elog_xml="";
    registerProcessorParameter("Elog_xml" ,"Elog xml",_Elog_xml,_Elog_xml);
}

AnalysisProcessor::~AnalysisProcessor() {}
void AnalysisProcessor::init()
{
    Intro();
    _maxRecord= Global::parameters->getIntVal("MaxRecordNumber")-1;
    _skip= Global::parameters->getIntVal("SkipNEvents");
    if(EstimateNoiseContamination==true)
    {
      names.push_back("NOISE_ESTIMATION_BEFORE");
      names.push_back("NOISE_ESTIMATION_AFTER");
    }
    std::vector<std::string>LCIOFiles;
    Global::parameters->getStringVals("LCIOInputFiles" ,LCIOFiles );

   for(unsigned int i=0;i!=LCIOFiles.size();++i)
   {
    LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
    std::vector<std::string>colee{"SDHCAL_HIT"};
    lcReader->open( LCIOFiles[i] ) ;
    lcReader->setReadCollectionNames( colee ) ;
    _GlobalEvents+=lcReader->getNumberOfEvents()-1;
    delete lcReader;
   }
   std::cout<<_GlobalEventsSc<<"  "<<_GlobalEvents<<std::endl;
    printParameters();
    Reader* myReader = readerFactory.CreateReader(_ReaderType);

    if(myReader)
    {
      myReader->Read(_FileNameGeometry,geom);
      geom.PrintGeom();
      std::map<int, Dif > Difs=geom.GetDifs();
      std::map<int,int> PlansType;
      for(std::map<int, Dif >::iterator it=Difs.begin(); it!=Difs.end(); ++it)
      {
        if(geom.GetDifType(it->first)==scintillator) names.push_back("SDHCAL_HIT_SC");
        if(geom.GetDifType(it->first)!=temporal&&geom.GetDifType(it->first)!=scintillator&&geom.GetDifType(it->first)!=tcherenkov)
        {
          PlansType.insert(std::pair<int,int>(geom.GetDifNbrPlate(it->first)-1,geom.GetDifType(it->first)));
        }
      }
      FillDelimiter(_Delimiters,PlansType.size(),Delimiter);
      int Xmax=geom.GetXMax();
      int Ymax=geom.GetYMax();
      for(std::map<int, int >::iterator it=PlansType.begin(); it!=PlansType.end(); ++it)
	    {
	      geometryplans[it->first]=geometryplan(it->first,it->second,geom,Delimiter);
	      testedPlanList.push_back(geometryplans[it->first]);
        std::string a="Distribution_hit_selectionner_par_analysis"+ std::to_string(it->first +1 );
        std::string b="Hit_expected_"+ std::to_string(it->first +1 );
        std::string AllTheHits="All_The_Hits_From_Trivent"+ std::to_string(it->first +1 );
        std::string c="Efficiency_of_the_voisinage_of_the_pad"+ std::to_string(it->first +1 );
        std::string d="Efficiency_of_the Asic_"+ std::to_string(it->first +1 );
        std::string e="Distance_to_the_expected_hit_in_X_"+ std::to_string(it->first +1 );
        std::string f="Distance_to_the_expected_hit_in_Y_"+ std::to_string(it->first +1 );
        std::string g="Multiplicity_of_the_voisinage_of_the_pad_"+ std::to_string(it->first +1 );
        std::string h="Threshold_"+ std::to_string(it->first +1 );
        std::string hh="gain_"+ std::to_string(it->first +1 );
        
        std::string estrate="Estimate_Rate_"+std::to_string(it->first +1 )+"_nbrplanetouched_";
        std::string estrate2="Estimate_Rate_"+std::to_string(it->first +1 )+"_nbrplanetouched_min";
        std::string nbrviewedd="Number_tracks_viewed_"+std::to_string(it->first +1 );
        std::string efff="sum_effi_"+std::to_string(it->first +1 );
        std::string effsq="sum_effi_squarred_"+std::to_string(it->first +1 );
        std::string estimateurrmin="estimateur_min_"+std::to_string(it->first +1 );
        std::string estimateurr="estimateur_"+std::to_string(it->first +1 );
        std::string estimateurrr="estimateur_suposed_more_robust_"+std::to_string(it->first +1 );
        std::string rateseen="Rate_Seen"+std::to_string(it->first +1 );
        std::string rateestim="Rate_estim"+std::to_string(it->first +1 );
        int X=geom.GetSizeX(it->first)+1;
        int Y=geom.GetSizeY(it->first)+1;
        std::string ndifxy="difxy_"+std::to_string(it->first +1 );
        std::string ndifr="difr_"+std::to_string( it->first +1 );
        estimateurmax[it->first]=new TH1D(estimateurrmin.c_str(),estimateurrmin.c_str(),PlansType.size()+1-_NbrPlaneUseForTrackingRate,double(_NbrPlaneUseForTrackingRate),double(PlansType.size()+1));
        estimateur[it->first]=new TH1D(estimateurr.c_str(),estimateurr.c_str(),PlansType.size()+1-_NbrPlaneUseForTrackingRate,double(_NbrPlaneUseForTrackingRate),double(PlansType.size()+1));
        estimateur2[it->first]=new TH1D(estimateurrr.c_str(),estimateurrr.c_str(),PlansType.size()+1-_NbrPlaneUseForTrackingRate,double(_NbrPlaneUseForTrackingRate),double(PlansType.size()+1));
        RateSeen[it->first]=new TH2F(rateseen.c_str(),rateseen.c_str(),X,0,X,Y,0,Y);
        RateSeenPon[it->first]=new TH2F(rateestim.c_str(),rateestim.c_str(),X,0,X,Y,0,Y);
        for(unsigned int Nbrplatetouched=_NbrPlaneUseForTrackingRate;Nbrplatetouched<=PlansType.size();++Nbrplatetouched)
        {
          std::string nbrplan=std::to_string(Nbrplatetouched);
          std::string name1=estrate+"_"+nbrplan;
          std::string name2=estrate2+"_"+nbrplan;
          std::string name3=nbrviewedd+"_"+nbrplan;
          std::string name4=efff+"_"+nbrplan;
          std::string name5=effsq+"_"+nbrplan;
          RealRateWithSelectedZone[it->first][Nbrplatetouched]=new TH2F(name1.c_str(),name1.c_str(),X,0,X,Y,0,Y);
          
          RealRateWithSelectedZonemin[it->first][Nbrplatetouched]=new TH2F(name2.c_str(),name2.c_str(),X,0,X,Y,0,Y);
          Nbrviewedtracks[it->first][Nbrplatetouched]=new TH2F(name3.c_str(),name3.c_str(),X,0,X,Y,0,Y);
          eff[it->first][Nbrplatetouched]=new TH2F(name4.c_str(),name4.c_str(),X,0,X,Y,0,Y);
          effsquared[it->first][Nbrplatetouched]=new TH2F(name5.c_str(),name5.c_str(),X,0,X,Y,0,Y);
        }
        
        int bin[3]={Xmax,Ymax,int(PlansType.size())*2};
        double xmin[3]={0.,0.,0.};
        double xmax[3]={double(Xmax),double(Ymax),double(PlansType.size())};
        hss2= new THnSparseD("Tracks", "Tracks", 3, bin, xmin, xmax);
        hss= new THnSparseD("TracksPon", "TracksPon", 3, bin, xmin, xmax);
        Gain[names[0]].push_back(new TH2F(hh.c_str(),hh.c_str(),X,0,X,Y,0,Y));
        if(_Config_xml!=""||_Elog_xml!="")
        {
          for(unsigned int j=0;j<ThresholdMap.size();++j)
          {
            ThresholdMap[j]["SDHCAL_HIT"].push_back(new TH2F((h+"_"+Thresholds_name[j]).c_str(),(h+"_"+Thresholds_name[j]).c_str(),X,0,X,Y,0,Y));
          }
        }
        for(unsigned int i=0;i<names.size();++i)
        {
          Distribution_exp_tgraph[names[i]]=new TGraph2D();
          Distribution_hits_tgraph[names[i]]=new TGraph2D();
          RateSeen3D[names[i]]=new TGraph2D();
          if(it->second==positional)
          {
            difxy[names[i]].push_back(new TH2D((ndifxy+names[i]).c_str(),(ndifxy+names[i]).c_str(),346,-50,50,1000,-50,50));
            difr[names[i]].push_back(new TH1D((ndifr+names[i]).c_str(),(ndifr+names[i]).c_str(),346,-50,50));
            Distribution_hits[names[i]].push_back(new TH2F((a+names[i]).c_str(),(a+names[i]).c_str(),X,0,X,Y,0,Y));
            Distribution_exp[names[i]].push_back(new TH2F((b+names[i]).c_str(),(b+names[i]).c_str(),X,0,X,Y,0,Y));
            for(unsigned int j=0;j<Efficiency_pads.size();++j)
            {
              All_The_Hits_From_Trivent[j][names[i]].push_back(new TH2F((AllTheHits+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(AllTheHits+names[i]).c_str(),X,0,X,Y,0,Y));
              Efficiency_pads[j][names[i]].push_back(new TH2F((c+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(c+names[i]).c_str(),X,0,X,Y,0,Y));
              Efficiency_asics[j][names[i]].push_back(new TH2F((d+names[i]+Thresholds_name[j]).c_str(),(d+names[i]).c_str(),X,0,X,Y,0,Y));
              Multiplicity_pads[j][names[i]].push_back(new TH2F((g+names[i]+Thresholds_name[j]).c_str(),(g+names[i]).c_str(),X,0,X,Y,0,Y));
            }

		        HowLongFromExpectedX[names[i]].push_back(new TH1F((e+names[i]).c_str(),(e+names[i]).c_str(),2*(_dlimforStrip),-_dlimforStrip,_dlimforStrip));
            HowLongFromExpectedY[names[i]].push_back(new TH1F((f+names[i]).c_str(),(f+names[i]).c_str(),2*(_dlimforStrip),-_dlimforStrip,_dlimforStrip));
          }
          else
          {
            difxy[names[i]].push_back(new TH2D((ndifxy+names[i]).c_str(),(ndifxy+names[i]).c_str(),346,-50,50,1000,-50,50));
            difr[names[i]].push_back(new TH1D((ndifr+names[i]).c_str(),(ndifr+names[i]).c_str(),346,-50,50));
            Distribution_hits[names[i]].push_back(new TH2F((a+names[i]).c_str(),(a+names[i]).c_str(),X,0,X,Y,0,Y));
            Distribution_exp[names[i]].push_back(new TH2F((b+names[i]).c_str(),(b+names[i]).c_str(),X,0,X,Y,0,Y));
            for(unsigned int j=0;j<Efficiency_pads.size();++j)
            {
            //All_The_Hits_From_Trivent[j][names[i]].push_back(new TH2F((AllTheHits+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(AllTheHits+names[i]).c_str(),X,0,X,Y,0,Y));
            Efficiency_pads[j][names[i]].push_back(new TH2F((c+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(c+names[i]).c_str(),X,0,X,Y,0,Y));
            Efficiency_asics[j][names[i]].push_back(new TH2F((d+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(d+names[i]).c_str(),X,0,X,Y,0,Y));
            Multiplicity_pads[j][names[i]].push_back(new TH2F((g+"_"+names[i]+"_"+Thresholds_name[j]).c_str(),(g+names[i]).c_str(),X,0,X,Y,0,Y));
            }
            HowLongFromExpectedX[names[i]].push_back(new TH1F((e+names[i]).c_str(),(e+names[i]).c_str(),2*(_dlimforPad),-_dlimforPad,_dlimforPad));
            HowLongFromExpectedY[names[i]].push_back(new TH1F((f+names[i]).c_str(),(f+names[i]).c_str(),2*(_dlimforPad),-_dlimforPad,_dlimforPad));

          }
        }
      }
      if(_ShortEfficiency!=0)
	    {
		    for(unsigned int i=0;i<testedPlanList.size();++i)
		    {
			    std::string name="Short_Efficiency"+patch::to_string(i+1);
			    std::string namee="Short_Multiplicity"+patch::to_string(i+1);
			    for(unsigned int i=0;i<names.size();++i)
          {
            name+=names[i];
            namee+=names[i];
            std::vector<TH1F*>vec;
            std::vector<TH1F*>vec2;
            for(unsigned int j=0;j!=Thresholds_name.size();++j)
            {
			        vec.push_back(new TH1F((name+Thresholds_name[j]).c_str(),(name+Thresholds_name[j]).c_str(),1000,0,_ShortEfficiency*1000));
			        vec2.push_back(new TH1F((namee+Thresholds_name[j]).c_str(),(namee+Thresholds_name[j]).c_str(),1000,0,_ShortEfficiency*1000));
			      }
			      Short_Efficiency[names[i]].push_back(vec);
			      Short_Multiplicity[names[i]].push_back(vec2);
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
   

}

void AnalysisProcessor::processEvent( LCEvent * evtP )
{
     
  _NbrRun=evtP->getRunNumber();

  if(isFirstEvent()==true)
    { 
	if(_Elog_xml!="")
    	{
      		Reader* Elog =readerFactory.CreateReader("XMLReaderElog");
      		if(Elog)
      		{
        		Elog->Read(_Elog_xml,conf,_NbrRun);
      		}
      		delete Elog;
    	}
	else
	{
		if(_Config_xml!="")
		{
			Reader* Conf =readerFactory.CreateReader("XMLReaderConfig");
      			if(Conf)
      			{
        			Conf->Read(_Config_xml,conf);
      			}
			delete Conf;
		}
        }
	std::map<unsigned int,DifInfo>ggg;
      	ggg=conf.ReturnMe();
    	for(std::map<unsigned int,DifInfo>::iterator it=ggg.begin();it!=ggg.end();++it)
    	{
      		int dif_id=it->first;
		std::cout<<dif_id<<std::endl;
      		std::map<unsigned int,AsicInfo> ppp=(it->second).ReturnMe();
      		for(std::map<unsigned int,AsicInfo>::iterator itt=ppp.begin();itt!=ppp.end();++itt)
      		{
        		int asic_id=itt->first;
        		std::vector<unsigned int> thee= (itt->second).getThresholds();
         		std::array<unsigned int,64>ooo=(itt->second).ReturnMe();
         		for(unsigned int i=0;i!=ooo.size();++i)
         		{
				if(geom.GetDifType(it->first)!=temporal&&geom.GetDifType(it->first)!=scintillator&&geom.GetDifType(it->first)!=tcherenkov)
				{
           			if(geom.GetDifNbrPlate(it->first)-1>=0)Gain[names[0]][geom.GetDifNbrPlate(it->first)-1]->Fill((1+MapILargeHR2[i]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id),(32-(MapJLargeHR2[i]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id),ooo[i]);
           			for(unsigned int hh=0;hh!=thee.size();++hh)
           				if(geom.GetDifNbrPlate(it->first)-1>=0)ThresholdMap[hh][names[0]][geom.GetDifNbrPlate(it->first)-1]->Fill((1+MapILargeHR2[i]+AsicShiftI[asic_id])+geom.GetDifPositionX(dif_id),(32-(MapJLargeHR2[i]+AsicShiftJ[asic_id]))+geom.GetDifPositionY(dif_id),thee[hh]);}
         		}
      		}
    	}
    }

  Planss.clear();

  if (evtP == nullptr) 
    {
      PrintStatShort();  // utile à faire ici ???    
      return;
    }

  for (std::vector<std::string>::iterator itcol=_hcalCollections.begin(); itcol!= _hcalCollections.end(); ++itcol)
    {
      std::string currentCollectionName=*itcol;
      LCCollection* currentCollection=nullptr;
      try { currentCollection= evtP ->getCollection(currentCollectionName); }
      catch (DataNotAvailableException &e) {continue;} 

      Numbers[currentCollectionName]++;
      Progress(_skip,_GlobalEvents,_maxRecord,Numbers[currentCollectionName]);


      CellIDDecoder<CalorimeterHit> cd(currentCollection);
      int numElements = currentCollection->getNumberOfElements();
      for (int ihit=0; ihit < numElements; ++ihit)
	{
	  CalorimeterHit *raw_hit = dynamic_cast<CalorimeterHit*>( currentCollection->getElementAt(ihit)) ;
	  if (raw_hit != nullptr)
	    {
	      int dif_id=cd(raw_hit)["Dif_id"];
	      //int I=cd(raw_hit)["I"];
	      //int J=cd(raw_hit)["J"];
	      RealNumberPlane[dif_id]++;
	      Planss[currentCollectionName][geom.GetDifNbrPlate(dif_id)-1].addHit(raw_hit);
	      Planss[currentCollectionName][geom.GetDifNbrPlate(dif_id)-1].SetType(geom.GetDifType(dif_id));
	    } // end if (raw_hit != nullptr)
	} // end for (int ihit=0; ihit < numElements; ++ihit)
    } // for(std::map<std::string,LCCollection*>::iterator it=Collections.begin(); it!=Collections.end(); ++it)

  // NB si Planss est vide, Tracks ne fait rien.
  Tracks(Planss,geom,geometryplans,useforrealrate);
  // NB si Planss est vide, testedPlan::testYou ne fait rien
  for (std::vector<testedPlan>::iterator iter=testedPlanList.begin(); iter != testedPlanList.end(); ++iter) iter->testYou(Planss,testedPlanList);
  PrintStatShort();     
 
}

void AnalysisProcessor::processRunHeader( LCRunHeader* run)
{
    LCTOOLS::dumpRunHeader(run);
}
void AnalysisProcessor::end()
{

    for(std::map<std::string,std::array<std::map<std::vector<int>,std::vector<int>>,6>>::iterator itt=Efficiency_per_pad.begin();itt!=Efficiency_per_pad.end();++itt)
    {
      for(unsigned int i=0;i!=(itt->second).size();++i)
      {
        for(std::map<std::vector<int>,std::vector<int>>::iterator it=Efficiency_per_pad[itt->first][i].begin();it!=Efficiency_per_pad[itt->first][i].end();++it)
        {
          unsigned int was_at_least_a_hit=0;
          unsigned int number_of_hits=0;
          for(unsigned int j =0;j!=(it->second).size();++j)
	        {
            if(it->second[j]>0)
            {
			        was_at_least_a_hit++;
			        number_of_hits+=it->second[j];
		        }
	        }
	        std::cout<<red<<i<<"  "<<itt->first<<"  "<<it->first[2]-1<<"  "<<it->first[0]<<"  "<<it->first[1]<<normal<<std::endl;
          Efficiency_pads[i][itt->first][it->first[2]-1]->Fill(it->first[0],it->first[1],was_at_least_a_hit*1.0/(it->second).size());
          if(was_at_least_a_hit!=0)Multiplicity_pads[i][itt->first][it->first[2]-1]->Fill(it->first[0],it->first[1],number_of_hits*1.0/was_at_least_a_hit);
      }
      }
    }
    
    
    
    for(unsigned int i=0;i!=_hcalCollections.size();++i)
    {
      std::cout<<magenta<<_hcalCollections[i]<<normal<<std::endl;
      static std::map<std::string,std::string>Mess;
      for(unsigned int i=0;i!=_hcalCollections.size();++i)Mess[_hcalCollections[i]]="List of counters "+ _hcalCollections[i]+" : ";
      std::cout<<white<<Mess[_hcalCollections[i]]<<normal<<std::endl;
      for (std::vector<testedPlan>::iterator it=testedPlanList.begin(); it != testedPlanList.end(); ++it) it->print(_hcalCollections[i]);
    }
    std::string b="Results_"+std::to_string( (long long int) _NbrRun)+".root";
    std::string c="Results_Analysis"+std::to_string( (long long int) _NbrRun)+".root";
    TFile *hfile = new TFile(c.c_str(),"RECREATE");
    TFile *hfile2 = new TFile(b.c_str(),"UPDATE");
    TH1F* ugly = (TH1F*) hfile2->Get("ugly");
    double total_time;
    if(ugly!=nullptr)total_time=ugly->GetBinContent(1);
    else total_time=1;
    std::cout<<red<<"hhhfhfhdhdhdhdh "<<total_time<<" kjjffjfj"<<normal<<std::endl; 
    hfile2->Close();
    
 
    hfile->cd("../");
    /*for(std::map<std::string,std::vector<std::vector<std::array<int,4>>>>::iterator it=useforrealrate.begin();it!=useforrealrate.end();++it)
    {*/
     for(unsigned int i=0;i!=useforrealrate/*[it->first]*/.size();++i)
     {
       double inbin=1.0;int II=0;int JJ=0;int ZZ=0;/*double x=0;double y=0;double z=0;*/unsigned int number_touched=0;
      for(unsigned int j=0;j!=useforrealrate/*[it->first]*/[i].size();++j)
      {
        static int numberrrr=0;
        
        II=int(useforrealrate/*[it->first]*/[i][j][0]);
        JJ=int(useforrealrate/*[it->first]*/[i][j][1]);
        ZZ=int(useforrealrate/*[it->first]*/[i][j][2]);
        //x=useforrealrate/*[it->first]*/[i][j][4];
        //y=useforrealrate/*[it->first]*/[i][j][5];
        //z=useforrealrate/*[it->first]*/[i][j][6];
        
        RateSeen3D["SDHCAL_HIT"]->SetPoint(numberrrr,II,JJ,ZZ);
        numberrrr++;
        int xmax=Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetNbinsX();
        int ymax=Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetNbinsY();
        if(int(useforrealrate/*[it->first]*/[i][j][3])==1)
        {
          number_touched++;
          
          std::cout<<II<<"  "<<JJ<<"  "<<ZZ<<std::endl;
         
          if(II>0&&II<xmax&&JJ>0&&JJ<ymax)
          {
             
            if(Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetBinContent(II,JJ)>0)inbin*=Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetBinContent(II,JJ);
          }
        }
        else
        {
          if(II>0&&II<xmax&&JJ>0&&JJ<ymax)
          {
            if(Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetBinContent(II,JJ)>=0&&Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetBinContent(II,JJ)!=1)inbin*=1-Efficiency_pads[5]["SDHCAL_HIT"][ZZ]->GetBinContent(II,JJ);
          }
        }
        
      }
      for(unsigned int j=0;j!=useforrealrate/*[it->first]*/[i].size();++j)
      {      
        II=int(useforrealrate/*[it->first]*/[i][j][0]);
        JJ=int(useforrealrate/*[it->first]*/[i][j][1]);
        ZZ=int(useforrealrate/*[it->first]*/[i][j][2]);
        double fillr[3]={double(II),double(JJ),double(ZZ)};
        hss2->Fill(fillr,1);
        hss->Fill(fillr,1.0/inbin);
        for(std::map<int,TH2F*>::iterator it=RealRateWithSelectedZonemin[ZZ].begin();it!=RealRateWithSelectedZonemin[ZZ].end();++it)
        {
          if((int)number_touched>=it->first)RealRateWithSelectedZonemin[ZZ][it->first]->Fill(II,JJ,1.0/inbin);
          RateSeen[ZZ]->Fill(II,JJ,1.0);
          RateSeenPon[ZZ]->Fill(II,JJ,1.0/inbin);
        }
        if((int)number_touched>=_NbrPlaneUseForTrackingRate)
        {
          RealRateWithSelectedZone[ZZ][number_touched]->Fill(II,JJ,1.0/inbin);
          Nbrviewedtracks[ZZ][number_touched]->Fill(II,JJ,1);
          eff[ZZ][number_touched]->Fill(II,JJ,inbin);
          effsquared[ZZ][number_touched]->Fill(II,JJ,inbin*inbin);
        }
      }
     }
    /*}*/
    std::string namei="";
    for(std::map<int,std::map<int,TH2F*>>::iterator it=Nbrviewedtracks.begin();it!=Nbrviewedtracks.end();++it)
    {
      for(std::map<int,TH2F*>::iterator itt=Nbrviewedtracks[it->first].begin();itt!=Nbrviewedtracks[it->first].end();++itt)
      {
         namei="rate_estimated_"+std::to_string(it->first)+"_"+std::to_string(itt->first)+"_planes_touched";
         estimatedrate[it->first][itt->first]=dynamic_cast<TH2F*>( (itt->second)->Clone(namei.c_str()) ) 		;
         estimatedrate[it->first][itt->first]->Multiply(eff[it->first][itt->first]);
         estimatedrate[it->first][itt->first]->Divide(effsquared[it->first][itt->first]);
      }
    }



    TF1 * tf = new TF1("TransferFunction", transfer_function);
    tt->Write();
    hss2->Write();
    hss->Write();
    TH3D* h2= hss2->Projection(2,1,0);
    TH3D* h1= hss->Projection(2,1,0);
    Make_good_TH3(h2);
    Make_good_TH3(h1);
    h2->GetListOfFunctions()->Add(tf);
    h2->Write();
    h1->GetListOfFunctions()->Add(tf);
    h1->Write();
    RateSeen3D["SDHCAL_HIT"]->Write("Rate_Seen_3D");
    std::string traces="Traces";
    Rate_Estimation->Fill (0.5,useforrealrate.size());
    Rate_Estimation->Fill(0.5,numbertracks);
    Rate_Estimation->Scale(1.0/total_time);
    Rate_Estimation->Write();
      hfile->mkdir(traces.c_str());
      hfile->cd(traces.c_str());
      
      for(unsigned int i=0;i!=xzaxis.size();++i)
      {
        xzaxis[i].Write((traces+"_xz_"+patch::to_string(i)).c_str());
        yzaxis[i].Write((traces+"_yz_"+patch::to_string(i)).c_str());
      }
      hfile->cd("../");
    
    std::string SlowControl="";
      if(_Config_xml!=""||_Elog_xml!="")
     {
      for(unsigned int o=0;o!=Gain[names[0]].size();++o)
      {
        SlowControl="SlowControl_"+patch::to_string(o+1);
        hfile->mkdir(SlowControl.c_str(),SlowControl.c_str());
        hfile->cd(SlowControl.c_str());
        Gain[names[0]][o]->Write();
        for(unsigned int j=0;j<ThresholdMap.size();++j)
        {
          ThresholdMap[j][names[0]][o]->Write();
        }
        hfile->cd("../");
      }
    }
    for(unsigned int naa=0;naa<names.size();++naa)
    {
      
      Distribution_exp_tgraph[names[naa]]->Write(("Distribution_exp"+names[naa]).c_str());
      Distribution_hits_tgraph[names[naa]]->Write(("Distribution_hits"+names[naa]).c_str());
      
     std::string plate="";
    
      for(unsigned int i =0 ;i!=HowLongFromExpectedX["SDHCAL_HIT"].size();++i)
      {
        
        plate="Plate "+ patch::to_string(i+1);

       
        if(naa==0)hfile->mkdir(plate.c_str(),plate.c_str());
    	  hfile->cd(plate.c_str());
    	  if(naa==0)hfile->mkdir((plate+"/Alignement").c_str());
        hfile->cd((plate+"/Alignement").c_str());
        difr[names[naa]][i]->Write();
    	  difxy[names[naa]][i]->Write();
    	  difxy[names[naa]][i]->ProjectionY()->Write();
    	  difxy[names[naa]][i]->ProjectionX()->Write();
        hfile->cd(plate.c_str());   	
    	  if(naa==0)
        {
        
     
        
        hfile->mkdir((plate+"/Rate_Estimation").c_str());
        hfile->cd((plate+"/Rate_Estimation").c_str());
        std::string name="Rate_estimated_nbr_plane_used_";
        std::string name2="Rate_estimated_nbr_plane_used_min_";
        //EfficacityVsRate[i]->Write();
        RateSeen[i]->Scale(1.0/total_time);
        RateSeen[i]->Write();
        RateSeenPon[i]->Scale(1.0/total_time);
        RateSeenPon[i]->Write();
        for(std::map<int,TH2F*>::iterator it=RealRateWithSelectedZonemin[i].begin();it!=RealRateWithSelectedZonemin[i].end();++it)
        {
            RealRateWithSelectedZonemin[i][it->first]->Write();
            RealRateWithSelectedZonemin[i][it->first]->Scale(1.0/(total_time*SumCombinaison(RealNumberPlane.size(),it->first)));
            estimateur[i]->Fill(it->first,RealRateWithSelectedZonemin[i][it->first]->Integral());
            RealRateWithSelectedZonemin[i][it->first]->Write((name2+std::to_string(it->first)).c_str());
            
        }
          for(std::map<int,TH2F*>::iterator it=RealRateWithSelectedZone[i].begin();it!=RealRateWithSelectedZone[i].end();++it)
        {
            RealRateWithSelectedZone[i][it->first]->Write();
            RealRateWithSelectedZone[i][it->first]->Scale(1.0/(total_time*TMath::Binomial(RealNumberPlane.size(),it->first)));
            RealRateWithSelectedZone[i][it->first]->Write((name+std::to_string(it->first)).c_str());
            estimateurmax[i]->Fill(it->first,RealRateWithSelectedZone[i][it->first]->Integral());
            Nbrviewedtracks[i][it->first]->Write();
            eff[i][it->first]->Write();
            effsquared[i][it->first]->Write();
            
            estimatedrate[i][it->first]->Scale(1.0/(total_time*TMath::Binomial(RealNumberPlane.size(),it->first)));
            estimatedrate[i][it->first]->Write();
            estimateur2[i]->Fill(it->first,estimatedrate[i][it->first]->Integral());
        }
        estimateur[i]->Write();
        estimateur2[i]->Write();
        estimateurmax[i]->Write();
        }
        hfile->cd(plate.c_str());
	      HowLongFromExpectedX[names[naa]][i]->Write();
    	  HowLongFromExpectedY[names[naa]][i]->Write();
        Distribution_hits[names[naa]][i]->Write();

        if(_ShortEfficiency!=0) 
        { 
          if(naa==0)hfile->mkdir((plate+"/Short_Efficiency").c_str());
          hfile->cd((plate+"/Short_Efficiency").c_str());
          for(unsigned int j=0;j!=Thresholds_name.size();++j)
          {
            Short_Efficiency[names[naa]][i][j]->Write();
            Short_Multiplicity[names[naa]][i][j]->Write();
          }
          hfile->cd(plate.c_str());
        }
        /*for(unsigned int j=0;j<All_The_Hits_From_Trivent.size();++j)
        {
          All_The_Hits_From_Trivent[j][names[naa]][i]->Scale(1.0/total_time);
          All_The_Hits_From_Trivent[j][names[naa]][i]->Write();
        }*/
        Distribution_exp[names[naa]][i]->Write();
        if(naa==0)hfile->mkdir((plate+"/Efficiency_Multiplicity_Maps").c_str());
        hfile->cd((plate+"/Efficiency_Multiplicity_Maps").c_str());
        for(unsigned int j=0;j<Efficiency_pads.size();++j)
        {
        
        Efficiency_pads[j][names[naa]][i]->Write();
        Multiplicity_pads[j][names[naa]][i]->Write();
        
        }
        hfile->cd(plate.c_str());
        hfile->cd(plate.c_str());
      }

      
        for(unsigned int i=0; i<Distribution_hits.size(); ++i)
        {
          delete Distribution_hits[names[naa]][i];
          delete Distribution_exp[names[naa]][i];
          delete HowLongFromExpectedX[names[naa]][i];
          delete HowLongFromExpectedY[names[naa]][i];
          for(unsigned int j=0;j<Efficiency_pads.size();++j)
          {
            delete Efficiency_pads[j][names[naa]][i];
            delete Multiplicity_pads[j][names[naa]][i];
          }
          if(_ShortEfficiency!=0) 
          {
            for(unsigned int j=0;j!=Thresholds_name.size();++j)
            {
              delete Short_Efficiency[names[naa]][i][j];
              delete Short_Multiplicity[names[naa]][i][j];
            }
          } 
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
    for(unsigned int f=0;f!=_hcalCollections.size();++f)
    {
      PrintStat(_hcalCollections[f]);    
   }
}
