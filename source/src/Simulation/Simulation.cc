#include<iostream>
#include<vector>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "TMath.h"
#include<map>
#include "TH1F.h"
#include<string>
#include "TFile.h"
#include "Colors.h"
#include "TRandom3.h"


double frand()
{
    return ( rand()/(double)RAND_MAX ) * 100;
}


int SumCombinaison(int n,int kmin)
{
  int sumcomb=0;
  for(unsigned int i=kmin;i<=n;++i) sumcomb+=TMath::Binomial(n,i);
  return sumcomb;
}


int main()
{
  
  srand(time(NULL));
  std::vector<double>Efficiency_Pads;
  unsigned int NbrPlansMin=4;
  unsigned int NbrPlansMax=9;
  unsigned int Nbriteration=1000;
  unsigned int NbrTracks=10000;
  unsigned int NbrPlansPourTraceMin=3;
  TFile *f=nullptr;
  //iteration sur le nombre de plans dans l'analyse
  for(unsigned int NbrPlans=NbrPlansMin;NbrPlans<=NbrPlansMax;++NbrPlans)
  {
    std::cout<<red<<"Cas avec "<<NbrPlans<<"  "<<normal<<std::endl;
    //Creation du fichier pour NbrPlans;
    std::string name="Result_Simulation_for_"+std::to_string(NbrPlans)+"_Planes.root"; 
    f = new TFile(name.c_str(),"recreate");
    //map qui contiendra les TH1;
    std::map<int,std::map<int,TH1F*>>TH1Fs;
    std::map<int,std::map<int,TH1F*>>TH1FsMax;
    //iteration sur l'efficacité des chambres
    for(unsigned int effi=1;effi<=100;++effi)
    {
    
      //map avec nombre de pad toucher = clef
      std::map<int,double>Sum;
      //map avec nombre de pas >=clef
      std::map<int,double>SumMax;
       //map avec nombre de pad toucher = clef
      std::map<int,std::vector<double>>ProdEffi;
      //map avec nombre de pas >=clef
      std::map<int,std::vector<double>>ProdEffiMax;
      //initialisation des maps
      for(unsigned int NbrplanMin=NbrPlansPourTraceMin;NbrplanMin<=NbrPlans;++NbrplanMin)
      {
        
        std::string name_th1="The_"+std::to_string(NbrPlans)+"_with_efficiency_"+std::to_string(effi)+"_trackused_"+std::to_string(NbrplanMin);
        std::cout<<yellow<<name_th1<<normal<<std::endl;
        std::string name_th1_sup="The_"+std::to_string(NbrPlans)+"_with_efficiency_"+std::to_string(effi)+"_trackused_supegal_"+std::to_string(NbrplanMin);
        TH1Fs[effi][NbrplanMin]=new TH1F(name_th1.c_str(),name_th1.c_str(),400,-200,200); 
        TH1FsMax[effi][NbrplanMin]=new TH1F(name_th1_sup.c_str(),name_th1_sup.c_str(),400,-200,200); 
        Sum[NbrplanMin]=0;
        SumMax[NbrplanMin]=0;
        ProdEffi[NbrplanMin]={};
        ProdEffiMax[NbrplanMin]={};
      }
      std::cout<<"Efficacite "<<effi<<std::endl;
      
      //Creation des plans avec leur efficacité 
      Efficiency_Pads.clear();
      for(unsigned int Plans=1;Plans<=NbrPlans;++Plans)
      {
        Efficiency_Pads.push_back(effi);
      }
     
      double nombre_aleatoire = 0;
      for(unsigned int nbriter=1;nbriter<=Nbriteration;++nbriter)
      {
        
        for(unsigned int NbrplanMin=NbrPlansPourTraceMin;NbrplanMin<=NbrPlans;++NbrplanMin)
        {
          Sum[NbrplanMin]=0;
          SumMax[NbrplanMin]=0;
          ProdEffi[NbrplanMin]={};
          ProdEffiMax[NbrplanMin]={};
        }
        //iteration sur le nombre de traces voulues
        std::cout<<"Traces numero ";
        for(unsigned int Track=1;Track<=NbrTracks;++Track)
        {
          if(Track%1000==0)std::cout<<Track<<" ";
          double prod=1;
          double number_plane_touched=0;
          //Generation des traces 
          for(unsigned int Pad=0;Pad!=Efficiency_Pads.size();++Pad)
          {
            //nombre aleatoire en 0 et 100
            nombre_aleatoire=frand();//Rndm()*100 	;
            //
            if(nombre_aleatoire<=Efficiency_Pads[Pad]) 
            {
              prod*=Efficiency_Pads[Pad]*1.0/100;
              number_plane_touched++;
            }
          
            else prod*=1-(1.0*Efficiency_Pads[Pad]/100);
          }
          for(unsigned int NbrplanMin=NbrPlansPourTraceMin;NbrplanMin<=NbrPlans;++NbrplanMin)
          {
            if(number_plane_touched>=NbrplanMin) ProdEffiMax[NbrplanMin].push_back(prod);
            if(number_plane_touched==NbrplanMin) ProdEffi[NbrplanMin].push_back(prod);
           // std::cout<<magenta<<ProdEffiMax[NbrplanMin].size()<<normal<<std::endl;
          }
        } 
        std::cout<<std::endl;
             
      for(std::map<int,std::vector<double>>::iterator it=ProdEffi.begin();it!=ProdEffi.end();++it)
      {
        for(unsigned int nbr=0;nbr!=(it->second).size();++nbr)
        {
          Sum[it->first]+=1.0/(it->second)[nbr];
        }
      }
      for(std::map<int,std::vector<double>>::iterator it=ProdEffiMax.begin();it!=ProdEffiMax.end();++it)
      {
        for(unsigned int nbr=0;nbr!=(it->second).size();++nbr)
        {
          SumMax[it->first]+=1.0/(it->second)[nbr];
        }
      }
      for(std::map<int,double>::iterator it=Sum.begin();it!=Sum.end();++it)
      {
          Sum[it->first]*=1.0/TMath::Binomial(NbrPlans,it->first);
          //std::cout<<Sum[it->first]<<normal<<std::endl;
          TH1Fs[effi][it->first]->Fill((Sum[it->first]-NbrTracks)*100.0/NbrTracks);
      }
      for(std::map<int,double>::iterator it=SumMax.begin();it!=SumMax.end();++it)
      {
          SumMax[it->first]*=1.0/SumCombinaison(NbrPlans,it->first);
          TH1FsMax[effi][it->first]->Fill((SumMax[it->first]-NbrTracks)*100.0/NbrTracks);
          std::cout<<(SumMax[it->first]-NbrTracks)*100.0/NbrTracks<<std::endl;
      }
      }

      
    }
  
  
    //Ecrire TH1F dans le fichier
    for(std::map<int,std::map<int,TH1F*>>::iterator it=TH1Fs.begin();it!=TH1Fs.end();++it)
    {
      for(std::map<int,TH1F*>::iterator itt=TH1Fs[it->first].begin();itt!=TH1Fs[it->first].end();++itt)
      {
        (itt->second)->Write();
        //std::cout<<"ecrit"<<std::endl;
      }
    }
    for(std::map<int,std::map<int,TH1F*>>::iterator it=TH1FsMax.begin();it!=TH1FsMax.end();++it)
    {
      for(std::map<int,TH1F*>::iterator itt=TH1FsMax[it->first].begin();itt!=TH1FsMax[it->first].end();++itt)
      {
        (itt->second)->Write();
      }
    }
    std::cout<<red<<"Fini"<<normal<<std::endl;
    f->Close();
  }
}
 
  
/*  
  unsigned int NbrPlates=7;
  unsigned int iter=100;F
  unsigned int Nbriteration=100;
  for(unsigned int y=1;y!=100;++y)
  {
    std::string name=std::to_string(y);
    t[y]=new TH1F(name.c_str(),name.c_str(),20000,0,2000000);
    for(unsigned int k=0;k!=iter;++k)
    {
     
      Efficiency_Pads.clear();
      for(unsigned int i=1;i!=NbrPlates;++i)
      {
        Efficiency_Pads.push_back(y);
      }
      double product=1.0;
      std::vector<double>ProdEffi;
      int nombre_aleatoire = 0;
      std::map<int,int>Nbr;
    
  
      for(unsigned int i=1;i!=Nbriteration;++i)
      {
        double prod=1;
        double number_plane_touched=0;
        for(unsigned int j=0;j!=Efficiency_Pads.size();++j)
        {
          nombre_aleatoire=frand();
          //std::cout<<"Nombre aleatoire : "<<nombre_aleatoire<<"  "<<std::endl;
          if(nombre_aleatoire<=Efficiency_Pads[j]) 
          {
            prod*=Efficiency_Pads[j]*1.0/100;
            number_plane_touched++;
          }
        
          else prod*=1-(1.0*Efficiency_Pads[j]/100);
        }
        if(number_plane_touched>=3)
        {
          Nbr[number_plane_touched]++;
          ProdEffi.push_back(prod);
          //std::cout<<prod<<std::endl;
        }
      }
      double sum=0;
      for(unsigned int i=0;i!=ProdEffi.size();++i)
      {
        sum+=1.0/ProdEffi[i];
      }
      std::cout<<SumCombinaison(NbrPlates,3)<<std::endl;
      double number_reconstructed=sum/SumCombinaison(NbrPlates,3);
      for(unsigned int i=0;i!=Efficiency_Pads.size();++i) std::cout<<Efficiency_Pads[i]<<"  ";
      std::cout<<std::endl;
      std::cout<<"Reconstructed tracks "<<number_reconstructed<<" Number of tracks real : "<<Nbriteration<<" % "<<std::fabs(number_reconstructed-Nbriteration)*1.0/Nbriteration*100<<std::endl;
      t[y]->Fill(number_reconstructed);
      }
    }
 for(std::map<int,TH1F*>::iterator it=t.begin();it!=t.end();++it) (it->second)->Write();
  f->Close();
}*/
