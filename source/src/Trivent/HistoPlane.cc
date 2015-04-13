#include "Trivent/HistoPlane.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include<string>
#include<vector>
#include<map>
#include<iostream>
#include "Colors.h"
#include "Patch.h"
#include "TROOT.h"

HistoPlane::HistoPlane(int NbrPlate,int SizeX, int SizeY, std::vector< std::string  >& vec_name_th1,std::vector< std::string >& vec_name_th2,std::vector< std::string >& vec_name_th2_Asic):NbrPlatee(NbrPlate),Means(0),Nbrof0Hits(0),local_max(-1),local_min(99999999),total_time(0),_SizeX(SizeX),_SizeY(SizeY)
{  
        std::string addnbr;
	for(std::vector< std::basic_string<char>  >::iterator it=vec_name_th1.begin();it!=vec_name_th1.end();++it)
	{
                addnbr=(*it)+patch::to_string(NbrPlate);
                if (gROOT->FindObject(addnbr.c_str()) != NULL) continue;
  		TH1Fs.insert(std::pair<std::string,TH1F*>((*it),new TH1F(addnbr.c_str(),(*it).c_str(),25000,0,25000)));
	}
	for(std::vector< std::basic_string<char>  >::iterator it=vec_name_th2.begin();it!=vec_name_th2.end();++it)
	{
                addnbr=(*it)+patch::to_string(NbrPlate);
                 if (gROOT->FindObject(addnbr.c_str()) != NULL) continue;
  		TH2Fs.insert(std::pair<std::string,TH2F*>((*it),new TH2F(addnbr.c_str(),(*it).c_str(),(int)SizeX+1,0,(int)SizeX+1,(int)SizeY+1,0,(int)SizeY+1)));
	}
        for(std::vector< std::basic_string<char>  >::iterator it=vec_name_th2_Asic.begin();it!=vec_name_th2_Asic.end();++it)
	{
                addnbr=(*it)+patch::to_string(NbrPlate);
                if (gROOT->FindObject(addnbr.c_str()) != NULL) continue;
  		TH2Fs.insert(std::pair<std::string,TH2F*>((*it),new TH2F(addnbr.c_str(),(*it).c_str(),(int)SizeX/8,0,(int)SizeX/8,(int)SizeY/8,1,(int)SizeY/8+1)));
	}
        if (TH1Fs.size()!=0&&TH2Fs.size()!=0) std::cout<<red<<"Creating "<<TH1Fs.size()<<" TH1 and "<<TH2Fs.size()<<" TH2F for plate "<<patch::to_string(NbrPlate+1)<<normal<<std::endl;
}

HistoPlane::~HistoPlane()
{
	/*for(std::map<std::string,TH1F*>::iterator it=TH1Fs.begin();it!=TH1Fs.end();++it)
	{
  		
		std::cout<<it->first<<std::endl;
                delete it->second;          
                TH1Fs.erase(it->first);
                std::cout<<green<<TH1Fs.size()<<normal<<std::endl;
  		
	}
	for(std::map<std::string,TH2F*>::iterator it=TH2Fs.begin();it!=TH2Fs.end();++it)
	{
  		
		std::cout<<red<<it->first<<normal<<std::endl;
                delete it->second;
                TH2Fs.erase(it->first);std::cout<<green<<TH2Fs.size()<<normal<<std::endl;
                
  		
	}*/
}

TH1F* HistoPlane::Return_TH1F(const char* name)
{
	if (TH1Fs.find(name)!=TH1Fs.end())return TH1Fs.find(name)->second;
        else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}
TH2F* HistoPlane::Return_TH2F(const char* name)
{
	if (TH2Fs.find(name)!=TH2Fs.end())return TH2Fs.find(name)->second;
        else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}


void HistoPlane::ScaleHisto(const char* name,float i)
{
	if (TH1Fs.find(name)!=TH1Fs.end())(TH1Fs.find(name))->second->Scale(i);
        if (TH2Fs.find(name)!=TH2Fs.end())(TH2Fs.find(name))->second->Scale(i);
}

void HistoPlane::WriteAll()
{
   
  
  for(std::map<std::string,TH1F*>::iterator it=TH1Fs.begin();it!=TH1Fs.end();++it)
  {
  		
  		(it->second)->Write();
  		
  }
  for(std::map<std::string, TH2F*>::iterator it=TH2Fs.begin();it!=TH2Fs.end();++it)
  {
  		
  		(it->second)->Write();
  		
  }
}

void HistoPlane::Save(TFile* file)
{
    std::string plate="Plate "+ std::to_string( (long long int) NbrPlatee+1);
    file->mkdir(plate.c_str(),plate.c_str());
    file->cd(plate.c_str());
    for(std::map<int,int>::iterator it = Time_Plates.begin(); it!=Time_Plates.end(); ++it) 
    {
      TH1Fs["Time_Distr"]->Fill(it->first,it->second);
      TH1Fs["Hits_Distr"]->Fill(it->second,1);
    }
    TH1Fs["Hits_Distr"]->Fill(0.0,(double)Nbrof0Hits);
    
    TH1Fs["Time_Distr"]->GetXaxis()->SetRange(0.0,Time_Plates.size()+10);
    TH1Fs["Hits_Distr"]->GetXaxis()->SetRange(0,150);
    Return_TH1F("Hits_Distr_Noise")->Fill(0.0,(double)Nbrof0Hits);
    WriteAll();
    std::string name="Noise_Flux_Hz";
    TH2Fs["Flux_Noise"]->Scale(1/(total_time*2e-7));
    TH2Fs["Flux_Noise"]->Write(name.c_str());
    name="Flux_Noise_Mean_Scaled";
    Means=TH2Fs["Flux_Noise"]->Integral()/(1.0*GetArea());
    TH2Fs["Flux_Noise"]->Scale(1/Means);
    TH2Fs["Flux_Noise"]->Write(name.c_str());
    name="Asic_Noise_Flux";
    TH2Fs["Flux_Noise_Asic"]->Scale(1/(total_time*2e-7));
    TH2Fs["Flux_Noise_Asic"]->Write(name.c_str());
    name="Asic_Noise_Flux_Mean_Scaled";
    TH2Fs["Flux_Noise_Asic"]->Scale(1/Means);
    TH2Fs["Flux_Noise_Asic"]->Write(name.c_str());
    name="Event_Flux_Hz";
    TH2Fs["Flux_Events"]->Scale(1/(total_time*2e-7));
    TH2Fs["Flux_Events"]->Write(name.c_str());
    name="Asic_Event_Flux";
    TH2Fs["Flux_Events_Asic"]->Scale(1/(total_time*2e-7));
    TH2Fs["Flux_Events_Asic"]->Write(name.c_str());
}
