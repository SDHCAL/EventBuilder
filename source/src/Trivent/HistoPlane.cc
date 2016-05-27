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
#include <array>
#include "TCanvas.h"
#include "TPDF.h"

TCanvas* canvas = new  TCanvas("canvas");
extern bool pdf;

HistoPlane::HistoPlane(bool Distr,int NbrPlate,std::vector<int>Difs_Names,int SizeX, int SizeY, std::vector< std::string  >& vec_name_th1,std::vector< std::string >& vec_name_th2,std::vector< std::string >& vec_name_th2_Asic):_Distr(Distr),NbrPlatee(NbrPlate),Means(0),Nbrof0Hits(0),local_max(-1),local_min(99999999),total_time(0),_SizeX(SizeX),_SizeY(SizeY),_Difs_Names(Difs_Names)
{       
        int Addone=0;
        std::string addnbr="";
        std::string addnbr2="";
        if(Distr==true)
	{
        Addone=1;
        addnbr="Distribution Nbr Hits in Plane : "+patch::to_string(NbrPlatee+1);
        addnbr2="Distribution timestamps in Plane : "+patch::to_string(NbrPlatee+1);
        std::string a = "Distribution Nbr Hits in Plane : "+patch::to_string(NbrPlatee+1);
        std::string a2 = "Distribution timestamps in Plane : "+patch::to_string(NbrPlatee+1);
        if (gROOT->FindObject(addnbr.c_str()) == NULL)Plate_Dist=new TH1D (a.c_str(),a.c_str(),20,0.,20.);
	if (gROOT->FindObject(addnbr2.c_str()) == NULL)Plate_Dist2=new TH1D (a2.c_str(),a2.c_str(),1000,0.,1000.);
        }
        for(int i=0;i<_Difs_Names.size();++i)
	{
		Calibration.insert(std::pair<int,std::array<std::array<double,64>, 48>>(_Difs_Names[i],std::array<std::array<double,64>, 48>{}));
                if(Distr==true)
		{
			addnbr="Distribution Nbr Hits in Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
                        addnbr2="Distribution timestamps in Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
			std::string a = "Distribution Nbr Hits in Dif : "+patch::to_string(_Difs_Names[i]);
                        std::string a2 = "Distribution timestamp in Dif : "+patch::to_string(_Difs_Names[i]);
			if (gROOT->FindObject(addnbr.c_str()) != NULL)continue;
		 	Difs_Distr[_Difs_Names[i]]=new TH1D (addnbr.c_str(),a.c_str(),10,0.,10.);
			if (gROOT->FindObject(addnbr2.c_str()) != NULL)continue;
                        Difs_Distr2[_Difs_Names[i]]=new TH1D (addnbr2.c_str(),a2.c_str(),1000,0.,1000.);
			for(int j=1;j<=48;++j)
			{
				std::vector<int>vec{_Difs_Names[i],j};
                        	addnbr="Distribution Nbr Hits in Asic : "+patch::to_string(j)+"Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
                                addnbr2="Distribution timespamp in Asic : "+patch::to_string(j)+"Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
				std::string b = "Distribution Nbr Hits in Asic : "+patch::to_string(j);
                                std::string b2 = "Distribution timestamps in Asic : "+patch::to_string(j);
				if (gROOT->FindObject(addnbr.c_str()) != NULL)continue;
				Asics_Distr[vec]= new TH1D (addnbr.c_str(),b.c_str(),10,0.,10.);
				if (gROOT->FindObject(addnbr2.c_str()) != NULL)continue;
				Asics_Distr2[vec]= new TH1D (addnbr2.c_str(),b2.c_str(),1000,0.,1000.);
				for(int k=0;k<=63;++k)
		        	{
  					std::vector<int>vec{_Difs_Names[i],j,k};
                        		addnbr="Distribution Nbr Hits in Pad : "+patch::to_string(k)+"Asic : "+patch::to_string(j)+"Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
					addnbr2="Distribution timestamp Hits in Pad : "+patch::to_string(k)+"Asic : "+patch::to_string(j)+"Dif : "+patch::to_string(_Difs_Names[i])+" Plane : "+patch::to_string(NbrPlatee+1);
					std::string c = "Distribution Nbr Hits in Pad : "+patch::to_string(k);
					std::string c2 = "Distribution timestamp Hits in Pad : "+patch::to_string(k);
					if (gROOT->FindObject(addnbr.c_str()) != NULL)continue;
					Pads_Distr[vec]= new TH1D (addnbr.c_str(),c.c_str(),5,0.,5.);
					if (gROOT->FindObject(addnbr2.c_str()) != NULL)continue;
					Pads_Distr2[vec]= new TH1D (addnbr2.c_str(),c2.c_str(),1000,0.,1000.);
				}
			}
		}
	}
	for(std::vector< std::basic_string<char>  >::iterator it=vec_name_th1.begin();it!=vec_name_th1.end();++it)
	{
                addnbr=(*it)+patch::to_string(NbrPlate);
                if (gROOT->FindObject(addnbr.c_str()) != NULL) continue;
  		TH1Fs.insert(std::pair<std::string,TH1D*>((*it),new TH1D(addnbr.c_str(),(*it).c_str(),25000,0,25000)));
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
       // if (TH1Fs.size()!=0&&TH2Fs.size()!=0) std::cout<<red<<"Creating "<<TH1Fs.size()+Addone+Difs_Distr.size()+Asics_Distr.size()+Pads_Distr.size()<<" TH1 and "<<TH2Fs.size()<<" TH2F for plate "<<patch::to_string(NbrPlate+1)<<normal<<std::endl;
}

HistoPlane::HistoPlane(const HistoPlane &source):_Distr(source._Distr),NbrPlatee(source.NbrPlatee),Means(source.Means),Nbrof0Hits(source.Nbrof0Hits),local_max(source.local_max),local_min(source.local_min),total_time(source.total_time),_SizeX(source._SizeX),_SizeY(source._SizeY),_Difs_Names(source._Difs_Names)
{   
        std::cout<<"Create copy"<<std::endl;
        Plate_Dist=new TH1D (*(source.Plate_Dist));
	Plate_Dist2=new TH1D (*(source.Plate_Dist2));
        for(std::map<int,TH1D*>::const_iterator it=(source.Difs_Distr2).begin();it!=(source.Difs_Distr2).end();++it)
	{
		Difs_Distr2[it->first]=new TH1D (*it->second);
        }
	for(std::map<std::vector<int>,TH1D*>::const_iterator it=(source.Asics_Distr2).begin();it!=(source.Asics_Distr2).end();++it)
        {
		Asics_Distr2[it->first]= new TH1D (*it->second);
        }
	for(std::map<std::vector<int>,TH1D*>::const_iterator it=(source.Pads_Distr2).begin();it!=(source.Pads_Distr2).end();++it)
	{
		Pads_Distr2[it->first]= new TH1D (*it->second);
	}
        for(std::map<int,TH1D*>::const_iterator it=(source.Difs_Distr2).begin();it!=(source.Difs_Distr2).end();++it)
	{
		Difs_Distr2[it->first]=new TH1D (*it->second);
        }
	for(std::map<std::vector<int>,TH1D*>::const_iterator it=(source.Asics_Distr).begin();it!=(source.Asics_Distr).end();++it)
        {
		Asics_Distr[it->first]= new TH1D (*it->second);
        }
	for(std::map<std::vector<int>,TH1D*>::const_iterator it=(source.Pads_Distr).begin();it!=(source.Pads_Distr).end();++it)
	{
		Pads_Distr[it->first]= new TH1D (*it->second);
	}
        for(std::map<std::string,TH1D*>::const_iterator it=(source.TH1Fs).begin();it!=(source.TH1Fs).end();++it)
	{
  		TH1Fs.insert(std::pair<std::string,TH1D*>(it->first,new TH1D (*it->second)));
	}
	for(std::map<std::string,TH2F*>::const_iterator it=(source.TH2Fs).begin();it!=(source.TH2Fs).end();++it)
	{
  		TH2Fs.insert(std::pair<std::string,TH2F*>(it->first,new TH2F (*it->second)));
	}
}   



HistoPlane::~HistoPlane()
{
        for(int i=0;i<_Difs_Names.size();++i)
	{
		delete Difs_Distr[_Difs_Names[i]];
		delete Difs_Distr2[_Difs_Names[i]];
		for(int j=1;j<=48;++j)
		{
			std::vector<int>vec{_Difs_Names[i],j};
			delete Asics_Distr[vec];
			delete Asics_Distr2[vec];
			for(int k=0;k<=63;++k)
		        {
  				std::vector<int>vec{_Difs_Names[i],j,k};
				delete Pads_Distr[vec];
				delete Pads_Distr2[vec];
			}
		}
	}
	for(std::map<std::string,TH1D*>::iterator it=TH1Fs.begin();it!=TH1Fs.end();++it)
	{
                delete it->second;          
	}
        TH1Fs.clear();
	for(std::map<std::string,TH2F*>::iterator it=TH2Fs.begin();it!=TH2Fs.end();++it)
	{
                delete it->second;
	}
        TH2Fs.clear();
}

TH1D* HistoPlane::Return_TH1F(const char* name)
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

void HistoPlane::WriteAll(int NbrPlatee,std::string namepdf)
{
   
  std::string name="Title:Plate "+ patch::to_string(NbrPlatee+1);
  for(std::map<std::string,TH1D*>::iterator it=TH1Fs.begin();it!=TH1Fs.end();++it)
  {
  		std::string nameth1= name+"TH1";
  		(it->second)->Write();
                if(pdf)(it->second)->Draw();
                if(pdf)canvas->Print(namepdf.c_str(),nameth1.c_str());
  		
  }
  for(std::map<std::string, TH2F*>::iterator it=TH2Fs.begin();it!=TH2Fs.end();++it)
  {
  		std::string nameth2= name+"TH2";
  		(it->second)->Write();
  		if(pdf)(it->second)->Draw("COLZ");
                if(pdf)canvas->Print(namepdf.c_str(),nameth2.c_str());
  }
}

void HistoPlane::Save(TFile* file,std::string namepdf)
{
    TH1F* ugly=new TH1F("ugly","ugly",1,0,1);
    ugly->Fill(0.5,global_total_time*200e-9);
    ugly->Write();
    delete ugly;
    std::string plate="Plate "+ patch::to_string(NbrPlatee+1);
    file->mkdir(plate.c_str(),plate.c_str());
    file->cd(plate.c_str());
    for(std::map<int,int>::iterator it = Time_Plates.begin(); it!=Time_Plates.end(); ++it) 
    {
      TH1Fs["Time_Distr"]->Fill(it->first,it->second);
      TH1Fs["Hits_Distr"]->Fill(it->second,1);
    }
    TH1Fs["Hits_Distr"]->Fill(0.0,(double)Nbrof0Hits);
    
    TH1Fs["Time_Distr"]->GetXaxis()->SetRange(0.0,Time_Plates.size()+10);
    TH1Fs["Hits_Distr"]->GetXaxis()->SetRange(0,2000);
    Return_TH1F("Hits_Distr_Noise")->Fill(0.0,(double)Nbrof0Hits);
    WriteAll(NbrPlatee,namepdf);
    std::string name="Noise_Flux_Hz";
    TH2Fs["Flux_Noise"]->Scale(1/(global_total_time*2e-7));
    TH2Fs["Flux_Noise"]->Write(name.c_str());
    name="Flux_Noise_Mean_Scaled";
    Means=TH2Fs["Flux_Noise"]->Integral()/(1.0*GetArea());
    TH2Fs["Flux_Noise"]->Scale(1/Means);
    TH2Fs["Flux_Noise"]->Write(name.c_str());
    name="Asic_Noise_Flux";
    TH2Fs["Flux_Noise_Asic"]->Scale(1/(global_total_time*2e-7));
    TH2Fs["Flux_Noise_Asic"]->Write(name.c_str());
    name="Asic_Noise_Flux_Mean_Scaled";
    TH2Fs["Flux_Noise_Asic"]->Scale(1/Means);
    TH2Fs["Flux_Noise_Asic"]->Write(name.c_str());
    name="Event_Flux_Hz";
    TH2Fs["Flux_Events"]->Scale(1/(global_total_time*2e-7));
    TH2Fs["Flux_Events"]->Write(name.c_str());
    name="Asic_Event_Flux";
    TH2Fs["Flux_Events_Asic"]->Scale(1/(global_total_time*2e-7));
    TH2Fs["Flux_Events_Asic"]->Write(name.c_str());
    
    if(_Distr==true)Write_TH1_Hit_In_Pad_Per_RamFull(file,plate);
    
}

 void HistoPlane::Fill_TH1_Timestamps_Distribution(unsigned int& Dif_Id ,int& Asic_Id,int& Channel_Id,double& timestamp)
  {
       std::vector<int>vec{Dif_Id,Asic_Id,Channel_Id};
       Pads_Distr2[vec]->Fill(timestamp);
       vec.pop_back();
                        Asics_Distr2[vec]->Fill(timestamp);
			vec.pop_back();
                        Difs_Distr2[vec[0]]->Fill(timestamp);
  }
