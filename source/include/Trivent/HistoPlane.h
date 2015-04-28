#ifndef HISTO_PLANE
#define HISTO_PLANE
#include <string>
#include "TH1.h"
#include "Patch.h"
#include "TFile.h"
#include "TH2.h"
#include <map>
#include <iostream>
#include "Colors.h"
#include <cmath>

class HistoPlane
{
 public:
  ~HistoPlane();
  HistoPlane(){};
  HistoPlane(bool Distr,int NbrPlate,std::vector<int>Difs_Names,int SizeX, int SizeY,std::vector< std::basic_string<char>  >& vec_name_th1,std::vector< std::basic_string<char>  >& vec_name_th2,std::vector< std::basic_string<char>  >& vec_name_th2_Asic);
  void inline Clear_Time_Plates(){Time_Plates.clear();};
  void inline Clear_Time_Plates_perRun(){Times_Plates_perRun.clear();};
  void Fill_Time_Plates(int timeStamp){Time_Plates[timeStamp]++;};
  void inline Fill_Time_Plates_perRun(int  timeStamp){Times_Plates_perRun[timeStamp]++;};
  int  inline Get_local_max(){return local_max;};
  void inline Init_local_min_max(){Set_local_max(0);Set_local_min(99999999);};
  void inline Init_Hit_In_Pad_Per_RamFull()
  {
        Hit_In_Pad_Per_RamFull.clear();
  }
  HistoPlane(const HistoPlane &source);
  int inline Get_local_min(){return local_min;};
  void inline Set_local_max( int lm){local_max=lm;};
  void inline Set_local_min( int lm){local_min=lm;};
  void inline Set_Total_Time(){if(local_max==0) return ;else if (local_max==local_min) total_time+=local_max;else total_time+=local_max-local_min;};
  void inline Set_Nbrof0Hits(){if(local_min==99999999)return;else Nbrof0Hits+=(local_max-local_min-Times_Plates_perRun.size());};
  unsigned long long int  inline Get_Total_Time(){return total_time;};
  double inline Get_Means(){return Means;};
  int long long inline Get_Nbrof0Hits(){return Nbrof0Hits;};
  int inline Get_NbrPlate(){return NbrPlatee;};
  void  Save(TFile* file);
  void inline Set_hit_other(){hit_other+=1;};
  int inline Get_hit_other(){return hit_other;};
  void inline Set_hit_trigger(){hit_trigger+=1;};
  void inline Fill_Hit_In_Pad_Per_RamFull(int &Dif_Id ,int &Asic_Id,int &Channel_Id)
  {
  	std::vector<int>vec{Dif_Id,Asic_Id,Channel_Id};Hit_In_Pad_Per_RamFull[vec]+=1.0;
  };
  void inline Fill_TH1_Hit_In_Pad_Per_RamFull()
  {
        for(std::map<std::vector<int>,double>::iterator it=Hit_In_Pad_Per_RamFull.begin();it!=Hit_In_Pad_Per_RamFull.end();++it)
        {
			//std::cout<<(it->first)[0]<<"  "<<(it->first)[1]<<" "<<(it->first)[2]<<"  "<<it->second<<std::endl;
			Pads_Distr[it->first]->Fill(it->second);
                        std::vector<int>a=(it->first);
			a.pop_back();
                        Asics_Distr[a]->Fill(it->second);
			a.pop_back();
                        Difs_Distr[a[0]]->Fill(it->second);
   			Plate_Dist->Fill(it->second);  
	}
  }
  void inline Write_TH1_Hit_In_Pad_Per_RamFull(TFile* file,std::string plate)
  {
        Plate_Dist->Write();
	for(unsigned int i=0;i<_Difs_Names.size();++i)
	{
                std::string Diff=plate+"/Dif"+ patch::to_string(_Difs_Names[i])+"Plane"+patch::to_string(NbrPlatee+1);
    		file->mkdir(Diff.c_str(),Diff.c_str());
    		file->cd(Diff.c_str());
                std::vector<int>a={_Difs_Names[i]};
	 	Difs_Distr[a[0]]->Write();
		for(unsigned int j=1;j<=24;++j)
		{
                        std::string Asicc=Diff+"/AsicinDif"+ patch::to_string(a[0]);
                        if(j==1)
			{
    	        	file->mkdir(Asicc.c_str(),Asicc.c_str());
			}
                        file->cd(Asicc.c_str());
                        a.push_back(j);
                        Asics_Distr[a]->Write();
			for(unsigned int k=0;k<=63;++k)
		        {
                          std::string Padd=Asicc+"/PadsinAsic"+ patch::to_string(j);
    	        	  if(k==0)
			  {
			  file->mkdir(Padd.c_str(),Padd.c_str());
                	  file->cd(Padd.c_str());
			  }
                          a.push_back(k);
                          Pads_Distr[a]->Write();
                          a.pop_back();
			}
                        a.pop_back();
		}
	 }
  }
  int inline Get_hit_trigger(){return hit_trigger;};
  void inline Fill_Calibration(int &Dif_Id ,int &Asic_Id,int &Channel_Id){ Calibration[Dif_Id][Asic_Id-1][Channel_Id]+=1;};
  int inline Get_Calibration(int &Dif_Id,int &Asic_Id,int &Channel_Id){return Calibration[Dif_Id][Asic_Id-1][Channel_Id] ;};
  void inline Get_Flux()
  {
    	double MEAN=0;
    	double RMS=0;
    	unsigned int  Number_Pads_Touched=0;
    	for(std::map<int,std::array<std::array<double,64>,24>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{   
      				if((it->second)[i][j]>=1)
				{
					MEAN+=(it->second)[i][j];
					Number_Pads_Touched++;
				}
    			}
		}
	}
    	MEAN/=Number_Pads_Touched;
    	for(std::map<int,std::array<std::array<double,64>,24>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{    
      				RMS+=(it->second[i][j]-MEAN)*(it->second[i][j]-MEAN);
    			}
		}
	}
    	RMS/=Number_Pads_Touched;
    	RMS=sqrt(RMS);
    	std::cout<<MEAN<<" "<<RMS<<std::endl;
    	for(std::map<int,std::array<std::array<double,64>,24>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{      
      				if(it->second[i][j]==0) (it->second)[i][j]=1;
      				else if(it->second[i][j]==10*RMS) it->second[i][j]=0;
      				else it->second[i][j]=254.*(1+std::exp( -MEAN/RMS))/(1+std::exp( ((it->second[i][j])-MEAN)/RMS))+1;
      
    			}
  		}
	}
  }
  void inline Print_Calibration(std::ostream& file)
  {
	for(std::map<int,std::array<std::array<double,64>,24>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{
		         	file<<"s.SetGain("<<it->first<<","<<i+1<<","<<j<<","<<(it->second)[i][j]<<")"<<std::endl;
			}
		}
	}
  };
  TH1F* Return_TH1F(const char* name);
  TH2F* Return_TH2F(const char* name);
  double inline GetArea(){return _SizeX*_SizeY;};
  void  ScaleHisto(const char* name,float i);
  void  WriteAll();
  double inline Efficiency(){std::cout<<hit_other<<"  "<<hit_trigger<<std::endl;return 1.0*hit_trigger/(hit_other+hit_trigger);};
  
 private:
  TH1F* Plate_Dist;
  std::map<int,TH1F*>Difs_Distr;
  std::map<std::vector<int>,TH1F*>Asics_Distr;
  std::map<std::vector<int>,TH1F*>Pads_Distr;
  unsigned int hit_other;
  unsigned int hit_trigger;
  std::map<int,int>Time_Plates;
  std::map<int,int>Times_Plates_perRun;
  std::map<std::string,TH1F*>TH1Fs;
  std::map<std::string,TH2F*>TH2Fs;
  double Means;
  int NbrPlatee;
  int long long Nbrof0Hits;
  int  long long local_max;
  int long long  local_min;
  unsigned long long int total_time;
  double _SizeX;
  double _SizeY;
  std::map<std::vector<int>,double>Hit_In_Pad_Per_RamFull;
  std::vector<int>_Difs_Names;
  std::map<int,std::array<std::array<double,64>, 24>>Calibration;
  bool _Distr;
};
#endif
