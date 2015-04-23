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
  HistoPlane()
  {
   
  };
  HistoPlane(int NbrPlate,int SizeX, int SizeY,std::vector< std::basic_string<char>  >& vec_name_th1,std::vector< std::basic_string<char>  >& vec_name_th2,std::vector< std::basic_string<char>  >& vec_name_th2_Asic);
  void inline Clear_Time_Plates(){Time_Plates.clear();};
  void inline Clear_Time_Plates_perRun(){Times_Plates_perRun.clear();};
  void Fill_Time_Plates(int timeStamp){Time_Plates[timeStamp]++;};
  void inline Fill_Time_Plates_perRun(int  timeStamp){Times_Plates_perRun[timeStamp]++;};
  int  inline Get_local_max(){return local_max;};
  void inline Init_local_min_max(){Set_local_max(0);Set_local_min(99999999);};
  void inline Init_Hit_In_Asic_Per_RamFull()
  {
	Hit_In_Asic_Per_RamFull.clear();
	for(unsigned int i=0;i<Hit_In_Pad_Per_RamFull.size();++i) Hit_In_Pad_Per_RamFull[i].clear();
        Hit_In_Pad_Per_RamFull.clear();
  }
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
  void inline Fill_Hit_In_Asic_Per_RamFull(int Asic_Id,int Channel_Id){Hit_In_Asic_Per_RamFull[Asic_Id]+=1;Hit_In_Pad_Per_RamFull[Asic_Id][Channel_Id]+=1;};
  void inline Fill_TH1_Hit_In_Asic_Per_RamFull()
  {
        int sum=0;
	for(std::map<int,int>::iterator it=Hit_In_Asic_Per_RamFull.begin();it!=Hit_In_Asic_Per_RamFull.end();++it)
        {
		//std::cout<<it->first<<"  "<<it->second<<std::endl;
		Asic[it->first-1]->Fill(it->second,1.0);
                sum += it->second;
	}
        for(std::map<int,std::map<int,int>>::iterator it=Hit_In_Pad_Per_RamFull.begin();it!=Hit_In_Pad_Per_RamFull.end();++it)
        {
		for(std::map<int,int>::iterator itt=(it->second).begin();itt!=(it->second).end();++itt)
		{
			//std::cout<<it->first<<"  "<<itt->first<<" "<<itt->second<<std::endl;
			Difs_Distr[it->first][itt->first]->Fill(itt->second,1.0);
		}
	}
        //std::cout<<sum<<std::endl;
        Dif_Dist->Fill(sum,1);
  }
  void inline Write_TH1_Hit_In_Asic_Per_RamFull(TFile* file,std::string plate)
  {
        for(unsigned int i=1;i<25;++i)
	{
		std::string Asicc=plate+"/Asic"+ patch::to_string(i)+"Plane"+patch::to_string(NbrPlatee+1);
    		file->mkdir(Asicc.c_str(),Asicc.c_str());
    		file->cd(Asicc.c_str());
                
	 	Asic[i-1]->Write();
                std::string Padss=Asicc+"/PadsinAsic"+ patch::to_string(i)+"Plane"+patch::to_string(NbrPlatee+1);
    	        file->mkdir(Padss.c_str(),Padss.c_str());
                file->cd(Padss.c_str());
			for(unsigned int j=0;j<Difs_Distr[i].size();++j)
			{
				Difs_Distr[i][j]->Write();
			}
	}
  }
  int inline Get_hit_trigger(){return hit_trigger;};
  void inline Fill_Calibration(int &a ,int &b,int &c){std::vector<int>vec{a,b,c};Calibration[vec]+=1.0;};
  double inline Get_Calibration(int &a,int &b,int & c){std::vector<int>vec{a,b,c};return Calibration[vec];};
  void inline Get_Flux()
  {
    double MEAN=0;
    double RMS=0;
    unsigned int  Number_Pads_Touched=Calibration.size();
    for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    {   
      MEAN+=it->second;
    }
    MEAN/=Number_Pads_Touched;
    for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    {   
      RMS+=(it->second-MEAN)*(it->second-MEAN);
    }
    RMS/=Number_Pads_Touched;
    RMS=sqrt(RMS);
    std::cout<<MEAN<<" "<<RMS<<std::endl;
    for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    {  
      if(it->second==0) continue;
      else if(it->second==10*RMS) it->second=0;
      else it->second=254.*(1+std::exp( -MEAN/RMS))/(1+std::exp( ((it->second)-MEAN)/RMS))+1;
      
    }
  }
  void inline Print_Calibration(std::ostream& file){for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it){/*std::cout<<"s.ChangeGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(it->second)<<")"<<std::endl;*/file<<"s.SetGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(int)(it->second)<<")"<<std::endl;}};
  TH1F* Return_TH1F(const char* name);
  TH2F* Return_TH2F(const char* name);
  double inline GetArea(){return _SizeX*_SizeY;};
  void  ScaleHisto(const char* name,float i);
  void  WriteAll();
  double inline Efficiency(){std::cout<<hit_other<<"  "<<hit_trigger<<std::endl;return 1.0*hit_trigger/(hit_other+hit_trigger);};
  
 private:
  std::vector<TH1F*>Asic;
  std::map<int,std::vector<TH1F*>>Difs_Distr;
  TH1F* Dif_Dist;
  unsigned int hit_other;
  unsigned int hit_trigger;
  std::map<int,int>Time_Plates;
  std::map<int,int>Hit_In_Asic_Per_RamFull;
  std::map<int,std::map<int,int>>Hit_In_Pad_Per_RamFull;
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
  std::map<std::vector<int>,double>Calibration;
};
#endif
