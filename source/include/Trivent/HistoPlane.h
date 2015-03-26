#ifndef HISTO_PLANE
#define HISTO_PLANE
#include<string>
#include "TH1.h"
#include "TFile.h"
#include "TH2.h"
#include <map>
#include<iostream>
#include "Colors.h"
#include <fstream>
extern std::ofstream file;

class HistoPlane
{
public:
~HistoPlane();
HistoPlane(){};
HistoPlane(int NbrPlate,int SizeX, int SizeY,std::vector< std::basic_string<char>  >& vec_name_th1,std::vector< std::basic_string<char>  >& vec_name_th2,std::vector< std::basic_string<char>  >& vec_name_th2_Asic);
void inline Clear_Time_Plates(){Time_Plates.clear();};
void inline Clear_Time_Plates_perRun(){Times_Plates_perRun.clear();};
void Fill_Time_Plates(int timeStamp){Time_Plates[timeStamp]++;};
void inline Fill_Time_Plates_perRun(int  timeStamp){Times_Plates_perRun[timeStamp]++;};
int  inline Get_local_max(){return local_max;};
void inline Init_local_min_max(){Set_local_max(0);Set_local_min(99999999);};
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
int inline Get_hit_trigger(){return hit_trigger;};
void inline Fill_Calibration(int &a ,int &b,int &c){std::vector<int>vec{a,b,c};Calibration[vec]+=1;};
double inline Get_Calibration(int &a,int &b,int & c){std::vector<int>vec{a,b,c};return Calibration[vec];};
void inline Get_Flux()
{
        double max=Calibration.begin()->second;
        double min=Calibration.begin()->second;
	for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it)
        {   
                        if(it->second>=1000*Means)it->second=-1;
                        else
  			{
				double val=it->second;
        			if(val>max)max=val;
                		if(val<min)min=val;
			}
	}
        for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it)
	{ 
                if (it->second==-1)it->second=1;
		else (it->second)=((1-(it->second-min)/(max-min))*254+1); 
	}
};
void inline Print_Calibration(){for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it){/*std::cout<<"s.ChangeGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(it->second)<<")"<<std::endl;*/file<<"s.SetGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(it->second)<<")"<<std::endl;}};
TH1F* Return_TH1F(const char* name);
TH2F* Return_TH2F(const char* name);
double inline GetArea(){return _SizeX*_SizeY;};
void  ScaleHisto(const char* name,float i);
void  WriteAll();
double inline Efficiency(){std::cout<<hit_other<<"  "<<hit_trigger<<std::endl;return 1.0*hit_trigger/(hit_other+hit_trigger);};

private:
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
std::map<std::vector<int>,double>Calibration;
};
#endif
