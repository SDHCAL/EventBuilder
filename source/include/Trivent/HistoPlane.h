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
#include <array>
#include<numeric>
#include<algorithm>
#include "Trivent/Mapping.h"
#include "Geometry/Geometry.h"
#include "TTree.h"
#define size_pad 10.4125
#define size_strip 2.5
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
  void inline Init_local_min_max(){Set_local_max(0);Set_local_min(99999999);};
  void inline Init_Hit_In_Pad_Per_RamFull()
  {
        Hit_In_Pad_Per_RamFull.clear();
  }
  HistoPlane(const HistoPlane &source);
  int inline Get_local_min(){return local_min;};
  int  inline Get_local_max(){return local_max;};
  void inline Set_local_max( int lm){local_max=lm;};
  void inline Set_local_min( int lm){local_min=lm;};
  void inline Set_Total_Time(){if(local_max==0) return ;else if (local_max==local_min) total_time+=local_max;else total_time+=local_max-local_min;};
  void inline Set_Global_Total_Time(unsigned long long int time){global_total_time+=time;};
  unsigned long long int inline Get_Global_Total_Time(){return global_total_time;};
  void inline Set_Nbrof0Hits(){if(local_min==99999999)return;else Nbrof0Hits+=(local_max-local_min-Times_Plates_perRun.size());};
  unsigned long long int  inline Get_Total_Time(){return total_time;};
  double inline Get_Means(){return Means;};
  int long long inline Get_Nbrof0Hits(){return Nbrof0Hits;};
  int inline Get_NbrPlate(){return NbrPlatee;};
  void  Save(TFile* file,std::string);
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
	Plate_Dist2->Write();
	for(unsigned int i=0;i<_Difs_Names.size();++i)
	{
                std::string Diff=plate+"/Dif"+ patch::to_string(_Difs_Names[i])+"Plane"+patch::to_string(NbrPlatee+1);
    		file->mkdir(Diff.c_str(),Diff.c_str());
    		file->cd(Diff.c_str());
                std::vector<int>a={_Difs_Names[i]};
	 	Difs_Distr[a[0]]->Write();
		Difs_Distr2[a[0]]->Write();
		for(unsigned int j=1;j<=48;++j)
		{
                        std::string Asicc=Diff+"/AsicinDif"+ patch::to_string(a[0]);
                        if(j==1)
			{
    	        	file->mkdir(Asicc.c_str(),Asicc.c_str());
			}
                        file->cd(Asicc.c_str());
                        a.push_back(j);
                        Asics_Distr[a]->Write();
			Asics_Distr2[a]->Write();
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
			  Pads_Distr2[a]->Write();
                          a.pop_back();
			}
                        a.pop_back();
		}
	 }
  }
  int inline Get_hit_trigger(){return hit_trigger;};
  //void inline Set_Calibration(int &Dif_Id ,int &Asic_Id,int &Channel_Id,double &v){Calibration[Dif_Id][Asic_Id-1][Channel_Id]+=v;};
  void inline Set_Calibration(int Dif_Id ,int Asic_Id,int Channel_Id,double v){Calibration[Dif_Id][Asic_Id-1][Channel_Id]+=v;};
  //int inline Get_Calibration(int &Dif_Id,int &Asic_Id,int &Channel_Id){return Calibration[Dif_Id][Asic_Id-1][Channel_Id] ;};
  void inline Fill_NumberHitsDistribution(int &a ,int &b,int &c){std::vector<int>vec{a,b,c};NumberHitsDistribution[vec]+=1;};
  void inline SaveCalibrationMap(std::ostream& file,std::map<int,std::vector<double>>& SinCos,Geometry& geom,int ok,int end,TTree& t1)
  {
        Int_t I, J, dif_id,asic_id,channel_id,Num;
	if(ok==0)
	{
   		t1.Branch("I",&I,"I/I");
   		t1.Branch("J",&J,"J/I");
   		t1.Branch("dif_id",&dif_id,"dif_id/I");
   		t1.Branch("asic_id",&asic_id,"asic_id/I");
   		t1.Branch("channel_id",&channel_id,"channel_id/I");
		t1.Branch("Num",&Num,"Num/I");
	}
	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{   

				double ca=SinCos[it->first][0];
	  			double sa=SinCos[it->first][1];
    				double cb=SinCos[it->first][2];
	  			double sb=SinCos[it->first][3];
    				double cg=SinCos[it->first][4];
	 			double sg=SinCos[it->first][5];
				float pos[3];
				dif_id=it->first;
                                asic_id=i+1;
				channel_id=j;
				Num=(it->second)[i][j];
				unsigned int NbrPlate =geom.GetDifNbrPlate(it->first)-1;
    				float Z= geom.GetPlatePositionZ(NbrPlate);
    				unsigned int K =geom.GetDifNbrPlate(it->first);
    				I=0;
    				J=0;
    				if(geom.GetDifType(it->first)==pad) 
    				{
      					I =(1+MapILargeHR2[j]+AsicShiftI[i+1])+geom.GetDifPositionX(dif_id);
      					J =(32-(MapJLargeHR2[j]+AsicShiftJ[i+1]))+geom.GetDifPositionY(dif_id);
      					pos[0] = cg*cb*I*size_pad+(-sg*ca+cg*sb*sa)*J*size_pad+(sg*sa+cg*sb*ca)*Z+geom.GetPlatePositionX(NbrPlate);
      					pos[1] = sg*cb*I*size_pad+(cg*ca+sg*sb*sa)*J*size_pad+(-cg*sa+sg*sb*ca)*Z+geom.GetPlatePositionY(NbrPlate);
      					pos[2] = -sb*I*size_pad+cb*sa*J*size_pad+cb*ca*Z;
    				}
      				file<<I<<";"<<J<<";"<<it->first<<";"<<i+1<<";"<<j<<";"<<(it->second)[i][j]<<std::endl;
				t1.Fill();
    			}
		}
	}
  }
  //double inline Get_Calibration(int &a,int &b,int & c){std::vector<int>vec{a,b,c};return Calibration[vec];};
 // void inline Set_Calibration(int &a,int &b,int & c,double &v){std::vector<int>vec{a,b,c};Calibration[vec]=v;};


/*void inline Get_Calibration(double minvalue=1.0,double maxvalue=254.0,double RMScutfactor=10.0,bool update=false)
  {
    //double MEAN=0;
    //double RMS=0;
    //double MPV=-1;
    double MEAN2=0;
    double RMS2=0;
    double MPV2=-1;
    //TH1F *C =  new TH1F("C","C",1000000,0,1000000);
    
    unsigned int  Number_Pads_Touched=NumberHitsDistribution.size();
    //std::cout << Number_Pads_Touched << std::endl;
    //for(std::map<std::vector<int>,long int>::iterator it=NumberHitsDistribution.begin();it!=NumberHitsDistribution.end();++it)
    //{   
      //MEAN+=it->second;
      //RMS+=it->second*it->second;
      //C->Fill(it->second);
      //if (it->second>MPV) MPV=it->second;
    //}
    for(std::map<std::vector<int>,long int>::iterator it=NumberHitsDistribution.begin();it!=NumberHitsDistribution.end();++it)
    {   
      MEAN2+=it->second;
      if(it->second>MPV2)MPV2=it->second;
    }
    MEAN2/=Number_Pads_Touched;
    for(std::map<std::vector<int>,long int>::iterator it=NumberHitsDistribution.begin();it!=NumberHitsDistribution.end();++it)
    {   
      RMS2+=(MEAN2-it->second)*(MEAN2-it->second);
    }
    RMS2/=Number_Pads_Touched;
    RMS2=std::sqrt(RMS2);
    //MEAN=C->GetMean();
    //MPV=C->GetMaximumBin();
    //RMS=C->GetRMS();
    //MEAN/=Number_Pads_Touched;
    //RMS=sqrt(RMS/Number_Pads_Touched-MEAN*MEAN);
    std::cout<<MEAN2<<" "<<RMS2<<" "<<MPV2<<std::endl;
    for(std::map<std::vector<int>,long int>::iterator it=NumberHitsDistribution.begin();it!=NumberHitsDistribution.end();++it)
    {  
      double facteur=0;
      if(it->second==0  || it->second>=RMScutfactor*RMS2) facteur=0;
      else facteur=maxvalue*(1+std::exp( -MPV2/RMS2))/(1+std::exp( (it->second-MPV2)/RMS2))+minvalue;
      if (update) Calibration[it->first]*=facteur; else Calibration[it->first]=facteur;      
    }
  }
*/

void inline Get_Calibration()
  {
    	double MEAN=0;
    	double RMS=0;
    	unsigned int  Number_Pads_Touched=0;
    	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
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
    	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
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
    	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
    	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{      
      				if(it->second[i][j]==0) (it->second)[i][j]=1;
      				else if(it->second[i][j]==4*RMS) it->second[i][j]=0;
      				else it->second[i][j]=254.*(1+std::exp( -MEAN/RMS))/(1+std::exp( ((it->second[i][j])-MEAN)/RMS))+1;
      
    			}
  		}
	}
  }
  void inline SaveCalibration(std::ostream& file)
  {
	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{
		         	file<<"s.SetGain("<<it->first<<","<<i+1<<","<<j<<","<<int((it->second)[i][j])<<")"<<std::endl;
			}
		}
	}
  };
    void inline Print_Calibration(std::ostream& file)
  {
	for(std::map<int,std::array<std::array<double,64>,48>>::iterator it=Calibration.begin();it!=Calibration.end();++it)
	{
		for(unsigned int i=0;i<(it->second).size();++i)
		{
                	for(unsigned int j=0;j<(it->second)[i].size();++j)
			{
		         	std::cout<<"s.SetGain("<<it->first<<","<<i+1<<","<<j<<","<<int((it->second)[i][j])<<")"<<std::endl;
			}
		}
	}
  };
  //void inline Print_Calibration(std::ostream& file){for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it){std::cout<<"s.ChangeGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(it->second)<<")"<<std::endl;*/file<<"s.SetGain("<<(it->first)[0]<<","<<(it->first)[1]<<","<<(it->first)[2]<<","<<(int)(it->second)<<")"<<std::endl;}};
 // void inline SaveCalibration(std::ostream& file){for(std::map<std::vector<int>,double>::iterator it=Calibration.begin();it!=Calibration.end();++it){file<<(it->first)[0]<<" "<<(it->first)[1]<<" "<<(it->first)[2]<<" "<<int(it->second)<<std::endl;}}
  TH1D* Return_TH1F(const char* name);
  TH2F* Return_TH2F(const char* name);
  double inline GetArea(){return _SizeX*_SizeY;};
  void  ScaleHisto(const char* name,float i);
  void  WriteAll(int NbrPlatee, std::string);
  double inline Efficiency(){std::cout<<hit_other<<"  "<<hit_trigger<<std::endl;return 1.0*hit_trigger/(hit_other+hit_trigger);};
  void Fill_TH1_Timestamps_Distribution(unsigned int& Dif_Id ,int& Asic_Id,int& Channel_Id,double& timestamp);
 private:
  TH1D* Plate_Dist;
  TH1D* Plate_Dist2;
  std::map<int,TH1D*>Difs_Distr;
  std::map<std::vector<int>,TH1D*>Asics_Distr;
  std::map<std::vector<int>,TH1D*>Pads_Distr;
  std::map<int,TH1D*>Difs_Distr2;
  std::map<std::vector<int>,TH1D*>Asics_Distr2;
  std::map<std::vector<int>,TH1D*>Pads_Distr2;
  unsigned int hit_other;
  unsigned int hit_trigger;
  std::map<int,int>Time_Plates;
  std::map<int,int>Times_Plates_perRun;
  std::map<std::string,TH1D*>TH1Fs;
  std::map<std::string,TH2F*>TH2Fs;
  double Means;
  int NbrPlatee;
  int long long Nbrof0Hits;
  int  long long local_max;
  int long long  local_min;
  unsigned long long int total_time;
  static unsigned long long int global_total_time;
  double _SizeX;
  double _SizeY;
  std::map<std::vector<int>,double>Hit_In_Pad_Per_RamFull;
  std::vector<int>_Difs_Names;
  std::map<int,std::array<std::array<double,64>, 48>>Calibration;
  //std::map<std::vector<int>,double>Calibration;
  std::map<std::vector<int>,long int>NumberHitsDistribution;
  bool _Distr;
};
#endif
