#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include<cstdlib> 
#include"TGraphErrors.h"
#include"TFile.h"
#include<string>
#include"TAxis.h"
#include"Patch.h"
#include<cstdlib>

int main (int argc,char* argv[]) 
{
  std::ifstream ifs ("Results.txt", std::ifstream::in);
  if(ifs.is_open())
  {
  TFile* file=new TFile("Plots.root","RECREATE");
  std::map<unsigned int,std::vector<double> >Efficiency;
  std::map<unsigned int,std::vector<double> >Multiplicity;
  std::map<unsigned int,std::vector<double> >Multiplicity_error;
  std::map<unsigned int,std::vector<double> >Efficiency_error;
  
  std::string line;
  while(std::getline(ifs,line))
  {
    std::stringstream iss(line,std::ios_base::in);
    std::string Value;
    unsigned int RunNumber=0;
    double XAxis=0;
    unsigned int index=1;
    while(std::getline(iss,Value,';'))
    {
      double Valuee=atof(Value.c_str());
      if(index==1)RunNumber=Valuee;
      if(index==2)XAxis=Valuee;
      if(index>2 && (index-2)%4==0)Multiplicity_error[XAxis].push_back(Valuee);
      else if(index>2 && (index-2)%4==3)Multiplicity[XAxis].push_back(Valuee);
      else if(index>2 && (index-2)%4==2)Efficiency_error[XAxis].push_back(Valuee);
      else if (index>2) Efficiency[XAxis].push_back(Valuee);
      ++index;
    }
  }
  std::map<unsigned int,std::vector<double> >::iterator it=Efficiency.begin();
  //std::cout<<"kkkkk"<<(it->second).size()<<"  "<<Efficiency.size()<<std::endl;
  std::vector<TGraphErrors>Efficiencies((it->second).size(),TGraphErrors());
  std::vector<TGraphErrors>Multiplicities((it->second).size(),TGraphErrors());
  std::vector<std::vector<std::vector<double> > >Eff((it->second).size());
  std::vector<std::vector<std::vector<double> > >Mul((it->second).size());
  double index=0;
  for(std::map<unsigned int,std::vector<double> >::iterator it=Efficiency.begin();it!=Efficiency.end();++it)
  {
    
    for(unsigned int i=0;i!=(it->second).size();++i)
    {
      Eff[i].push_back({index,double(it->first),(it->second)[i],0.0,Efficiency_error[it->first][i]});
      Mul[i].push_back({index,double(it->first),Multiplicity[it->first][i],0.0,Multiplicity_error[it->first][i]});
      
    }
    ++index;
  }
    
    for(unsigned int i=0;i!=Eff.size();++i)
    
    {
      for(unsigned int j=0;j!=Eff[0].size();++j)
      {  
     int index=int(Eff[i][j][0]);
     double X=Eff[i][j][1];
     double Y=Eff[i][j][2];
     double x_error=Eff[i][j][3];
     double y_error=Eff[i][j][4];
     double X1=Mul[i][j][1];
     double Y1=Mul[i][j][2];
     double x1_error=Mul[i][j][3];
     double y1_error=Mul[i][j][4];
     Efficiencies[i].SetPoint(index,X,Y);
     Efficiencies[i].SetTitle(argv[1]);
     Efficiencies[i].SetLineColor(i+1);
     Efficiencies[i].SetLineWidth(2);
     Efficiencies[i].GetXaxis()->SetTitle(argv[3]);
     Efficiencies[i].GetYaxis()->SetTitle("Efficiency");
     Efficiencies[i].SetPointError(index,x_error,y_error);
     Multiplicities[i].SetPoint(index,X1,Y1);
     Multiplicities[i].SetTitle(argv[2]);
     Multiplicities[i].SetLineColor(i+1);
     Multiplicities[i].GetXaxis()->SetTitle(argv[3]);
     Multiplicities[i].SetLineWidth(2);
     Multiplicities[i].GetYaxis()->SetTitle("Multiplicity");
     Multiplicities[i].SetPointError(index,x1_error,y1_error);
    }
  }
  
  std::string nameeff="Efficiency_Plate_";
  std::string namemul="Multiplicity_Plate_";
  for(unsigned int i=0;i!=Efficiencies.size();++i)
  {
    
    Multiplicities[i].Write((nameeff+patch::to_string(i)).c_str());
    Efficiencies[i].Write((namemul+patch::to_string(i)).c_str());
  }
  file->Close();
   
  }
  else
  {
  std::cout<<"Impossible to find the file"<<std::endl;
   std::exit(1); 
  }
  
  return 0;
}
