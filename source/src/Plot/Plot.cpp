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
    std::map<std::string,std::map<unsigned int,std::vector<double> > >Efficiency;
    std::map<std::string,std::map<unsigned int,std::vector<double> > >Multiplicity;
    std::map<std::string,std::map<unsigned int,std::vector<double> > >Multiplicity_error;
    std::map<std::string,std::map<unsigned int,std::vector<double> > >Efficiency_error;
    std::string ThresholdReal="";
    std::string line;
    while(std::getline(ifs,line))
    {
      std::stringstream iss(line,std::ios_base::in);
      std::stringstream iss2(line,std::ios_base::in);
      std::string Value;
      unsigned int RunNumber=0;
      double XAxis=0;
      unsigned int index=1;
      unsigned int index2=1;
      std::string Threshold="";
      while(std::getline(iss2,Value,';'))
      {
        Threshold=Value;
        ++index2;
      }
      if(index2==2) 
      {
        ThresholdReal=Threshold;
        //std::cout<<ThresholdReal<<std::endl;
        continue;
      }
      else
      while(std::getline(iss,Value,';'))
      {
        double Valuee=atof(Value.c_str());
        if(index==1) RunNumber=(unsigned int)(Valuee);
        else if(index==2) XAxis=Valuee;
        else if(index>2 && (index-2)%4==0)Multiplicity_error[ThresholdReal][XAxis].push_back(Valuee);
        else if(index>2 && (index-2)%4==3)Multiplicity[ThresholdReal][XAxis].push_back(Valuee);
        else if(index>2 && (index-2)%4==2)Efficiency_error[ThresholdReal][XAxis].push_back(Valuee);
        else if(index>2) Efficiency[ThresholdReal][XAxis].push_back(Valuee);
        ++index;
      }
    }
    std::map<std::string,std::vector<TGraphErrors>>Efficiencies;
    std::map<std::string,std::vector<TGraphErrors>>Multiplicities;
    std::map<std::string,std::vector<std::vector<std::vector<double> > > >Eff;
    std::map<std::string,std::vector<std::vector<std::vector<double> > > >Mul;
    for(std::map<std::string,std::map<unsigned int,std::vector<double> > >::iterator it=Efficiency.begin();it!=Efficiency.end();++it)
    {
      std::map<unsigned int,std::vector<double> >::iterator itt=(it->second).begin();
      unsigned int size=(itt->second).size();
      Efficiencies[it->first]=std::vector<TGraphErrors>(size,TGraphErrors());
      Multiplicities[it->first]=std::vector<TGraphErrors>(size,TGraphErrors());
      Eff[it->first]=std::vector<std::vector<std::vector<double>>>(size);
      Mul[it->first]=std::vector<std::vector<std::vector<double>>>(size);
    }
 
    for(std::map<std::string,std::map<unsigned int,std::vector<double> > >::iterator it=Efficiency.begin();it!=Efficiency.end();++it)
    {
      double index=0;
      for(std::map<unsigned int,std::vector<double> >::iterator itt=Efficiency[it->first].begin();itt!=Efficiency[it->first].end();++itt)
      {
        for(unsigned int i=0;i!=(itt->second).size();++i)
        {
          Eff[it->first][i].push_back({index,double(itt->first),(itt->second)[i],0.0,Efficiency_error[it->first][itt->first][i]});
          Mul[it->first][i].push_back({index,double(itt->first),Multiplicity[it->first][itt->first][i],0.0,Multiplicity_error[it->first][itt->first][i]});
        }
        ++index;
      }
    }
    
    for(std::map<std::string,std::vector<std::vector<std::vector<double> > > >::iterator it=Eff.begin();it!=Eff.end();++it)
    {
      for(unsigned int i=0;i!=Eff[it->first].size();++i)
      {
        for(unsigned int j=0;j!=Eff[it->first][i].size();++j)
        {  
      int index=int(Eff[it->first][i][j][0]);
      double X=Eff[it->first][i][j][1];
     double Y=Eff[it->first][i][j][2];
     double x_error=Eff[it->first][i][j][3];
     double y_error=Eff[it->first][i][j][4];
     double X1=Mul[it->first][i][j][1];
     double Y1=Mul[it->first][i][j][2];
     double x1_error=Mul[it->first][i][j][3];
     double y1_error=Mul[it->first][i][j][4];
     Efficiencies[it->first][i].SetPoint(index,X,Y);
     Efficiencies[it->first][i].SetTitle(argv[1]);
     Efficiencies[it->first][i].SetLineColor(i+1);
     Efficiencies[it->first][i].SetLineWidth(2);
     Efficiencies[it->first][i].GetXaxis()->SetTitle(argv[3]);
     Efficiencies[it->first][i].GetYaxis()->SetTitle("Efficiency");
     Efficiencies[it->first][i].SetPointError(index,x_error,y_error);
     Multiplicities[it->first][i].SetPoint(index,X1,Y1);
     Multiplicities[it->first][i].SetTitle(argv[2]);
     Multiplicities[it->first][i].SetLineColor(i+1);
     Multiplicities[it->first][i].GetXaxis()->SetTitle(argv[3]);
     Multiplicities[it->first][i].SetLineWidth(2);
     Multiplicities[it->first][i].GetYaxis()->SetTitle("Multiplicity");
     Multiplicities[it->first][i].SetPointError(index,x1_error,y1_error);
    }
  }
  }
  for(std::map<std::string,std::vector<std::vector<std::vector<double> > > >::iterator it=Eff.begin();it!=Eff.end();++it)
    {
      std::string nameeff="Efficiency_"+(it->first)+"_Plate_";
      std::string namemul="Multiplicity_"+(it->first)+"_Plate_";
  for(unsigned int i=0;i!=Efficiencies[it->first].size();++i)
  {
    
    Multiplicities[it->first][i].Write((namemul+patch::to_string(i)).c_str());
    Efficiencies[it->first][i].Write((nameeff+patch::to_string(i)).c_str());
  }
  
   
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
