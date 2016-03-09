#include "HistoHandler/HistoHandler.h"
#include<vector>
#include<string>
#include<iterator>
#include<iostream>
#include<sstream> 
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TGraph.h"
#include"TGraph2D.h"
#include<array>
#include"Colors.h"
#define  NbrChamber 2
#define  NbrAsic 3
#define  NbrPad 4

TH1* HistoHandler::ReturnTH1(std::string& name)
{
  if (TH1s.find(name)!=TH1s.end())return TH1s.find(name)->second;
  else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}
TH2* HistoHandler::ReturnTH2(std::string& name)
{
  if (TH2s.find(name)!=TH2s.end())return TH2s.find(name)->second;
  else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}
TH3* HistoHandler::ReturnTH3(std::string& name)
{
  if (TH3s.find(name)!=TH3s.end())return TH3s.find(name)->second;
  else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}
TGraph* HistoHandler::ReturnTGraph(std::string& name)
{
  if (TGraphs.find(name)!=TGraphs.end())return TGraphs.find(name)->second;
  else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}
TGraph2D* HistoHandler::ReturnTGraph2D(std::string& name )
{
  if (TGraph2Ds.find(name)!=TGraph2Ds.end())return TGraph2Ds.find(name)->second;
  else {std::cout<<"Impossible to find "<<name<<normal<<std::endl;return NULL;}
}


HistoHandler& HistoHandler::getInstance( )
{
    static HistoHandler hists;
    return hists;
}


void HistoHandler::Delete()
{
  for(std::map<std::string,TH1*>::iterator it=TH1s.begin();it!=TH1s.end();++it) delete it->second;
  for(std::map<std::string,TH2*>::iterator it=TH2s.begin();it!=TH2s.end();++it) delete it->second;
  for(std::map<std::string,TH3*>::iterator it=TH3s.begin();it!=TH3s.end();++it) delete it->second;
  for(std::map<std::string,TGraph*>::iterator it=TGraphs.begin();it!=TGraphs.end();++it) delete it->second;
  for(std::map<std::string,TGraph2D*>::iterator it=TGraph2Ds.begin();it!=TGraph2Ds.end();++it) delete it->second;
  std::exit(1);
}

std::vector<std::string> HistoHandler::Parse(std::string& ToParse)
{
  std::vector<std::string>tokens;
  std::stringstream iss(ToParse);
  std::string item;
  while(std::getline(iss, item,';'))
  {
        tokens.push_back(item);
        std::cout<<green<<item<<normal<<std::endl;
  }
  return tokens;
}

void HistoHandler::CreateIt(std::string& ToParse,unsigned int& who)
{
  std::vector<std::string> tokens=Parse(ToParse);
  std::string Original_Title=tokens[Title];
  std::cout<<yellow<<who<<normal<<std::endl;
  if(who>=Chamber)
  {
    for(unsigned int i=0;i!=NbrChamber;++i)
    {
      std::string tit=Original_Title+"_"+std::to_string(i);
      tokens[Title]=Original_Title+"_"+std::to_string(i);
      if(who>=Chamber)
      {
        for(unsigned int j=0;j!=NbrAsic;++j)
        {
          std::string ti=tit+"_"+std::to_string(j);
          tokens[Title]=tit+"_"+std::to_string(j);
          if(who>=Pad)
          {
            for(unsigned int k=0;k!=NbrPad;++k)
            {
              tokens[Title]=ti+"_"+std::to_string(k);
              std::cout<<tokens[Title]<<std::endl;
              CreateIt(tokens);
            }
          }
          else CreateIt(tokens);
        }
      }
      else CreateIt(tokens);
    }
  }
  else CreateIt(tokens);
}
void HistoHandler::CreateIt(std::vector<std::string>& tokens)
{
  //std::vector<std::string> tokens=Parse(ToParse);

  //TH1
  if(tokens[Type]=="TH1C"||tokens[Type]=="TH1S"||tokens[Type]=="TH1I"||tokens[Type]=="TH1F"||tokens[Type]=="TH1D")
  {
    if(tokens.size()!=6)
    {
      std::cout<<red<<"Missing argument(s) or too many argument(s) for "<<tokens[Title]<<normal<<std::endl;
      Delete();
    }
    else
    {
      unsigned int GenerateItForThresholds=Thresholds_name.size()-1;
      if(tokens[Threshold]=="ThrOn")GenerateItForThresholds=0;
      for(unsigned int i=GenerateItForThresholds;i!=Thresholds_name.size();++i)
      {
        std::string Titl=tokens[Title]+"_"+Thresholds_name[i];
        if(tokens[Type]=="TH1C") TH1s[Titl]=new TH1C(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]));
        else if (tokens[Type]=="TH1S") TH1s[Titl]=new TH1S(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]));
        else if (tokens[Type]=="TH1I") TH1s[Titl]=new TH1I(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]));
        else if (tokens[Type]=="TH1F") TH1s[Titl]=new TH1F(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]));
        else if (tokens[Type]=="TH1D"){TH1s[Titl]=new TH1D(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]));
        std::cout<<red<<Titl<<normal<<std::endl;
        }
        else Delete();
        Types[Titl]=tokens[Type];
      }
    }
  }
  //TH2
  if(tokens[Type]=="TH2C"||tokens[Type]=="TH2S"||tokens[Type]=="TH2I"||tokens[Type]=="TH2F"||tokens[Type]=="TH2D")
  {
    if(tokens.size()!=9)
    {
      std::cout<<red<<"Missing argument(s) or too many argument(s) for "<<tokens[Title]<<normal<<std::endl;
      Delete();
    }
    else
    {
      unsigned int GenerateItForThresholds=Thresholds_name.size()-1;
      if(tokens[Threshold]=="ThrOn")GenerateItForThresholds=0;
      for(unsigned int i=GenerateItForThresholds;i!=Thresholds_name.size();++i)
      {
        std::string Titl=tokens[Title]+"_"+Thresholds_name[i];
        if(tokens[Type]=="TH2C") TH2s[Titl]=new TH2C(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]));
        else if (tokens[Type]=="TH2S") TH2s[Titl]=new TH2S(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]));
        else if (tokens[Type]=="TH2I") TH2s[Titl]=new TH2I(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]));
        else if (tokens[Type]=="TH2F") TH2s[Titl]=new TH2F(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]));
        else if (tokens[Type]=="TH2D") TH2s[Titl]=new TH2D(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]));
        else Delete();
        Types[Titl]=tokens[Type];
      }
    }
  }
  //TH3
  if(tokens[Type]=="TH3C"||tokens[Type]=="TH3S"||tokens[Type]=="TH3I"||tokens[Type]=="TH3F"||tokens[Type]=="TH3D")
  {
    if(tokens.size()!=12)
    {
      std::cout<<red<<"Missing argument(s) or too many argument(s) for "<<tokens[Title]<<normal<<std::endl;
      Delete();
    }
    else
    {
      unsigned int GenerateItForThresholds=Thresholds_name.size()-1;
      if(tokens[Threshold]=="ThrOn")GenerateItForThresholds=0;
      for(unsigned int i=GenerateItForThresholds;i!=Thresholds_name.size();++i)
      {
        std::string Titl=tokens[Title]+"_"+Thresholds_name[i];
        if(tokens[Type]=="TH3C") TH3s[Titl]=new TH3C(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]),std::stoi(tokens[NBinsZ]),std::stof(tokens[Zmin]),std::stof(tokens[Zmax]));
        else if (tokens[Type]=="TH3S") TH3s[Titl]=new TH3S(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]),std::stoi(tokens[NBinsZ]),std::stof(tokens[Zmin]),std::stof(tokens[Zmax]));
        else if (tokens[Type]=="TH3I") TH3s[Titl]=new TH3I(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]),std::stoi(tokens[NBinsZ]),std::stof(tokens[Zmin]),std::stof(tokens[Zmax]));
        else if (tokens[Type]=="TH3F") TH3s[Titl]=new TH3F(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]),std::stoi(tokens[NBinsZ]),std::stof(tokens[Zmin]),std::stof(tokens[Zmax]));
        else if (tokens[Type]=="TH3D") TH3s[Titl]=new TH3D(Titl.c_str(),Titl.c_str(),std::stoi(tokens[NBinsX]),std::stof(tokens[Xmin]),std::stof(tokens[Xmax]),std::stoi(tokens[NBinsY]),std::stof(tokens[Ymin]),std::stof(tokens[Ymax]),std::stoi(tokens[NBinsZ]),std::stof(tokens[Zmin]),std::stof(tokens[Zmax]));
        else Delete();
        Types[Titl]=tokens[Type];
      }
    }
  }
  //TGraph
  if(tokens[Type]=="TGraph"||tokens[Type]=="TGraph2D")
  {
    if(tokens.size()==3)
    {
      unsigned int GenerateItForThresholds=Thresholds_name.size()-1;
      if(tokens[Threshold]=="ThrOn")GenerateItForThresholds=0;
      for(unsigned int i=GenerateItForThresholds;i!=Thresholds_name.size();++i)
      {
        std::string Titl=tokens[Title]+"_"+Thresholds_name[i];
        if(tokens[Type]=="TGraph") TGraphs[Titl]=new TGraph(std::stoi(tokens[NPoints]));
        else if(tokens[Type]=="TGraph2D") TGraph2Ds[Titl]=new TGraph2D(std::stoi(tokens[NPoints]));
        else Delete();
        Types[Titl]=tokens[Type];
      }
    }
    else if (tokens.size()==2)
    {
      unsigned int GenerateItForThresholds=Thresholds_name.size()-1;
      if(tokens[Threshold]=="ThrOn")GenerateItForThresholds=0;
      for(unsigned int i=GenerateItForThresholds;i!=Thresholds_name.size();++i)
      {
        std::string Titl=tokens[Title]+"_"+Thresholds_name[i];
        if(tokens[Type]=="TGraph") TGraphs[Titl]=new TGraph();
        else if(tokens[Type]=="TGraph2D") TGraph2Ds[Titl]=new TGraph2D();
        else Delete();
        Types[Titl]=tokens[Type];
      }
    }
  }
}

void HistoHandler::Iterate(std::vector<std::string>& ToIterate)
{
  static unsigned int  who=0;
  if(ToIterate.size()!=0)
  {
    for(unsigned int i=0;i!=ToIterate.size();++i)
    {
      CreateIt(ToIterate[i],who);
    }
  }
  who++;
}


void HistoHandler::RegisterHistos( std::vector<std::string>& HistosDetector,  std::vector<std::string>& HistosChambers,  std::vector<std::string>& HistosAsics,  std::vector<std::string>& HistosPads)
{
  Iterate(HistosDetector);
  Iterate(HistosChambers);
  Iterate(HistosAsics);
  Iterate(HistosPads);
}
