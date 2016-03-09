#ifndef HISTOHANDLER_H
#define HISTOHANDLER_H
#include<map>
#include<string>
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TGraph.h"
#include"TGraph2D.h"
#include<iostream>
#include<vector>
#include<array>
//Singleton which will own all the Histos


class HistoHandler
{
 public:
  static HistoHandler& getInstance( );
  void RegisterHistos(std::vector<std::string>& HistosDetector,std::vector<std::string>& HistosChambers,std::vector<std::string>& HistosAsics, std::vector<std::string>& HistosPads);
 private : 
  enum Who{Detector,Chamber,Asic,Pad};
  enum What{Type,Threshold,Title,NBinsX,Xmin,Xmax,NBinsY,Ymin,Ymax,NBinsZ,Zmin,Zmax,NPoints=2};
  void Iterate(std::vector<std::string>& ToIterate);
  std::vector<std::string> Parse(std::string& ToParse);
  void CreateIt(std::vector<std::string>& tokens);
  void CreateIt(std::string& tokens,unsigned int& who);
  void Delete();
  TH1* ReturnTH1(std::string&);
  TH2* ReturnTH2(std::string&);
  TH3* ReturnTH3(std::string&);
  TGraph* ReturnTGraph(std::string&);
  TGraph2D* ReturnTGraph2D(std::string&);
  HistoHandler(){};
  ~HistoHandler()
  {
    for(std::map<std::string,TH1*>::iterator it=TH1s.begin();it!=TH1s.end();++it) delete it->second;
  for(std::map<std::string,TH2*>::iterator it=TH2s.begin();it!=TH2s.end();++it) delete it->second;
  for(std::map<std::string,TH3*>::iterator it=TH3s.begin();it!=TH3s.end();++it) delete it->second;
  for(std::map<std::string,TGraph*>::iterator it=TGraphs.begin();it!=TGraphs.end();++it) delete it->second;
  for(std::map<std::string,TGraph2D*>::iterator it=TGraph2Ds.begin();it!=TGraph2Ds.end();++it) delete it->second;
  };
  HistoHandler& operator= (const HistoHandler&)=delete;
  HistoHandler (const HistoHandler&)=delete;
  std::map<std::string,std::string>Types;
  std::map<std::string,TH1*>TH1s;
  std::map<std::string,TH2*>TH2s;
  std::map<std::string,TH3*>TH3s;
  std::map<std::string,TGraph*>TGraphs;
  std::map<std::string,TGraph2D*>TGraph2Ds;
  std::array<std::string,6>Thresholds_name{{"Threshold1","Threshold2","Threshold3","Threshold12","Threshold23","Thresholdall"}};
};
#endif
