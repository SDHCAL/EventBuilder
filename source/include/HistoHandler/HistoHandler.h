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
  //For Detector Plot and Detector_Threshold Plot
  TH1* ReturnTH1(std::string,unsigned int Threshold);
  TH2* ReturnTH2(std::string,unsigned int Threshold);
  TH3* ReturnTH3(std::string,unsigned int Threshold);
  TGraph* ReturnTGraph(std::string,unsigned int Threshold);
  TGraph2D* ReturnTGraph2D(std::string,unsigned int Threshold);
  //For Chamber Plot and Chamber_Threshold Plot
  TH1* ReturnTH1(std::string,unsigned int Chamber,unsigned int Threshold);
  TH2* ReturnTH2(std::string,unsigned int Chamber,unsigned int Threshold);
  TH3* ReturnTH3(std::string,unsigned int Chamber,unsigned int Threshold);
  TGraph* ReturnTGraph(std::string,unsigned int Chamber,unsigned int Threshold);
  TGraph2D* ReturnTGraph2D(std::string,unsigned int Chamber,unsigned int Threshold);
  //For Asic Plot and Asic_Threshold Plot
  TH1* ReturnTH1(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Threshold);
  TH2* ReturnTH2(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Threshold);
  TH3* ReturnTH3(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Threshold);
  TGraph* ReturnTGraph(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Threshold);
  TGraph2D* ReturnTGraph2D(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Threshold);
  //For Pad PLot and Pad_Threshold Plot
  TH1* ReturnTH1(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Pad,unsigned int Threshold);
  TH2* ReturnTH2(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Pad,unsigned int Threshold);
  TH3* ReturnTH3(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Pad,unsigned int Threshold);
  TGraph* ReturnTGraph(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Pad,unsigned int Threshold);
  TGraph2D* ReturnTGraph2D(std::string,unsigned int Chamber,unsigned int Asic,unsigned int Pad,unsigned int Threshold);
  static HistoHandler& getInstance( );
  void RegisterHistos(std::vector<std::string>& HistosDetector,std::vector<std::string>& HistosChambers,std::vector<std::string>& HistosAsics, std::vector<std::string>& HistosPads);
 private : 
  //For Internal Return
  TH1* InternReturnTH1(std::string&);
  TH2* InternReturnTH2(std::string&);
  TH3* InternReturnTH3(std::string&);
  TGraph* InternReturnTGraph(std::string&);
  TGraph2D* InternReturnTGraph2D(std::string&);
  enum Who{Detector,Chamber,Asic,Pad};
  enum What{Type,Threshold,Title,NBinsX,Xmin,Xmax,NBinsY,Ymin,Ymax,NBinsZ,Zmin,Zmax,NPoints=2};
  enum Threshold{Threshold1,Threshold2,Threshold3,Threshold12,Threshold23,Thresholdall};
  void Iterate(std::vector<std::string>& ToIterate);
  std::vector<std::string> Parse(std::string& ToParse);
  void CreateIt(std::vector<std::string>& tokens);
  void CreateIt(std::string& tokens,unsigned int& who);
  void Delete();
  std::string GiveName(std::string,std::string,std::string,std::string,unsigned int);
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
  std::array<std::string,7>Thresholds_name{{"Threshold1","Threshold2","Threshold3","Threshold12","Threshold23","Thresholdall",""}};
};
#endif
