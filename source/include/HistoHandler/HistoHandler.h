#ifndef HISTOHANDLER_H
#define HISTOHANDLER_H
#include<map>
#include<string>
#include"TH1.h"
#include<iostream>
//Singleton which will own all the Histos

class HistoHandler
{
 public:
  static HistoHandler& getInstance( );
  void RegisterHistos(std::string HistosDetector, std::string HistosChambers, std::string HistosAsics, std::string HistosPads);
  void RegisterHistos(std::string a){Histos[a]=new TH1F();
      std::cout<<Histos.size()<<std::endl;}
 private : 
  HistoHandler(){};
  ~HistoHandler(){};
  HistoHandler& operator= (const HistoHandler&){};
  HistoHandler (const HistoHandler&){};
  std::map<std::string,TH1*>Histos;
};
#endif
