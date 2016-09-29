#ifndef CONFIG_H
#define CONFIG_H
#include <iostream>
#include "Colors.h"
#include <string>
#include <sstream>
#include <utility>
#include <map>
#include <array>


class AsicInfo
{
   public:
   AsicInfo():ASIC_TYPE(""),B0(0),B1(0),B2(0),DIF_ID(0),ID(0),HEADER(0),NUM(0)
   {
    PAGAIN.fill(0);
   };
   AsicInfo(std::string& asic_type,std::array<unsigned int,64> pagain,unsigned int& b0,unsigned int& b1,unsigned int& b2,unsigned int& dif_id,unsigned int& id, unsigned int& header,unsigned int& num)
   {
     ASIC_TYPE=asic_type;
     PAGAIN=pagain;
     B0=b0;
     B1=b1;
     B2=b2;
     DIF_ID=dif_id;
     ID=id;
     HEADER=header;
     NUM=num;
   };
   std::array<unsigned int,64> ReturnMe(){return PAGAIN;};
   ~AsicInfo(){};
   std::string getAsicType(){return ASIC_TYPE;};
   std::array<unsigned int,64> getGain(){return PAGAIN;};
   unsigned int getGain(unsigned int& index_pad){return PAGAIN[index_pad];};
   std::vector<unsigned int> getThresholds(){return std::vector<unsigned int>{B0,B1,B2};}
   unsigned int getThreshold(int& threshold)
   {
     if(threshold==0) return getThreshold0();
     if(threshold==1) return getThreshold1();
     if(threshold==2) return getThreshold2();
     else return -1;
   }
   
   unsigned int getThreshold0(){return B0;};
   unsigned int getThreshold1(){return B1;};
   unsigned int getThreshold2(){return B2;};
   unsigned int getDif_ID(){return DIF_ID;};
   unsigned int getID(){return ID;};
   unsigned int getHeader(){return HEADER;};
   unsigned int getNUM(){return NUM;};
   private:
   std::string ASIC_TYPE;
   std::array<unsigned int,64>PAGAIN;
   unsigned int B0,B1,B2;
   unsigned int DIF_ID, ID, HEADER,NUM;
};



class DifDatabaseInfo
{
   public:
   DifDatabaseInfo():NAME(""),DIF_TYPE(""),ID(-1),ENABLED(false){};
   DifDatabaseInfo(std::string& name,std::string& dif_type,unsigned int& id,bool& enabled):NAME(name),DIF_TYPE(dif_type),ID(id),ENABLED(enabled){};
   std::string getDifType(){return DIF_TYPE;};
   std::string getName(){return NAME;};
   bool getEnabled(){return ENABLED;};
   unsigned int getID(){return ID;};
   void AddAsic(AsicInfo zz)
   {
    AsicInfos.insert(std::pair<unsigned int,AsicInfo>(zz.getHeader(),zz));
   };
   AsicInfo getAsicInfo(unsigned int asic_id)
   {
      if(AsicInfos.find(asic_id)!=AsicInfos.end())return AsicInfos.find(asic_id)->second;
      else return AsicInfo();
   }
   std::map<unsigned int,AsicInfo> ReturnMe(){return AsicInfos;};
   private:
   std::string NAME, DIF_TYPE;
   unsigned int ID;
   bool ENABLED;
   std::map<unsigned int,AsicInfo> AsicInfos;
};


class ConfigInfos
{
  public:
  ConfigInfos(){};
  ~ConfigInfos(){};
  void AddElogParam(std::string key ,std::string value)
  {
    ElogParams.insert(std::pair<std::string,std::string>(key,value));
  }
  std::string ElogParam(std::string key)
  {
    if(ElogParams.find(key)!=ElogParams.end()) return ElogParams[key];
    else return "";
  }
  void AddDif(DifDatabaseInfo zz)
  {
    DifInfos.insert(std::pair<unsigned int,DifDatabaseInfo>(zz.getID(),zz)); 
  }
  void AddAsic(AsicInfo zz)
  {
    DifInfos[zz.getDif_ID()].AddAsic(zz);
  } 
  std::vector<unsigned int> getThresholds(unsigned int& dif_id,unsigned int& asic_id)
  {
    return (DifInfos.find(dif_id)->second).getAsicInfo(asic_id).getThresholds();
  }
  unsigned int getGain(unsigned int& dif_id,unsigned int& asic_id, unsigned& pad_id)
  {
     return (DifInfos.find(dif_id)->second).getAsicInfo(asic_id).getGain(pad_id);
  }
  std::map<unsigned int,DifDatabaseInfo> ReturnMe(){return DifInfos;};
  private:
  std::map<unsigned int,DifDatabaseInfo>DifInfos;
  std::map<std::string,std::string>ElogParams;

};
#endif
