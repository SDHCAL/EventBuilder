#ifndef TESTEDPLAN_H
#define TESTEDPLAN_H
#include "GeometryPlan.h"
#include "Geometry/Geometry.h"
#include "HitsInPlan.h"
#include "Colors.h"
#include <string>
#include <map>
#include <vector>
#include <array>
#include "HistoHandler/HistoHandler.h"
class testedPlan
{

   HistoHandler& a =HistoHandler::getInstance( );
public:
    void print(std::string name);
    testedPlan(geometryplan geo):geomplan(geo)
    {
        Counts.clear();
        nombreTests.clear();
        nombreTestsShort.clear();
        nombreTestsOK.clear();
        nombreTestsOKShort.clear();
        sommeNombreHits.clear();
        sommeNombreHitsShort.clear();
    }
    inline float GetZexp(const double & pxz0,const double & pyz0,const double & pxz1,const double& pyz1)
    {
        return geomplan.GetZexp(pxz0,pyz0,pxz1,pyz1);
    }
    inline float GetProjectioni(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return geomplan.GetProjectioni(Xexp,Yexp,Zexp);
    }
    inline float GetProjectionj(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return geomplan.GetProjectionj(Xexp,Yexp,Zexp);
    }
    inline void ClearShort(std::string name)
    {
      nombreTestsShort.clear();
      nombreTestsOKShort.clear();
      sommeNombreHitsShort.clear();
    }
        inline int NbrPlate()
    {
        return geomplan.NbrPlate();
    }
    inline float GetXY()
    {
        return geomplan.GetXY();
    }
    inline float GetXZ()
    {
        return geomplan.GetXZ();
    }
    inline float GetYZ()
    {
        return geomplan.GetYZ();
    }
    inline float GetX0()
    {
        return geomplan.GetX0();
    }
    inline float GetY0()
    {
        return geomplan.GetY0();
    }
    inline float GetZ0()
    {
        return geomplan.GetZ0();
    }
    inline int GetType()
    {
        return geomplan.GetType();
    }
    inline double GetIp()
    {
        return geomplan.GetIp();
    }
    inline double GetIm()
    {
        return geomplan.GetIm();
    }
    inline double GetJp()
    {
        return geomplan.GetJp();
    }
    inline double GetJm()
    {
        return geomplan.GetJm();
    }
    inline void clear()
    {
        nombreTests.clear();
        nombreTestsOK.clear();
        sommeNombreHits.clear(); 
    }
   
    inline double multiplicityShort(int i, std::string name)
    {
      if(sommeNombreHitsShort.find(name)==sommeNombreHitsShort.end()||nombreTestsOKShort.find(name)==nombreTestsOKShort.end())return -1;
      else return  sommeNombreHitsShort[name][i]/nombreTestsOKShort[name][i];
    }
    inline double GetNumberOKShort(int i,std::string name )
    {
      if(nombreTestsOKShort.find(name)==nombreTestsOKShort.end())return -1;
      else return nombreTestsOKShort[name][i];
    }
    inline double efficiencyShort(int i, std::string name)
    { 
      if(nombreTestsOKShort.find(name)==nombreTestsOKShort.end()||nombreTestsShort.find(name)==nombreTestsShort.end())return -1;
      else return nombreTestsOKShort[name][i]/nombreTestsShort[name];
    }
    inline double efficiency(int i,std::string name)
    {
        if(nombreTestsOK.find(name)==nombreTestsOK.end()||nombreTests.find(name)==nombreTests.end())return -1;
        else return nombreTestsOK[name][i]/nombreTests[name];
    }
    inline double GetNumberOK(int i,std::string name)
    {
       if(nombreTestsOK.find(name)==nombreTestsOK.end())return -1;
       else return nombreTestsOK[name][i];
    }
    inline double multiplicity(int i,std::string name)
    {
        if(sommeNombreHits.find(name)==sommeNombreHits.end()||sommeNombreHits.find(name)==sommeNombreHits.end())return -1; 
        else return sommeNombreHits[name][i]/nombreTestsOK[name][i];
    }
     inline double GetNombreHits(int i,std::string name)
    {
        if(sommeNombreHits.find(name)==sommeNombreHits.end())
        {
          std::cout<<red<<name<<" unknown"<<std::endl;
          return -1; 
        }
        else return sommeNombreHits[name][i];
    }
    
    
    inline float get_ca(){return geomplan.get_ca();};
    inline float get_sa(){return geomplan.get_sa();};
    inline float get_cb(){return geomplan.get_cb();};
    inline float get_sb(){return geomplan.get_sb();};
    inline float get_cg(){return geomplan.get_cg();};
    inline float get_sg(){return geomplan.get_sg();};
    void testYou(std::map<std::string,std::map<int,hitsInPlan>>&mapDIFplan,std::vector<testedPlan>& tested);
    
    private:
    
    std::map<std::string,std::array<double,5>>Counts;
    std::map<std::string,double>nombreTests;
    std::map<std::string,double>nombreTestsShort;
    std::map<std::string,std::array<double,6>>nombreTestsOK;
    std::map<std::string,std::array<double,6>>nombreTestsOKShort;
    std::map<std::string,std::array<double,6>>sommeNombreHits;
    std::map<std::string,std::array<double,6>>sommeNombreHitsShort;
    geometryplan geomplan;
    enum {TESTYOUCALLED,NOTOOMUCHHITSINPLAN,XZTRACKFITPASSED,YZTRACKFITPASSED,NOHITINPLAN, NCOUNTERS};
};
#endif
