#ifndef ANALYSIS_H
#define ANALYSIS_H
// marlin
#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "Geometry/Geometry.h"
#include <cmath>
#include <assert.h>
#include "UTIL/CellIDDecoder.h"
#include <map>
#include <array>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <vector>
#include "Config/Config.h"
#define degtorad 0.0174532925

enum Threshold{Threshold_2=1,Threshold_1,Threshold_3};
class hitsInPlan;


class geometryplan
{
  public:
  geometryplan(){};
  ~geometryplan(){};
  geometryplan(int numeroPlan,float x ,float y,float z,float xy, float xz, float yz, int type,double _Ip,double _Im, double _Jp, double _Jm) : Nbr(numeroPlan),X0(x),Y0(y),Z0(z),XY(xy),XZ(xz),YZ(yz),Type(type),Ip(_Ip),Im(_Im),Jp(_Jp),Jm(_Jm)
  {
        ca=cos(xy*degtorad);
        sa=sin(xy*degtorad);
        cb=cos(xz*degtorad);
        sb=sin(xz*degtorad);
        cg=cos(yz*degtorad);
        sg=sin(yz*degtorad);
          /*****************************/
        /* Xexp=pxz0+pxz1*Zexp
           Yexp = pyz0+pyz1*Zexp
           Zexp=Zexp
           (xnorm,ynorm,znorm).(Xexp-X0,Yexp-Y0,Zexp-z0)=0
        */
        xnorm=sg*sa+cg*sb*ca;
        ynorm=-cg*sa+sg*sb*ca;
        znorm=cb*ca;
        xi=cg*cb;
        yi=sg*cb;
        zi=-sb;
        xj=-sg*ca+cg*sb*sa;
        yj=cg*ca+sg*sb*sa;
        zj=cb*sa;
    }
    geometryplan(const int& numeroPlan,const int& type,Geometry& geom,std::map<int,std::vector<double>>& Delimiter)
    {
      Nbr=numeroPlan;
      X0=geom.GetPlatePositionX(numeroPlan);
	    Y0=geom.GetPlatePositionY(numeroPlan);
	    Z0=geom.GetPlatePositionZ(numeroPlan);
	    XY=geom.GetDifPlateAlpha(numeroPlan);
	    XZ=geom.GetDifPlateBeta(numeroPlan);
	    YZ=geom.GetDifPlateGamma(numeroPlan);
	    Type=type;
	    ca=cos(XY*degtorad);
      sa=sin(XY*degtorad);
      cb=cos(XZ*degtorad);
      sb=sin(XZ*degtorad);
      cg=cos(YZ*degtorad);
      sg=sin(YZ*degtorad);
      /*****************************/
      /* Xexp=pxz0+pxz1*Zexp
           Yexp = pyz0+pyz1*Zexp
           Zexp=Zexp
           (xnorm,ynorm,znorm).(Xexp-X0,Yexp-Y0,Zexp-z0)=0
      */
      xnorm=sg*sa+cg*sb*ca;
      ynorm=-cg*sa+sg*sb*ca;
      znorm=cb*ca;
      xi=cg*cb;
      yi=sg*cb;
      zi=-sb;
      xj=-sg*ca+cg*sb*sa;
      yj=cg*ca+sg*sb*sa;
      zj=cb*sa;
	    if(Delimiter.size()==0)
	    {
	      Ip=0;
	      Im=0;
	      Jp=0;
	      Jm=0;
	    }
	    else
	    {
	      Ip=Delimiter[numeroPlan+1][1];
	      Im=Delimiter[numeroPlan+1][0];
	      Jp=Delimiter[numeroPlan+1][3];
	      Jm=Delimiter[numeroPlan+1][2];
	    } 
    }
    inline float GetZexp(const double & pxz0,const double & pyz0,const double & pxz1,const double& pyz1)
    {
        return (xnorm*(pxz0-X0)+ynorm*(pyz0-Y0)+znorm*Z0)/(xnorm*pxz1+ynorm*pyz1+znorm);
    }
    inline float GetProjectioni(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return (Xexp-X0)*xi+(Yexp-Y0)*yi+(Zexp-Z0)*zi;
    }
    inline float GetProjectionj(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return (Xexp-X0)*xj+(Yexp-Y0)*yj+(Zexp-Z0)*zj;
    }
    inline float get_ca(){return ca;};
    inline float get_sa(){return sa;};
    inline float get_cb(){return cb;};
    inline float get_sb(){return sb;};
    inline float get_cg(){return cg;};
    inline float get_sg(){return sg;};
     inline int NbrPlate()
    {
        return Nbr;
    }
    inline float GetXY()
    {
        return XY;
    }
    inline float GetXZ()
    {
        return XZ;
    }
    inline float GetYZ()
    {
        return YZ;
    }
    inline float GetX0()
    {
        return X0;
    }
    inline float GetY0()
    {
        return Y0;
    }
    inline float GetZ0()
    {
        return Z0;
    }
    inline int GetType()
    {
        return Type;
    }
    inline double GetIp()
    {
        return Ip;
    }
    inline double GetIm()
    {
        return Im;
    }
    inline double GetJp()
    {
        return Jp;
    }
    inline double GetJm()
    {
        return Jm;
    }
    private:
    int  Nbr;
    float X0,Y0,Z0,XY,XZ,YZ,xnorm,ynorm,znorm,xi,yi,zi,xj,yj,zj,ca,sa,cb,sb,cg,sg;
    int Type;
    double Ip;
    double Im;
    double Jp;
    double Jm;

};




class trackFitter
{
 public:
 trackFitter() : _XZpassed(false),_YZpassed(false)  {}
  bool Find(std::vector<hitsInPlan*>& hitsByPlan,double MaxChi2,int planType,std::string collectionName);
  bool Find(TGraphErrors &grxz,TGraphErrors &gryz,double MaxChi2,int planType);
  double getXZParameter(int i) {assert(_XZpassed); return _xzFit->Parameter(i);}
  double getYZParameter(int i) {assert(_YZpassed); return _yzFit->Parameter(i);}
  double getXZChisquare() {assert(_XZpassed); return _xzFit->Chi2();}
  double getYZChisquare() {assert(_YZpassed); return _yzFit->Chi2();}

  bool foundTracksXZ() {return _XZpassed;}
  bool foundTracksYZ() {return _YZpassed;}

 private:
  bool _XZpassed,_YZpassed;
  TFitResultPtr _xzFit,_yzFit;
};


class testedPlan
{

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



class hitsInPlan
{
 public:
 hitsInPlan() : _oldDataBarycenter(false) ,_oldDataExtrema(false) {;}

  //modifier
  inline void addHit(CalorimeterHit* a) { hits.push_back(a);_oldDataBarycenter=_oldDataExtrema=false;}
  inline void Clear() { hits.clear(); _oldDataBarycenter=_oldDataExtrema=false;}
  inline void SetType(int i ) { _type=i; }
  void computeBarycentre();
  void computeMaxima();


  //getters
  inline unsigned int nHits() { return hits.size();}
  inline int GetType() { return _type; }
  inline double barycentreX() { return barycentre[0];}
  inline double barycentreY() { return barycentre[1];}
  inline double barycentreZ() { return barycentre[2];}
  inline double minX() { return min[0];}
  inline double minY() { return min[1];}
  inline double minZ() { return min[2];}
  inline double maxX() { return max[0];}
  inline double maxY() { return max[1];}
  inline double maxZ() { return max[2];}

  std::vector<CalorimeterHit*>& GetHits() { return hits; }

  //analysis
  std::array<double,6> countHitAt(double& x, double& y, double dlim,int Xexpected,int Yexpected,int Kexpected,double Imin,double Imax,double Jmin,double Jmax,std::string collectionName);
  std::map<std::string,int> countHitAtStrip(double& x, double dlim,std::string collectionName);



  //Error functions
  inline double ErrorX()
  {
        //double e=fabs(maxX()-minX());
        //return (e<1 ? 1 : e);
        return 10;
  }
  inline double ErrorY()
  {
    //double e=fabs(maxY()-minY());
    //if(this->GetType()==pad)return (e<1 ? 1 : e);
    //else return 0;
    //return 2.5;
    return 10;
    if(this->GetType()==pad)return 0;
    
  }
  inline double ErrorZ()
  {
    return 10;
  }
  
  
 private:
  int _type;
  std::vector<CalorimeterHit*> hits;
  double barycentre[3];
  double min[3];
  double max[3];
  
  bool _oldDataBarycenter,_oldDataExtrema;

};




std::map<std::string,std::array<std::map<std::vector<int>,std::vector<int>>,6>>Efficiency_per_pad;
double _Chi2;
double _Chi2Rate;
int _ShortEfficiency;
int _NbrHitPerPlaneMax ;
int _NbrPlaneUseForTracking ;
int _NbrPlaneUseForTrackingRate ;
double _dlimforPad;
int _NbrRun;
double _dlimforStrip;
bool EstimateNoiseContamination;
std::map<int ,std::vector<double> >Delimiter;
std::string _Delimiters;
std::vector<std::string> _hcalCollections;
class AnalysisProcessor : public marlin::Processor
{
public:
    AnalysisProcessor();
    ~AnalysisProcessor();
    void PrintStat(std::string);
    void PrintStatShort();
    void PrintStatShort(std::string);
    void operator=(const AnalysisProcessor&);
    virtual marlin::Processor* newProcessor() {
      return new AnalysisProcessor;
    }
    void init();
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
protected:
 
    LCWriter* _EventWriter;
    std::string _FileNameGeometry;
    unsigned int _eventNr;
    unsigned int _skip;
    unsigned int _maxRecord;
    unsigned int _GlobalEvents;
    unsigned int _GlobalEventsSc;
    
    std::string _Config_xml;   

    ConfigInfos conf;
    Geometry geom;
    std::string _ReaderType;
    std::map<std::string,std::map<int,hitsInPlan>> Planss;
    std::map<int,geometryplan> geometryplans;
    std::vector<testedPlan> testedPlanList;

};
#endif
