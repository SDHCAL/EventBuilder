#ifndef ANALYSIS_H
#define ANALYSIS_H
// marlin
#include "marlin/Processor.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "Geometry/Geometry.h"
#include <cmath>
#include "UTIL/CellIDDecoder.h"
#include <map>
#define degtorad 0.0174532925
class plan;
class testedPlan
{
public:
    testedPlan(int numeroPlan,float x ,float y,float z,float xy, float xz, float yz, int type,double _Ip,double _Im, double _Jp, double _Jm) : Nbr(numeroPlan),X0(x),Y0(y),Z0(z),XY(xy),XZ(xz),YZ(yz),Type(type),Ip(_Ip),Im(_Im),Jp(_Jp),Jm(_Jm), nombreTests(0), nombreTestsOK(0), sommeNombreHits(0)
    {
        for (int i=0; i<NCOUNTERS; i++) counts[i]=0;
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
    inline double multiplicityShort() {return sommeNombreHitsShort/nombreTestsOKShort;}
    inline double GetNumberOKShort(){return nombreTestsOKShort;}
    inline void ClearShort(){nombreTestsShort=0;nombreTestsOKShort=0;sommeNombreHitsShort=0;}
    inline double efficiencyShort() {return nombreTestsOKShort/nombreTestsShort;}
    inline double efficiency()
    {
        return nombreTestsOK/nombreTests;
    }
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
    inline double GetNumberOK()
    {
        return nombreTestsOK;
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
    inline double multiplicity()
    {
        return sommeNombreHits/nombreTestsOK;
    }
    inline void clear()
    {
        nombreTests=nombreTestsOK=sommeNombreHits=0;
    }
    void testYou(std::map<int,plan>& mapDIFplan);
    void print();
private:
    int Nbr;
    float X0,Y0,Z0,XY,XZ,YZ,xnorm,ynorm,znorm,xi,yi,zi,xj,yj,zj,ca,sa,cb,sb,cg,sg;
    int Type;
    double Ip;
    double Im;
    double Jp;
    double Jm;
    double nombreTests;
    double nombreTestsOK;
    double sommeNombreHits;
    double nombreTestsShort;
    double nombreTestsOKShort;
    double sommeNombreHitsShort;
    enum {TESTYOUCALLED,NOTOOMUCHHITSINPLAN,XZTRACKFITPASSED,YZTRACKFITPASSED,NOHITINPLAN, NCOUNTERS};
    double counts[NCOUNTERS];
};





class plan
{
public:
    plan()
    {
      ;
    }
    inline void addHit(CalorimeterHit* a)
    {
        hits.push_back(a);
    }
    inline int nHits()
    {
        return hits.size();
    }
    inline void computeBarycentre();
    inline double barycentreX()
    {
        return barycentre[0];
    }
    inline double barycentreY()
    {
        return barycentre[1];
    }
    inline double barycentreZ()
    {
        return barycentre[2];
    }
    inline void computeMaxima();
    inline double minX()
    {
        return min[0];
    }
    inline void SetType(int i )
    {
        _type=i;
    }
    inline double minY()
    {
        return min[1];
    }
    inline double minZ()
    {
        return min[2];
    }
    inline double maxX()
    {
        return max[0];
    }
    inline double maxY()
    {
        return max[1];
    }
    inline double maxZ()
    {
        return max[2];
    }
    inline double ErrorX()
    {
        double e=fabs(maxX()-minX());
        return (e<1 ? 1 : e);
    }
    inline double ErrorY()
    {
        double e=fabs(maxY()-minY());
        if(this->GetType()==pad)return (e<1 ? 1 : e);
        else return 0;
    }
    inline double ErrorZ()
    {
        return 2;
    }
    inline bool operator==(plan b);
    inline bool operator!=(plan b);
    void printBarycentre();
    void printMaxima();
    inline void Clear()
    {
        hits.clear();
    }
    inline int GetType()
    {
        return _type;
    }
    inline int countHitAt(double& x, double& y, double dlim);
    inline int countHitAtStrip(double& x, double dlim);
private:
    int _type;
    std::vector<CalorimeterHit*> hits;
    double barycentre[3];
    double min[3];
    double max[3];
};



double _Chi2;
unsigned int _eventNr;
int _ShortEfficiency;
int _NbrHitPerPlaneMax ;
int _NbrPlaneUseForTracking ;
double _dlimforPad;
double _dlimforStrip;
std::map<int ,std::vector<double> >Delimiter;
std::string _Delimiters;

class AnalysisProcessor : public marlin::Processor
{
public:
    AnalysisProcessor();
    ~AnalysisProcessor();
    void PrintStat();
    void PrintStatShort();
    marlin::Processor *newProcessor(){return new AnalysisProcessor();}
    void init();
    void processEvent(EVENT::LCEvent *evtP);
    void processRunHeader( LCRunHeader* run);
    void end();
protected:
    std::vector<std::string> _hcalCollections;
    LCWriter* _EventWriter;
    std::string _FileNameGeometry;
    int _eventNr;
    int _NbrRun;

    Geometry geom;
    std::string _ReaderType;
    std::map<int,plan>Plans;
    std::vector<testedPlan> testedPlanList;
};
#endif
