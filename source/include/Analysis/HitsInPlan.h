#ifndef HITSINPLAN_H
#define HITSINPLAN_H
#include "EVENT/CalorimeterHit.h"
#include <vector>
#include <string>
#include <array>
#include "Types.h"
class hitsInPlan
{
 public:
 hitsInPlan() : _oldDataBarycenter(false) ,_oldDataExtrema(false) {;}

  //modifier
  inline void addHit(EVENT::CalorimeterHit* a) { hits.push_back(a);_oldDataBarycenter=_oldDataExtrema=false;}
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

  std::vector<EVENT::CalorimeterHit*>& GetHits() { return hits; }

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
  std::vector<EVENT::CalorimeterHit*> hits;
  double barycentre[3];
  double min[3];
  double max[3];
  
  bool _oldDataBarycenter,_oldDataExtrema;

};
#endif
