#ifndef TRACKFITTER_H
#define TRACKFITTER_H
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "HitsInPlan.h"
#include "TGraphErrors.h"
#include <cassert>
#include <string>
#include <vector>
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

#endif
