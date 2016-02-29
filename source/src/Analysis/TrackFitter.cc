#include "Analysis/TrackFitter.h"
#include "Types.h"
bool trackFitter:: Find(std::vector<hitsInPlan*>& hitsByPlan,double MaxChi2,int planType,std::string collectionName)
{
  TGraphErrors grxz(hitsByPlan.size());
  TGraphErrors gryz(hitsByPlan.size());
  for (unsigned int i=0; i < hitsByPlan.size(); ++i)
      {
	hitsInPlan &hp=*(hitsByPlan[i]);
        hp.computeBarycentre();
        hp.computeMaxima(); // NB : computation results not used after
        grxz.SetPoint(i,hp.barycentreZ(),hp.barycentreX());
        if(hp.GetType()==pad)
        {
          gryz.SetPoint(i,hp.barycentreZ(),hp.barycentreY());
          gryz.SetPointError(i,hp.ErrorZ(),hp.ErrorY());
        }
        grxz.SetPointError(i,hp.ErrorZ(),hp.ErrorX());
      }
  return Find(grxz,gryz,MaxChi2,planType);

}

bool trackFitter::Find(TGraphErrors &grxz,TGraphErrors &gryz,double MaxChi2,int planType)
{
  _XZpassed=_YZpassed=false;
  _xzFit=grxz.Fit("pol1","SQRO","",-50000.0,50000.0);
  if (_xzFit->Chi2() > MaxChi2) return false;
  _XZpassed=true;
  _yzFit=gryz.Fit("pol1","SQRO","",-50000.0,50000.0);
  if (_yzFit->Chi2() > MaxChi2 && planType!=positional) return false;
  _YZpassed=true;
  return true;
}
