#include "Analysis/HitsInPlan.h"
void hitsInPlan::computeBarycentre()
{
  if (_oldDataBarycenter) return;
    for (int i=0; i<3; i++) barycentre[i]=0;
    for (std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it) 
    {
        for (int i=0; i<3; i++) 
        {
            barycentre[i]+=(*it)->getPosition()[i];
        }
    }
    if (nHits() != 0) for (int i=0; i<3; i++) barycentre[i]/=nHits();
    _oldDataBarycenter=true; 
}


void hitsInPlan::computeMaxima()
{
  if (_oldDataExtrema) return;
    for (int i=0; i<3; i++) 
    {
        min[i]=10000000;
        max[i]=-10000000;
    }
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it) 
    {
      for (int i=0; i<3; i++) 
      {
        if((*it)->getPosition()[i]<min[i])min[i]=(*it)->getPosition()[i];
        if((*it)->getPosition()[i]>max[i])max[i]=(*it)->getPosition()[i];
      }
    }
    _oldDataExtrema=true;
}
