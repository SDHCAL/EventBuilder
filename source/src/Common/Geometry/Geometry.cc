#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
#include"Progress.h"
void Geometry::PrintGeom()
{ 
    for(int i=0;i<GetNumberPlates();++i)
    {
    std::cout<<"Plane : "<<Shift(i+1)<<"  "<<green<<"X : "<<Shift(GetPlatePositionX(i))<<"  Y : "<<Shift(GetPlatePositionY(i))<<"  Z : "<<Shift(GetPlatePositionZ(i))<<" Alpha : "<<Shift(GetDifPlateAlpha(i))<<" Beta : "<<Shift(GetDifPlateBeta(i))<<" Gamma : "<<Shift(GetDifPlateGamma(i))<<" SizeX : "<<Shift(GetSizeX(i))<<" SizeY : "<<Shift(GetSizeY(i))<<normal<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
     if(GetDifNbrPlate((it->first))-1==i)
     std::cout<<green<<" DifId : "<<Shift((it->second).GetDifId())<<" I : "<<Shift((it->second).GetPositionX())<<"  J : "<<Shift((it->second).GetPositionY())<<" DifType : "<<GetDifTypeName((it->second).GetDifId())<<normal<<std::endl;
    }
    }
   
}
