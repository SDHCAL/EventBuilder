#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
#include"Utilities.h"
void Geometry::PrintGeom()
{ 
    for(int i=0;i<GetNumberPlates();++i)
    {
    std::cout<<"Plane : "<<Shift(i+1)<<i+1<<"  "<<green<<"X : "<<Shift(GetPlatePositionX(i))<<GetPlatePositionX(i)<<"  Y : "<<Shift(GetPlatePositionY(i))<<GetPlatePositionY(i)<<"  Z : "<<Shift(GetPlatePositionZ(i))<<GetPlatePositionZ(i)<<" Alpha : "<<Shift(GetDifPlateAlpha(i))<<GetDifPlateAlpha(i)<<" Beta : "<<Shift(GetDifPlateBeta(i))<<GetDifPlateBeta(i)<<" Gamma : "<<Shift(GetDifPlateGamma(i))<<GetDifPlateGamma(i)<<" SizeX : "<<Shift(GetSizeX(i))<<GetSizeX(i)<<" SizeY : "<<Shift(GetSizeY(i))<<GetSizeY(i)<<normal<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
     if(GetDifNbrPlate((it->first))-1==i)
     std::cout<<green<<" DifId : "<<Shift((it->second).GetDifId())<<(it->second).GetDifId()<<" I : "<<Shift((it->second).GetPositionX())<<(it->second).GetPositionX()<<"  J : "<<Shift((it->second).GetPositionY())<<(it->second).GetPositionY()<<" DifType : "<<GetDifTypeName((it->second).GetDifId())<<normal<<std::endl;
    }
    }
   
}
