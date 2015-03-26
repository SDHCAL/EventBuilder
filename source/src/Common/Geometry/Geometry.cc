#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
void Geometry::PrintGeom()
{ 
    for(int i=0;i<GetNumberPlates();++i)
    {
    std::cout<<"Plane : "<<i+1<<"  "<<green<<"X : "<<GetPlatePositionX(i)<<"  Y : "<<GetPlatePositionY(i)<<"  Z : "<<GetPlatePositionZ(i)<<" Alpha : "<<GetDifPlateAlpha(i)<<" Beta : "<<GetDifPlateBeta(i)<<" Gamma : "<<GetDifPlateGamma(i)<<" SizeX : "<<GetSizeX(i)<<" SizeY : "<<GetSizeY(i)<<normal<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
     if(GetDifNbrPlate((it->first))-1==i)
     std::cout<<green<<" DifId : "<<(it->second).GetDifId()<<" I : "<<(it->second).GetPositionX()<<"  J : "<<(it->second).GetPositionY()<<" DifType : "<<GetDifTypeName((it->second).GetDifId())<<normal<<std::endl;
    }
    }
   
}
