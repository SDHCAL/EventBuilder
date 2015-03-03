#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
void Geometry::PrintGeom()
{ 
    for(int i=0;i<GetNumberPlates();++i)
    {
    std::cout<<"Plans : "<<i<<"  "<<red<<"X:"<<GetPlatePositionX(i)<<"  Y:"<<GetPlatePositionY(i)<<"  Z:"<<GetPlatePositionZ(i)<<normal<<std::endl;
    std::cout<<"Difs : "<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
     if(GetDifNbrPlate((it->first))-1==i)
     std::cout<<red<<"I:"<<(it->second).GetPositionX()<<"  J:"<<(it->second).GetPositionY()<<"  DifId:"<<(it->second).GetDifId()<<" DifType : "<<GetDifTypeName((it->second).GetDifId())<<normal<<std::endl;
    }
    }
   
}
