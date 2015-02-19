#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
void Geometry::PrintGeom()
{ 
    std::cout<<"Plans : "<<std::endl;
    for(unsigned int i=0;i<GetNumberPlates();++i)
    {
   
   std::cout<<red<<"X:"<<GetPlatePositionX(i)<<"  Y:"<<GetPlatePositionY(i)<<"  Z:"<<GetPlatePositionZ(i)<<normal<<std::endl;
    }
    std::cout<<"Difs : "<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
    std::cout<<red<<"I:"<<(it->second).GetPositionX()<<"  J:"<<(it->second).GetPositionY()<<"  DifId:"<<(it->second).GetDifId()<<normal<<std::endl;
    }
}
