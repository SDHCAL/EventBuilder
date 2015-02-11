#include"../include/Geometry.h"
#include"../../Colors.h"
#include<iostream>
void Geometry::PrintGeom()
{ 
    std::cout<<"Plans : "<<std::endl;
    for(unsigned int i=0;i<GetNumberPlates();++i)
    {
   
    std::cout<<red<<GetPlatePositionX(i)<<"  "<<GetPlatePositionY(i)<<"  "<<GetPlatePositionZ(i)<<normal<<std::endl;
    }
    std::cout<<"Difs : "<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
    std::cout<<red<<(it->second).GetPositionX()<<"  "<<(it->second).GetPositionY()<<"  "<<(it->second).GetDifId()<<normal<<std::endl;
    }
}
