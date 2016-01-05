#include"Geometry/Geometry.h"
#include"Colors.h"
#include<iostream>
#include<vector>
#include<string>
#include"Progress.h"
std::vector<std::string>glass_types{"standard","chinese"};
void Geometry::PrintGeom()
{ 
    for(int i=0;i<GetNumberPlates();++i)
    {
      std::cout<<"Plane : "<<Shift(i+1)<<"  "<<green<<"X : "<<Shift(GetPlatePositionX(i))<<"  Y : "<<Shift(GetPlatePositionY(i))<<"  Z : "<<Shift(GetPlatePositionZ(i))<<" Alpha : "<<Shift(GetDifPlateAlpha(i))<<" Beta : "<<Shift(GetDifPlateBeta(i))<<" Gamma : "<<Shift(GetDifPlateGamma(i))<<" SizeX : "<<Shift(GetSizeX(i))<<" SizeY : "<<Shift(GetSizeY(i))<<normal<<std::endl;
    if(GetGlassType(i)!=-1)std::cout<<green<<" Glass type : "<<glass_types[GetGlassType(i)]<<" "<<normal;
    else std::cout<<red<<" Glass type : UNKNOW "<<normal;
    if(GetHVChannel(i)!="")std::cout<<green<<" HV_Channel : "<<GetHVChannel(i)<<" "<<normal;
    if(GetGazChannel(i)!="")std::cout<<green<<" Gaz_Channel : "<<GetGazChannel(i)<<" "<<normal;
    if(GetGazNumber(i)!=-1)std::cout<<green<<" ,Position in the gaz circuit : "<<GetGazNumber(i)<<" "<<normal;
    std::cout<<std::endl;
    for(std::map<int,Dif>::iterator it=Difs.begin();it!=Difs.end();++it)
    {
     if(GetDifNbrPlate((it->first))-1==i)
     std::cout<<green<<" DifId : "<<Shift((it->second).GetDifId())<<" I : "<<Shift((it->second).GetPositionX())<<"  J : "<<Shift((it->second).GetPositionY())<<" DifType : "<<GetDifTypeName((it->second).GetDifId())<<normal<<std::endl;
    }
    }
   
}
