#ifndef GEOMETRY_H
#define GEOMETRY_H
#include<vector>
#include<map>
#include<string>
#include"Plate.h"
#include"Dif.h"
#include"Patch.h"
#include"Types.h"
//enum Types {pad,positional,temporal,tcherenkov,tricot,scintillator};


class Geometry
{
 public:
  void AddDif(const double& x,const double& y, int DifId,const double& xy ,const double& xz,const double& yz,const int& nbr, const int& up_down,const int& DifType)
  {
    Difs[DifId]=Dif(x,y,DifId,xy,xz,yz,nbr,up_down,DifType);
  }
  void AddPlate(const double& x,const double& y,const double& z,const double& xy ,const double& xz ,const double& yz , std::vector<int>&DifInPlate,const double& sizeX, const double& sizeY,const double& glass_type, const double& Gaz_number,const std::string& HV_channel,const std::string& Gaz_channel)
  {
    Plates.push_back(Plate(x,y,z,xy,xz,yz,DifInPlate,sizeX,sizeY,glass_type,Gaz_number,HV_channel,Gaz_channel));
    
  }
  inline const double GetDifPositionX( int& i){ return ((Difs.find(i))->second).GetPositionX();};
  inline const double GetDifPositionY( int& i){ return ((Difs.find(i))->second).GetPositionY();};
  inline const double GetDifAlpha( const int& i){ return ((Difs.find(i))->second).GetAngleXY();};
  inline const double GetDifBeta( const int& i){ return ((Difs.find(i))->second).GetAngleXZ();};
  inline const double GetDifGamma( const int& i){ return ((Difs.find(i))->second).GetAngleYZ();};
  inline const double GetDifPositionXMax(const int& i )
  {
    unsigned currentMax = 0;
    for(std::map<int, Dif >::iterator it = Difs.begin(); it != Difs.end(); ++it ) 
    {
      if ((it->second).GetNbrPlate()==i && (it ->second).GetPositionX() > currentMax) 
      {
        currentMax = (it ->second).GetPositionX();
      }
    } 
    return currentMax;
  };
  inline const double GetDifPositionYMax(const int& i )
  {
    unsigned currentMax = 0;
    for(std::map<int, Dif >::iterator it = Difs.begin(); it != Difs.end(); ++it ) 
    {
      if ((it->second).GetNbrPlate()==i && (it ->second).GetPositionY() > currentMax) 
      {
        currentMax = (it ->second).GetPositionY();
      }
    } 
    return currentMax;
  };
  inline const double GetDifPlateAlpha( const int& i){ return Plates[i].GetAngleXY();};
  inline const double GetDifPlateBeta( const int& i){ return Plates[i].GetAngleXZ();};
  inline const double GetDifPlateGamma( const int& i){ return Plates[i].GetAngleYZ();};
  inline const int GetDifId( int& i){ return ((Difs.find(i))->second).GetDifId();};
  inline const int GetNbrDifInPlate(int& i){return  Plates[i].GetNbrDifInPlate();};
  inline const double GetPlatePositionX(const unsigned int& i){ return Plates[i].GetPositionX();};
  inline const double GetSizeX(const unsigned int& i){ return Plates[i].GetSizeX();};
  inline const double GetSizeY(const unsigned int& i){ return Plates[i].GetSizeY();};
  inline const double GetPlatePositionY(const unsigned int& i){ return Plates[i].GetPositionY();};
  inline const double GetPlatePositionZ(const unsigned int& i){ return Plates[i].GetPositionZ();};
  inline const int GetDifNbrPlate( const int& i){ if(Difs.find(i)==Difs.end())return -1; else return ((Difs.find(i))->second).GetNbrPlate();};
  inline const int GetDifType( const int& i){ return ((Difs.find(i))->second).GetDifType();};
  inline const std::string GetDifTypeName( const int& i){std::vector<std::string>Types_names{"Pad","Positional","Temporal","Tcherenkov","Tricot","Scintillator"};if(Difs.find(i)==Difs.end())return "Unknown Difs";else return Types_names[((Difs.find(i))->second).GetDifType()];};
  inline const int GetDifUpDown( int& i){ return ((Difs.find(i))->second).GetDifUpDown();};
  inline const unsigned int GetNumberDifs(){return Difs.size();};
  inline const unsigned int GetNumberPlates(){return Plates.size();};
  inline const std::vector<Plate> GetPlates(){return Plates;};
  inline const Plate GetPlate(const int& i){return Plates[i];};
  inline const std::map<int,Dif> GetDifs(){return Difs;};
  inline const std::vector<int> GetDifsInPlane(unsigned int& i){return Plates[i].GetDifInPlate();};
  inline double GetGlassType(const unsigned int& i){return Plates[i].GetGlassType();};
  inline double GetGazNumber(const unsigned int& i){return Plates[i].GetGazNumber();};
  inline std::string GetHVChannel(const unsigned int& i){return Plates[i].GetHVChannel();};
  inline std::string GetGazChannel(const unsigned int& i){return Plates[i].GetGazChannel();};
  inline int GetXMax()
  {
    int max=-999999;
    for(unsigned int i=0;i!=Plates.size();++i)
    {
      if(Plates[i].GetSizeX()>max) max=Plates[i].GetSizeX();
    }
    return max;
  }
  inline int GetYMax()
  {
    int max=-999999;
    for(unsigned int i=0;i!=Plates.size();++i)
    {
      if(Plates[i].GetSizeY()>max) max=Plates[i].GetSizeY();
    }
    return max;
  }
  void PrintGeom(); 
 private:
  std::vector<Plate> Plates;
  std::map<int, Dif > Difs;
};
#endif
