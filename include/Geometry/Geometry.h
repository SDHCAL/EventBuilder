#ifndef GEOMETRY_H
#define GEOMETRY_H
#include<vector>
#include<map>
#include"Plate.h"
#include"Dif.h"
enum Types {pad,positional,temporal,tcherenkov};
class Geometry
{
public:
void AddDif(const double& x,const double& y, int DifId,const double& xy ,const double& xz,const double yz,const int& nbr, const int& up_down,const int& DifType)
{
Difs[DifId]=Dif(x,y,DifId,xy,xz,yz,nbr,up_down,DifType);
}
void AddPlate(const double& x,const double& y,const double& z,const double& xy ,const double& xz ,const double& yz , std::vector<int>DifInPlate)
{
Plates.push_back(Plate(x,y,z,xy,xz,yz,DifInPlate));
}
inline double GetDifPositionX( int& i){ return ((Difs.find(i))->second).GetPositionX();};
inline double GetDifPositionY( int& i){ return ((Difs.find(i))->second).GetPositionY();};
inline double GetDifAlpha( const int& i){ return ((Difs.find(i))->second).GetAngleXY();};
inline double GetDifBeta( const int& i){ return ((Difs.find(i))->second).GetAngleXZ();};
inline double GetDifGamma( const int& i){ return ((Difs.find(i))->second).GetAngleYZ();};
 inline double GetDifPositionXMax(const int& i ){unsigned currentMax = 0;for(std::map<int, Dif >::iterator it = Difs.begin(); it != Difs.end(); ++it ) {if ((it->second).GetNbrPlate()==i && (it ->second).GetPositionX() > currentMax) {currentMax = (it ->second).GetPositionX();}} return currentMax;};
inline double GetDifPositionYMax(const int& i ){unsigned currentMax = 0;for(std::map<int, Dif >::iterator it = Difs.begin(); it != Difs.end(); ++it ) {if ((it->second).GetNbrPlate()==i && (it ->second).GetPositionY() > currentMax) {currentMax = (it ->second).GetPositionY();}} return currentMax;};
inline double GetDifPlateAlpha( const int& i){ return Plates[i].GetAngleXY();};
inline double GetDifPlateBeta( const int& i){ return Plates[i].GetAngleXZ();};
inline double GetDifPlateGamma( const int& i){ return Plates[i].GetAngleYZ();};
inline double GetDifId( int& i){ return ((Difs.find(i))->second).GetDifId();};
inline int GetNbrDifInPlate(int& i){return  Plates[i].GetNbrDifInPlate();};
inline double GetPlatePositionX(const unsigned int& i){ return Plates[i].GetPositionX();};
inline double GetPlatePositionY(const unsigned int& i){ return Plates[i].GetPositionY();};
inline double GetPlatePositionZ(const unsigned int& i){ return Plates[i].GetPositionZ();};
inline int GetDifNbrPlate( const int& i){ return ((Difs.find(i))->second).GetNbrPlate();};
inline int GetDifType( const int& i){ return ((Difs.find(i))->second).GetDifType();};
inline int GetDifUpDown( int& i){ return ((Difs.find(i))->second).GetDifUpDown();};
inline unsigned int GetNumberDifs(){return Difs.size();};
inline unsigned int GetNumberPlates(){return Plates.size();};
inline std::vector<Plate> GetPlates(){return Plates;};
inline Plate GetPlate(int i){return Plates[i];};
inline std::map<int,Dif> GetDifs(){return Difs;};
void PrintGeom(); 
private:
std::vector<Plate> Plates;
std::map<int, Dif > Difs;
};
#endif
