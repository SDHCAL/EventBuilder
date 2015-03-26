#ifndef PLATE_H
#define PLATE_H
#include<vector>
class Plate
{
public:
Plate(){};
Plate(const Plate& _Plate):_PositionX(_Plate._PositionX),_PositionY(_Plate._PositionY),_PositionZ(_Plate._PositionZ),_AngleXY(_Plate._AngleXY),_AngleXZ(_Plate._AngleXZ),_AngleYZ(_Plate._AngleYZ),_DifInPlate(_Plate._DifInPlate),_SizeX(_Plate._SizeX),_SizeY(_Plate._SizeY){};
Plate(const double& x,const double& y,const double& z, const double& xy,const double& xz,const double& yz,std::vector<int>DifInPlate,const double& sizeX, const double& sizeY):_PositionX(x),_PositionY(y),_PositionZ(z),_AngleXY(xy),_AngleXZ(xz),_AngleYZ(yz),_DifInPlate(DifInPlate),_SizeX(sizeX),_SizeY(sizeY){};
//void SetType(const std::string& Type);
//void AddDif(const int& NbrDif);
inline double GetPositionX(){return _PositionX;};
inline double GetPositionY(){return _PositionY;};
inline double GetPositionZ(){return _PositionZ;};
inline double GetAngleXY(){return _AngleXY;};
inline double GetAngleXZ(){return _AngleXZ;};
inline double GetAngleYZ(){return _AngleYZ;};
inline int GetNbrDifInPlate(){return _DifInPlate.size();};
inline std::vector<int> GetDifInPlate(){return _DifInPlate;};
inline double GetSizeX(){return _SizeX;};
inline double GetSizeY(){return _SizeY;};
//double GetType(const std::string& Type);
//Dif SetType(const int& NbrDif);

private:
double _PositionX;
double _PositionY;
double _PositionZ;
double _AngleXY;
double _AngleXZ;
double _AngleYZ;
double _SizeX;
double _SizeY;

//std::string _type;
//int _NbrDifs;
std::vector<int> _DifInPlate;
};
#endif
