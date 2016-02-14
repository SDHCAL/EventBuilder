#ifndef PLATE_H
#define PLATE_H
#include<vector>
class Plate
{
 public:
  Plate()
  {
  };
  Plate(const Plate& _Plate):_PositionX(_Plate._PositionX),_PositionY(_Plate._PositionY),_PositionZ(_Plate._PositionZ),_AngleXY(_Plate._AngleXY),_AngleXZ(_Plate._AngleXZ),_AngleYZ(_Plate._AngleYZ),_DifInPlate(_Plate._DifInPlate),_SizeX(_Plate._SizeX),_SizeY(_Plate._SizeY),m_glass_type(_Plate.m_glass_type),m_Gaz_number(_Plate.m_Gaz_number),m_HV_channel(_Plate.m_HV_channel),m_Gaz_channel(_Plate.m_Gaz_channel){};
 Plate(const double& x,const double& y,const double& z, const double& xy,const double& xz,const double& yz,std::vector<int>DifInPlate,const double& sizeX, const double& sizeY,const double& glass_type, const double& Gaz_number,const std::string& HV_channel,const std::string& Gaz_channel):_PositionX(x),_PositionY(y),_PositionZ(z),_AngleXY(xy),_AngleXZ(xz),_AngleYZ(yz),_DifInPlate(DifInPlate),_SizeX(sizeX),_SizeY(sizeY),m_glass_type(glass_type),m_Gaz_number(Gaz_number),m_HV_channel(HV_channel),m_Gaz_channel(Gaz_channel){};
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
  inline double GetGlassType(){return m_glass_type;};
  inline double GetGazNumber(){return m_Gaz_number;};
  inline std::string GetHVChannel(){return m_HV_channel;};
  inline std::string GetGazChannel(){return m_Gaz_channel;};

  
 private:
  double _PositionX;
  double _PositionY;
  double _PositionZ;
  double _AngleXY;
  double _AngleXZ;
  double _AngleYZ;
  std::vector<int> _DifInPlate;
  double _SizeX;
  double _SizeY;
  double m_glass_type;
  double m_Gaz_number;
  std::string m_HV_channel;
  std::string m_Gaz_channel;
};
#endif
