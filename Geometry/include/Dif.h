#ifndef DIF_H
#define DIF_H
class Dif
{
 public:
  Dif(){};
  Dif(const double& x,const double& y, int DifId, const double& xy, const double& xz,const double& yz, const int& nbr, const int& up_down, const int& DifType):_PositionX(x),_PositionY(y),_DifId(DifId),_AngleXY(xy),_AngleXZ(xz),_AngleYZ(yz),_nbrPlate(nbr),_up_down(up_down),_DifType(DifType){};
  //void SetType(const std::string& Type);
  inline double GetPositionX(){return _PositionX;};
  inline double GetPositionY(){return _PositionY;};
  inline int GetDifId(){return _DifId;};
  inline double GetAngleXY(){return _AngleXY;};
  inline double GetAngleXZ(){return _AngleXZ;};
  inline double GetAngleYZ(){return _AngleYZ;};
  inline int GetNbrPlate(){return _nbrPlate;};
  inline int GetDifType(){return _DifType;};
  inline int GetDifUpDown(){return _up_down;};
  //double GetType(const std::string& Type);
 private:
  double _PositionX;
  double _PositionY;
  int _DifId;
  double _AngleXY;
  double _AngleXZ;
  double _AngleYZ;
  int _nbrPlate;
  int _up_down;
  int _DifType;
};
#endif
