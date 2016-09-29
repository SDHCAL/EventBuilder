#ifndef GEOMETRYPLAN_H
#define GEOMETRYPLAN_H
#include <cmath>
#include <map>
#include <vector>
#include "Geometry.h"
#define degtorad 0.0174532925

class geometryplan
{
  public:
  geometryplan(){};
  ~geometryplan(){};
  geometryplan(int numeroPlan,float x ,float y,float z,float xy, float xz, float yz, int type,double _Ip,double _Im, double _Jp, double _Jm) : Nbr(numeroPlan),X0(x),Y0(y),Z0(z),XY(xy),XZ(xz),YZ(yz),Type(type),Ip(_Ip),Im(_Im),Jp(_Jp),Jm(_Jm)
  {
        ca=cos(xy*degtorad);
        sa=sin(xy*degtorad);
        cb=cos(xz*degtorad);
        sb=sin(xz*degtorad);
        cg=cos(yz*degtorad);
        sg=sin(yz*degtorad);
          /*****************************/
        /* Xexp=pxz0+pxz1*Zexp
           Yexp = pyz0+pyz1*Zexp
           Zexp=Zexp
           (xnorm,ynorm,znorm).(Xexp-X0,Yexp-Y0,Zexp-z0)=0
        */
        xnorm=sg*sa+cg*sb*ca;
        ynorm=-cg*sa+sg*sb*ca;
        znorm=cb*ca;
        xi=cg*cb;
        yi=sg*cb;
        zi=-sb;
        xj=-sg*ca+cg*sb*sa;
        yj=cg*ca+sg*sb*sa;
        zj=cb*sa;
    }
    geometryplan(const int& numeroPlan,const int& type,Geometry& geom,std::map<int,std::vector<double>>& Delimiter)
    {
      Nbr=numeroPlan;
      X0=geom.GetPlatePositionX(numeroPlan);
	    Y0=geom.GetPlatePositionY(numeroPlan);
	    Z0=geom.GetPlatePositionZ(numeroPlan);
	    XY=geom.GetDifPlateAlpha(numeroPlan);
	    XZ=geom.GetDifPlateBeta(numeroPlan);
	    YZ=geom.GetDifPlateGamma(numeroPlan);
	    Type=type;
	    ca=cos(XY*degtorad);
      sa=sin(XY*degtorad);
      cb=cos(XZ*degtorad);
      sb=sin(XZ*degtorad);
      cg=cos(YZ*degtorad);
      sg=sin(YZ*degtorad);
      /*****************************/
      /* Xexp=pxz0+pxz1*Zexp
           Yexp = pyz0+pyz1*Zexp
           Zexp=Zexp
           (xnorm,ynorm,znorm).(Xexp-X0,Yexp-Y0,Zexp-z0)=0
      */
      xnorm=sg*sa+cg*sb*ca;
      ynorm=-cg*sa+sg*sb*ca;
      znorm=cb*ca;
      xi=cg*cb;
      yi=sg*cb;
      zi=-sb;
      xj=-sg*ca+cg*sb*sa;
      yj=cg*ca+sg*sb*sa;
      zj=cb*sa;
	    if(Delimiter.size()==0)
	    {
	      Ip=0;
	      Im=0;
	      Jp=0;
	      Jm=0;
	    }
	    else
	    {
	      Ip=Delimiter[numeroPlan+1][1];
	      Im=Delimiter[numeroPlan+1][0];
	      Jp=Delimiter[numeroPlan+1][3];
	      Jm=Delimiter[numeroPlan+1][2];
	    } 
    }
    inline float GetZexp(const double & pxz0,const double & pyz0,const double & pxz1,const double& pyz1)
    {
        return (xnorm*(pxz0-X0)+ynorm*(pyz0-Y0)+znorm*Z0)/(xnorm*pxz1+ynorm*pyz1+znorm);
    }
    inline float GetProjectioni(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return (Xexp-X0)*xi+(Yexp-Y0)*yi+(Zexp-Z0)*zi;
    }
    inline float GetProjectionj(const double & Xexp,const double & Yexp,const double & Zexp)
    {
        return (Xexp-X0)*xj+(Yexp-Y0)*yj+(Zexp-Z0)*zj;
    }
    inline float get_ca(){return ca;};
    inline float get_sa(){return sa;};
    inline float get_cb(){return cb;};
    inline float get_sb(){return sb;};
    inline float get_cg(){return cg;};
    inline float get_sg(){return sg;};
     inline int NbrPlate()
    {
        return Nbr;
    }
    inline float GetXY()
    {
        return XY;
    }
    inline float GetXZ()
    {
        return XZ;
    }
    inline float GetYZ()
    {
        return YZ;
    }
    inline float GetX0()
    {
        return X0;
    }
    inline float GetY0()
    {
        return Y0;
    }
    inline float GetZ0()
    {
        return Z0;
    }
    inline int GetType()
    {
        return Type;
    }
    inline double GetIp()
    {
        return Ip;
    }
    inline double GetIm()
    {
        return Im;
    }
    inline double GetJp()
    {
        return Jp;
    }
    inline double GetJm()
    {
        return Jm;
    }
    private:
    int  Nbr;
    float X0,Y0,Z0,XY,XZ,YZ,xnorm,ynorm,znorm,xi,yi,zi,xj,yj,zj,ca,sa,cb,sb,cg,sg;
    int Type;
    double Ip;
    double Im;
    double Jp;
    double Jm;

};

#endif
