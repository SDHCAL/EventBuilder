#ifndef _INTVAL_H_
#define _INTVAL_H_

#include "DBValue.h"

using namespace std;


class IntVal : public DBValue{
 private:
  int val;
 public:
  IntVal(){val=0;isNull=true;}
  IntVal(int i){val=i;isNull=false;}
  string getType(){
      return "int";
  }
  int getInt(){
    if(isNull)
      return 0;
    return val;
  }
  void setInt(int i){val=i;isNull=false;}
  string toString(){
    ostringstream oss;
    if(isNull)
      oss<<"";
    else
      oss<<val;
    return oss.str();
  }
};
#endif
