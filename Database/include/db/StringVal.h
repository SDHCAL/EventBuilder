#ifndef _STRINGVAL_H_
#define _STRINGVAL_H_

#include "DBValue.h"

using namespace std;


class StringVal : public DBValue{
 private:
  string val;
 public:
  StringVal(){val="";isNull=true;};
  StringVal(string s){val=s;isNull=false;};
  string getType(){
    return "string";
  }
  string getString(){return val;}
  void setString(string s){val=s;isNull=false;}
  string toString(){
    if(isNull)
      return "";
    return val;
  }
};
#endif
