#ifndef _INTVECTORVAL_H_
#define _INTVECTORVAL_H_

#include "DBValue.h"

using namespace std;


class IntVectorVal : public DBValue{
 private:
  vector<int> val;
 public:
  IntVectorVal(){isNull=true;};
  IntVectorVal(vector<int> v){val=v;isNull=false;};
  string getType(){
    return "vector<int>";
  }
  
  vector<int> getIntVector(){
    if(isNull){
      return vector<int>();
    }
    return val;
  }
  
  void setIntVector(vector<int> v){val=v;isNull=false;};

  string toString(){
    ostringstream oss;
    if(isNull)
      oss<<"";
    else{
      for(vector<int>::const_iterator itr=val.begin();itr!=val.end();itr++){
	oss<<(*itr)<<" ";
      }
    }
    return oss.str();
  }
};
#endif
