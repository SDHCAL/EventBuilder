#ifndef _DBVALUE_H_
#define _DBVALUE_H_

#include <iostream>
#include <sstream>
#include <vector>
#include "Exceptions/Exception.h"

using namespace std;
using namespace ILCException;

class DBValue{
 protected:
  bool isNull;
  
 public:

  virtual string getType()=0;
  virtual bool null();
  virtual string toString()=0;
  virtual ~DBValue(){};

  virtual string getString();
  virtual int getInt();
  virtual vector<int> getIntVector();

  virtual void setString(string s);
  virtual void setInt(int i);
  virtual void setIntVector(vector<int> v);
  virtual void setNull();

  
};
#endif
