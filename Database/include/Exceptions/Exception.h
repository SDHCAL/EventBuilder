#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <iostream>
#include <string>

using namespace std;

namespace ILCException{
/**
Class describing a default exception
**/
class Exception{
 private:
  string message;
 public:
  /**
     Constructor
  **/
  Exception(string);
  /**
     Returns the message contained in the exception
     @return A string containing the message
  **/
  string getMessage();
};
}
#endif
