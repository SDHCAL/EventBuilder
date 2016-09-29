#include "Exceptions/Exception.h"

namespace ILCException{
Exception::Exception(string c){
  cout<<c<<endl;
  message = c;
}

string Exception::getMessage(){
  return message;
}
}
