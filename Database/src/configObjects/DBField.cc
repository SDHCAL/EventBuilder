#include "configObjects/DBField.h"

DBField::DBField(string n, string t, int l){
  name = n;
  type = t;
  length = l;
}

string DBField::getName(){
  return name;
}

string DBField::getType(){
  return type;
}

int DBField::getLength(){
  return length;
}
