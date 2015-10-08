#ifndef Progress_H
#define Progress_H
#include <iostream> 
#include"Colors.h"
#include <vector>
#include <map>
#include <string>
#include <cstdlib>

std::string Shift(double val)
{
        std::string ret="";
	if(val<10) return ret="  ";
	if(val>=10&&val<100) return ret=" ";
        else return ret;
}
std::string Shift(int val)
{
        std::string ret="";
	if(val<10) return ret="  ";
	if(val>=10&&val<100) return ret=" ";
        else return ret;
}
unsigned int Every(unsigned int & _maxRecord)
{
  if(_maxRecord<=0) return 1000;
  else if(_maxRecord<=10) return 1;
  else if(_maxRecord<=100) return 10;
  else if(_maxRecord<=1000) return 100;
  else  return 1000;
}


void Progress(unsigned int& _skip,unsigned int& _GlobalEvents, unsigned int& _maxRecord, unsigned int& _eventNr)
{
  unsigned int _rolling=Every(_maxRecord);
  unsigned int skip=0;
  if(_skip!=0)skip=_skip+1;
  int maxRecordplusskip=0;
  if(_maxRecord+skip>=_GlobalEvents) 
  {
     maxRecordplusskip=_GlobalEvents;  
  }
  else maxRecordplusskip=_maxRecord+skip;
  if(_maxRecord>=_GlobalEvents)_maxRecord=_GlobalEvents ;
  if(_eventNr %_rolling ==0 || _eventNr==_GlobalEvents || _eventNr==maxRecordplusskip )
  {
    	if(_maxRecord==-1)
	{
		std::cout<<red<<"[";
                int percent=int((_eventNr-skip)*100.0/(_GlobalEvents-skip));
                if(percent<10)std::cout<<"  ";
   		if(percent>=10&&percent!=100)std::cout<<" ";
		std::cout<<percent<<"%]"<<normal<<" Event Number : "<<_eventNr<<"/"<<_GlobalEvents<<std::endl;
	}
        else 
	{
		std::cout<<red<<"[";
		int percent=int((_eventNr-skip)*100.0/(_maxRecord));
		if(percent<10)std::cout<<"  ";
    if(percent>=10&&percent!=100)std::cout<<" ";
		std::cout<<percent<<"%]"<<normal<<" Event Number : "<<_eventNr<<"/"<<maxRecordplusskip<<" Total : "<<_GlobalEvents<<std::endl;
	}
  }





}


#endif
