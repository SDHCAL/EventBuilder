#ifndef Progress_H
#define Progress_H
#include <iostream> 
#include "Colors.h"
#include <cstdlib>
#include "Patch.h"
#include <string>

inline std::string Shift(double val)
{
        std::string ret="";
	if(val<10) return ret="  "+patch::to_string(val);
	if(val>=10&&val<1000) return ret=" "+patch::to_string(val);
	if(val>=1000&&val<10000) return ret=patch::to_string(val);
        else return ret+patch::to_string(val);
}

inline unsigned int Every(unsigned int & _maxRecord)
{
  if(_maxRecord<=0) return 1000;
  else if(_maxRecord<=10) return 1;
  else if(_maxRecord<=100) return 10;
  else if(_maxRecord<=1000) return 100;
  else  return 1000;
}


inline void Progress(unsigned int& _skip,unsigned int& _GlobalEvents, unsigned int& _maxRecord, unsigned int& _eventNr,std::string Perso="")
{
  unsigned int _rolling=Every(_maxRecord);
  unsigned int skip=0;
  if(_skip!=0)skip=_skip+1;
  unsigned int maxRecordplusskip=0;
  if(_maxRecord+skip>=_GlobalEvents) 
  {
     maxRecordplusskip=_GlobalEvents;  
  }
  else maxRecordplusskip=_maxRecord+skip;
  if(_maxRecord>=_GlobalEvents)_maxRecord=_GlobalEvents ;
  if(_eventNr %_rolling ==0 || _eventNr==_GlobalEvents || _eventNr==maxRecordplusskip )
  {
    if(_maxRecord==0)  
	{
    int percent=int((_eventNr-skip)*100.0/(_GlobalEvents-skip));
    if(percent<100) 
		{
		  if(Perso!="")std::cout<<Perso<<std::endl;
		  std::cout<<red<<"["<<Shift(percent)<<"%]"<<normal<<" Event Number : "<<Shift(_eventNr)<<"/"<<_GlobalEvents<<std::endl;
		}
	}
        else 
	{
	  int percent=int((_eventNr-skip)*100.0/(_maxRecord));
	  if(percent<100) 
		{
		if(Perso!="")std::cout<<Perso<<std::endl;
		std::cout<<red<<"["<<Shift(percent)<<"%]"<<normal<<" Event Number : "<<Shift(_eventNr)<<"/"<<maxRecordplusskip<<" Total : "<<_GlobalEvents<<std::endl;
		}
	}
  }





}


#endif
