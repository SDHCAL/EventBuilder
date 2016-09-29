#ifndef _RUNINFO_H_
#define _RUNINFO_H_

#include <iostream>
#include <sstream>
#include "db/DBConnection.h"
#include "Daq.h"

using namespace std;


/**
   Informations about a Run
**/
class RunInfo{

 private :

  int number;
  int status;
  int type;
  string startDate;
  string stopDate;
  int daqID;
  string comment;

  RunInfo();

 public :

  /**
     Constructor
     @param t The type of the run
     @param d The name of the Daq used for the run
  **/
  RunInfo(int t, string d) throw (ILCException::Exception);
  /**
     Change the status of the Run
     @param s The new status : 0 INITIALIZED, 1 CONFIGURED, 2 RUNNING, 3 STOPPED, 4 END
  **/  
  //void setStatus(int s) throw (ILCException::Exception);
 /**
     Change the status of the run to RUNNING (2) and record the startDate
  **/
  //void start() throw (ILCException::Exception);
 /**
     Change the status of the run to STOPPED (3) and record the stopDate
  **/
  //void stop() throw (ILCException::Exception);
  //void setDescription(string d) throw (ILCException::Exception);

  int getRunNumber();
  int getStatus();
  int getType();
  string getStartTime();
  string getStopTime();
  /**
     Get a pointer on the DAQ used for the run : you have to delete this pointer.
     @return A pointer on the used DAQ
  **/
  Daq* getDaq() throw (ILCException::Exception);
  string getDescription();

  /**
     Download the RunInfo object from the database
     @param n The run number
     @return A pointer on the object containing all informations. You have to delete this pointer.
  **/
  static RunInfo* getRunInfo(int n) throw (ILCException::Exception);

 /**
     Get the definition of a run status code
  **/
  static string getStatusDefinition(int val) throw (ILCException::Exception);

 /**
     Get the definition of a run type code
  **/
  static string getTypeDefinition(int val) throw (ILCException::Exception);

  /**
     Get the list of RunInfo using a given DAQ
     @param dname The name of the DAQ
     @return A vector of RunInfo pointers. You have to delete these pointers yourself!
  **/
  static vector<RunInfo*> getRunsFromDaqName(string dname) throw (ILCException::Exception);

};
#endif
