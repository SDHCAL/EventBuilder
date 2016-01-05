#ifndef _JOB_H_
#define _JOB_H_

#include <iostream>
#include <sstream>
#include "../db/DBConnection.h"

using namespace std;


/**
   Description of a job from the DAQ
**/
class Job{

 private :

  int id;
  int status;
  int jobID;
  int daqID;
  string host;

 public :

  Job(string h, int jid);
  /**
     Change the status of the Job
     @param s The new status : 0 CREATED, 1 DESTROY
  **/    
  void setStatus(int s);
  /**
     Used for database purpose, do not use this method!
  **/
 void setInternalId(int i);
  /**
     You should not have to use this method (will be called from the DAQ class)
  **/
  void setDaqId(int d);
  
  int getID();
  int getStatus();
  string getHost();

  /**
     Send the current object to the database. You should not have to use this method (it will be called while uploading a Daq object)
  **/
  //void upload() throw (ILCException::Exception);

  /**
     Get the definition of a job status code
  **/
  static string getStatusDefinition(int val) throw (ILCException::Exception);

};
#endif
