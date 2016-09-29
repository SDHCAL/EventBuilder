#ifndef _DAQ_H_
#define _DAQ_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include "db/DBConnection.h"
#include "Job.h"

using namespace std;

/**
   Informations about the DAQ system
**/
class Daq{

 private :

  int id;
  int status;
  string configName;
  string host;
  string xml;
  string comment;
  string name;
  string creationDate;
  vector<Job*> jobList;

 public :

  Daq(string n);
  ~Daq();
  
 // void setConfigName(string s);
  /**
     Change the status of the DAQ
     @param s The new status : 0 CREATED, 1 DESTROY
  **/
  //void setStatus(int s);
 /**
     Set the name of the host running the DAQ
     @param h The name of the host
  **/
  //void setHost(string h);
  /**
     Set the XML content of the configuration
     @param fileName The name of the XML file containing the configuration
  **/
  void setXML(string fileName) throw (ILCException::Exception);
  //void setComment(string c);
  
  string getName();
  string getTime();
  string getConfigName();
  int getStatus();
  string getHost();
  /**
     Get the XML configuration
     @return The XML listing of the configuration
  **/
  string getXML();
  string getComment();

  /**
     Download an existing DAQ from the database. You have to delete the object after usage!
  **/
  static Daq* getDaq(string n) throw (ILCException::Exception);

  static string getNameFromId(int i) throw (ILCException::Exception);

  /**
     Send the current object to the database. If the object is already in the database, it is updated with the new values.
  **/
  //void uploadToDatabase() throw (ILCException::Exception);

  /**
     Add a job to the DAQ jobs list. The job will be deleted along with the Daq object (you should not do it yourself)
     @param j The job to add
  **/
  void addJob(Job* j);

  /**
     Get the list of jobs associated to this DAQ
     @return A vector containing pointers on job objects (This is not a copy, do not delete them!).
  **/
  vector<Job*> getJobs();

  /**
     Get the definition of a daq status code
  **/
  static string getStatusDefinition(int val) throw (ILCException::Exception);

  /**
     Get the names of existing DAQ
  **/
  static vector<string> getListOfDaqNames() throw (ILCException::Exception);

};
#endif
