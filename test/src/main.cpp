#include "Global.h"
#include "XMLParser.h"
#include "ProcessorLoader.h"

#include "lcio.h"
#include <stdio.h>
#include "Colors.h"
#include "Intro.h"
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCRunHeader.h" 

#include "EVENT/SimCalorimeterHit.h" 
#include "EVENT/CalorimeterHit.h" 
#include "EVENT/RawCalorimeterHit.h" 

#include "UTIL/CellIDDecoder.h"

#include <cstdlib>

using namespace std ;
using namespace lcio ;
using namespace marlin;

//#include "TH1F.h"
//#include "TFile.h"
#include <iostream>
#include <string>
#include <vector>
//Boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;




int main(int argc, char** argv )
{
  Intro();
  std::cout<<red<<"ttetetet"<<normal<<std::endl;
  //------ load shared libraries with processors ------
  StringVec libs ;
  LCTokenizer t( libs, ':' ) ;
  std::string marlinProcs("") ;
  char * var =  getenv("MARLINO_DLL" ) ;
  if( var != nullptr ) marlinProcs = var ;
  else std::cout<<red<<"<!-- You have no MARLINO_DLL variable in your environment - so no processors will be loaded. ! -->"<<normal<< std::endl ;
  std::for_each( marlinProcs.begin(), marlinProcs.end(), t ) ;
  ProcessorLoader loader( libs.begin() , libs.end()  ) ;
  if( loader.failedLoading() )return(1);
  //------- end processor libs -------------------------
  XMLParser* parser =nullptr;
  std::string filen="";
  if(argc==2) filen=argv[1];
  else 
  {
    std::cout<<red<<"Please provide a .xml file "<<normal<<std::endl;
    std::exit(1);
  }
  if( filen.rfind(".xml") == std::string::npos ||!(  filen.rfind(".xml")+ strlen(".xml") == filen.length() ) ) 
  {
    std::cout<<red<<"Please provide a .xml file "<<normal<<std::endl;  
    std::exit(1);
  } 
  else 
  {
    std::cout<<red<<"Here"<<normal<<std::endl;
    parser = new XMLParser(filen.c_str()) ;
    std::cout<<green<<"Here"<<normal<<std::endl;
    //parser->setCmdLineParameters( cmdlineparams ) ;
  }
  parser->parse() ;
  std::cout<<red<<"Here"<<normal<<std::endl;
  std::vector<std::string> DirectoryNecessary{"RawData","XML","DetectorGeometry"};
  std::vector<std::string> DirectoryToCheck{"Streamout","Trivent","Analysis","Results"};
  fs::path data_dir(fs::current_path());
  for(unsigned int i=0;i!=DirectoryNecessary.size();++i)
  {
    std::string pa=data_dir.string()+"/"+DirectoryNecessary[i];
    std::cout<<pa<<std::endl;
    if(!fs::is_directory(pa))
    {
      std::cout<<red<<"Folder "<<DirectoryNecessary[i]<<" is missing. Please add it with the right files in it ! "<<normal<<std::endl;
      std::exit(1); 
    }
  }
  for(unsigned int i=0;i!=DirectoryToCheck.size();++i)
  {
    std::string pa=data_dir.string()+"/"+DirectoryToCheck[i];
    if(!fs::is_directory(pa))
    {
      if (boost::filesystem::create_directory(pa))
      std::cout <<green<< "Creating "<<DirectoryToCheck[i]<<normal << "\n";
    }
    else std::cout<<green <<DirectoryToCheck[i]<<" founded" <<normal<< "\n";
  }
  std::string w = "test"; // word to look for
  for(auto& p: fs::recursive_directory_iterator(data_dir.string()+"/RawData"))
  {
    if(p.path().extension() == ".slcio") 
    {
     std::cout<<p.path()<<std::endl;    
    }
  }

  char* FILEN ;
  int runNumber=0 ;
  int evtNumber=0 ;
  int nthEvent=1 ;

  // read file name from command line (only argument) 
  if( argc < 3 ) {

    cout << " usage: dumpevent filename runNum evtNum " << endl ;
    cout << "    or: dumpevent filename n      " << endl ;
    cout << "  where the first dumps the event with the specified run and event number" << endl ;
    cout << "  and the second simply dumps the n-th event in the file" << endl << endl ;
    cout << "  set the environment variable LCIO_READ_COL_NAMES to a space separated list" << endl ;
    cout << "  of collection names that you would like to dump (all are dumped if not set)" << endl ;

    exit(1) ;
  }
  
  FILEN = argv[1] ;

  bool dumpNthEvent( argc == 3 ) ;
 

  if( dumpNthEvent ) {

    nthEvent  = atoi( argv[2] ) ;

    if( nthEvent < 1 ) {

      cout << " usage: dumpevent filename n   -   whith  n > 0 !  " << endl ;
      
      exit(1) ;
    }

  }else{

    runNumber = atoi( argv[2] ) ;
    evtNumber = atoi( argv[3] ) ;
  }
  

  // set the default encoding for cellid's according to the old Mokka convention
  CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
  CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;
  CellIDDecoder<RawCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

  LCReader* lcReader ;
  if( dumpNthEvent ) 
    lcReader = LCFactory::getInstance()->createLCReader() ;
  else
    lcReader = LCFactory::getInstance()->createLCReader(LCReader::directAccess) ;
  

  // ------ check if LCIO_READ_COL_NAMES is set -------------
  
  char* rColChar = getenv ("LCIO_READ_COL_NAMES");
  
  if ( rColChar != 0 ) {

    std::vector< std::string > colSubset ;
    std::stringstream sts( rColChar ) ;
    std::string colName;

    while( sts >> colName) {
    
      colSubset.push_back( colName ) ;
    }

    lcReader->setReadCollectionNames(  colSubset ) ;
  }
  //-----------------------------------------------------------



  LCEvent* evt(0) ;

  try{
    
     lcReader->open( FILEN ) ;
     
     if( dumpNthEvent ) {
       
       if( nthEvent > 1 )
	 lcReader->skipNEvents(  nthEvent - 1 ) ;

       evt = lcReader->readNextEvent() ; 
       
     }else{
       
       evt = lcReader->readEvent(runNumber,  evtNumber) ; 
     }
  
     
     //   } catch( EndOfDataException& e) {
     //     cout << " couldn't find event " << evtNumber << " - run " << runNumber 
     // 	 << " in file " << FILEN << endl ;    
     //     exit(1) ;
   
     if( !evt  ){

       if(dumpNthEvent){

	 cout << " less than " << nthEvent << "  events in  file " << FILEN << endl ;    
	 
       }else{

	 cout << " couldn't find event " << evtNumber << " - run " << runNumber 
	      << " in file " << FILEN << endl ;    
       } 
       
       exit(1) ;
     }

     LCTOOLS::dumpEventDetailed( evt ) ;
     
     
     lcReader->close() ;
     
   }
   catch( IOException& e) {
     cout << e.what() << endl ;
     exit(1) ;
   }
   return 0 ;
}

