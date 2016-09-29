#include "lcio.h"

#include "Intro.h"

#include "LCIOSTLTypes.h"

#include "ProcessorMgr.h"
#include "Processor.h"
#include "Exceptions.h"
#include "IO/LCReader.h"

#include "XMLParser.h"

#include "Global.h"

#include <sstream>
#include <fstream>
#include <string>
#include <assert.h>


#include <cstring>
#include <algorithm>


#include "ProcessorLoader.h"

#include "VerbosityLevels.h"
#include "streamlog.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/** Helper class for Parser
 */
class LCTokenizer{

  std::vector< std::string >& _tokens ;
  char _del ;
  char _last ;
 public:

  LCTokenizer( std::vector< std::string >& tokens, char del ) : _tokens(tokens) , _del(del), _last(del) {
  }


  void operator()(const char& c) { 

    if( c != _del  ) {

	if( _last == _del  ) {
	  _tokens.push_back("") ; 
	}
      _tokens.back() += c ;
      result() ;
    }
    _last = c ;

  } 

  ~LCTokenizer(){
  }
  
  std::vector<std::string> & result()  { 
    
    return _tokens ; 
    
  }
};
// void createProcessors( XMLParser&  parser ) ;
void  createProcessors( const IParser&  parser) ;

void listAvailableProcessors() ;
void listAvailableProcessorsXML() ;
int printUsage() ;


/** LCIO framework that can be used to analyse LCIO data files
 *  in a modular way. All tasks have to be implemented in Subclasses
 *  of Processor. They will be called in the order specified in the steering file.
 *
 */
int main(int argc, char** argv ){
  Intro();
  // ---- catch all uncaught exceptions in the end ...
  try{
 
    
    if( argc > 1 ){
        if( std::string(argv[1]) == "-x" ){
            std::cout  << "<?xml version=\"1.0\" encoding=\"us-ascii\"?>" << std::endl
                << "<!-- ?xml-stylesheet type=\"text/xsl\" href=\"http://ilcsoft.desy.de/marlin/marlin.xsl\"? -->" << std::endl
                << "<!-- ?xml-stylesheet type=\"text/xsl\" href=\"marlin.xsl\"? -->" << std::endl << std::endl;
        }
    }

    //#ifndef MARLIN_NO_DLL

    //------ load shared libraries with processors ------

    StringVec libs ;
    LCTokenizer t( libs, ':' ) ;

    std::string marlinProcs("") ;

    char * var =  getenv("MARLINO_DLL" ) ;

    if( var != 0 ) {
        marlinProcs = var ;
    } else {
        std::cout << std::endl << "<!-- You have no MARLINO_DLL variable in your environment "
            " - so no processors will be loaded. ! --> " << std::endl << std::endl ;
    }

    std::for_each( marlinProcs.begin(), marlinProcs.end(), t ) ;

    ProcessorLoader loader( libs.begin() , libs.end()  ) ;
    if( loader.failedLoading() ){
        return(1);
    }

    //------- end processor libs -------------------------

    //#endif


    const char* steeringFileName = "none"  ;

    //map<string, map<string,string> > cmdlineparams; 
    CommandLineParametersMap cmdlineparams; 

    // check for dynamic command line arguments
    for( int i = 1 ; i < argc ; i++ ) {
        // cout << "argv[" << i << "]:\t" << argv[i] << endl ;

        if( string( argv[i] ).substr( 0, 2 ) == "--" ){
            // cout << "dynamic opt:\t" << string( argv[i] ).substr( 2 ) << endl ;

            // split dynamic argument by '=', i.e.
            // --global.LCIOInputFiles="1.slcio 2.slcio 3.slcio" --> global.LCIOInputFiles , "1.slcio 2.slcio 3.slcio"
            StringVec cmdlinearg, cmdlinekey ;
            LCTokenizer t( cmdlinearg, '=' ) ;

            string param( argv[i] ) ;

            string s( param.substr( 2 ) ) ;
            for_each( s.begin(), s.end(), t ) ;

            // cout << "split opt:\tkey: " << cmdlinearg[0] << ", value: " << cmdlinearg[1] << endl ;

            if( cmdlinearg.size() != 2 ){
                cerr << endl << "*** invalid command line option: " << argv[i] << endl << endl;
                return printUsage();
            }

            // split first arg by '.'
            // --global.LCIOInputFiles --> global , LCIOInputFiles
            s = cmdlinearg[0] ;
            LCTokenizer t2( cmdlinekey, '.' ) ;

            for_each( s.begin(), s.end(), t2 ) ;

            if( cmdlinekey.size() != 2 ){
                cerr << endl << "*** invalid command line option: " << argv[i] << endl << endl;
                return printUsage();
            }

            // save dynamic options into map
            cmdlineparams[ cmdlinekey[0] ][ cmdlinekey[1] ] = cmdlinearg[1] ;

            cout << "<!-- steering file parameter: [ " << cmdlinekey[0] << "." << cmdlinekey[1] << " ] will be OVERWRITTEN with value: [\"" << cmdlinearg[1] << "\"] -->" << endl;

            // erase dynamic options from **argv
            for( int j = i ; j < argc-1 ; j++ ){
                argv[j] = argv[j+1];
            }
            argc--;
            i--;
        }
    }

    cout << endl ;

    // read file name from command line
    if( argc > 1 ){

        if( std::string(argv[1]) == "-l" ){
            listAvailableProcessors() ;
            return(0) ;
        }
        else if( std::string(argv[1]) == "-x" ){
            listAvailableProcessorsXML() ;
            return(0) ;
        }
        else if( std::string(argv[1]) == "-h"  || std::string(argv[1]) == "-?" ){

            return printUsage() ;
        }


        // one argument given: the steering file for normal running :
        steeringFileName = argv[1] ;

    } else {

        return printUsage() ;
    }


    //###### init streamlog ######
    std::string binname = argv[0]  ;
    binname = binname.substr( binname.find_last_of("/") + 1 , binname.size() ) ;
    streamlog::out.init( std::cout , binname ) ;



    IParser* parser ;

    // for now allow xml and old steering
    std::string filen(  steeringFileName ) ;

    if( filen.rfind(".xml") == std::string::npos ||  // .xml not found at all
            !(  filen.rfind(".xml")
                + strlen(".xml") == filen.length() ) ) {  
        std::exit(1);

    } else {

        parser = new XMLParser( steeringFileName ) ;

        // tell parser to take into account any options defined on the command line
        parser->setCmdLineParameters( cmdlineparams ) ;

    }

    parser->parse() ;

    Global::parameters = parser->getParameters("Global")  ;

    //fg: can't use assert, as this generates no code when compiled with NDEBUG
    if( Global::parameters  == 0 ) {
        std::cout << "  Could not get global parameters from steering file ! " << std::endl  
            << "   The program has to exit - sorry ! " 
            << std::endl ;
        return(1) ;
    }


    // //-----  register log level names with the logstream ---------
    streamlog::out.addLevelName<DEBUG>() ;
    streamlog::out.addLevelName<DEBUG0>() ;
    streamlog::out.addLevelName<DEBUG1>() ;
    streamlog::out.addLevelName<DEBUG2>() ;
    streamlog::out.addLevelName<DEBUG3>() ;
    streamlog::out.addLevelName<DEBUG4>() ;
    streamlog::out.addLevelName<DEBUG5>() ;
    streamlog::out.addLevelName<DEBUG6>() ;
    streamlog::out.addLevelName<DEBUG7>() ;
    streamlog::out.addLevelName<DEBUG8>() ;
    streamlog::out.addLevelName<DEBUG9>() ;
    streamlog::out.addLevelName<MESSAGE>() ;
    streamlog::out.addLevelName<MESSAGE0>() ;
    streamlog::out.addLevelName<MESSAGE1>() ;
    streamlog::out.addLevelName<MESSAGE2>() ;
    streamlog::out.addLevelName<MESSAGE3>() ;
    streamlog::out.addLevelName<MESSAGE4>() ;
    streamlog::out.addLevelName<MESSAGE5>() ;
    streamlog::out.addLevelName<MESSAGE6>() ;
    streamlog::out.addLevelName<MESSAGE7>() ;
    streamlog::out.addLevelName<MESSAGE8>() ;
    streamlog::out.addLevelName<MESSAGE9>() ;
    streamlog::out.addLevelName<WARNING>() ;
    streamlog::out.addLevelName<WARNING0>() ;
    streamlog::out.addLevelName<WARNING1>() ;
    streamlog::out.addLevelName<WARNING2>() ;
    streamlog::out.addLevelName<WARNING3>() ;
    streamlog::out.addLevelName<WARNING4>() ;
    streamlog::out.addLevelName<WARNING5>() ;
    streamlog::out.addLevelName<WARNING6>() ;
    streamlog::out.addLevelName<WARNING7>() ;
    streamlog::out.addLevelName<WARNING8>() ;
    streamlog::out.addLevelName<WARNING9>() ;
    streamlog::out.addLevelName<ERROR>() ;
    streamlog::out.addLevelName<ERROR0>() ;
    streamlog::out.addLevelName<ERROR1>() ;
    streamlog::out.addLevelName<ERROR2>() ;
    streamlog::out.addLevelName<ERROR3>() ;
    streamlog::out.addLevelName<ERROR4>() ;
    streamlog::out.addLevelName<ERROR5>() ;
    streamlog::out.addLevelName<ERROR6>() ;
    streamlog::out.addLevelName<ERROR7>() ;
    streamlog::out.addLevelName<ERROR8>() ;
    streamlog::out.addLevelName<ERROR9>() ;
    streamlog::out.addLevelName<SILENT>() ;


    //-------- init logging level ------------
    std::string verbosity = Global::parameters->getStringVal("Verbosity" ) ;
    streamlog::logscope scope( streamlog::out ) ;

    scope.setLevel( verbosity ) ;

    createProcessors( *parser ) ;




    StringVec lcioInputFiles ; 



        int maxRecord = Global::parameters->getIntVal("MaxRecordNumber") ;
        int skipNEvents = Global::parameters->getIntVal("SkipNEvents");

        bool modify = ( Global::parameters->getStringVal("AllowToModifyEvent") == "true" ) ;

	if( modify ) {

	  streamlog_out( WARNING )  << " ******************************************************************************* \n" 
				    << " *    AllowToModifyEvent is set to 'true'                                      * \n"
				    << " *    => all processors can modify the input event in processEvent() !!        * \n"
				    << " *        consider setting this flag to 'false'                                * \n"
				    << " *        unless you really need it...                                         * \n"
				    << " *    - if you need a processor that modifies the input event                  * \n"
				    << " *      please implement the EventModifier interface and use the modifyEvent() * \n"
				    << " *      method for this                                                        * \n"
 				    << " ******************************************************************************* \n" 
				    << std::endl ;
	}

        // create lcio reader 
        LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;

	StringVec readColNames ; 
	if( (Global::parameters->getStringVals("LCIOReadCollectionNames" , readColNames ) ).size() != 0 ){
	  
	  streamlog_out( WARNING )  << " *********** Parameter LCIOReadCollectionNames given - will only read the following collections: **** " 
				    << std::endl ;

	  for( unsigned i=0,N=readColNames.size() ; i<N ; ++i ) {
	    streamlog_out( WARNING )  << "     " << readColNames[i] << std::endl ;
	  } 
	  streamlog_out( WARNING )  << " *************************************************************************************************** " << std::endl ;

#if  LCIO_PATCHVERSION_GE( 2,4,0 )

	  lcReader->setReadCollectionNames( readColNames ) ;
#endif
	} 

        lcReader->registerLCRunListener( ProcessorMgr::instance() ) ; 
        lcReader->registerLCEventListener( ProcessorMgr::instance() ) ; 

        ProcessorMgr::instance()->init() ; 

        bool rewind = true ;

        while( rewind ) {

            rewind = false ;

            // process the data
            lcReader->open( lcioInputFiles  ) ; 


            if( skipNEvents > 0 ){

                streamlog_out( WARNING ) << " --- Marlin.cc - will skip first " << skipNEvents << " event(s)" 
                    << std::endl << std::endl ;

                lcReader->skipNEvents(  skipNEvents ) ;
            }

            
                if( maxRecord > 0 ){

                    
                        lcReader->readStream( maxRecord ) ;
                    
               

                } else {

                    lcReader->readStream() ;
                }


            


            lcReader->close() ;

            if( !rewind ) {

                ProcessorMgr::instance()->end() ; 

                delete lcReader ;
            }

        } // end rewind



    return 0 ;

  } catch( std::exception& e) {

    std::cerr << " ***********************************************\n" 
	      << " A runtime error occured - (uncaught exception):\n" 
	      << "      " <<    e.what() << "\n"
	      << " Marlin will have to be terminated, sorry.\n"  
	      << " ***********************************************\n" 
	      << std:: endl ; 

    return 1 ;

  }

}

//   void  createProcessors(XMLParser&  parser) {
void createProcessors( const IParser&  parser) {

    StringVec activeProcessors ;
    Global::parameters->getStringVals("ActiveProcessors" , activeProcessors ) ;

    StringVec procConds ;
    Global::parameters->getStringVals("ProcessorConditions" , procConds ) ;

    bool useConditions = ( activeProcessors.size() == procConds.size() ) ;

    //     for( StringVec::iterator m = activeProcessors.begin() ; m != activeProcessors.end() ; m++){
    for(unsigned int i=0 ; i<  activeProcessors.size() ; i++ ) {

        StringParameters* p = parser.getParameters( activeProcessors[i] )  ;

        if( p!=0 ){
            std::string type = p->getStringVal("ProcessorType") ;

            if( useConditions ) 
                ProcessorMgr::instance()->addActiveProcessor( type , activeProcessors[i] , p , procConds[i] )  ;
            else
                ProcessorMgr::instance()->addActiveProcessor( type , activeProcessors[i] , p )  ;

        } else{
            std::stringstream sstr ;
            sstr << "Undefined processor : " << activeProcessors[i] << std::endl ;
            streamlog_out( ERROR )  << sstr.str() ;	
            throw Exception( sstr.str() );
        }
    }
}



void listAvailableProcessors() {

    ProcessorMgr::instance()->dumpRegisteredProcessors() ;
}

void listAvailableProcessorsXML() {

    ProcessorMgr::instance()->dumpRegisteredProcessorsXML() ;
}


int printUsage() {

  std::cout << " Usage: Marlin [OPTION] [FILE]..." << std::endl 
	    << "   runs a Marlin application " << std::endl 
	    << std::endl 
	    << " Running the application with a given steering file:" << std::endl 
	    << "   Marlin steer.xml   " << std::endl 
	    << std::endl 
	    << "   Marlin [-h/-?]             \t print this help information" << std::endl 
	    << "   Marlin -x                  \t print an example steering file to stdout" << std::endl 
	    << "   Marlin -l                  \t [deprecated: old format steering file example]" << std::endl 
	    << std::endl 
	    << " Example: " << std::endl 
	    << " To create a new default steering file from any Marlin application, run" << std::endl 
	    << "     Marlin -x > mysteer.xml" << std::endl 
	    << " and then use either an editor or the MarlinGUI to modify the created steering file " << std::endl 
	    << " to configure your application and then run it. e.g. : " << std::endl 
	    << "     Marlin mysteer.xml > marlin.out 2>&1 &" << std::endl << std::endl
	    << " Dynamic command line options may be specified in order to overwrite individual steering file parameters, e.g.:" << std::endl 
	    << "     Marlin --global.LCIOInputFiles=\"input1.slcio input2.slcio\" --global.GearXMLFile=mydetector.xml" << std::endl 
	    << "            --MyLCIOOutputProcessor.LCIOWriteMode=WRITE_APPEND --MyLCIOOutputProcessor.LCIOOutputFile=out.slcio steer.xml" << std::endl << std::endl
	    << "     NOTE: Dynamic options do NOT work together with Marlin options (-x, -f, -l) nor with the MarlinGUI or old steering files" << std::endl 
	    << std::endl ;
  
  return(0) ;

}

