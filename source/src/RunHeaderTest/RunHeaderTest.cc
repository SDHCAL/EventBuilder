#include <iostream>
#include <iomanip>
#include <string>

#include "lcio.h"

#include "IO/LCWriter.h"
#include "EVENT/LCIO.h"
#include "DATA/LCIntVec.h"

#include "IMPL/LCEventImpl.h" 
#include "IMPL/LCRunHeaderImpl.h" 
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCParametersImpl.h"

using namespace lcio ;

void createLCIOfile(std::string filename)
{
  int runnumber=12;
  std::string DetectorName="Albator";

  // create sio writer
  LCWriter* lcWrt = LCFactory::getInstance()->createLCWriter()  ; 
  lcWrt->open( filename , LCIO::WRITE_NEW )  ;
  LCRunHeaderImpl* runHdr = new LCRunHeaderImpl ; 
  runHdr->setRunNumber( runnumber ) ;
  runHdr->setDetectorName( DetectorName ) ;
  runHdr->setDescription( "Test paramaters in run headers" );
  
  runHdr->parameters().setValue("UnEntier",33);
  runHdr->parameters().setValue("UnFloat",float(3.14159));
  runHdr->parameters().setValue("UneString","Le camion de pompier est rouge");

  lcWrt->writeRunHeader( runHdr ) ;

  LCEventImpl*  evt = new LCEventImpl() ;
  evt->setRunNumber( runnumber );
  evt->setEventNumber( 7 );
  evt->setDetectorName( DetectorName ) ;
  evt->parameters().setValue("LePetitChaperon","Rouge");
 
  LCCollectionVec* vecIntCol=new LCCollectionVec( LCIO::LCINTVEC )  ;
  LCIntVec* vecInt=new LCIntVec;
  for (int i=10;i>=0;--i) vecInt->push_back(i);
  vecIntCol->push_back(vecInt);
  evt->addCollection( vecIntCol, "CompteARebours" );

  lcWrt->writeEvent( evt ) ;
  
  delete evt;
  delete runHdr;
  lcWrt->close() ;
  delete lcWrt ;
}


int main(int argc, char** argv)
{
  for (int i=0; i<argc; ++i)
    std::cout << std::setw(4) << i << " : "  << argv[i] << std::endl;

  createLCIOfile("test.slcio");

  return 0;
}
