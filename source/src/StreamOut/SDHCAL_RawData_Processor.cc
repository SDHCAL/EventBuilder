#include "Streamout/SDHCAL_RawData_Processor.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <assert.h>
#include "Colors.h"
#include <EVENT/LCCollection.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <UTIL/LCTOOLS.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "TH1F.h"
#include "TFile.h"
#include <bitset>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "Reader/ReaderFactory.h"
#include "Reader/Reader.h"
#include "Streamout/DIFUnpacker.h"
#include <sstream>
namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}
int _NbrRun=0;
std::map<int,TH1F*>HistoTimeAsic1;
std::map<int,TH1F*>HistoTimeAsic2;
using namespace lcio ;
using namespace marlin ;
float CalibT0_0 = 7.14;
float CalibDeltaT_0 =0.053;
float CalibT0_1 = 7.142;
float CalibDeltaT_1 =0.0996;
void SDHCAL_buffer::printBuffer(unsigned int start, unsigned int stop,std::ostream& flux)
{
  flux << std::hex;
  for (unsigned int k=start; k<stop; k++) flux << (unsigned int)(first[k]) << " - ";
  flux << std::dec <<  std::endl;
}

SDHCAL_RawBuffer_Navigator::SDHCAL_RawBuffer_Navigator(SDHCAL_buffer b) :_buffer(b),_SCbuffer(0,0)
{
  _DIFstartIndex=DIFUnpacker::getStartOfDIF(_buffer.buffer(),_buffer.getsize(),24); //92 was here
  _theDIFPtr=NULL;
  _badSCdata=false;
}

DIFPtr* SDHCAL_RawBuffer_Navigator::getDIFPtr()
{
  if (NULL==_theDIFPtr) _theDIFPtr=new DIFPtr(getDIFBufferStart(),getDIFBufferSize());
  return _theDIFPtr;
}

uint32_t SDHCAL_RawBuffer_Navigator::getDIF_CRC() 
{
  uint32_t i=getEndOfDIFData();
  uint32_t ret=0;
  ret |= ( (_buffer.buffer()[i-2])<<8 );
  ret |= _buffer.buffer()[i-1];
  return ret;
}

void SDHCAL_RawBuffer_Navigator::setSCBuffer()
{
  if (! hasSlowControlData() ) return;
  if (_SCbuffer.getsize()!=0 ) return; //deja fait
  if (_badSCdata) return;
  _SCbuffer.first=&(getDIFBufferStart()[getEndOfDIFData()]);
  //compute Slow Control size
  uint32_t maxsize=_buffer.getsize()-_DIFstartIndex-getEndOfDIFData()+1; // should I +1 here ?
  uint32_t k=1; //SC Header
  uint32_t dif_ID=_SCbuffer.first[1];
  uint32_t chipSize=_SCbuffer.first[3];
  while ((dif_ID != 0xa1 && _SCbuffer.first[k] != 0xa1 && k <maxsize) ||(dif_ID == 0xa1 && _SCbuffer.first[k+2]==chipSize && k<maxsize))
  {
    k+=2; //DIF ID + ASIC Header
    uint32_t scsize=_SCbuffer.first[k];
    if (scsize != 74 && scsize != 109) 
	  {
	    //std::cout << "PROBLEM WITH SC SIZE " << scsize << std::endl;
	    k=0;
	    _badSCdata=true;
	    break;
	  }
    k++; //skip size bit
    k+=scsize; // skip the data
  }
  if (_SCbuffer.first[k] == 0xa1 && !_badSCdata ) _SCbuffer.second=k+1; //add the trailer
  else 
  {
    _badSCdata=true;
    //std::cout << "PROBLEM SC TRAILER NOT FOUND " << std::endl;
  }
}

SDHCAL_buffer SDHCAL_RawBuffer_Navigator::getEndOfAllData()
{
  setSCBuffer();
  if (hasSlowControlData() && !_badSCdata)
  {
    return SDHCAL_buffer( &(_SCbuffer.buffer()[_SCbuffer.getsize()]), getSizeAfterDIFPtr()-3-_SCbuffer.getsize() );
  }
  else return SDHCAL_buffer( &(getDIFBufferStart()[getEndOfDIFData()]), getSizeAfterDIFPtr()-3 ); //remove the 2 bytes for CRC and the DIF trailer
}

SDHCAL_RawData_Processor aSDHCAL_RawData_Processor;

SDHCAL_RawData_Processor::SDHCAL_RawData_Processor() : Processor("SDHCAL_RawData_Processor") 
{
  // modify processor description
  _description = "SDHCAL_RawData_Processor converts SDHCAL Raw Data into collection of RawCalorimeterHit and RawCalorimeterTime and TcherenkovSignal" ;
  // register steering parameters: name, description, class-variable, default value
  bool debugDefault=false;
  _debugMode=false;
  registerProcessorParameter("DebugMode","Turn ON/OFF debug mode : Warning Debug mode uses assert and may crash the application",_debugMode,debugDefault);

  registerInputCollection( LCIO::LCGENERICOBJECT,"XDAQCollectionName","XDAQ produced collection name",_XDAQCollectionNames,std::string("RU_XDAQ"));
  registerOutputCollection( LCIO::RAWCALORIMETERHIT,"OutputRawCaloHitCollectionName","Name of output collection containing raw calorimeter hits",_RawHitCollectionName,std::string("DHCALRawHits"));
  registerOutputCollection( LCIO::RAWCALORIMETERHIT,"OutputCaloHitTimeCollectionName","Name of output collection containing Times",_RawHitCollectionNameTime,std::string("DHCALRawTimes"));
  _TcherenkovSignal="";
  registerOutputCollection( LCIO::LCGENERICOBJECT,"OutputTcherenkovCollectionName","Name of output collection containing Tcherenkov signal",_TcherenkovSignal,std::string("Tcherenkov")); 
  _nevt=_nWrongObj=_nProcessedObject=_hasSlowControl=_hasBadSlowControl=0;
  _ReaderType="";
  registerProcessorParameter("ReaderType","Type of the Reader needed to read InFileName",_ReaderType ,_ReaderType);
  _FileNameGeometry="";
  registerProcessorParameter("FileNameGeometry","Name of the Geometry File",_FileNameGeometry,_FileNameGeometry);
}

void SDHCAL_RawData_Processor::init() 
{ 
  printParameters() ;
  ReaderFactory readerFactory;
  Reader* myReader = readerFactory.CreateReader(_ReaderType);
  if(myReader)
  {
    myReader->Read(_FileNameGeometry,geom);
    geom.PrintGeom();
    std::map<int, Dif > Difs=geom.GetDifs();
    for(std::map<int, Dif >::iterator it=Difs.begin();it!=Difs.end();++it)
    {
      if(geom.GetDifType(it->first)==temporal)
      {
       std::string name1="Times_given_by_TDC_Asics_1"+patch::to_string(it->first);
       std::string name2="Times_given_by_TDC_Asics_2"+patch::to_string(it->first);
       HistoTimeAsic1.insert ( std::pair<int,TH1F*>(it->first,new TH1F(name1.c_str(),name1.c_str(),501,-250,250)) );
       HistoTimeAsic2.insert ( std::pair<int,TH1F*>(it->first,new TH1F(name2.c_str(),name2.c_str(),501,-250,250)) );
      }
    }
  }
  else
  {
    std::cout << "Reader type n'existe pas !!" << std::endl;
    exit(1);
  }
  delete myReader;
  //Prepare a flag to tag data type in RawVec (dit les types de data qu'on va enregistrer?)
  
}

void SDHCAL_RawData_Processor::processRunHeader( LCRunHeader* run) 
{ 
  LCTOOLS::dumpRunHeader(run);
  
} 

void SDHCAL_RawData_Processor::processEvent( LCEvent * evt ) 
{ 

  IMPL::LCFlagImpl chFlag(0) ;
  EVENT::LCIO bitinfo;
  chFlag.setBit(bitinfo.RCHBIT_LONG );// raw calorimeter data -> format long //(sert a qq chose?)
  chFlag.setBit(bitinfo.RCHBIT_BARREL );// barrel
  chFlag.setBit(bitinfo.RCHBIT_ID1 );// cell ID 
  chFlag.setBit(bitinfo.RCHBIT_TIME );//timestamp
  chFlag.setBit(bitinfo.RCHBIT_ENERGY_ERROR);  
  lcio::IntVec trig(8);
  IMPL::LCCollectionVec *RawVec=new IMPL::LCCollectionVec(LCIO::RAWCALORIMETERHIT) ;
  IMPL::LCCollectionVec *RawVecTime=new IMPL::LCCollectionVec(LCIO::CALORIMETERHIT) ;
  IMPL::LCCollectionVec *RawVecTcherenkov=new IMPL::LCCollectionVec(LCIO::LCGENERICOBJECT) ;
  RawVec->setFlag(chFlag.getFlag());   
  RawVecTime->setFlag(chFlag.getFlag());
  _eventNr=evt->getEventNumber();
  if(_eventNr %1000 ==0)std::cout<<"Event Number : "<<_eventNr<<std::endl;
  _NbrRun=evt->getRunNumber();
  _nevt++;
  
  
  
  
  try
  {
    LCCollection* col = evt->getCollection(_XDAQCollectionNames);
    
    int nElement=col->getNumberOfElements();
    _CollectionSizeCounter[nElement]++;
    for (int iel=0; iel<nElement; iel++)
    {
    	LCGenericObject* obj=dynamic_cast<LCGenericObject*>(col->getElementAt(iel));
        LMGeneric *lmobj=(LMGeneric *) obj;
	if (obj==NULL)
	{
	_nWrongObj++;
	continue;
	}
	_nProcessedObject++;
	//LMGeneric *lmobj=(LMGeneric *) obj;
	//unsigned char* debug_variable_1=lmobj->getSDHCALBuffer().endOfBuffer();
	SDHCAL_RawBuffer_Navigator bufferNavigator(lmobj->getSDHCALBuffer());
	//unsigned char* debug_variable_2=bufferNavigator.getDIFBuffer().endOfBuffer();
	//streamlog_out(DEBUG)<<"DIF BUFFER END "<<(unsigned int *) debug_variable_1<<" "<< (unsigned int *) debug_variable_2 << std::endl;
	//if (_debugMode) assert (debug_variable_1 == debug_variable_2);
	uint32_t idstart=bufferNavigator.getStartOfDIF();
	//streamlog_out(DEBUG)<<"DIF header starts at "<<idstart<<" and is equal to "<<std::hex<<(unsigned int)lmobj->getCharBuffer()[idstart] << std::dec << std::endl;
	if(idstart==0 && streamlog::out.write< streamlog::DEBUG >() ) lmobj->getSDHCALBuffer().printBuffer( 0 , streamlog::out() );
	_DIFStarter[idstart]++;	
	if (!  bufferNavigator.validBuffer() ) continue;
        DIFPtr *d=bufferNavigator.getDIFPtr();
	//if (_debugMode) assert(d!=NULL);
	if (d!=NULL)
	{
		_DIFPtrValueAtReturnedPos[bufferNavigator.getDIFBufferStart()[d->getGetFramePtrReturn()] ]++;
	      	//if (_debugMode) assert( bufferNavigator.getDIFBufferStart()[d->getGetFramePtrReturn()]==0xa0 );
	}
	_SizeAfterDIFPtr[bufferNavigator.getSizeAfterDIFPtr()]++;

	//Make something with the Tcherenkov signal
	if(geom.GetDifType(d->getID())==tcherenkov)
	{
		std::cout<<"I'm tchrenkov's Signal !!!!! Dif : "<<d->getID()<<std::endl;
		IMPL::LCGenericObjectImpl* Tcherenkov= new IMPL::LCGenericObjectImpl(1,0,0) ;
		Tcherenkov->setIntVal(0,(unsigned long int)(d->getFrameTimeToTrigger(0)));
		RawVecTcherenkov->addElement(Tcherenkov);
	}

	//if Difs are temporal ones -> createCalorimeterHit
	
        else if(geom.GetDifType(d->getID())==temporal)
	{
	  
  	for (uint32_t i=0;i<d->getNumberOfFrames();i++)
	{
	    unsigned long int ID0;
	  unsigned long int ID1;
	  unsigned long int TTT;
	  float Time0_6=0;
    	double t=d->getFrameTimeToTrigger(i)*2E-7;
		  int32_t iasic=(d->getFrameData(i,uint32_t(0)));
		  //std::cout<<red<<iasic<<normal<<std::endl;
		  if (t>3.8) 
			{
	  		printf("Wrong Time %f %x \n",t,d->getFrameTimeToTrigger(i));
	 		 	continue;
			}
   		int32_t i1=(d->getFrameData(i,uint32_t(18))) |(d->getFrameData(i,uint32_t(17))<<8) ;
		  int32_t i0=(d->getFrameData(i,uint32_t(16))) |(d->getFrameData(i,uint32_t(15))<<8) ;
		  //int32_t c0=(d->getFrameData(i,uint32_t(19))) &0x7;
			//float Time0_6=0;
		  //for (int iw=0;iw<50;iw++) printf("%.2x ",d->getFrameData(i,iw));
      if (iasic==1)
			{
				Time0_6 =(i0 -i1)*CalibT0_0+i1*CalibDeltaT_0;
			  //printf("=> %d %d %d-> %f %f ns\n",c0,i0,i1,Time0_6,c0*25+Time0_6);
			 	//std::cout<<"iasic : "<<iasic<<"Time :"<<Time0_6<<std::endl;
			  HistoTimeAsic1[d->getID()]->Fill(Time0_6);			  
			}
		  else
			{
				Time0_6=(i0 -i1)*CalibT0_1+i1*CalibDeltaT_1;
				//printf("=> %d %d %d-> %f %f ns\n",c0,i0,i1,Time0_6,c0*25+Time0_6);	
				//std::cout<<"iasic : "<<iasic<<"Time :"<<Time0_6<<std::endl;	
				HistoTimeAsic2[d->getID()]->Fill(Time0_6);	  
			}
			
			
		  ID0=(unsigned long int)(((unsigned short)d->getID())&0xFF);			//8 firsts bits: DIF Id
		  ID0+=(unsigned long int)(((unsigned short)d->getFrameAsicHeader(i)<<8)&0xFF00);	//8 next bits:   Asic Id
		  bitset<6> Channel(0
			//j
			);														
		  ID0+=(unsigned long int)((Channel.to_ulong()<<16)&0x3F0000);			//6 next bits:   Asic's Channel
		  unsigned long BarrelEndcapModule=0;  //(40 barrel + 24 endcap) modules to be coded here  0 for testbeam (over 6 bits)
		  ID0+=(unsigned long int)((BarrelEndcapModule<<22)&0xFC00000);	
		  ID1 = (unsigned long int)(d->getFrameBCID(i));
		  bitset<2> ThStatus;
		  ThStatus[0]=d->getFrameLevel(i,0
			//j
			,0);
		  ThStatus[1]=d->getFrameLevel(i,0
			//j
			,1);
		  //ThStatus[2]=isSynchronised; //I'm not computing the synchronisation here.
	    
		      	TTT = (unsigned long int)(d->getFrameTimeToTrigger(i));
		      	
			      //Use setTime to stock TimeStamp !!!!!!!	      	
			     IMPL::CalorimeterHitImpl *hit2=new IMPL::CalorimeterHitImpl() ;
		      	hit2->setCellID0((unsigned long int)ID0);               
		      	hit2->setCellID1(ID1);
		      	//Use setEnergyError to stock Time from the TDC !!!!!!!
            hit2->setEnergyError((float)Time0_6);
            hit2->setTime((float)TTT);//Time stamp of this event from Run Begining
		      	RawVecTime->addElement(hit2);
        	}
        	  
      }
      	else
	{
      		//create RawCalorimeterHit
	    	for (uint32_t i=0;i<d->getNumberOfFrames();i++)
	    	{
        		for (uint32_t j=0;j<64;j++)
	      		{
		      		if (!(d->getFrameLevel(i,j,0) || d->getFrameLevel(i,j,1))) continue; // skip empty pads
		      		//  std::cout <<" New hit "<<std::endl;
		      		unsigned long int ID0;
		      		ID0=(unsigned long int)(((unsigned short)d->getID())&0xFF);			//8 firsts bits: DIF Id
		      		ID0+=(unsigned long int)(((unsigned short)d->getFrameAsicHeader(i)<<8)&0xFF00);	//8 next bits:   Asic Id
		      		bitset<6> Channel(j);														
		      		ID0+=(unsigned long int)((Channel.to_ulong()<<16)&0x3F0000);			//6 next bits:   Asic's Channel
		      		unsigned long BarrelEndcapModule=0;  //(40 barrel + 24 endcap) modules to be coded here  0 for testbeam (over 6 bits)
		      		ID0+=(unsigned long int)((BarrelEndcapModule<<22)&0xFC00000);	
		      		unsigned long int ID1 = (unsigned long int)(d->getFrameBCID(i));
		      		bitset<2> ThStatus;
		      		ThStatus[0]=d->getFrameLevel(i,j,0);
		      		ThStatus[1]=d->getFrameLevel(i,j,1);
		      		//ThStatus[2]=isSynchronised; //I'm not computing the synchronisation here.
	        		IMPL::RawCalorimeterHitImpl *hit = new IMPL::RawCalorimeterHitImpl() ;
		      		hit->setCellID0((unsigned long int)ID0);               
		      		hit->setCellID1(ID1);
		      		hit->setAmplitude(ThStatus.to_ulong());
		      		unsigned long int TTT = (unsigned long int)(d->getFrameTimeToTrigger(i));
		      		hit->setTimeStamp(TTT);			      		//Time stamp of this event from Run Begining
		      		RawVec->addElement(hit);
	      		}//for (uint32_t j=0;j<64;j++)
      		}//for (uint32_t i=0;i<d->getNumberOfFrames();i++)
      		//register Triggers'time : lots of values here ?
	}
	
	trig[0] = d->getDTC();
	trig[1] = d->getGTC();
	trig[2] = d->getBCID();
	trig[3] = d->getAbsoluteBCID()&0xFFFFFF;
	trig[4] = (d->getAbsoluteBCID()/(0xFFFFFF+1))&0xFFFFFF;
	trig[5] = d->getTASU1(); //what happen if no temperature ?
	trig[6] = d->getTASU2();
	trig[7] = d->getTDIF();
	std::stringstream ss("");
	ss<<"DIF"<<d->getID()<<"_Triggers";
	RawVec->parameters().setValues(ss.str(),trig);
	RawVecTime->parameters().setValues(ss.str(),trig);
	if (bufferNavigator.hasSlowControlData()) _hasSlowControl++;
	if (bufferNavigator.badSCData())  _hasBadSlowControl++;
	SDHCAL_buffer eod=bufferNavigator.getEndOfAllData();
	_SizeAfterAllData[eod.getsize()]++;
	//unsigned char* debug_variable_3=eod.endOfBuffer();
	//streamlog_out(DEBUG)<<"END DATA BUFFER END "<<(unsigned int *)debug_variable_1<<" "<<(unsigned int *)debug_variable_3<<std::endl;
	//if (_debugMode) assert (debug_variable_1 == debug_variable_3);
	//streamlog_out(DEBUG) << "End of Data remaining stuff : ";
	if(streamlog::out.write< streamlog::DEBUG >() ) eod.printBuffer( 0 , streamlog::out() );
	int nonzeroCount=0;
	for (unsigned char* it=eod.buffer(); it != eod.endOfBuffer(); it++) if (int(*it) !=0) nonzeroCount++;
	_NonZeroValusAtEndOfData[nonzeroCount]++;
	} //for (int iel=0; iel<nElement; iel++)
}
catch(DataNotAvailableException &e)
{
	//streamlog_out(WARNING) << _XDAQCollectionNames << " collection not available" << std::endl;
}
  	evt->addCollection(RawVec,_RawHitCollectionName);
  	evt->addCollection(RawVecTime,_RawHitCollectionNameTime);
  	evt->addCollection(RawVecTcherenkov,_TcherenkovSignal);
}

void SDHCAL_RawData_Processor::printCounter(std::string description, std::map<int,int> &m)
{
  streamlog_out(MESSAGE) << " statistics for " << description << " : ";
  for (std::map<int,int>::iterator it=m.begin(); it != m.end(); it++)
  {
    if (it != m.begin() ) streamlog_out(MESSAGE) << ",";
    streamlog_out(MESSAGE) << " [" << it->first << "]=" << it->second;
  }
  streamlog_out(MESSAGE) << std::endl;
}

void SDHCAL_RawData_Processor::end()
{ 
  if(HistoTimeAsic1.size()!=0)
  {
  std::string name="TimeGivenByTDC_"+ patch::to_string(_NbrRun)+".root";
  TFile *hfile = new TFile(name.c_str(),"RECREATE","Results");	
  for(std::map<int, TH1F* >::iterator it=HistoTimeAsic1.begin();it!=HistoTimeAsic1.end();++it) it->second->Write();
  for(std::map<int, TH1F* >::iterator it=HistoTimeAsic2.begin();it!=HistoTimeAsic2.end();++it) it->second->Write();
  hfile->Close();
  delete hfile;
  }
  streamlog_out(MESSAGE) << "FINAL STATISTICS : " << std::endl;
  streamlog_out(MESSAGE) << " runs for " << _nevt << " events." << std::endl;
  printCounter("Size of GenericObject collections",_CollectionSizeCounter);
  streamlog_out(MESSAGE) << "Number of badly extracted LCGenericObject : " <<  _nWrongObj << std::endl;
  streamlog_out(MESSAGE) << "Number of processed LCGenericObject : " <<  _nProcessedObject << std::endl;
  printCounter("Start of DIF header",_DIFStarter);
  printCounter("Value after DIF data are processed",_DIFPtrValueAtReturnedPos);
  printCounter("Size remaining in buffer after end of DIF data",  _SizeAfterDIFPtr );
  streamlog_out(MESSAGE) << "Number of Slow Control found " << _hasSlowControl << " out of which " << _hasBadSlowControl << " are bad" << std::endl;
  printCounter("Size remaining after all of data have been processed",_SizeAfterAllData );
  printCounter("Number on non zero values in end of data buffer",_NonZeroValusAtEndOfData);
}
