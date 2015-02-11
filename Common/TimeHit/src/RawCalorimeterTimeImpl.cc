#include "RawCalorimeterTimeImpl.h"

namespace IMPL{
  
  RawCalorimeterTimeImpl::RawCalorimeterTimeImpl() :
    _cellID0(0),
    _cellID1(0),
    _amplitude(0),
    _timeStamp(0){
  }
  
  
  RawCalorimeterTimeImpl::~RawCalorimeterTimeImpl(){
  }
  
  int RawCalorimeterTimeImpl::getCellID0() const {
    return _cellID0 ;
  }
  
  int RawCalorimeterTimeImpl::getCellID1() const {
    return _cellID1 ;
  }
  
  double RawCalorimeterTimeImpl::getTime() const {
    return _amplitude ;
  }

  int RawCalorimeterTimeImpl::getTimeStamp() const {
    return _timeStamp ;
  }
  
  
  
  void RawCalorimeterTimeImpl::setCellID0(int id0){
    checkAccess("RawCalorimeterTimeImpl::setCellID0") ;
    _cellID0 = id0 ;
  }
  
  void RawCalorimeterTimeImpl::setCellID1(int id1){
    checkAccess("RawCalorimeterTimeImpl::setCellID1") ;
    _cellID1 = id1 ;
  }
  
  void RawCalorimeterTimeImpl::setTime(double amplitude){
    checkAccess("RawCalorimeterTimeImpl::setTime") ;
    _amplitude = amplitude ;
  }

  void RawCalorimeterTimeImpl::setTimeStamp(int timeStamp){
    checkAccess("RawCalorimeterTimeImpl::setTimeStamp") ;
    _timeStamp = timeStamp ;
  }
  

  
 
} // namespace IMPL
