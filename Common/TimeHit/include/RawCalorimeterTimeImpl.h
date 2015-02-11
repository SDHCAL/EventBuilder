#ifndef EVENT_RAWCALORIMETERTIMEIMPL_H
#define EVENT_RAWCALORIMETERTIMEIMPL_H 1


#include "EVENT/MCParticle.h"
#include "RawCalorimeterTime.h"
#include <vector>
#include "IMPL/AccessChecked.h"

namespace IMPL {
  
  
/** Implementation of the real data RawCalorimeterHit.
 * 
 * @author gaede
 * @version $Id: RawCalorimeterHitImpl.h,v 1.3 2006-08-04 16:52:46 gaede Exp $
 */

  class RawCalorimeterTimeImpl : public EVENT::RawCalorimeterTime , public IMPL::AccessChecked {
    
  public: 
    
    /** Default constructor, initializes values to 0.
     */
    RawCalorimeterTimeImpl() ;

    /// Destructor.
    virtual ~RawCalorimeterTimeImpl() ;
    
    virtual int id() const { return simpleUID() ; }

    /** Returns the detector specific (geometrical) cell id.
     */
    virtual int getCellID0() const ;

    /** Returns the second detector specific (geometrical) cell id. Optional, check/set 
     *  flag(LCIO::RCHBIT_ID1)==1.
     */
    virtual int getCellID1() const ;

    /** Returns the amplitude of the hit in ADC counts.
     */
    virtual double getTime() const ;

    /** Returns a time stamp for the hit. Optional, check/set 
     *  flag(LCIO::RCHBIT_TIME)==1.
     */
    virtual int getTimeStamp() const ;

    
    void setCellID0(int id0) ;

    /** Sets the second cell id;
     *  Only store if the flag word (bit CHBIT_ID1) 
     *  of the collection is set. Default is to store only cellid0.
     */
    void setCellID1(int id1) ;
    
    /** Sets the amplitude.
     */
    void setTime(double amplitude) ;
    
    /** Set the time stamp
     */
    void setTimeStamp(int time) ;

  protected:

    int _cellID0 ;
    int _cellID1 ;
    int _amplitude ;
    int _timeStamp ;
    
  }; // class
} // namespace IMPL
#endif /* ifndef EVENT_RAWCALORIMETERHITIMPL_H */
