/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 00:36:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PeerNGAMotion.h,v $
                                                                        
                                                                        
#ifndef PeerNGAMotion_h
#define PeerNGAMotion_h

// Written: fmk 
// Created: 10/06
//
// Description: This file contains the class definition for PeerNGAMotion.
// PeerNGAMotion is a TimeSeries class which obtains the reference points from
// the Peer Strong Motion Database, interpolates the load factor using 
// these specified control points and the time. The dT for the control
// points is obtained from the Peer Database.

// object. The control points are obtained from the peer strong motion database
//
// What: "@(#) PeerNGAMotion.h, revA"

#include <TimeSeries.h>

class Vector;

class PeerNGAMotion : public TimeSeries
{
  public:
    // constructors  
  PeerNGAMotion(int tag,
		const char *earthquake,
		const char *station,
		const char *type,
		double cfactor = 1.0);

  PeerNGAMotion(int tag,
		const char *earthquakeStation,
		const char *station,
		double cfactor = 1.0);

  PeerNGAMotion();    
    
  // destructor    
  ~PeerNGAMotion();
  
  TimeSeries *getCopy(void);


  // method to get factor
  double getFactor(double pseudoTime);
  double getDuration ();
  double getPeakFactor ();
  double getTimeIncr (double pseudoTime);
  double getDt();
  int getNPts();
  
  // methods for output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag =0);    
  
 protected:
  PeerNGAMotion(int tag,
		Vector *thePath,
		double dT, 
		double cFactor);
  
 private:
  Vector *thePath;      // vector containing the data points
  double dT;
  int currentTimeLoc;   // current location in time
  double cFactor;       // additional factor on the returned load factor
  int dbTag1, dbTag2;   // additional database tags needed for vector objects
  int lastSendCommitTag;
  int otherDbTag;       // a database tag needed for the vector object   
  
  Channel *lastChannel;
};

#endif

