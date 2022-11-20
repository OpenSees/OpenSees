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
                                                                        
// $Revision: 1.10 $
// $Date: 2006-08-04 22:33:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DriftRecorder.h,v $
                                                                        
#ifndef DriftRecorder_h
#define DriftRecorder_h

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for 
// DriftRecorder. 

#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

class Domain;
class Node;

class DriftRecorder: public Recorder
{
 public:
  DriftRecorder();
  DriftRecorder(int ndI, int ndJ, int dof, int perpDirn,
		Domain &theDomain, 
		OPS_Stream &theHandler,
		bool echoTime = false,
		double dT = 0.0,
		double relDeltaTTol = 0.00001);

  DriftRecorder(const ID &ndI, const ID &ndJ, int dof, int perpDirn,
		Domain &theDomain, 
		OPS_Stream &theHandler,
		bool echoTime = false,
		double dT = 0.0,
		double relDeltaTTol = 0.00001);
  
  ~DriftRecorder();

  int record(int commitTag, double timeStamp);
  int restart(void);    
  int flush(void);    

  int setDomain(Domain &theDomain);
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
  virtual double getRecordedValue(int clmnId, int rowOffset, bool reset); //added by SAJalali

 protected:
  
 private:	
  int initialize(void);

  ID *ndI;
  ID *ndJ;
  Node **theNodes; // i & j nodes
  int dof;
  int perpDirn;
  Vector *oneOverL;
  Vector *data;

  Domain *theDomain;
  OPS_Stream *theOutputHandler;

  bool initializationDone;
  int numNodes;
  bool echoTimeFlag;   // flag indicating whether time to be included in o/p
  
  double deltaT, nextTimeStampToRecord;
  double relDeltaTTol;
};

#endif
