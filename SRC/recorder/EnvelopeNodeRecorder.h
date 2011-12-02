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
                                                                        

// $Revision: 1.8 $
// $Date: 2006-08-04 22:33:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeNodeRecorder.h,v $
                                                                        

                                                                        
#ifndef EnvelopeNodeRecorder_h
#define EnvelopeNodeRecorder_h

// Written: fmk 
//
// Description: This file contains the class definition for EnvelopeRecorder.
// A EnvelopeRecorder is used to record the envelop of specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) EnvelopeNodeRecorder.h, revA"


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Domain;
class FE_Datastore;
class Node;

class EnvelopeNodeRecorder: public Recorder
{
  public:
    EnvelopeNodeRecorder();
    EnvelopeNodeRecorder(const ID &theDof, 
			 const ID &theNodes, 
			 const char *dataToStore,
			 Domain &theDomain,
			 OPS_Stream &theOutputHandler,
			 double deltaT = 0.0,
			 bool echoTimeFlag = false); 
    
    ~EnvelopeNodeRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:	
    int initialize(void);

    ID *theDofs;
    ID *theNodalTags;
    Node **theNodes;

    Vector *currentData;
    Matrix *data;

    Domain *theDomain;
    OPS_Stream *theHandler;

    int dataFlag; // flag indicating what it is to be stored in recorder

    double deltaT;
    double nextTimeStampToRecord;

    bool first;
    bool initializationDone;
    int numValidNodes;

    bool echoTimeFlag;   // flag indicating whether time to be included in o/p
};

#endif
