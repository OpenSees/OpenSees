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
                                                                        
// $Revision: 1.14 $
// $Date: 2009-04-14 21:14:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.h,v $
                                                                        
                                                                        
#ifndef ElementRecorder_h
#define ElementRecorder_h


// Written: fmk 
// Created: 09/99
// Revision: A
//
// Description: This file contains the class definition for ElementRecorder.
// A ElementRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) ElementRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <ID.h>

class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class ElementRecorder: public Recorder
{
  public:
    ElementRecorder();
    ElementRecorder(const ID *eleID, 
		    const char **argv, 
		    int argc,
		    bool echoTime, 
		    Domain &theDomain, 
		    OPS_Stream &theOutputHandler,
		    double deltaT = 0.0,
		    const ID *dof = 0);

    ~ElementRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:

    
  private:	
    int initialize(void);

    int numEle;
    int numDOF;

    ID *eleID;
    ID *dof;

    Response **theResponses;

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    bool echoTimeFlag;             // flag indicating if pseudo time also printed

    double deltaT;
    double nextTimeStampToRecord;

    Vector *data;
    bool initializationDone;
    char **responseArgs;
    int numArgs;

    int addColumnInfo;
};


#endif
