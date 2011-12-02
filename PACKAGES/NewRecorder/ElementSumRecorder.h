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
                                                                        
// $Revision: 1.1 $
// $Date: 2009-03-20 22:44:17 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewRecorder/ElementSumRecorder.h,v $
                                                                        
#ifndef ElementSumRecorder_h
#define ElementSumRecorder_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for ElementSumRecorder.
// A ElementSumRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) ElementSumRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <ID.h>

class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class ElementSumRecorder: public Recorder
{
  public:
    ElementSumRecorder();
    ElementSumRecorder(const ID eleID, 
		       char **argv, 
		       int argc,
		       bool echoTime, 
		       OPS_Stream *theOutputHandler);

    ~ElementSumRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

  protected:
    
  private:	
    int numEle;
    ID *eleID;

    Response **theResponses;

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    bool echoTimeFlag;             // flag indicating if pseudo time also printed

    Vector *data;
    char **responseArgs;
    int numArgs;
};


#endif
