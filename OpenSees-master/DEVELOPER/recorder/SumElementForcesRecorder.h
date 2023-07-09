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
// $Date: 2009/03/20 22:44:17 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewRecorder/SumElementForcesRecorder.h,v $
                                                                        
#ifndef SumElementForcesRecorder_h
#define SumElementForcesRecorder_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for SumElementForcesRecorder.
// A SumElementForcesRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) SumElementForcesRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <ID.h>

class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class SumElementForcesRecorder: public Recorder
{
  public:
    SumElementForcesRecorder();
    SumElementForcesRecorder(const ID eleID, 
			     bool echoTime, 
			     OPS_Stream *theOutputHandler);

    ~SumElementForcesRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    
    int domainChanged(void);
    int setDomain(Domain &theDomain);

    const char *getClassType(void) const;
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:
    
  private:	
    int numEle;           // the number of elements
    Element **theElements;// pointer to array of element pointers
    ID eleID;            // ID (integer list) of element tags to record

    Domain *theDomain;    // pointer to domain holding elements
    OPS_Stream *theOutput;// pointer to output location
    bool echoTimeFlag;    // flag indicating if pseudo time to be printed
    Vector *data;         // Vector (double array) to store sum of element forces
};


#endif
