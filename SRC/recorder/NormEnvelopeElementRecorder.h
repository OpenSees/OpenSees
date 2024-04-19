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
// $Date: 2009-04-14 21:14:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NormEnvelopeElementRecorder.h,v $
                                                                        
#ifndef NormEnvelopeElementRecorder_h
#define NormEnvelopeElementRecorder_h

// Written: fmk 
//
// What: "@(#) NormEnvelopeElementRecorder.h, revA"

#include <Recorder.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class NormEnvelopeElementRecorder: public Recorder
{
  public:
    NormEnvelopeElementRecorder();
    NormEnvelopeElementRecorder(const ID *eleID, 
				const char **argv, 
				int argc,
				Domain &theDomain, 
				OPS_Stream &theOutputHandler,
				double deltaT = 0.0,
				double relDeltaTTol = 0.00001,
				bool echoTimeFlag = true,
				const ID *dof =0); 


    ~NormEnvelopeElementRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);    
    int flush(void);    

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
    OPS_Stream *theHandler;

    double deltaT;
    double relDeltaTTol;
    double nextTimeStampToRecord;

    Matrix *data;
    Vector *currentData;
    bool first;

    bool initializationDone;
    char **responseArgs;
    int numArgs;

    bool echoTimeFlag; 

    int addColumnInfo;
};


#endif
