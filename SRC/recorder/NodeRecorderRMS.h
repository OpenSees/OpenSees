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
                                                                        
#ifndef NodeRecorderRMS_h
#define NodeRecorderRMS_h

// Written: fmk 


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <TimeSeries.h>

class Domain;
class FE_Datastore;
class Node;

class NodeRecorderRMS: public Recorder
{
  public:
    NodeRecorderRMS();
    NodeRecorderRMS(const ID &theDof, 
			 const ID *theNodes, 
			 const char *dataToStore,
			 Domain &theDomain,
			 OPS_Stream &theOutputHandler,
			 double deltaT = 0.0,
			 TimeSeries **theTimeSeries =0); 
    
    ~NodeRecorderRMS();

    int record(int commitTag, double timeStamp);
    int restart(void);    

    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
    virtual double getRecordedValue(int clmnId, int rowOffset, bool reset); //added by SAJalali

  protected:
    
  private:	
    int initialize(void);

    ID *theDofs;
    ID *theNodalTags;
    Node **theNodes;

    Vector *currentData;
    Vector *runningTotal;
    int count;

    Domain *theDomain;
    OPS_Stream *theHandler;

    int dataFlag; // flag indicating what it is to be stored in recorder

    double deltaT;
    double nextTimeStampToRecord;

    bool initializationDone;
    int numValidNodes;

    int addColumnInfo;
    TimeSeries **theTimeSeries;
    double *timeSeriesValues;

};

#endif
