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
                                                                        

// $Revision: 1.4 $
// $Date: 2004-01-29 23:30:30 $
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

#include <fstream>
using std::ofstream;

class Domain;
class FE_Datastore;

class EnvelopeNodeRecorder: public Recorder
{
  public:
    EnvelopeNodeRecorder(const ID &theDof, 
		 const ID &theNodes, 
		 Domain &theDomain,
		 const char *fileName,
		 const char *dataToStore,
		 double deltaT = 0.0,
		 int startFlag = 0); 

    EnvelopeNodeRecorder(const ID &theDof, 
		 const ID &theNodes, 
		 Domain &theDomain,
		 FE_Datastore *database,
		 const char *dbTable,
		 const char *dataToStore,
		 double deltaT = 0.0,
		 int startFlag = 0); 
    
    ~EnvelopeNodeRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    void restart(void);    
    
  protected:
    
  private:	
    ID *theDofs;
    ID *theNodes;
    Vector *currentData;
    Matrix *data;
    Domain *theDomain;

    int flag;   // flag indicating whether time, load factor or nothing printed
	        // at start of each line in file
    char *fileName;
    ofstream theFile;     

    int dataFlag; // flag indicating what it is to be stored in recorder

    double deltaT;
    double nextTimeStampToRecord;

    FE_Datastore *db;
    char **dbColumns;
    int numDbColumns;

    bool first;
};

#endif
