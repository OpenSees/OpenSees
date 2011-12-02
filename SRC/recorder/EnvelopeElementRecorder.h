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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-25 23:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeElementRecorder.h,v $
                                                                        
#ifndef EnvelopeElementRecorder_h
#define EnvelopeElementRecorder_h

// Written: fmk 
//
// What: "@(#) EnvelopeElementRecorder.h, revA"

#include <Recorder.h>

#include <fstream>
using std::ofstream;

#include <Information.h>
#include <OPS_Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;

class EnvelopeElementRecorder: public Recorder
{
  public:
    EnvelopeElementRecorder(const ID &eleID, 
			    Domain &theDomain, 
			    const char **argv, 
			    int argc,
			    double deltaT = 0.0, 
			    const char *fileName =0);

    EnvelopeElementRecorder(const ID &eleID, 
			    Domain &theDomain, 
			    const char **argv, 
			    int argc,
			    FE_Datastore *db, 
			    const char *tableName, 
			    double deltaT = 0.0);

    ~EnvelopeElementRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);

    void restart(void);    
    
  protected:
    
  private:	
    int numEle;
    ID responseID;                 // integer element returns in setResponse
    Response **theResponses;

    Domain *theDomain;
    char *fileName;                // file name  
    ofstream theFile; 	           // output stream

    double deltaT;
    double nextTimeStampToRecord;

    FE_Datastore *db;
    char **dbColumns;
    int numDbColumns;

    Matrix *data;
    Vector *currentData;
    bool first;
};


#endif
