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
// $Date: 2004-11-24 22:45:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DamageRecorder.h,v $
                                                                        
#ifndef DamageRecorder_h
#define DamageRecorder_h


// Written: Arash Altoontash, Gregory Deierlein, 04/04
// Created: 04/04
// Revision: Arash Altoontash
//
// Description: This file contains the class definition for DamageRecorder.
// A DamageRecorder is used to obtain a response from an element section/material during 
// the analysis and apply the information to the damage model and record the damage index.
//
// What: "@(#) DamageRecorder.h, revA"

#include <Recorder.h>
#include <fstream>
using std::ofstream;

#include <Information.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;
class FE_Datastore;
class DamageModel;

class DamageRecorder: public Recorder
{
  public:
    DamageRecorder( int elemid, ID &secIDs, int dofid, DamageModel *dmgPtr, Domain &theDomainPtr,
		    bool echotimeflag, double deltat , const char *filename );

    ~DamageRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);

    int restart(void);    
    
  protected:
    
  private:	

    int eleID, numSec, dofID;
    ID responseID;                 // integer element returns in setResponse
	ID sectionTags;

    Response **theResponses;
	DamageModel **theDamageModels;

    Domain *theDomain;
    bool echoTimeFlag;             // flag indicating if pseudo time also printed
    char *fileName;                // file name  
    ofstream theFile; 	           // output stream

    double deltaT;
    double nextTimeStampToRecord;

};


#endif
