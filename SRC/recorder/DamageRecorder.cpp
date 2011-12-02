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
// $Date: 2006-08-04 22:33:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DamageRecorder.cpp,v $
                                                                        
// Written: Arash Altoontash, Gregory Deierlein, 
// Created: 04/04
// Revision: A
//
// Description: This file contains the class implementatation of ElementRecorder.
//
// What: "@(#) ElementRecorder.C, revA"

#include <DamageRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <DamageModel.h>

#include <OPS_Globals.h>

#include <iomanip>
using std::ios;

DamageRecorder::DamageRecorder( int elemid, ID &secIDs, int dofid, DamageModel *dmgPtr,
				Domain &theDomainPtr, bool echotimeflag, double deltat , OPS_Stream &output)
  :Recorder(RECORDER_TAGS_DamageRecorder),
   eleID(elemid) , numSec(secIDs.Size()), dofID(dofid),
   responseID(secIDs.Size()), sectionTags(secIDs.Size()),theDomain(&theDomainPtr),
   echoTimeFlag(echotimeflag), deltaT(deltat), nextTimeStampToRecord(0.0),
   theOutput(&output), data(0)
{
  // make copy of the damage model
  if ( dmgPtr == NULL ) {
    opserr << "DamageRecorder::DamageRecorder - no damage pointer associated with the damge recorder" <<endln;
    exit(-1);
  }

  theOutput->tag("OpenSeesOutput");

  int numDbColumns = 0;
  if (echoTimeFlag == true) {
    theOutput->tag("TimeOutput");
    theOutput->attr("ResponseType", "time");
    theOutput->endTag();
    numDbColumns += 1;
  }

  theDamageModels = new DamageModel *[numSec];

  int j;
  for (j= 0; j<numSec; j++) {
    theDamageModels[j] = dmgPtr->getCopy();
    if ( theDamageModels [j] == NULL ) {
      opserr << "DamageRecorder::DamageRecorder - out of memory copying damage models ";
      exit(-1);
    }
    theDamageModels[j]->revertToStart();
  }
  
  // Get the element	
  Information eleInfo(1.0);
  Element *theEle = theDomainPtr.getElement(eleID);
  if ( theEle == NULL ) {
    opserr << "WARNING DamageRecorder::DamageRecorder() - no element with tag: "
	   << eleID << " exists in Domain\n";
    exit(-1);
  }
  
  // allocate pointers to element responses and damage models
  theResponses = new Response *[3*numSec];
  for ( j=0; j<3*numSec; j++)
    theResponses[j] = NULL;
  
  // establish response for the element sections or materials
  const int argc = 3;
  char *argv[argc];
  for ( j = 0; j<argc ; j++ )
    argv[j] = new char[20];
  
  strcpy( argv[0] , "-section" );
  strcpy( argv[2] , "deformation" );
  for ( j=0 ; j<numSec ; j++) {
    sectionTags(j) = secIDs(j);
    sprintf(argv[1],"%d",sectionTags(j));
    // itoa( sectionTags(j) , argv[1] , 10);
    
    theResponses[j] = theEle->setResponse( ( const char**)argv, argc, eleInfo, *theOutput);
    if (theResponses[j] == 0) {
      opserr << "DamageRecorder::DamageRecorder - out of memory creating deformation response ";
      exit(-1);
    }
  }
  
  
  strcpy( argv[2] , "force" );
  for ( j=0 ; j<numSec ; j++) {
    sectionTags(j) = secIDs(j);
    sprintf(argv[1],"%d",sectionTags(j));
    //_itoa( sectionTags(j) , argv[1] , 10);
    theResponses[j+numSec] = theEle->setResponse( ( const char**) argv, argc, eleInfo, *theOutput);
    if (theResponses[j+numSec] == 0) {
      opserr << "DamageRecorder::DamageRecorder - out of memory creating force response ";
      exit(-1);
    }
  }
  
  strcpy( argv[2] , "stiffness" );
  for ( j=0 ; j<numSec ; j++) {
    sectionTags(j) = secIDs(j);
    sprintf(argv[1],"%d",sectionTags(j));
    // _itoa( sectionTags(j) , argv[1] , 10);
    theResponses[j+2*numSec] = theEle->setResponse( ( const char**) argv, argc, eleInfo, *theOutput);
    if (theResponses[j+2*numSec] == 0) {
      opserr << "DamageRecorder::DamageRecorder - out of memory creating tanegnt response ";
      exit(-1);
    }
  }
  
  for ( j = 0; j<argc ; j++ )
    delete( argv[j] );


  numDbColumns += numSec;

  // create the vector to hold the data
  data = new Vector(numDbColumns);

  theOutput->tag("Data");
}


DamageRecorder::~DamageRecorder()
{
  if (data != 0)
    delete data;
  
  // 
  // invoke destructor on response args
  //
  
  if (theResponses != 0) {
    for (int i = 0; i < 3*numSec; i++)
      delete theResponses[i];
    delete [] theResponses;
  }
  
  if (theDamageModels!=0)	{
    for ( int i=0; i<numSec; i++)
      delete theDamageModels[i];
    delete []theDamageModels;
  }

  theOutput->endTag(); // Data
  theOutput->endTag(); // OpenSeesOutput

  if (theOutput != 0)
    delete theOutput;
}


int 
DamageRecorder::record(int commitTag, double timeStamp)
{
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    
    // print out the pseudo time if requested
    int counter = 0;
    if (echoTimeFlag == true) {
      (*data)(counter++) = timeStamp;
    }
    
    Vector DamageInformation(3);
    // get the responses and write to file if file or opserr specified
    // for each element do a getResponse() & print the result
    for (int i=0; i< numSec; i++) {
      DamageInformation.Zero();
      for ( int j=0 ; j<2 ; j++) {
	if ( theResponses[i+numSec*j] == 0) {
	  DamageInformation(j) = 0.0;
	} else {
	  if ( theResponses[i+numSec*j]->getResponse() < 0) {
	    DamageInformation(j) = 0.0;
	  } else {
	    // ask the element for the reponse
	    Information &eleinfo = theResponses[i+numSec*j]->getInformation();
	    const Vector &infovector = eleinfo.getData();
	    DamageInformation(j) = infovector(dofID);
	  }
	}
      }
      DamageInformation(2) = 0.0;
      theDamageModels[i]->setTrial(DamageInformation);
      theDamageModels[i]->commitState();
      double Damageindex = theDamageModels[i]->getDamage();
      
      // print results to file or stderr depending on whether
      // a file was opened

      (*data)(counter++) = Damageindex;      
    }
  }

  theOutput->write(*data);

  // succesfull completion - return 0
  return result;
}


int 
DamageRecorder::playback(int commitTag)
{
    return 0;
}

int
DamageRecorder::restart(void)
{
  return 0;
}

