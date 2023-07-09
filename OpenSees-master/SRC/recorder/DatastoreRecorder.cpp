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
                                                                        
// $Revision: 1.5 $
// $Date: 2004-11-25 00:53:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DatastoreRecorder.cpp,v $
                                                                        
                                                                        
// File: ~/DatastoreRecorder/DatastoreRecorder.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for DatastoreRecorder.
// A DatastoreRecorder object is used in the program to store/restore the
// Domain information in a FE_Datastore object.
//
// What: "@(#) DatastoreRecorder.h, revA"

#include <DatastoreRecorder.h>
#include <Domain.h>
#include <FE_Datastore.h>

DatastoreRecorder::DatastoreRecorder(FE_Datastore &theDb)
:Recorder(RECORDER_TAGS_DatastoreRecorder), theDatastore(&theDb)
{
    
}


DatastoreRecorder::~DatastoreRecorder()
{

}


int 
DatastoreRecorder::record(int commitTag, double timeStamp)
{
    return theDatastore->commitState(commitTag);
}


int 
DatastoreRecorder::playback(int commitTag)
{
    return theDatastore->restoreState(commitTag);
}

int
DatastoreRecorder::restart(void)
{
	return 0;
}










