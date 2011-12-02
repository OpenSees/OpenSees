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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DatastoreRecorder.h,v $
                                                                        
                                                                        
// File: ~/DatastoreRecorder/DatastoreRecorder.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for DatastoreRecorder.
// A DatastoreRecorder object is used in the program to invoke commitState()/
// on an FE\_datastore object when commit() is invoked on a Domain. This means 
// we do not have to add additional methods in the Domain for FE_datastore objects.
//
// What: "@(#) DatastoreRecorder.h, revA"

#ifndef DatastoreRecorder_h
#define DatastoreRecorder_h

#include <Recorder.h>

class Domain;
class FE_Datastore;

class DatastoreRecorder: public Recorder
{
  public:
    DatastoreRecorder(FE_Datastore &theDatastore);
    ~DatastoreRecorder();
    int record(int commitTag);
    int playback(int commitTag);
    void restart(void);
    
  protected:
    
  private:	
    FE_Datastore *theDatastore;
};


#endif

