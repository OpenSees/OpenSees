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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-04-23 22:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/Recorder.h,v $
                                                                        
                                                                        
#ifndef Recorder_h
#define Recorder_h

// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for Recorder.
// Recorder is an abstract base class. An Recorder object is used
// in the program to store/restore information at each commit().
//
// What: "@(#) Recorder.h, revA"

class Domain;
#include <MovableObject.h>
#include <TaggedObject.h>


class Recorder: public MovableObject, public TaggedObject
{
  public:
    Recorder(int classTag);
    virtual ~Recorder();

    virtual int record(int commitTag, double timeStamp) =0;
    
    virtual int restart(void);
    virtual int domainChanged(void);
    virtual int setDomain(Domain &theDomain);
    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    virtual void Print(OPS_Stream &s, int flag); 

  protected:
    
  private:	
    static int lastRecorderTag;
};


#endif

