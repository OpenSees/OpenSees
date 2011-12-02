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
// $Date: 2004-11-24 22:41:22 $
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

class Recorder: public MovableObject
{
  public:
    Recorder(int classTag);
    virtual ~Recorder();

    virtual int record(int commitTag, double timeStamp) =0;

    
    virtual int restart(void);
    virtual int setDomain(Domain &theDomain);
    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
  protected:
    
  private:	
};


#endif

