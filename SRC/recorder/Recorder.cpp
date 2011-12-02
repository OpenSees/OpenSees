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
// $Date: 2005-01-10 21:58:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/Recorder.cpp,v $
                                                                        

// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class implementation for Recorder.
//
// What: "@(#) Recorder.cpp, revA"

#include <Recorder.h>
#include <OPS_Globals.h>

Recorder::Recorder(int classTag)
  :MovableObject(classTag)
{

}

Recorder::~Recorder() 
{

}

int 
Recorder::restart(void)
{
  return 0;
}

int 
Recorder::setDomain(Domain &theDomain)
{
  return 0;
}

int 
Recorder::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "Recorder::sendSelf() - not yet implemented\n";
  return 0;
}

int 
Recorder::recvSelf(int commitTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker)
{
  opserr << "Recorder::recvSelf() - not yet implemented\n";
  return 0;
}
