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
// $Date: 2000-09-15 08:23:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/actor/MovableObject.h,v $
                                                                        
                                                                        
#ifndef MovableObject_h
#define MovableObject_h

// File: ~/actor/MovableObject.h
//
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for MovableObject.
// MovableObject is meant to be an abstract base class and thus no objects 
// of it's type can be instantiated. A movable object is an object which
// can send/receive itself to/from a Channel object.
//
// What: "@(#) MovableObject.h, revA"

#include <classTags.h>

class Channel;
class FEM_ObjectBroker;

class MovableObject
{
  public:
    MovableObject(int classTag, int dbTag);        
    MovableObject(int classTag);    
    virtual ~MovableObject();

    int getClassTag(void) const;
    int getDbTag(void) const;
    void setDbTag(int dbTag);

    virtual int sendSelf(int commitTag, Channel &theChannel) =0;  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker) =0;
    
  protected:
    
  private:
    int classTag;
    int dbTag;
};

#endif
