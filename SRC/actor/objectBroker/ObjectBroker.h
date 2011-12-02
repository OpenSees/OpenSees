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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/ObjectBroker.h,v $
                                                                        
                                                                        
// File: ~/actor/broker/ObjectBroker.h
//
// Written: fmk
// Created: Fri Sept 20 12:27:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for ObjectBroker.
// ObjectBroker is meant to be an abstract base class and thus no objects of it's type
// can be instantiated. 
//
// What: "@(#) ObjectBroker.h, revA"

#ifndef ObjectBroker_h
#define ObjectBroker_h

class Channel;
class MovableObject;
class GraphNumberer;

class ObjectBroker
{
  public:
    ObjectBroker();
    virtual ~ObjectBroker();

    
  protected:
    
  private:

};

#endif
