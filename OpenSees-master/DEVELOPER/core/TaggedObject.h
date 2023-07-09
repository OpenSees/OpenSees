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
// $Date: 2003-02-14 23:02:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/TaggedObject.h,v $
                                                                        
                                                                        
#ifndef TaggedObject_h
#define TaggedObject_h

// File: ~/tagged/TaggedObject.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for TaggedObject.
// A TaggedObject is an object with an integer identifier. It is used as
// a base class by DomainComponent, Graph and other classes in the framework.
//
// What: "@(#) TaggedObject.h, revA"

#include <OPS_Stream.h>

class Domain;

class TaggedObject 
{
  public:
    TaggedObject(int tag);
    virtual ~TaggedObject();

    inline int getTag(void) const;

    virtual void Print(OPS_Stream &s, int flag =0) =0;       
    friend OPS_Stream &operator<<(OPS_Stream &s, TaggedObject &m);        

  protected:
    void setTag(int newTag);  // CAUTION: this is a dangerous method to call
    
  private:    
    int theTag;    
};

// INLINED TAGGED_OBJECT FUNCTIONS

inline int 
TaggedObject::getTag(void) const
{
    return theTag;
}

#endif


