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
                                                                        
// $Revision: 1.4 $
// $Date: 2009-08-25 22:11:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/TaggedObjectStorage.h,v $
                                                                        
                                                                        
#ifndef TaggedObjectStorage_h
#define TaggedObjectStorage_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for 
// TaggedObjectStorage. A TaggedObjectStorage object a container object
// used to hold objects of type TaggedObject; each object of which has
// some UNIQUE identifier. 
//
// What: "@(#) TaggedObjectStorage.h, revA"

class TaggedObject;
class TaggedObjectIter;

#include <OPS_Globals.h>

class TaggedObjectStorage 
{
  public:
    TaggedObjectStorage() {};
    virtual  ~TaggedObjectStorage() {};    

    // public methods to populate the container
    virtual  int   setSize(int newSize) =0;
    virtual  bool  addComponent(TaggedObject *newComponent) =0;
//				bool allowMultipleTags = false) =0;
    virtual  TaggedObject  *removeComponent(int tag) =0;    
    virtual  int   getNumComponents(void) const =0;
    
    virtual  TaggedObject *getComponentPtr(int tag) =0;
    virtual  TaggedObjectIter  &getComponents(void) =0;

    virtual  TaggedObjectStorage *getEmptyCopy(void) =0;
    virtual  void clearAll(bool invokeDestructors = true) =0;
    
    virtual  void Print(OPS_Stream &s, int flag =0) =0;
    
  protected:    
    
  private:
};

#endif


