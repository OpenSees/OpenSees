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
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/ArrayOfTaggedObjects.h,v $
                                                                        
                                                                        
#ifndef ArrayOfTaggedObjects_h
#define ArrayOfTaggedObjects_h

// File: ~/tagged/storage/ArrayOfTaggedObjects.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for 
// ArrayOfTaggedObjects. ArrayOfTaggedObjects is a storage class. The class 
// is responsible for holding and providing access to objects of type 
// TaggedObject. The data structure used to hold the objects is a simple 
// array of pointers. As a one dimensional array is used certain ideas are tried 
// to improve performance: (1) if the array needs to be larger to hold more 
// components, the array size is doubled and (2) when adding/retrieving components,
// the array location given by the components tag is first checked. 
//
// What: "@(#) ArrayOfTaggedObjects.h, revA"


#include <TaggedObjectStorage.h>
#include <ArrayOfTaggedObjectsIter.h>

class ArrayOfTaggedObjects : public TaggedObjectStorage
{
  public:
    ArrayOfTaggedObjects(int size);
    ~ArrayOfTaggedObjects();    

    // public methods to populate a domain
    int  setSize(int newSize);
    bool addComponent(TaggedObject *newComponent);
//		       bool allowMutltipleTags = false);
    TaggedObject *removeComponent(int tag);    
    int  getNumComponents(void) const;
    
    TaggedObject     *getComponentPtr(int tag);
    TaggedObjectIter &getComponents();

    ArrayOfTaggedObjectsIter  getIter();
    
    virtual TaggedObjectStorage *getEmptyCopy(void);
    virtual void clearAll(bool invokeDestructor = true);
    
    void Print(OPS_Stream &s, int flag =0);
    friend class ArrayOfTaggedObjectsIter;
    
  protected:    
    
  private:
    int numComponents;          // num of components added
    int sizeComponentArray;     // size of the array
    int positionLastEntry;      // marker of last position used in the array
    int positionLastNoFitEntry; // marker of place array filled up to
    bool fitFlag;               // flag indicating if all components in nicely
    TaggedObject **theComponents; // the array
    ArrayOfTaggedObjectsIter  myIter; // an iter for accessing the objects
};

#endif


