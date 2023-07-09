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
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/MapOfTaggedObjects.h,v $
                                                                        
                                                                        
#ifndef MapOfTaggedObjects_h
#define MapOfTaggedObjects_h

// File: ~/tagged/storage/MapOfTaggedObjects.h
// 
// Written: fmk 
// Created: 02/00
// Revision: A
//
// Description: This file contains the class definition for 
// MapOfTaggedObjects. MapOfTaggedObjects is a storage class. The class 
// is responsible for holding and providing access to objects of type 
// TaggedObject. A map template of the standard template class is used to store
// the pointers to these objects.
//
// What: "@(#) MapOfTaggedObjects.h, revA"


#include <TaggedObjectStorage.h>
#include <MapOfTaggedObjectsIter.h>

class MapOfTaggedObjects : public TaggedObjectStorage
{
  public:
    MapOfTaggedObjects();
    ~MapOfTaggedObjects();    

    // public methods to populate a domain
    int  setSize(int newSize);
    bool addComponent(TaggedObject *newComponent);
//		      bool allowMutltipleTags = false);
    TaggedObject *removeComponent(int tag);    
    int getNumComponents(void) const;
    
    TaggedObject     *getComponentPtr(int tag);
    TaggedObjectIter &getComponents();

    MapOfTaggedObjectsIter getIter();
    
    TaggedObjectStorage *getEmptyCopy(void);
    void clearAll(bool invokeDestructor = true);
    
    void Print(OPS_Stream &s, int flag =0);
    friend class MapOfTaggedObjectsIter;
    
  protected:    
    
  private:
    map<int, TaggedObject *> theMap; // the map container for storing the pointers
    MapOfTaggedObjectsIter  myIter;  // the iter for this object
};

#endif



