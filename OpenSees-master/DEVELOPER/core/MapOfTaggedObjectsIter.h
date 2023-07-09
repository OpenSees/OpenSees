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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/MapOfTaggedObjectsIter.h,v $
                                                                        
                                                                        
#ifndef MapOfTaggedObjectsIter_h
#define MapOfTaggedObjectsIter_h

// File: ~/tagged/storage/ArrayComponentIter
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// MapOfTaggedObjectsIter. A MapOfTaggedObjectsIter is an iter for 
// returning the TaggedObjects of a storage objects of type 
// MapOfTaggedComponents.

#include <TaggedObjectIter.h>

#include <map>
using namespace std;

class MapOfTaggedObjects;

class MapOfTaggedObjectsIter: public TaggedObjectIter
{
  public:
    MapOfTaggedObjectsIter(MapOfTaggedObjects &theComponents);
    virtual ~MapOfTaggedObjectsIter();
    
    virtual void reset(void);
    virtual TaggedObject *operator()(void);
    
  private:
    map<int, TaggedObject *> *theMap;
    map<int, TaggedObject *>::iterator currentComponent;
};

#endif





