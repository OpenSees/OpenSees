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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/ID.cpp,v $
                                                                        
                                                                        
// File: ~/matrix/ID.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for ID.
//
// What: "@(#) ID.C, revA"

#include "ID.h"
#include <stdlib.h>


int ID::ID_NOT_VALID_ENTRY = 0;

// ID():
//	Standard constructor, sets size = 0;

ID::ID()
:sz(0), data(0), arraySize(0)
{

}


// ID(int size):
//	Constructor used to allocate a ID of size size.

ID::ID(int size)
  :sz(size), data(0), arraySize(size)
{

#ifdef _G3DEBUG
  if (sz <= 0) {
    g3ErrorHandler->warning("ID::ID(int) - size %d specified <= 0\n",size);
    sz = 1;
    arraySize = 1;
  }
#endif    

  // create the space for the data & check space was available
  data = (int *)malloc(size*sizeof(int));
  if (data == 0) {
    g3ErrorHandler->fatal("ID::ID(int): ran out of memory with size %d\n",
			  size);
    exit(-1);
  }

  // zero the data
  for (int i=0; i<size; i++)
    data[i] = 0;
}


// ID(int size):
//	Constructor used to allocate a ID of size size.

ID::ID(int size, int arraySz)
  :sz(size), data(0), arraySize(arraySz)
{
#ifdef _G3DEBUG
  if (sz < 0) {
    g3ErrorHandler->warning("ID::ID(size, arraySize) - size %d specified < 0\n",size);
    sz = 0;
  }
  if (arraySz <= 0) {
    g3ErrorHandler->warning("ID::ID(size, arraySize) - arraySize %d specified < 0\n",arraySz);
    if (sz != 0) 
      arraySz = sz;
    else
      arraySz = 1;
  }
  if (arraySz < sz) {
    g3ErrorHandler->warning("ID::ID(size, arraySize) - arraySize %d specified <  size\n",
			    arraySz, size);
    arraySz = sz;
  }
#endif    

  // create the space
  data = (int *)malloc(arraySize*sizeof(int));
  if (data == 0) {
    g3ErrorHandler->fatal("ID::ID(int, int): ran out of memory with arraySize %d\n",
			  arraySize);
    exit(-1);
  }

  // zero the data
  for (int i=0; i<arraySize; i++)
    data[i] = 0;
}



// ID(const ID&):
//	Constructor to init a ID from another.

ID::ID(const ID &other)
{
    sz = other.sz;
    arraySize = other.arraySize;
    data = 0;

  // create the space
  data = (int *)malloc(arraySize*sizeof(int));
  if (data == 0) {
    g3ErrorHandler->fatal("ID::ID(ID): ran out of memory with arraySize %d\n",
			  arraySize);
    exit(-1);
  }
  
  // copy the data 
  for (int i=0; i<sz; i++)
    data[i] = other.data[i];
}	



// ~ID():
// 	destructor, deletes the [] data

ID::~ID()
{
  if (data != 0) 
    free((void *)data);
}


void
ID::Zero(void)
{
  for (int i=0; i<sz; i++)
    data[i] =0;
}

int
ID::getLocation(int value) const
{
  // search through ID for the value
  for (int i=0; i<sz; i++)
    if (data[i] == value)
      return i;

  // if we get here the value is not in the array
  return -1;
}

int
ID::removeValue(int value)
{
  int place = -1;
  for (int i=0; i<sz; i++)
    if (data[i] == value) {
      place = i;
      // copy the rest of the components down one in ID
      for (int j=i; j<sz-1; j++)
	data[j] = data[j+1];		
      sz--;
    }
  return place;
}    


int &
ID::operator[](int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0) {
    g3ErrorHandler->warning("ID::[] - location x %d < 0\n",x);
    return ID_NOT_VALID_ENTRY;
  }
#endif

  // see if quick return
  if (x < sz)
    return data[x];
    
  /*
   * otherwise we have to enlarge the order of the ID
   */
    
  // see if we can just enlarge the array
  // without having to go get more space

  if (x < arraySize) {
    for (int i=sz; i<x; i++)
      data[i] = 0;
    sz = x+1;
    return data[x];
  }

  // otherwise we go get more space
  if (x >= arraySize) {
    int newArraySize = arraySize * 2;
    if (newArraySize < x) 
      newArraySize = x;
    int *newData = (int *)malloc(newArraySize*sizeof(int));    
    
    if (newData != 0) {
      // copy the old
      for (int i=0; i<sz; i++)
	newData[i] = data[i];
      // zero the new
      for (int j=sz; j<arraySize; j++)
	newData[j] = 0;
      
      sz = x+1;
      // release the memory held by the old
      free((void *)data);	    
      data = newData;
      arraySize = newArraySize;
      
      return newData[x];
    }
    else {
      // we could not allocate more mem .. leave the current size
      g3ErrorHandler->warning("ID::[]): ran out of memory with arraySize %d\n",
			      arraySize);
      return ID_NOT_VALID_ENTRY;
    }
  }
  
  // we should never get here, but some compilers need this line
  return ID_NOT_VALID_ENTRY;	
}
    


// ID &operator=(const ID  &V):
//	the assignment operator, This is assigned to be a copy of V. if sizes
//	are not compatable this.data [] is deleted. The data pointers will not
//	point to the same area in mem after the assignment.
//

ID &
ID::operator=(const ID &V) 
{
    // first check we are not trying v = v
    if (this != &V) {
	
	// check size compatability, if different delete
	// old and make room for new.
	if (sz != V.sz) {
	    if (arraySize < V.sz) {
		arraySize = V.sz;
		if (data != 0)
		    free((void *)data);
		sz = V.sz;
		data = (int *)malloc(arraySize*sizeof(int));		
		// check we got the memory requested
		if (data == 0) {
		    cerr << "WARNING ID::=(ID) - ran out of memory ";
		    cerr << "for new array of size" << arraySize << endl;
		    sz = 0;
		    arraySize = 0;
		}
	    }
	}
	
	// copy the data
	for (int i=0; i<sz; i++)
	    data[i] = V(i);
    }

    return *this;
}





// friend ostream &operator<<(ostream &s, const ID &V)
//	A function is defined to allow user to print the IDs using ostreams.

ostream &operator<<(ostream &s, const ID &V)
{
    for (int i=0; i<V.Size(); i++) 
    {
	s << V(i) << " ";
    }
    return s << "\n";
}

// friend istream &operator>>(istream &s, ID &V)
//	A function is defined to allow user to input the data into a ID which has already
//	been constructed with data, i.e. ID(int) or ID(const ID &) constructors.

istream &operator>>(istream &s, ID &V)
{
    for (int i=0; i<V.Size(); i++) 
	s >> V(i);

    return s;
}





