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
                                                                        
// $Revision: 1.18 $
// $Date: 2010-02-16 18:54:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/ID.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class implementation for ID.
//
// What: "@(#) ID.C, revA"

#include "ID.h"
#include <stdlib.h>
#include <map>
#include <list>

#include <iostream>
using std::nothrow;

int ID::ID_NOT_VALID_ENTRY = 0;

// ID():
//	Standard constructor, sets size = 0;

ID::ID()
  :sz(0), data(0), arraySize(0), fromFree(0)
{

}


// ID(int size):
//	Constructor used to allocate a ID of size size.

ID::ID(int size)
  :sz(size), data(0), arraySize(size), fromFree(0)
{

#ifdef _G3DEBUG
  if (sz < 0) {
    opserr << "ID::ID(int) - size " << size << " specified <= 0\n";
    sz = 0;
    arraySize = 0;
  }
#endif    

  // create the space for the data & check space was available
  //  data = (int *)malloc(size*sizeof(int));
  if (size > 0) {
    data = new (nothrow) int[size]; 
    if (data == 0) {
      opserr << "ID::ID(int): ran out of memory with size " << size << endln;
      exit(-1);
    }
    
    // zero the data
    for (int i=0; i<size; i++)
      data[i] = 0;
  }
}


// ID(int size):
//	Constructor used to allocate a ID of size size.

ID::ID(int size, int arraySz)
  :sz(size), data(0), arraySize(arraySz), fromFree(0)
{
#ifdef _G3DEBUG
  if (sz < 0) {
    opserr << "ID::ID(size, arraySize) - size " << size << " specified < 0\n";
    sz = 0;
  }
  if (arraySz <= 0) {
    opserr << "ID::ID(size, arraySize) - arraySize " << arraySz << " specified < 0\n";
    if (sz != 0) 
      arraySz = sz;
    else
      arraySz = 1;
  }
  if (arraySz < sz) {
    opserr << "ID::ID(size, arraySize) - arraySize " << arraySz  << " specified < " << size << endln;
    arraySz = sz;
  }
#endif    

  // create the space
  //  data = (int *)malloc(arraySize*sizeof(int));
  data = new (nothrow) int[arraySize];
  if (data == 0) {
    opserr << "ID::ID(int, int): ran out of memory with arraySize: " << arraySize << endln;
    exit(-1);
  }

  // zero the data
  for (int i=0; i<arraySize; i++)
    data[i] = 0;
}

ID::ID(int *d, int size, bool cleanIt)
  :sz(size), data(d), arraySize(size), fromFree(1)
{
  if (d == 0) { // OOPS must have been other constructor we wanted
    sz = 0;
    data = 0;
    arraySize = size;
    fromFree = 0;

    // create the space
    if (arraySize != 0) {
      data = (int *)malloc(arraySize*sizeof(int));
      if (data == 0) {
	opserr << "ID::ID(int, int): ran out of memory with arraySize " << arraySize << endln;
	exit(-1);
      }
    }

    // zero the data
    for (int i=0; i<arraySize; i++)
      data[i] = 0;
  }
  
  if (cleanIt == true)
    fromFree = 0;
}

// ID(const ID&):
//	Constructor to init a ID from another.

ID::ID(const ID &other)
  :sz(other.sz), data(0), arraySize(other.arraySize), fromFree(0)
{
  // create the space
  //  data = (int *)malloc(arraySize*sizeof(int));
  data = new (nothrow) int[arraySize]; 
  if (data == 0) {
    opserr << "ID::ID(ID): ran out of memory with arraySize " << arraySize << endln,
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
  if (data != 0 && fromFree == 0) 
    //    free((void *)data);
    delete [] data;
}

int 
ID::setData(int *newData, int size, bool cleanIt){
	
  if (data != 0 && fromFree == 0) 
    //    free((void *)data);
    delete [] data;

  sz = size;
  data = newData;
  
  if (cleanIt == false)
    fromFree = 1;
  else
    fromFree = 0;

  if (sz <= 0) {
    opserr << "ID::ID(int *, size) - size " << size << " specified <= 0\n";
    sz = 0;
  }

  return 0;
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
ID::getLocationOrdered(int value) const
{
  int middle = 0;
  int left = 0;
  int right = sz-1;
  if (sz != 0) {
    while (left <= right) {
      middle = (left + right)/2;
      double dataMiddle = data[middle];
      if (value == dataMiddle)
	return middle;   // already there
      else if (value > dataMiddle)
	left = middle + 1;
      else 
	right = middle-1;
    }
  }

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


int
ID::unique(void)
{
    // preserve order method, runs in O(nlogn)
    // order is preserved via list O(1)
    // map serves as dictionary O(logn)
    std::map<int,int> uniquesm;
    std::list<int> uniquesl;
    int count = 0;
    for (int i=0; i<sz; i++)  {
        int tmp = data[i];
        if (uniquesm.find(tmp) == uniquesm.end()) {
            uniquesm[tmp] = tmp;
            uniquesl.push_back(tmp);
        }
    }

    sz = uniquesl.size();
    int* newdata = new (nothrow) int[sz];
    for (std::list<int>::iterator pos=uniquesl.begin(); pos!=uniquesl.end(); pos++)
        newdata[count++] = *pos;

    delete [] data;
    arraySize = sz;
    data = newdata;

    return sz;
}


int &
ID::operator[](int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0) {
    opserr << "ID::[] - location " << x << " < 0\n";
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
    if (newArraySize <= x) 
      newArraySize = x+1;
    //    int *newData = (int *)malloc(newArraySize*sizeof(int));    
    int *newData = new (nothrow) int[newArraySize];

    if (newData != 0) {

      // copy the old
      for (int i=0; i<sz; i++)
	newData[i] = data[i];

      // zero the new
      for (int j=sz; j<newArraySize; j++)
	newData[j] = 0;
      
      sz = x+1;

      // release the memory held by the old
      //      free((void *)data);	    

      if (fromFree == 0 && data != 0)
	delete [] data;

      data = newData;
      arraySize = newArraySize;
      
      return newData[x];
    }
    else {
      // we could not allocate more mem .. leave the current size
      opserr << "ID::[]): ran out of memory with arraySize " << arraySize << endln;
      return ID_NOT_VALID_ENTRY;
    }
  }
  
  // we should never get here, but some compilers need this line
  return ID_NOT_VALID_ENTRY;	
}
    

int 
ID::resize(int newSize){

  // first check that newSize is valid
  if (newSize <= 0) {
    opserr << "ID::resize() - size specified " << newSize << " <= 0\n";
    return -1;
  } 
  

  if (sz >= newSize) {

    // is size smaller than current, simply reset sz
    sz = newSize;

  } else if (newSize <= arraySize) {

    // see if we can just enlarge the array
    // without having to go get more space
    
    for (int i=sz; i<newSize; i++)
      data[i] = 0;
    sz = newSize;

  } else if (newSize > arraySize) {

    // otherwise we go get more space
    
    int *newData = new (nothrow) int[newSize];
    if (newData != 0) {
      // copy the old
      for (int i=0; i<sz; i++)
	newData[i] = data[i];
      // zero the new
      for (int j=sz; j<newSize; j++)
	newData[j] = 0;
      
      sz = newSize;
      // release the memory held by the old
      //      free((void *)data);	    
      delete [] data;
      data = newData;
      arraySize = newSize;

    } else {
      opserr << "ID::resize() - out of memory creating ID of size " << newSize << "\n";
      return -1;      
    }
  }

  return 0;
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
		  //free((void *)data);
		  delete [] data;
		//		data = (int *)malloc(arraySize*sizeof(int));		
		data = new (nothrow) int[arraySize];
		// check we got the memory requested
		if (data == 0) {
		    opserr << "WARNING ID::=(ID) - ran out of memory ";
		    opserr << "for new array of size" << arraySize << endln;
		    sz = 0;
		    arraySize = 0;
		}
	    }
	    sz = V.sz;
	}
	
	// copy the data
	for (int i=0; i<sz; i++)
	    data[i] = V(i);
    }

    return *this;
}


// ID operator==(const ID &V):
//	The == operator checks the two IDs are of the same size.
// 	Then returns 1 if all the components of the two IDs are equal and 0 otherwise.

int 
ID::operator==(const ID &V) const
{
  if (sz != V.sz)
    return 0;

  int *dataThis = data;
  int *dataV = V.data;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 0;

  return 1;
}

int 
ID::operator==(int value) const
{
  int *dataThis = data;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 0;

  return 1;
}


// ID operator!=(const ID &V):
//	The != operator checks the two IDs are of the same size.
// 	Then returns 1 if any of the components of the two IDs are unequal and 0 otherwise.

int 
ID::operator!=(const ID &V) const
{
  if (sz != V.sz)
    return 1;

  int *dataThis = data;
  int *dataV = V.data;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 1;

  return 0;
}

int 
ID::operator!=(int value) const
{
  int *dataThis = data;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 1;

  return 0;
}

// ID operator!=(const ID &V):
//	The != operator checks the two IDs are of the same size.
// 	Then returns 1 if first component of the first ID is less than the second ID

int
ID::operator<(const ID &V) const
{
    if (sz != V.sz) {
	return sz < V.sz;
    }

    for (int i=0; i<sz; i++) {
	if (data[i] < V.data[i]) {
	    return 1;
	} else if (data[i] > V.data[i]) {
	    return 0;
	}
    }

    return 0;
}

// friend OPS_Stream &operator<<(OPS_Stream &s, const ID &V)
//	A function is defined to allow user to print the IDs using OPS_Streams.

OPS_Stream &operator<<(OPS_Stream &s, const ID &V)
{
    for (int i=0; i<V.Size(); i++) 
    {
	s << V(i) << " ";
    }
    return s << endln;
}

// friend istream &operator>>(istream &s, ID &V)
//	A function is defined to allow user to input the data into a ID which has already
//	been constructed with data, i.e. ID(int) or ID(const ID &) constructors.

/*
istream &operator>>(istream &s, ID &V)
{
    for (int i=0; i<V.Size(); i++) 
	s >> V(i);

    return s;
}
*/



int 
ID::insert(int x) 
{
  int middle = 0;
  int left = 0;
  int right = sz-1;
  if (sz != 0) {
    while (left <= right) {
      middle = (left + right)/2;
      double dataMiddle = data[middle];
      if (x == dataMiddle)
	return 1;   // already there
      else if (x > dataMiddle)
	left = middle + 1;
      else 
	right = middle-1;
    }
  }

  // we need to enlarge the array .. see if we can do it
  // without having to go get more space

  middle = left;

  /*
  for (int i=middle; i<sz; i++)
  (*this)[i+1] = (*this)[i];
  (*this)[i]=middle;
  return middle;
  */
  
  if (sz < arraySize) {

    int i = sz;
    while (i > middle) {
      data[i] = data[i-1];
      i--;
    }
    sz++;
    data[i] = x;
    return 0;
  } else {
    int newArraySize = (sz+1) * 2;
    int *newData = new (nothrow) int[newArraySize];
    if (newData != 0) {
      
      // copy the old
      for (int ii=0; ii<middle; ii++)
	newData[ii] = data[ii];
      newData[middle] = x;
      
      for (int jj=middle; jj<sz; jj++)
	newData[jj+1] = data[jj];
      
      sz++;
      
      if (data != 0 && fromFree == 0)
	delete [] data;
      data = newData;
      arraySize = newArraySize;
      
      return 0;
      
    }
  }
  return -1; 
}

