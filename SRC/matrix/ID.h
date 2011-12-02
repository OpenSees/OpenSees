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
// $Date: 2000-12-12 07:14:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/ID.h,v $
                                                                        
                                                                        
// File: ~/matrix/ID.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for ID.
// ID is a concrete class implementing the integer array abstraction.
// ID objects are Vectors of integers which only need a few
// operators defined on them.
//
// What: "@(#) ID.h, revA"


#ifndef ID_h
#define ID_h

#include <iostream.h>
#include <G3Globals.h>

class ID
{
  public:
    // constructors and destructor
    ID();
    ID(int);
    ID(int size, int arraySize);    
    ID(const ID &);    
    ~ID();
 
    // utility methods
    int Size(void) const;
    void Zero(void);
    
    // overloaded operators
    inline int &operator()(int x);
    inline int operator()(int x) const;
    int &operator[](int);    	    
    
    ID &operator=(const ID  &V);
    
    int getLocation(int value) const;
    int removeValue(int value);

    friend ostream &operator<<(ostream &s, const ID &V);
    friend istream &operator>>(istream &s, ID &V);    

    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketNoDelay;
    friend class MPI_Channel;
    
  private:
    static int ID_NOT_VALID_ENTRY;
    int sz;
    int *data;
    int arraySize;
};


inline int 
ID::Size(void) const {return sz;}

inline int &
ID::operator()(int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
      g3ErrorHandler->warning("ID::(loc) - loc %d outside range [0, %d]\n",x,sz-1);
      return ID_NOT_VALID_ENTRY;
  }
#endif

  return data[x];
}

inline int
ID::operator()(int x) const 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
      g3ErrorHandler->warning("ID::(loc) - loc %d outside range [0, %d]\n",x,sz-1);
      return ID_NOT_VALID_ENTRY;
  }
#endif

  return data[x];
}

#endif


