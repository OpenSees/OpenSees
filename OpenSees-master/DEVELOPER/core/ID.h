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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-09-23 22:49:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/ID.h,v $
                                                                        
                                                                        
// Written: fmk 
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

#include <OPS_Globals.h>

class ID
{
  public:
    // constructors and destructor
    ID();
    ID(int);
    ID(int size, int arraySize);
    ID(int *data, int size, bool cleanIt = false);
    ID(const ID &);    
    ~ID();
 
    // utility methods
    int Size(void) const;
    void Zero(void);
    int setData(int *newData, int size, bool cleanIt = false);
    int resize(int newSize);
    
    // overloaded operators
    inline int &operator()(int x);
    inline int operator()(int x) const;
    int &operator[](int);    	    
    
    ID &operator=(const ID  &V);

    int operator==(const ID &V) const;
    int operator==(int) const;
    int operator!=(const ID &V) const;
    int operator!=(int) const;
    int operator<(const ID &V) const;

    int insert(int value);  // differs from using [] in that inserted in order
    int getLocation(int value) const;
    int getLocationOrdered(int value) const; // for when insert was used to add elements
    int removeValue(int value);
    int unique(void);

    friend OPS_Stream &operator<<(OPS_Stream &s, const ID &V);
    //    friend istream &operator>>(istream &s, ID &V);    

    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;
    friend class MPI_Channel;
    friend class MySqlDatastore;
    friend class BerkeleyDbDatastore;
    
  private:
    static int ID_NOT_VALID_ENTRY;
    int sz;
    int *data;
    int arraySize;
    int fromFree;
};


inline int 
ID::Size(void) const {return sz;}

inline int &
ID::operator()(int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
    opserr << "ID::(loc) - loc " << x << " outside range 0 - " <<  sz-1 << endln;
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
    opserr << "ID::(loc) - loc " << x << " outside range 0 - " <<  sz-1 << endln;
    return ID_NOT_VALID_ENTRY;
  }
#endif

  return data[x];
}

#endif


