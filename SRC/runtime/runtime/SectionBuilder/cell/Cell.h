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
//
// File: Cell.h
//
// Written by Remo M. de Souza
// December 1998


#ifndef Cell_h 
#define Cell_h 

#include <OPS_Globals.h>

class Vector;

class Cell
{
  public:

    Cell();
    virtual ~Cell();
    
    // edition functions

    // reinforcing bar inquiring functions
    
    virtual        double getArea              (void) const = 0;
    virtual        double getdValue              (void) const = 0;
    virtual const  Vector &getCentroidPosition (void) = 0;
 
    virtual void   Print(OPS_Stream &s, int flag =0) const = 0;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const Cell &Cell);    
    
  protected:
    
  private:
};


#endif

