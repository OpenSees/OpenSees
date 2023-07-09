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
                                                                        
// $Revision: 1.3 $
// $Date: 2009-05-11 20:52:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/LinearSOESolver.h,v $
                                                                        
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for LinearSOESolver.
// LinearSOESolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of LinearSOESolver are used
// to solve a system of equations.
//
// What: "@(#) LinearSOESolver.h, revA"

#ifndef LinearSOESolver_h
#define LinearSOESolver_h

#include <MovableObject.h>
class LinearSOE;

class LinearSOESolver : public MovableObject
{
  public:
    LinearSOESolver(int classTag);    
    virtual ~LinearSOESolver();

    virtual int solve(void) = 0;
    virtual int setSize(void) = 0;
    virtual double getDeterminant(void) {return 1.0;};
    
  protected:
    
  private:

};

#endif
    
