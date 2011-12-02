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
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/SystemOfEqn.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/SystemOfEqn.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for SystemOfEqn.
// SystemOfEqn is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) SystemOfEqn.h, revA"

#ifndef SystemOfEqn_h
#define SystemOfEqn_h

#include <MovableObject.h>
#include <stdlib.h>

class Graph;
class AnalysisModel;
class FEM_ObjectBroker;

class SystemOfEqn : public MovableObject
{
  public:
    SystemOfEqn(int classTag);    
    virtual ~SystemOfEqn();

    virtual int solve(void) = 0;

  protected:
    
  private:

};

#endif


