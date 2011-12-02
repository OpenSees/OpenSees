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
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/QuadRule1d.h,v $
                                                                        
                                                                        
// File: ~/QuadRule/QuadRule1d.h
//
// Written: rms
// Created: 12/98
// Revision: 
//
// Description: This file contains the class definition for 
// QuadRule1d (Quadrature Rule 1dim). QuadRule1d is an abstract base class 
// and thus no objects of  it's type can be instatiated. It has pure 
// virtual functions which  must be implemented in it's derived classes.
//
// What: "@(#) QuadRule1d.h, revA"


#ifndef QuadRule1d_h
#define QuadRule1d_h


#include <QuadRule.h>

class Vector;
class Matrix;

class QuadRule1d: public QuadRule
{
  public:
    QuadRule1d ();
    virtual ~QuadRule1d();
    
  protected:
    
  private:
};


#endif

