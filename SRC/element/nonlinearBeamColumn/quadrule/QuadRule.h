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
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/QuadRule.h,v $
                                                                        
                                                                        
// File: ~/QuadRule/QuadRule.h
//
// Written: rms
// Created: 12/98
// Revision: 
//
// Description: This file contains the class definition for 
// QuadRule (Quadrature Rule). QuadRule is an abstract base class and 
// thus no objects of  it's type can be instatiated. It has pure 
// virtual functions which  must be implemented in it's derived classes.
//
// What: "@(#) QuadRule.h, revA"


#ifndef QuadRule_h
#define QuadRule_h

class Vector;
class Matrix;

class QuadRule
{
  public:
    QuadRule ();
    virtual ~QuadRule();

    virtual int            setOrder           (int quadOrder) = 0;
    virtual int            getOrder              (void) const = 0;
    virtual int            getNumIntegrPoints    (void) const = 0;
    virtual const Matrix & getIntegrPointCoords  (void) const = 0;
    virtual const Vector & getIntegrPointWeights (void) const = 0; 
    
  protected:
    
  private:
};


#endif

