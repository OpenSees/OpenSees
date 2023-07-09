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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/Pressure_ConstraintIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/Pressure_ConstraintIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// Pressure_ConstraintIter. Pressure_ConstraintIter is an abstract base class.
// An Pressure_ConstraintIter object is an iter for returning  the
// single point constraints  of an object of class Domain. 
// Pressure_ConstraintIters must be written for each subclass of Domain.

#ifndef Pressure_ConstraintIter_h
#define Pressure_ConstraintIter_h

class Pressure_Constraint;

class Pressure_ConstraintIter
{
  public:
    Pressure_ConstraintIter() {};
    virtual ~Pressure_ConstraintIter() {};
    
    virtual Pressure_Constraint *operator()(void) =0;
    
  protected:

  private:

};

#endif

