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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/SP_ConstraintIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/SP_ConstraintIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SP_ConstraintIter. SP_ConstraintIter is an abstract base class.
// An SP_ConstraintIter object is an iter for returning  the
// single point constraints  of an object of class Domain. 
// SP_ConstraintIters must be written for each subclass of Domain.

#ifndef SP_ConstraintIter_h
#define SP_ConstraintIter_h

class SP_Constraint;

class SP_ConstraintIter
{
  public:
    SP_ConstraintIter() {};
    virtual ~SP_ConstraintIter() {};
    
    virtual SP_Constraint *operator()(void) =0;
    
  protected:

  private:

};

#endif

