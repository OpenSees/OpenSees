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
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidRod.h,v $
                                                                        
// Written: fmk 12/99
// Revised:
//
// Purpose: This file contains the class definition for RigidRod.
// RigidRod is a class which constructs MP_Constraint objects
// for a rigid rod, all translational dof are constrained to be equal
// at the retained and constarined nodes.

#ifndef RigidRod_h
#define RigidRod_h

class Domain;
class ID;

class RigidRod
{
  public:
    RigidRod(Domain &theDomain, int nodeR, int nodeC);
    virtual ~RigidRod();
    
  protected:
    
  private:
};

#endif
