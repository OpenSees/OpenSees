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
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidDiaphragm.h,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidDiaphragm.h
//
// Written: fmk 1/99
// Revised:
//
// Purpose: This file contains the class definition for RigidDiaphragm.
// RigidDiaphragm is a class which constructs MP_Constraint objects
// for a 3d Frame with a rigid diaphragm .. suitable for small
// displacement problems only.

#ifndef RigidDiaphragm_h
#define RigidDiaphragm_h

class Domain;
class ID;

class RigidDiaphragm
{
  public:
    RigidDiaphragm(Domain &theDomain, int nodeR, ID &nodeC, 
		   int perpDirnToPlaneConstrained);
    virtual ~RigidDiaphragm();
    
  protected:
    
  private:
};

#endif
