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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidBeam.h,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidBeam.h
//
// Written: fmk 12/99
// Revised:
//
// Purpose: This file contains the class definition for RigidBeam.
// RigidBeam is a class which constructs an MP_Constraint object
// between two nodes which is similar to rigid beam

#ifndef RigidBeam_h
#define RigidBeam_h

class Domain;
class ID;

class RigidBeam
{
  public:
    RigidBeam(Domain &theDomain, int nodeR, int nodeC, int startMPtag);
    virtual ~RigidBeam();
    
  protected:
    
  private:
};

#endif
