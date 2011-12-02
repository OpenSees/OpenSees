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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/LagrangeDOF_Group.h,v $
                                                                        
                                                                        
#ifndef LagrangeDOF_Group_h
#define LagrangeDOF_Group_h

// File: ~/analysis/dof_grp/LagrangeDOF_Group.h
// 
// Written: fmk 
// Created: 02/99
// Revision: A
//
// Description: This file contains the class definition for LagrangeDOF_Group.
// A LagrangeDOF_Group object is instantiated by a LagrangeConstraintHandler for 
// every constrained node in the domain. 
//
// What: "@(#) LagrangeDOF_Group.h, revA"

#include <DOF_Group.h>
class SP_Constraint;
class MP_Constraint;

class LagrangeDOF_Group: public DOF_Group
{
  public:
    LagrangeDOF_Group(int tag, SP_Constraint &spPtr);    
    LagrangeDOF_Group(int tag, MP_Constraint &mpPtr);        
    virtual ~LagrangeDOF_Group();    

    // methods to form the tangent
    virtual const Matrix &getTangent(Integrator *theIntegrator);

    // methods to form the unbalance
    virtual const Vector &getUnbalance(Integrator *theIntegrator);

    // methods to obtain committed responses .. always 0
    virtual const Vector &getCommittedDisp(void);
    virtual const Vector &getCommittedVel(void);
    virtual const Vector &getCommittedAccel(void);
    
    // methods to update the trial response at the nodes
    virtual void setNodeDisp(const Vector &u);
    virtual void setNodeVel(const Vector &udot);
    virtual void setNodeAccel(const Vector &udotdot);

    virtual void incrNodeDisp(const Vector &u);
    virtual void incrNodeVel(const Vector &udot);
    virtual void incrNodeAccel(const Vector &udotdot);

    virtual void  zeroTangent(void);
    virtual void  addMtoTang(double fact = 1.0);    
    virtual void  zeroUnbalance(void);
    virtual void  addPtoUnbalance(double fact = 1.0);
    virtual void  addPIncInertiaToUnbalance(double fact = 1.0);    
    virtual void  addM_Force(const Vector &Udotdot, double fact = 1.0);        
    
  protected:
    
  private:

};

#endif

