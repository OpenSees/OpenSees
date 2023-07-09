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
// $Date: 2004-10-12 21:52:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/LagrangeDOF_Group.cpp,v $
                                                                        
// Written: fmk 
// Created: 02/99
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the LagrangeDOF_Group class interface.
//
// What: "@(#) LagrangeDOF_Group.C, revA"


#include <LagrangeDOF_Group.h>
#include <stdlib.h>

#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <TransientIntegrator.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>

LagrangeDOF_Group::LagrangeDOF_Group(int tag, SP_Constraint &spPtr)
:DOF_Group(tag, 1)
{

}

LagrangeDOF_Group::LagrangeDOF_Group(int tag, MP_Constraint &mpPtr)
:DOF_Group(tag, (mpPtr.getRetainedDOFs()).Size())
{

}


// ~LagrangeDOF_Group();    
//	destructor.

LagrangeDOF_Group::~LagrangeDOF_Group()
{

}    

const Matrix &
LagrangeDOF_Group::getTangent(Integrator *theIntegrator)
{
    // does nothing - the Lagrange FE_Elements provide coeffs to tangent
    if (tangent == 0) {
	int numDOF = this->getNumDOF();
	tangent = new Matrix(numDOF,numDOF);
	
	if (tangent == 0) {
	    opserr << "FATAL LagrangeDOF_Group::getTangent() ";
	    opserr << " ranout of memory\n";
	    exit(-1);
	}
    }
    
    tangent->Zero();
    return *tangent;
    
}

const Vector &
LagrangeDOF_Group::getUnbalance(Integrator *theIntegrator)
{
    // does nothing - the Lagrange FE_Elements provide residual 
    unbalance->Zero();
    return *unbalance;
}

// void setNodeDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
LagrangeDOF_Group::setNodeDisp(const Vector &u)
{
    return;
}
	
	
// void setNodeVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
LagrangeDOF_Group::setNodeVel(const Vector &udot)
{
    return;
}



// void setNodeAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;

void
LagrangeDOF_Group::setNodeAccel(const Vector &udotdot)
{
    return;
}


// void setNodeIncrDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
LagrangeDOF_Group::incrNodeDisp(const Vector &u)
{
    return;
}
	
	
// void setNodeincrVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
LagrangeDOF_Group::incrNodeVel(const Vector &udot)
{
    return;
}



// void setNodeIncrAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;

void
LagrangeDOF_Group::incrNodeAccel(const Vector &udotdot)
{
    return;
}


const Vector &
LagrangeDOF_Group::getCommittedDisp(void)
{
    unbalance->Zero();
    return *unbalance;
}

const Vector &
LagrangeDOF_Group::getCommittedVel(void)
{
    unbalance->Zero();
    return *unbalance;
}

const Vector &
LagrangeDOF_Group::getCommittedAccel(void)
{
    unbalance->Zero();
    return *unbalance;
}

void  
LagrangeDOF_Group::addMtoTang(double fact)
{
}

void  
LagrangeDOF_Group::zeroUnbalance(void)
{
}

void  
LagrangeDOF_Group::zeroTangent(void)
{
}

void  
LagrangeDOF_Group::addPtoUnbalance(double fact)
{

}

void  
LagrangeDOF_Group::addPIncInertiaToUnbalance(double fact)
{

}

void  
LagrangeDOF_Group::addM_Force(const Vector &Udotdot, double fact)
{

}

const Vector &
LagrangeDOF_Group::getTangForce(const Vector &disp, double fact)
{
  opserr << "WARNING LagrangeDOF_Group::getTangForce() - not yet implemented\n";
  unbalance->Zero();
  return *unbalance;
}

const Vector &
LagrangeDOF_Group::getC_Force(const Vector &disp, double fact)
{
  unbalance->Zero();
  return *unbalance;
}

const Vector &
LagrangeDOF_Group::getM_Force(const Vector &disp, double fact)
{
  unbalance->Zero();
  return *unbalance;
}


