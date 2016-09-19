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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-08-21 12:59:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/Pressure_Constraint.h,v $

// Written: Minjie
// Created: 08/12
// Revision: A
//
// Purpose: This file contains the class definition for Pressure_Constraint.
// Pressure_Constraint is a class which stores the information for a pressure
// constraint. A pressure constraint constrains the model with an incompressible
// or compressible continuity equation. In this way, the fluid is treated
// as solid with pressure constraints using Lagrangian formulation. 
//
// The Pressure_Constraint class is part of the implementation of PFEM.
// Nodel pressures are stores in the class. In each iteration,
// the pressure values are updated. Pressure values can be obtained from
// the Pressure_Constraint class.
//
// This class is different to other constraints classes. It does not contrain
// nodes to values or other nodes. It applies the continuity equations as
// constraints. So it has its own states like nodes and can update itself.
// For incompressible flows, or near incompressible flows, the pressure acts
// as a constraint upon the velocities to make the velocity field divergence-free.


// What: "@(#) Pressure_Constraint, revA"
                                                                        
                                                                        
#ifndef Pressure_Constraint_h
#define Pressure_Constraint_h

#include <DomainComponent.h>
#include <ID.h>

class Node;

class Pressure_Constraint : public DomainComponent
{
public:
    // constructors
    explicit Pressure_Constraint(int classTag);
    Pressure_Constraint(int classTag, int nodeId, int ptag, int ndf);
    Pressure_Constraint(int nodeId, int ptag, int ndf);
    Pressure_Constraint(int nodeId, int ndf);

    // destructor
    virtual ~Pressure_Constraint();

    // method to get information about the constraint
    virtual void setDomain(Domain* theDomain);
    virtual Node* getPressureNode();
    virtual double getPressure(int last=1);
    virtual void setPressure(double p);
    virtual const ID& getFluidElements();
    virtual const ID& getOtherElements();
    virtual void connect(int eleId, bool fluid=true);
    virtual void disconnect(int eleId);
    virtual void disconnect();
    virtual bool isFluid() const;
    virtual bool isInterface() const;
    virtual bool isStructure() const;
    virtual bool isIsolated() const;

    // methods for output
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
    virtual void Print(OPS_Stream &s, int flag =0);

private:

    static int findNodeTag(Domain* theDomain);

    int pTag;
    ID fluidEleTags;
    ID otherEleTags;
    int pndf;
};

#endif
