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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/SP_Constraint.h,v $
                                                                        
                                                                        
#ifndef SP_Constraint_h
#define SP_Constraint_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for SP_Constraint.
// SP_Constraint is a class which stores the information for a single
// point constraint. Each single point constraint specifies a particular
// degree-of-freedom response (displacement, rotation) at a node.
// The constraint may be time-varying .. time varying constarints however 
// must be implemented using subclasses.
//
// What: "@(#) SP_Constraint, revA"

#include <DomainComponent.h>

class SP_Constraint : public DomainComponent
{
  public:
    // constructors    
    SP_Constraint(int classTag);        
    SP_Constraint(int nodeTag, int ndof, int classTag);    
    SP_Constraint(int nodeTag, int ndof, double value, bool isConstant);

    // destructor
    virtual ~SP_Constraint();

    virtual int getNodeTag(void) const;
    virtual int getDOF_Number(void) const;
    virtual int applyConstraint(double loadFactor);    
    virtual double getValue(void);
    virtual bool isHomogeneous(void) const;
    virtual void setLoadPatternTag(int loadPaternTag);
    virtual int  getLoadPatternTag(void) const;
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    virtual void Print(OPS_Stream &s, int flag =0);

  protected:
    int nodeTag;     // to identify the node in the model
    int dofNumber;   // identifies which of the nodes dof is constrrained 
    double valueR;   // the reference value
    double valueC;   // if constant = the reference value, if not contant =
	             // the reference value * load factor
    bool isConstant; // flag indicating if constant
    int  loadPatternTag;    
};

#endif


