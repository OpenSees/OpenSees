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
// $Date: 2004-09-01 04:01:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/MP_Joint3D.h,v $
                                                                        
#ifndef MP_Joint3D_h
#define MP_Joint3D_h

// Written: Arash Altoontash, Gregory Deierlein
// Created: 04/03
// Revision: Arash

// Purpose: This file contains the class definition for MP_Joint3D.
// It is a sub class for MP_Constraint specialized to be used for simple joint 
// connection element. MP_Joint3D defines a nonlinear, time dependent multi 
// point constraint.
//

#include <DomainComponent.h>
#include <bool.h>
#include <MP_Constraint.h>
#include <Node.h>
#include <Domain.h>

class Matrix;
class ID;


class MP_Joint3D : public MP_Constraint
{
  public:
    // constructors        
    MP_Joint3D();

    MP_Joint3D( Domain *theDomain, int tag, int nodeRetain, int nodeConstr,
		int nodeRot, int Rotdof, int nodeDisp, int Dispdof, int LrgDsp = 0 );

    // destructor    
    ~MP_Joint3D();

    // method to get information about the constraint
    int getNodeRetained(void) const;
    int getNodeConstrained(void) const;    
    const ID &getConstrainedDOFs(void) const;        
    const ID &getRetainedDOFs(void) const;            
    int applyConstraint(double pseudoTime);
    bool isTimeVarying(void) const;
    const Matrix &getConstraint(void);    
	void setDomain(Domain *theDomain);

    // methods for output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);


  protected:
    
  private:
    int nodeRetained;
    int nodeConstrained;
    int nodeRotation;	       // tag for the node to define the rotation vector
    // for shear rotation
    int RotDOF;		       // tag for the shear mode that results in rotation
    int nodeDisplacement;      // tag for the node to define the rotation vector
    // for shear displacement
    int DispDOF;	       // tag for the shear mode that results in displacement
    int LargeDisplacement;     // flag for large displacements
    // 0 for constant constraint matrix(small deformations)
    // 1 for time varying constraint matrix(large deformations)
    // 2 for large deformations with length correction
    ID *constrDOF;	       // ID of constrained DOF at constrained node
    ID *retainDOF;	       // ID of related DOF at retained node
    Node *RetainedNode;	       // to identify the retained node
    Node *ConstrainedNode;     // to identify  the constrained node
    Node *RotationNode;
    Node *DisplacementNode;
    
    Vector RotNormVect;
    Vector DspNormVect;
    
    int dbTag1, dbTag2, dbTag3;	// need a dbTag for the two ID's
    
    double Length0;
    Matrix *constraint;	        // pointer to the constraint matrix
    Domain *thisDomain;		// pointer to domain the MP is defined on
};

#endif

