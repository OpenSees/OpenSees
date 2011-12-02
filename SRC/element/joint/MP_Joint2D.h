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
                                                                        
// $Revision: 1.4 $
// $Date: 2004-09-01 04:01:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/MP_Joint2D.h,v $
                                                                        
#ifndef MP_Joint2D_h
#define MP_Joint2D_h

// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/01
// Revision: Arash

// Purpose: This file contains the class definition for MP_Joint2D.
// It is a sub class for MP_Constraint specialized to be used for simple joint 
// connection element. MP_Joint2D defines a nonlinear, time dependent multi 
// point constraint.
//

#include <DomainComponent.h>
#include <bool.h>
#include <MP_Constraint.h>
#include <Node.h>
#include <Domain.h>

class Matrix;
class ID;


class MP_Joint2D : public MP_Constraint
{
  public:
    // constructors        
    MP_Joint2D();

    MP_Joint2D( Domain *theDomain, int tag, int nodeRetain, int nodeConstr,
		int Maindof, int fixedend , int LrgDsp = 0 );	//LrgDsp=0 means large displacement is not enabled

    // destructor    
    ~MP_Joint2D();

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
    int nodeRetained;			// to identify the retained node
    int nodeConstrained;		// to identify  the constrained node
	int MainDOF;				// main degree of freedom for rotation
	int AuxDOF;					// Auxilary degree of freedom for shear
	int FixedEnd;				// fixed rotational degree of freedom at the end 
                                // released = 0 , fixed = 1
    
    ID *constrDOF;				// ID of constrained DOF at constrained node
    ID *retainDOF;				// ID of related DOF at retained node
    Node *RetainedNode;			// to identify the retained node
    Node *ConstrainedNode;		// to identify  the constrained node

    int dbTag1, dbTag2, dbTag3;			// need a dbTag for the two ID's
	int LargeDisplacement;		// flag for large displacements enabled
	double Length0;
	Matrix *constraint;			// pointer to the constraint matrix
	Domain *thisDomain;			// pointer to domain the MP is defined on
};

#endif

