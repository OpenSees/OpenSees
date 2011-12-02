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
                                                                        
// $Revision: 1.1 $
// $Date: 2000-12-12 07:26:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/ImposedMotionSP.h,v $
                                                                        
#ifndef ImposedMotionSP_h
#define ImposedMotionSP_h

// File: ~/domain/constraints/ImposedMotionSP.h
//
// Written: fmk 
// Created: 11/00
// Revision: A
//
// Purpose: This file contains the class definition for ImposedMotionSP.
// ImposedMotionSP is a class which imposes the ground motion response
// values for a GroundMotion at a particular dof at a node.
//
// What: "@(#) ImposedMotionSP, revA"

#include <SP_Constraint.h>
#include <Vector.h>
class GroundMotion;
class Node;

class ImposedMotionSP : public SP_Constraint
{
  public:
    // constructors    
    ImposedMotionSP();        
    ImposedMotionSP(int spTag, int nodeTag, int ndof, 
		    GroundMotion &theGroundMotion, bool destroyMotion);

    // destructor
    ~ImposedMotionSP();

    int applyConstraint(double loadFactor);    
    double getValue(void);
    bool isHomogeneous(void) const;
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);

  protected:
    
  private:
    GroundMotion *theGroundMotion;  // pointer to ground motion
    Node *theNode;                  // pointer to node being constrained
    Vector *theNodeResponse;        // vector for setting nodal response
    Vector theGroundMotionResponse; // the GMotions response
    int destroyMotion;              // flag indiacting if destructor to be invoked
	                            // on GMotion when done
};

#endif


