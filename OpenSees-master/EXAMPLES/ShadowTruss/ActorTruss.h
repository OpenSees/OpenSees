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
// $Date: 2003-10-15 00:38:07 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ActorTruss.h,v $
                                                                        
#ifndef ActorTruss_h
#define ActorTruss_h

// Written: fmk 
// Created: 08/03
//
// Description: This file contains the interface for the ActorTruss class.
//
// What: "@(#) ActorTruss.h, revA"

#include <Actor.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class UniaxialMaterial;
class ElementalLoad;

class ActorTruss : public Actor
{
  public:
    // constructors
    ActorTruss(Channel &theChannel, 
	       FEM_ObjectBroker &theObjectBroker);

    // destructor
    ~ActorTruss();
    
    int run(void);
    int setMaterial(void);
    int setDomain(void);
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    int getTangentStiff(void);
    int getInitialStiff(void);
    int getResistingForce(void);

  protected:
    
  private:
    ID msgData;         // for receiving command

    Matrix trans;       // hold the transformation matrix, could use a Vector
	                // if ever bother to change the Vector interface for
			// x-product.

    double L;		// length of truss (undeformed configuration) - set in setDomain()
    double A; 		// area of truss

    UniaxialMaterial *theMaterial;       // pointer to a material

    // static data - single copy for all objects of the class
    static Matrix trussK;   // class wide matrix for returning stiffness
    static Vector trussR;   // class wide vector for returning residual
};
#endif
