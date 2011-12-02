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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidBeam.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidBeam.C
//
// Written: fmk 12/99
// Revised:
//
// Purpose: This file contains the class implementation for RigidBeam.

#include <G3Globals.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <RigidBeam.h>


RigidBeam::RigidBeam(Domain &theDomain, int nR, int nC, int startMPtag) {

    
    // get a pointer to the retained and constrained nodes - make sure they exist
    Node *nodeR = theDomain.getNode(nR);
    if (nodeR == 0) {
      g3ErrorHandler->warning("RigidBeam::RigidBeam - %s %d %s\n",
			      "retained Node", nR, "not in domain");
      return;
    }
    Node *nodeC = theDomain.getNode(nC);
    if (nodeR == 0) {
      g3ErrorHandler->warning("RigidBeam::RigidBeam - %s %d %s\n",
			      "constrained Node", nC, "not in domain");
      return;
    }

    // get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
    const Vector &crdR = nodeR->getCrds();
    const Vector &crdC = nodeC->getCrds();
    int dimR = crdR.Size();
    int dimC = crdC.Size();
    if (dimR != dimC) {
      g3ErrorHandler->warning("RigidBeam::RigidBeam - mismatch in dimension %s %d %s %d\n",
			      "between constrained Node", nC, "and Retained node",nR);
      return;
    }
    
    // check the number of dof at each node is the same
    int numDOF = nodeR->getNumberDOF();
    if (numDOF != nodeC->getNumberDOF()){ 
      g3ErrorHandler->warning("RigidBeam::RigidBeam - mismatch in numDOF %s %d %s %d\n",
			      "between constrained Node", nC, "and Retained node",nR);
      return;
    }

    // check the number of dof at the nodes >= dimension of problem
    if(numDOF < dimR){    
      g3ErrorHandler->warning("RigidBeam::RigidBeam - numDOF at nodes %d %d %s\n",
			      nR, nC, "must be >= dimension of problem");
      return;
    }

    
    // create the ID to identify the constrained dof 
    ID id(numDOF);

    // construct the tranformation matrix Ccr, where  Uc = Ccr Ur & set the diag, Ccr = I
    Matrix mat(numDOF,numDOF);
    mat.Zero();

    // set the values
    for (int i=0; i<numDOF; i++) {
      mat(i,i) = 1.0;
      id(i) = i;
    }

    // if there are rotational dof - we must modify Ccr DONE ASSUMING SMALL ROTATIONS
    if (dimR != numDOF) {
      if (dimR == 2 && numDOF == 3) {
	double deltaX = crdC(0) - crdR(0);
	double deltaY = crdC(1) - crdR(1);	    
	mat(0,2) = -deltaY;
	mat(1,2) = deltaX;
      } else if (dimR == 3 && numDOF == 6) {
	double deltaX = crdC(0) - crdR(0);
	double deltaY = crdC(1) - crdR(1);	    
	double deltaZ = crdC(2) - crdR(2);
	// rotation about z/3 axis
	mat(0,5) = -deltaY;
	mat(1,5) = deltaX;	

	// rotation about y/2 axis
	mat(0,4) = deltaZ;
	mat(2,4) = -deltaX;

	// rotation about x/1 axis
	mat(1,3) = -deltaZ;
	mat(2,3) = deltaY;
      } else { // not valid
	g3ErrorHandler->warning("RigidBeam::RigidBeam -  for nodes %d %d %s\n",
				nR, nC, "nodes do not have valid numDOF for their dimension");
	return;
      }
	
    }	
    
    // create the MP_Constraint
    MP_Constraint *newC = new MP_Constraint(startMPtag+1, nR, nC, 
					    mat, id, id);
    if (newC == 0) {
      g3ErrorHandler->warning("RigidBeam::RigidBeam - for nodes %d %d, out of memory\n",
			      nC, nR);
    } else {
      // add the constraint to the domain
      if (theDomain.addMP_Constraint(newC) == false) {
	g3ErrorHandler->warning("RigidBeam::RigidBeam - for nodes %d %d, could not add to domain\n",
				nC, nR);
	delete newC;
      }
    }
}
	
    
RigidBeam::~RigidBeam()
{
    // does nothing
}

 

