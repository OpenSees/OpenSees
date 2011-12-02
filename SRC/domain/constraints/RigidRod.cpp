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
// $Date: 2005-01-08 01:22:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidRod.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidRod.C
//
// Written: fmk 12/99
// Revised:
//
// Purpose: This file contains the class implementation for RigidRod.
#include <stdlib.h>

#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <RigidRod.h>


RigidRod::RigidRod(Domain &theDomain, int nR, int nC, int mPtag) {

    
    // get a pointer to the retained node and constrained nodes - ensure these exist
    Node *nodeR = theDomain.getNode(nR);
    if (nodeR == 0) {
      opserr << "RigidRod::RigidRod - retained Node" <<  nR <<  "not in domain\n";
      return;
    }
    Node *nodeC = theDomain.getNode(nC);
    if (nodeR == 0) {
      opserr << "RigidRod::RigidRod - constrained Node" <<  nC <<  "not in domain\n";
      return;
    }

    // get the coordinates of the two nodes - check dimensions are the same
    const Vector &crdR = nodeR->getCrds();
    const Vector &crdC = nodeC->getCrds();
    int dimR = crdR.Size();
    int dimC = crdC.Size();
    if (dimR != dimC) {
      opserr << "RigidRod::RigidRod - mismatch in dimension " <<
	"between constrained Node " <<  nC <<  " and Retained node " << nR << endln;
      return;
    }
    
    // check the number of dof at each node is the same 
    int numDOF = nodeR->getNumberDOF();
    if (numDOF != nodeC->getNumberDOF()){ 
      opserr << "RigidRod::RigidRod - mismatch in numDOF " <<
	"between constrained Node " <<  nC <<  " and Retained node " << nR << endln;
      return;
    }

    // check the number of dof at the nodes >= dimension of problem
    if(numDOF < dimR){    
      opserr << "RigidRod::RigidRod - numDOF at nodes " << nR << " and " << nC <<
	"must be >= dimension of problem\n";

      return;
    }

    
    // create the ID to identify the constrained dof 
    ID id(dimR);

    // construct the tranformation matrix Ccr, where  Uc = Ccr Ur & set the diag
    Matrix mat(dimR,dimR);
    mat.Zero();

    // set the values
    for (int i=0; i<dimR; i++) {
      mat(i,i) = 1.0;
      id(i) = i;
    }

    // create the MP_Constraint
    MP_Constraint *newC = new MP_Constraint(mPtag, nR, nC, mat, id, id);
					    
    if (newC == 0) {
      opserr << "RigidRod::RigidRod - for nodes " << nR << " and " << nC << " out of memory\n";
      exit(-1);
    } else {
      // add the constraint to the domain
      if (theDomain.addMP_Constraint(newC) == false) {
	opserr << "RigidRod::RigidRod - for nodes " << nC << " and " << nR << " could not add to domain\n",
	delete newC;
      }
    }
}
	
    
RigidRod::~RigidRod()
{
    // does nothing
}

 

