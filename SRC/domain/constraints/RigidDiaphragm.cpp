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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidDiaphragm.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidDiaphragm.C
//
// Written: fmk 1/99
// Revised:
//
// Purpose: This file contains the class implementation for RigidDiaphragm.

#include <G3Globals.h>
#include <stdlib.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <RigidDiaphragm.h>


RigidDiaphragm::RigidDiaphragm(Domain &theDomain, int nR, ID &nC, 
			       int perpPlaneConstrained, int startMPtag) {

    // check plane is valid, i.e. perpPlaneConstrained must be 0, 1 or 2
    if (perpPlaneConstrained < 0 || perpPlaneConstrained > 2) {
      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - %s %d %s\n",
			      "the dirn of perpendicular to constrained plane",
			      perpPlaneConstrained, "not valid");
      return;
    }

    // check constrainedNodes ID does not contain the retained node
    if (nC.getLocation(nR) >= 0) {
      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - %s %d %s\n",
			      "retained node", 
			      nR, 
			      "is in constrained node list");
      return;
    }	
    
    // get a pointer to the retained node and check node in 3d with 6 dof
    Node *nodeR = theDomain.getNode(nR);
    if (nodeR == 0) {
      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - %s %d %s\n",
			      "retained Node", nR, "not in domain",nR);
      return;
    }

    const Vector &crdR = nodeR->getCrds();
    if ((nodeR->getNumberDOF() != 6) || (crdR.Size() != 3)){
      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - %s %d %s\n",
			      "retained Node",
			      nR,
			      "not in 3d space with 6 dof");
      return;
    }	

    // 
    // create some objects which will be passed to the MP_Constraint 
    // constructor, elements of objects are filled in later
    //
    
    // create the ID to identify the constrained dof 
    ID id(3);

    // construct the tranformation matrix Ccr, where  Uc = Ccr Ur & set the diag
    Matrix mat(3,3);
    mat.Zero();
    mat(0,0) = 1.0; mat(1,1) = 1.0; mat(2,2) = 1.0;


    // now for each of the specified constrained dof we:
    // 1. check it's in the plane of the retained node, 
    // 2. set the ID and transformation matrix,
    // 3. create the MP_Constrainet and add it to the domain

    for (int i=0; i<nC.Size(); i++) {

      // get the constrained node
      int ndC = nC(i);
      Node *nodeC = theDomain.getNode(ndC);

      // ensure node exists
      if (nodeC != 0) {

	// get node coordinates
	const Vector &crdC = nodeC->getCrds();

	// check constrained node has correct dim and number of dof
	if ((nodeR->getNumberDOF() == 6) && (crdR.Size() == 3)){

	  // determine delta Coordintaes
	  double deltaX = crdC(0) - crdR(0);
	  double deltaY = crdC(1) - crdR(1);	    
	  double deltaZ = crdC(2) - crdR(2);
	  
	  // rigid diaphragm in xy plane
	  if (perpPlaneConstrained == 2) { 

	    // check constrained node in xy plane with retained node
	    if (deltaZ == 0.0) {

	      // dof corresponding to dX, dY and theta Z (0,1,5)
	      id(0) = 0; id(1) = 1; id(2) = 5;

	      // set up transformation matrix
	      mat(0,2) = - deltaY;
	      mat(1,2) = deltaX;

	    } else 
	      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, not in xy plane\n",ndC);

	  // rigid diaphragm in xz plane
	  } else if (perpPlaneConstrained == 1) { 

	    // check constrained node in xy plane with retained node
	    if (deltaY == 0.0) {

	      // dof corresponding to dX, dZ and theta Y (0,2,4)
	      id(0) = 0; id(1) = 2; id(2) = 4;

	      // set up transformation matrix
	      mat(0,2) = deltaZ;
	      mat(1,2) = -deltaX;

	    } else
	      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, not in xz plane\n",ndC);

	  // rigid diaphragm in yz plane
	  } else {	  

	    // check constrained node in xy plane with retained node
	    if (deltaX == 0.0) {

	      // dof corresponding to dY, dZ and theta X (1,2,3)
	      id(0) = 1; id(1) = 2; id(2) = 3;

	      // set up transformation matrix
	      mat(0,2) = -deltaZ;
	      mat(1,2) = deltaY;

	    } else
	      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, not in xz plane\n",ndC);
	  }
	      
	  // create the MP_Constraint
	  MP_Constraint *newC = new MP_Constraint(startMPtag+i, nR, ndC, 
						  mat, id, id);
	  if (newC == 0) {
	    g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, out of memory\n",ndC);
	  } else {
	    // add the constraint to the domain
	    if (theDomain.addMP_Constraint(newC) == false) {
	      g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, failed to add\n",ndC);
	      delete newC;
	    }
	  }

	} else  // node not in 3d space
	  g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d, not 3d node\n",ndC);
	
      } else // node does not exist
	g3ErrorHandler->warning("RigidDiaphragm::RigidDiaphragm - ignoring constrained Node %d as no node in domain\n",ndC);

    } // for each node in constrained nodes
}



	
	
    
RigidDiaphragm::~RigidDiaphragm()
{
    // does nothing
}

 

