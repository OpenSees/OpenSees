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
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidDiaphragm.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidDiaphragm.C
//
// Written: fmk 1/99
// Revised:
//
// Purpose: This file contains the class implementation for RigidDiaphragm.

#include <OPS_Globals.h>
#include <stdlib.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <RigidDiaphragm.h>
#include <elementAPI.h>

int OPS_RigidDiaphragm(Domain* theDomain)
{
    if (theDomain==0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    int num = OPS_GetNumRemainingInputArgs();
    if(num < 2) {
	opserr<<"WARNING: invalid # of args: rigidDiaphragm perpDirn rNode cNode1 ...\n";
	return -1;
    }

    // all data
    ID data(num);
    if(OPS_GetIntInput(&num, &data(0)) < 0) return -1;

    // constrained ndoes
    ID cNodes(num-2);
    for(int i=0; i<cNodes.Size(); i++) {
	cNodes(i) = data(i+2);
    }

    RigidDiaphragm theLink(*theDomain,data(1),cNodes,data(0)-1);
    return 0;
}

RigidDiaphragm::RigidDiaphragm(Domain &theDomain, int nR, ID &nC, 
			       int perpPlaneConstrained) {

    // check plane is valid, i.e. perpPlaneConstrained must be 0, 1 or 2
    if (perpPlaneConstrained < 0 || perpPlaneConstrained > 2) {
      opserr << "RigidDiaphragm::RigidDiaphragm - " << 
	"the dirn of perpendicular to constrained plane " << perpPlaneConstrained <<  " not valid\n";
      return;
    }

    // check constrainedNodes ID does not contain the retained node
    if (nC.getLocation(nR) >= 0) {
      opserr << "RigidDiaphragm::RigidDiaphragm - " << 
	"retained node " << nR << " is in constrained node list\n";
      return;
    }	
    
    // get a pointer to the retained node and check node in 3d with 6 dof
    Node *nodeR = theDomain.getNode(nR);
    if (nodeR == 0) {
      opserr << "RigidDiaphragm::RigidDiaphragm - " << 
	"retained Node " <<  nR <<  " not in domain\n";
      return;
    }

    const Vector &crdR = nodeR->getCrds();
    if ((nodeR->getNumberDOF() != 6) || (crdR.Size() != 3)){
      opserr << "RigidDiaphragm::RigidDiaphragm - " << 
	"retained Node " << nR << " not in 3d space with 6 dof\n";
			      
			      
      return;
    }	

    // 
    // create some objects which will be passed to the MP_Constraint 
    // constructor, elements of objects are filled in later
    //
    
    // create the ID to identify the constrained dof 
    ID id(3);

    // construct the transformation matrix Ccr, where  Uc = Ccr Ur & set the diag
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
	      opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << ", not in xy plane\n";

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
	      opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << ", not in xz plane\n";

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
	      opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << 
		", not in xz plane\n";
	  }
	      
	  // create the MP_Constraint
	  MP_Constraint *newC = new MP_Constraint(nR, ndC, mat, id, id);
						  
	  if (newC == 0) {
	    opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << 
	      ", out of memory\n";
	  } else {
	    // add the constraint to the domain
	    if (theDomain.addMP_Constraint(newC) == false) {
	      opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << 
		", failed to add\n";
	      delete newC;
	    }
	  }

	} else  // node not in 3d space
	  opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node  " << ndC << 
	    ", not 3d node\n";
	
      } else // node does not exist
      opserr << "RigidDiaphragm::RigidDiaphragm - ignoring constrained Node " << ndC << 
	" as no node in domain\n";

    } // for each node in constrained nodes
}



	
	
    
RigidDiaphragm::~RigidDiaphragm()
{
    // does nothing
}

 

