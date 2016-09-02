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
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/RigidBeam.cpp,v $
                                                                        
                                                                        
// File: ~/model/constraints/RigidBeam.C
//
// Written: fmk 12/99
// Revised:
//
// Purpose: This file contains the class implementation for RigidBeam.

#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <Matrix.h>
#include <ID.h>
#include <RigidBeam.h>
#include <RigidRod.h>
#include <elementAPI.h>
#include <string.h>

int OPS_RigidLink(Domain* theDomain)
{
    if (theDomain==0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }

    int num = OPS_GetNumRemainingInputArgs();
    if (num < 3) {
	opserr<<"WARNING: invalid # of args: rigidLink type rNode cNode\n";
	return -1;
    }

    // get type
    const char* type = OPS_GetString();

    // get all data
    num = 2;
    ID data(num);
    if(OPS_GetIntInput(&num,&data(0)) < 0) return -1;

    // construct a rigid rod or beam depending on type
    if (strcmp(type,"-bar")==0 || strcmp(type,"bar")==0) {
	RigidRod theLink(*theDomain,data(0),data(1));
    } else if (strcmp(type,"-beam")==0 || strcmp(type,"beam")==0) {
	RigidBeam theLink(*theDomain,data(0),data(1));
    } else {
	opserr<<"WARNING: unrecognised link type (-bar,-beam)\n";
	return -1;
    }

    return 0;
}

RigidBeam::RigidBeam(Domain &theDomain, int nR, int nC) {

    
    // get a pointer to the retained and constrained nodes - make sure they exist
    Node *nodeR = theDomain.getNode(nR);
    if (nodeR == 0) {
      opserr << "RigidBeam::RigidBeam - retained Node" <<  nR <<  "not in domain\n";
      return;
    }
    Node *nodeC = theDomain.getNode(nC);
    if (nodeR == 0) {
      opserr << "RigidBeam::RigidBeam - constrained Node" <<  nC <<  "not in domain\n";
      return;
    }

    // get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
    const Vector &crdR = nodeR->getCrds();
    const Vector &crdC = nodeC->getCrds();
    int dimR = crdR.Size();
    int dimC = crdC.Size();
    if (dimR != dimC) {
      opserr << "RigidBeam::RigidBeam - mismatch in dimension "  <<
	"between constrained Node " <<  nC <<  " and Retained node" << nR << endln;
      return;
    }
    
    // check the number of dof at each node is the same
    int numDOF = nodeR->getNumberDOF();
    if (numDOF != nodeC->getNumberDOF()){ 
      opserr << "RigidBeam::RigidBeam - mismatch in numDOF "  <<
	"between constrained Node " <<  nC <<  " and Retained node" << nR << endln;
      return;
    }

    // check the number of dof at the nodes >= dimension of problem
    if(numDOF < dimR){    
      opserr << "RigidBeam::RigidBeam - numDOF at nodes " << 
	nR << " and " <<  nC <<  "must be >= dimension of problem\n";
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
	opserr << "RigidBeam::RigidBeam -  for nodes " << 
	  nR << "and " << nC <<  "nodes do not have valid numDOF for their dimension\n";
	return;
      }
	
    }	
    
    // create the MP_Constraint
    MP_Constraint *newC = new MP_Constraint(nR, nC, mat, id, id);
					    
    if (newC == 0) {
      opserr << "RigidBeam::RigidBeam - for nodes " << nC << " and " << nR << ", out of memory\n";
    } else {
      // add the constraint to the domain
      if (theDomain.addMP_Constraint(newC) == false) {
	opserr << "RigidBeam::RigidBeam - for nodes " << nC << " and " << nR << ", could not add to domain\n";
	delete newC;
      }
    }
}
	
    
RigidBeam::~RigidBeam()
{
    // does nothing
}

 

