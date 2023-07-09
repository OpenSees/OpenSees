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
// $Date: 2010-04-23 22:53:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/MP_Joint3D.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 04/03
// Revision: Arash

// Purpose: This file contains the implementation of class MP_Joint3D.


#include <MP_Joint3D.h>

#include <stdlib.h>
#include <math.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


// constructor for FEM_ObjectBroker
MP_Joint3D::MP_Joint3D()
 :MP_Constraint(CNSTRNT_TAG_MP_Joint3D ), thisDomain(0),
nodeRetained(0), nodeConstrained(0), nodeRotation(0), RotDOF(0),
nodeDisplacement(0), DispDOF(0), LargeDisplacement(0), Length0(0.0),
constraint(0), constrDOF(0), retainDOF(0), RotNormVect(0), DspNormVect(0),
dbTag1(0), dbTag2(0), dbTag3(0), RetainedNode(0), ConstrainedNode(0),
RotationNode(0), DisplacementNode(0)
{
    
}


// general constructor for ModelBuilder
MP_Joint3D::MP_Joint3D( Domain *theDomain, int nodeRetain, int nodeConstr,
		int nodeRot, int Rotdof, int nodeDisp, int Dispdof, int LrgDsp )
:MP_Constraint(CNSTRNT_TAG_MP_Joint3D ), thisDomain(theDomain),
nodeRetained(nodeRetain), nodeConstrained(nodeConstr), nodeRotation(nodeRot),
RotDOF(Rotdof), nodeDisplacement(nodeDisp), DispDOF(Dispdof), 
LargeDisplacement(LrgDsp), Length0(0.0), constraint(0), constrDOF(0),
retainDOF(0), RotNormVect(3), DspNormVect(3), 
dbTag1(0), dbTag2(0), dbTag3(0), RetainedNode(0), ConstrainedNode(0),
RotationNode(0), DisplacementNode(0)
{

  if( thisDomain == NULL ) {
    opserr << "WARNING MP_Joint3D(): Specified domain does not exist";
    opserr << "Domain = 0\n";
    return;
  }

  // get node pointers of constrainted, retained, rotation and displacement nodes
  ConstrainedNode = theDomain->getNode(nodeConstrained);
  if (ConstrainedNode == NULL) {
    
    opserr << "MP_Joint3D::MP_Joint3D: nodeConstrained: ";
    opserr << nodeConstrained << "does not exist in model\n";
    exit(0);
  }
  
  RetainedNode = theDomain->getNode(nodeRetained);
  if (RetainedNode == NULL) {
    
    opserr << "MP_Joint3D::MP_Joint3D: nodeRetained: ";
    opserr << nodeRetained << "does not exist in model\n";
    exit(0);
  }
  
  RotationNode = theDomain->getNode(nodeRotation);
  if (RotationNode == NULL) {
    
    opserr << "MP_Joint3D::MP_Joint3D: nodeRotation: ";
    opserr << nodeRotation << "does not exist in model\n";
    exit(0);
  }
  
  DisplacementNode = theDomain->getNode(nodeDisplacement);
  if (DisplacementNode == NULL)
    {
      opserr << "MP_Joint3D::MP_Joint3D: nodeDisplacement: ";
      opserr << nodeDisplacement << "does not exist in model\n";
      exit(0);
    }
  
  
  // check for proper degrees of freedom
  int RnumDOF = RetainedNode->getNumberDOF();
  int CnumDOF = ConstrainedNode->getNumberDOF();
  if (RnumDOF != 9 || CnumDOF != 6 ){
    opserr << "MP_Joint3D::MP_Joint3D - mismatch in numDOF\n DOF not supported by this type of constraint";
    return;
  }
  
  
  // check the main degree of freedom. Assign auxiliary DOF 
  if ( RotDOF<6 || RotDOF>8 || DispDOF<6 || DispDOF>8 || RotDOF==DispDOF ) {
    opserr << "MP_Joint3D::MP_Joint3D - Wrong degrees of freedom" ;
    return;
  }
  
  // check for proper dimensions of coordinate space
  const Vector &crdRet = RetainedNode->getCrds();
  int dimRet = crdRet.Size();
  const Vector &crdCon = ConstrainedNode->getCrds();
  int dimCon = crdCon.Size();
  const Vector &crdRot = RotationNode->getCrds();
  int dimRot = crdRot.Size();
  const Vector &crdDsp = DisplacementNode->getCrds();
  int dimDsp = crdDsp.Size();
  
  if ( dimRet != 3 || dimCon != 3 || dimRot != 3 || dimDsp != 3 ){
    opserr << "MP_Joint3D::MP_Joint3D - mismatch in dimnesion\n dimension not supported by this type of constraint";
    return;
  }
  
  
  // calculate the initial length of the rigid link
  double deltaX = crdCon(0) - crdRet(0);
  double deltaY = crdCon(1) - crdRet(1);
  double deltaZ = crdCon(2) - crdRet(2);
  
  Length0 = sqrt( deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ );
  if ( Length0 <= 1.0e-12 ) {
    opserr << "MP_Joint3D::MP_Joint3D - The constraint length is zero\n";
  }
  
  // calculate the normal vectors for the rotation mode and displacement mode
  for (int i = 0 ; i<3 ; i++ ) {
    RotNormVect(i)= crdRot(i) - crdRet(i);
    DspNormVect(i)= crdDsp(i) - crdRet(i);	
  }
  
  if ( RotNormVect.Norm() <= 1.0e-12 || DspNormVect.Norm() <= 1.0e-12 ) {		
    opserr << "MP_Joint3D::MP_Joint3D - the normal vector for the rotation mode or the displacement mode is zero\n";
  }
  RotNormVect = RotNormVect / RotNormVect.Norm();
  DspNormVect = DspNormVect / DspNormVect.Norm();
  
  
  // allocate and set up the constranted and retained id's
  
  // the end is released
  constrDOF = new ID(6);
  retainDOF = new ID(8);
  for (int j = 0 ; j<6 ; j++ ) {
    (*constrDOF)(j) = j;
    (*retainDOF)(j) = j;	
  }
  (*retainDOF)(6) = RotDOF;
  (*retainDOF)(7) = DispDOF;
  
  
  if (constrDOF == NULL || retainDOF == NULL ) { 
    opserr << "MP_Joint3D::MP_Joint3D - ran out of memory \ncan not generate ID for nodes\n";
    exit(-1);
  }
  
  
  // allocate and set up the constraint matrix		
  constraint = new Matrix( constrDOF->Size() , retainDOF->Size() );
  
  (*constraint) (0,0) = 1.0 ;
  (*constraint) (1,1) = 1.0 ;
  (*constraint) (2,2) = 1.0 ;
  (*constraint) (1,3) = -deltaZ;
  (*constraint) (2,3) = deltaY;
  (*constraint) (3,3) = 1.0 ;
  (*constraint) (0,4) = deltaZ;
  (*constraint) (2,4) = -deltaX;
  (*constraint) (4,4) = 1.0 ;
  (*constraint) (0,5) = -deltaY;
  (*constraint) (1,5) = deltaX;
  (*constraint) (5,5) = 1.0 ;
  (*constraint) (3,6) = RotNormVect(0);
  (*constraint) (4,6) = RotNormVect(1);
  (*constraint) (5,6) = RotNormVect(2);
  (*constraint) (0,7) = deltaZ*DspNormVect(1) - deltaY*DspNormVect(2);
  (*constraint) (1,7) = deltaX*DspNormVect(2) - deltaZ*DspNormVect(0) ;
  (*constraint) (1,7) = deltaY*DspNormVect(0) - deltaX*DspNormVect(1) ;
  
  
  if (constraint == NULL ) {
    opserr << "MP_Joint3D::MP_Joint3D - ran out of memory \ncan not generate the constraint matrix";
    exit(-1);
  }
}



MP_Joint3D::~MP_Joint3D()
{
  // invoke the destructor on the matrix and the two ID objects
  if (constraint != NULL)
    delete constraint;
  if (constrDOF != NULL)
    delete constrDOF;
  if (retainDOF != NULL)
    delete retainDOF;    
}


int
MP_Joint3D::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Joint3D::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Joint3D::getConstrainedDOFs(void) const
{
  if (constrDOF == NULL) {
    opserr << "MP_Joint3D::getConstrainedDOF - no ID was set, ";
    opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
    exit(-1);
  }
  
  // return the ID corresponding to constrained DOF of Ccr
  return (*constrDOF);    
}


const ID &
MP_Joint3D::getRetainedDOFs(void) const
{
  if (retainDOF == NULL) {
    opserr << "MP_Joint3D::getRetainedDOFs - no ID was set\n ";
    opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
    exit(-1);
  }
  
  // return the ID corresponding to retained DOF of Ccr
  return (*retainDOF);    
}


int 
MP_Joint3D::applyConstraint(double timeStamp)
{
  if ( LargeDisplacement != 0 )
    {
      // calculate the constraint at this moment
      
      // get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
      const Vector &crdRet = RetainedNode->getCrds();
      const Vector &crdCon = ConstrainedNode->getCrds();
      const Vector &crdRot = RotationNode->getCrds();
      const Vector &crdDsp = DisplacementNode->getCrds();
      
      // get committed displacements of nodes to get updated coordinates
      const Vector &dispRet = RetainedNode->getDisp();
      const Vector &dispCon = ConstrainedNode->getDisp();
      const Vector &dispRot = RotationNode->getDisp();
      const Vector &dispDsp = DisplacementNode->getDisp();
      
      double deltaX = dispCon(0) + crdCon(0) - dispRet(0) - crdRet(0);
      double deltaY = dispCon(1) + crdCon(1) - dispRet(1) - crdRet(1);
      double deltaZ = dispCon(2) + crdCon(2) - dispRet(2) - crdRet(2);
      
      for ( int i = 0 ; i<3 ; i++ )
	{
	  RotNormVect(i)= dispRot(i) + crdRot(i) - dispRet(i) - crdRet(i);
	  DspNormVect(i)= dispDsp(i) + crdDsp(i) - dispRet(i) - crdRet(i);	
	}
      
      RotNormVect = RotNormVect / RotNormVect.Norm();
      DspNormVect = DspNormVect / DspNormVect.Norm();
      
      
      constraint->Zero();
      
      (*constraint) (0,0) = 1.0 ;
      (*constraint) (1,1) = 1.0 ;
      (*constraint) (2,2) = 1.0 ;
      (*constraint) (1,3) = -deltaZ;
      (*constraint) (2,3) = deltaY;
      (*constraint) (3,3) = 1.0 ;
      (*constraint) (0,4) = deltaZ;
      (*constraint) (2,4) = -deltaX;
      (*constraint) (4,4) = 1.0 ;
      (*constraint) (0,5) = -deltaY;
      (*constraint) (1,5) = deltaX;
      (*constraint) (5,5) = 1.0 ;
      (*constraint) (3,6) = RotNormVect(0);
      (*constraint) (4,6) = RotNormVect(1);
      (*constraint) (5,6) = RotNormVect(2);
      (*constraint) (0,7) = deltaZ*DspNormVect(1) - deltaY*DspNormVect(2);
      (*constraint) (1,7) = deltaX*DspNormVect(2) - deltaZ*DspNormVect(0) ;
      (*constraint) (2,7) = deltaY*DspNormVect(0) - deltaX*DspNormVect(1) ;
    }
  return 0;
}


bool
MP_Joint3D::isTimeVarying(void) const
{
  if ( LargeDisplacement != 0 ) return true;
  
  return false;
}


int MP_Joint3D::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}


int MP_Joint3D::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return 0;
}


const Matrix &MP_Joint3D::getConstraint(void)
{
    if (constraint == 0) {
      opserr << "MP_Joint3D::getConstraint - no Matrix was set\n";
      exit(-1);
    }    
    
    // Length correction
    // to correct the trial displacement
    if ( LargeDisplacement == 2 )
      {
	// get the coordinates of the two nodes - check dimensions are the same FOR THE MOMENT
	const Vector &crdR = RetainedNode->getCrds();
	const Vector &crdC = ConstrainedNode->getCrds();
	
	// get committed displacements of nodes to get updated coordinates
	const Vector &dispR = RetainedNode->getTrialDisp();
	const Vector &dispC = ConstrainedNode->getTrialDisp();
	
	double deltaX = dispC(0) + crdC(0) - dispR(0) - crdR(0);
	double deltaY = dispC(1) + crdC(1) - dispR(1) - crdR(1);
	double deltaZ = dispC(2) + crdC(2) - dispR(2) - crdR(2);
	
	
	Vector Direction(3);
	Direction(0) = deltaX;
	Direction(1) = deltaY;
	Direction(2) = deltaZ;
	double NewLength = Direction.Norm();
	if ( NewLength < 1e-12 ) opserr << "MP_Joint3D::applyConstraint : length of rigid link is too small or zero"; 
	Direction = Direction * (Length0/NewLength);		// correct the length
	// find new displacements of the constrainted node
	
	Vector NewLocation(6);
	NewLocation(0) = Direction(0) + dispR(0) + crdR(0) - crdC(0);
	NewLocation(1) = Direction(1) + dispR(1) + crdR(1) - crdC(1);
	NewLocation(2) = Direction(2) + dispR(2) + crdR(2) - crdC(2);
	NewLocation(3) = dispC(3);
	NewLocation(4) = dispC(4);
	NewLocation(5) = dispC(5);
	
	int dummy = ConstrainedNode->setTrialDisp( NewLocation );
      }
    // end of length correction procedure
    
    // return the constraint matrix Ccr
    return (*constraint);
}

void MP_Joint3D::Print(OPS_Stream &s, int flag )
{
  s << "MP_Joint3D: " << this->getTag() << "\n";
  s << "\tConstrained Node: " << nodeConstrained;
  s << " Retained Node: " << nodeRetained ;
  s << " Large Disp: " << LargeDisplacement;
  if (constrDOF != 0)
    s << " constrained dof: " << *constrDOF;    
  if (retainDOF != 0)
    s << " retained dof: " << *retainDOF;        
  if (constraint != 0)
    s << " constraint matrix: " << *constraint << "\n";
  
}


void
MP_Joint3D::setDomain(Domain *theDomain)
{
  this->DomainComponent::setDomain(theDomain);
  thisDomain = theDomain;
  
  RetainedNode = thisDomain->getNode(nodeRetained);
  ConstrainedNode = thisDomain->getNode(nodeConstrained);
  RotationNode = thisDomain->getNode(nodeRotation);
  DisplacementNode = thisDomain->getNode(nodeDisplacement);
}
