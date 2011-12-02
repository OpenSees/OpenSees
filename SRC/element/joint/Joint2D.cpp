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

// $Revision: 1.9 $
// $Date: 2006-08-04 22:22:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint2D.cpp,v $

// Written: Arash & GGD
// Created: 03/02
// Revision: Arash

// Joint2D.cpp: implementation of the Joint2D class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MP_Constraint.h>
#include <MP_Joint2D.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>
#include <Joint2D.h>
#include <DamageModel.h>


Matrix Joint2D::K(16,16);
Vector Joint2D::V(16);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Joint2D::Joint2D()
  :Element(0, ELE_TAG_Joint2D ), 
  ExternalNodes(5), InternalConstraints(4), 
  TheDomain(0), numDof(0), nodeDbTag(0), dofDbTag(0)
{
	for ( int i=0 ; i<5 ; i++ )
	{
		theSprings[i] = NULL;;
		fixedEnd[i] = 1;
		theNodes[i] = NULL;
	}
}


Joint2D::Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
			     UniaxialMaterial &spring1,	UniaxialMaterial &spring2,
			     UniaxialMaterial &spring3, UniaxialMaterial &spring4,
			     UniaxialMaterial &springC, Domain *theDomain, int LrgDisp)
  :Element(tag, ELE_TAG_Joint2D ), 
  ExternalNodes(5), InternalConstraints(4), 
  TheDomain(0), numDof(0), nodeDbTag(0), dofDbTag(0),theLoadSens(0)
{
  int i;
  numDof  = 16;

  K.Zero();
  V.Zero();

  TheDomain = theDomain;
  if( TheDomain==NULL ) {
    opserr << "WARNING Joint2D(): Specified domain does not exist , Domain = 0\n";
    return;
  }

  // Save external node id's
  ExternalNodes(0) = nd1;
  ExternalNodes(1) = nd2;
  ExternalNodes(2) = nd3;
  ExternalNodes(3) = nd4;
  ExternalNodes(4) = IntNodeTag;
  

  // get  the external nodes
  for ( i=0 ; i<4 ; i++)
  {
    theNodes[i] = NULL;
    theNodes[i] = TheDomain->getNode( ExternalNodes(i) );
    if (theNodes[i] == NULL) {
      opserr << "WARNING Joint2D::setDomain(): Nd" <<(i+1) <<": ";
      opserr << ExternalNodes(i) << "does not exist in model for element \n" << *this;
      return;
    }
  }
  
  // check for a two dimensional domain, since this element supports only two dimensions 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  const Vector &end3Crd = theNodes[2]->getCrds();
  const Vector &end4Crd = theNodes[3]->getCrds();
  
  int dimNd1 = end1Crd.Size();
  int dimNd2 = end2Crd.Size();
  int dimNd3 = end3Crd.Size();
  int dimNd4 = end4Crd.Size();

  if (dimNd1 != 2 || dimNd2 != 2 || dimNd3 != 2 || dimNd4 != 2 ) {
    opserr << "WARNING Joint2D::setDomain(): has incorrect space dimension \n";
    opserr << "                                    space dimension not supported by Joint2D";
    return;
  }
	
  // now verify the number of dof at node ends
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  int dofNd3 = theNodes[2]->getNumberDOF();
  int dofNd4 = theNodes[3]->getNumberDOF();

  if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3 ) {
    opserr << "WARNING Joint2D::Joint2D: has incorrect degrees of freedom \n";
    opserr << "                                    DOF not supported by Joint2D";
    return;
  }
  
  // check the joint size. The joint size must be non-zero
  Vector Center1(end1Crd);
  Vector Center2(end2Crd);
  Center1 = Center1 - end3Crd;
  Center2 = Center2 - end4Crd;
	
  double L1 = Center1.Norm();
  double L2 = Center2.Norm();
  
  if( Center1.Norm()<1e-12  || Center2.Norm()<1e-12 ) {
    opserr << "WARNING Joint2D::(): zero length\n";
    return;	
  }
	
  // check if nodes are not located on each other and they can construct
  // a parallelogram
  Center1 = end1Crd + end3Crd;
  Center2 = end2Crd + end4Crd;
  
  Center1 = 0.5 * Center1;
  Center2 = 0.5 * Center2;
  
  Vector Center3(Center2);
  Center3 = Center3 - Center1;

  if ( Center3.Norm() > 1e-6 ) {
    opserr << "WARNING Joint2D::(): can not construct a paralelogram over external nodes\n";
    return;	
  }
	
  // Generate internal node and add it up to domain
  theNodes[4]  = new Node ( IntNodeTag , 4, Center1(0) , Center1(1) );
  if ( theNodes[4] == NULL ) {
    opserr << "Joint2D::Joint2D - Unable to generate new nodes , out of memory\n" ;
  } else {
    if( TheDomain->addNode( theNodes[4] ) == false )		// add intenal nodes to domain
      opserr << "Joint2D::Joint2D - unable to add internal nodeto domain\n";
  }
  
  // make copy of the uniaxial materials for the element
  
  if ( &spring1 == NULL ) { fixedEnd[0] = 1;  theSprings[0] = NULL; } else { fixedEnd[0] = 0; theSprings[0] = spring1.getCopy(); }
  if ( &spring2 == NULL ) { fixedEnd[1] = 1;  theSprings[1] = NULL; } else { fixedEnd[1] = 0; theSprings[1] = spring2.getCopy(); }  
  if ( &spring3 == NULL ) { fixedEnd[2] = 1;  theSprings[2] = NULL; } else { fixedEnd[2] = 0; theSprings[2] = spring3.getCopy(); }
  if ( &spring4 == NULL ) { fixedEnd[3] = 1;  theSprings[3] = NULL; } else { fixedEnd[3] = 0; theSprings[3] = spring4.getCopy(); }
  if ( &springC == NULL ) { opserr << "ERROR Joint2D::Joint2D(): The central node does not exist "; exit(-1); } else { fixedEnd[4] = 0; theSprings[4] = springC.getCopy(); }
  
  
  for ( i=0 ; i<5 ; i++ )
    {
      if ( fixedEnd[i] == 0  && theSprings[i] == NULL ) {
	opserr << "ERROR Joint2D::Joint2D(): Can not make copy of uniaxial materials, out of memory ";
	exit(-1);
      }
    }
  
  // Generate and add constraints to domain
  
  // get the constraint numbers
  int startMPtag = theDomain->getNumMPs();
  for ( i=0 ; i<4 ; i++ ) InternalConstraints(i) = startMPtag + i ;
  
  // create MP_Joint constraint node 1
  if ( addMP_Joint( TheDomain, InternalConstraints(0), ExternalNodes(4), ExternalNodes(0), 2, fixedEnd[0], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 1\n";
    return;
  }
  
  // create MP_Joint constraint node 2
  if ( addMP_Joint( TheDomain, InternalConstraints(1), ExternalNodes(4), ExternalNodes(1), 3, fixedEnd[1], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 2\n";
		return;
  }
  
  // create MP_Joint constraint node 3
  if ( addMP_Joint( TheDomain, InternalConstraints(2), ExternalNodes(4), ExternalNodes(2), 2, fixedEnd[2], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 3\n";
    return;
  }
  
  // create MP_Joint constraint node 4
  if ( addMP_Joint( TheDomain, InternalConstraints(3), ExternalNodes(4), ExternalNodes(3), 3, fixedEnd[3], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 4\n";
    return;
  }

  // Zero the damage models
  for ( i = 0 ; i < 5 ; i++ ) theDamages[i] = NULL;
}


Joint2D::Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
			     UniaxialMaterial &spring1,	UniaxialMaterial &spring2,
			     UniaxialMaterial &spring3, UniaxialMaterial &spring4,
			     UniaxialMaterial &springC, Domain *theDomain, int LrgDisp,
				 DamageModel &dmg1, DamageModel &dmg2, DamageModel &dmg3,
				 DamageModel &dmg4, DamageModel &dmgC)
  :Element(tag, ELE_TAG_Joint2D ), 
  ExternalNodes(5), InternalConstraints(4), 
  TheDomain(0), numDof(0), nodeDbTag(0), dofDbTag(0),theLoadSens(0)
{
  int i;
  numDof  = 16;

  K.Zero();
  V.Zero();

  TheDomain = theDomain;
  if( TheDomain==NULL ) {
    opserr << "WARNING Joint2D(): Specified domain does not exist , Domain = 0\n";
    return;
  }

  // Save external node id's
  ExternalNodes(0) = nd1;
  ExternalNodes(1) = nd2;
  ExternalNodes(2) = nd3;
  ExternalNodes(3) = nd4;
  ExternalNodes(4) = IntNodeTag;
  

  // get  the external nodes
  for ( i=0 ; i<4 ; i++)
  {
    theNodes[i] = NULL;
    theNodes[i] = TheDomain->getNode( ExternalNodes(i) );
    if (theNodes[i] == NULL) {
      opserr << "WARNING Joint2D::setDomain(): Nd" <<(i+1) <<": ";
      opserr << ExternalNodes(i) << "does not exist in model for element \n" << *this;
      return;
    }
  }
  
  // check for a two dimensional domain, since this element supports only two dimensions 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  const Vector &end3Crd = theNodes[2]->getCrds();
  const Vector &end4Crd = theNodes[3]->getCrds();
  
  int dimNd1 = end1Crd.Size();
  int dimNd2 = end2Crd.Size();
  int dimNd3 = end3Crd.Size();
  int dimNd4 = end4Crd.Size();

  if (dimNd1 != 2 || dimNd2 != 2 || dimNd3 != 2 || dimNd4 != 2 ) {
    opserr << "WARNING Joint2D::setDomain(): has incorrect space dimension \n";
    opserr << "                                    space dimension not supported by Joint2D";
    return;
  }
	
  // now verify the number of dof at node ends
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  int dofNd3 = theNodes[2]->getNumberDOF();
  int dofNd4 = theNodes[3]->getNumberDOF();

  if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3 ) {
    opserr << "WARNING Joint2D::Joint2D: has incorrect degrees of freedom \n";
    opserr << "                                    DOF not supported by Joint2D";
    return;
  }
  
  // check the joint size. The joint size must be non-zero
  Vector Center1(end1Crd);
  Vector Center2(end2Crd);
  Center1 = Center1 - end3Crd;
  Center2 = Center2 - end4Crd;
	
  double L1 = Center1.Norm();
  double L2 = Center2.Norm();
  
  if( Center1.Norm()<1e-12  || Center2.Norm()<1e-12 ) {
    opserr << "WARNING Joint2D::(): zero length\n";
    return;	
  }
	
  // check if nodes are not located on each other and they can construct
  // a parallelogram
  Center1 = end1Crd + end3Crd;
  Center2 = end2Crd + end4Crd;
  
  Center1 = 0.5 * Center1;
  Center2 = 0.5 * Center2;
  
  Vector Center3(Center2);
  Center3 = Center3 - Center1;

  if ( Center3.Norm() > 1e-6 ) {
    opserr << "WARNING Joint2D::(): can not construct a paralelogram over external nodes\n";
    return;	
  }
	
  // Generate internal node and add it up to domain
  theNodes[4]  = new Node ( IntNodeTag , 4, Center1(0) , Center1(1) );
  if ( theNodes[4] == NULL ) {
    opserr << "Joint2D::Joint2D - Unable to generate new nodes , out of memory\n" ;
  } else {
    if( TheDomain->addNode( theNodes[4] ) == false )		// add intenal nodes to domain
      opserr << "Joint2D::Joint2D - unable to add internal nodeto domain\n";
  }
  
  // make copy of the uniaxial materials for the element
  
  if ( &spring1 == NULL )	{ fixedEnd[0] = 1;  theSprings[0] = NULL; }	else { fixedEnd[0] = 0; theSprings[0] = spring1.getCopy(); }
  if ( &spring2 == NULL ) { fixedEnd[1] = 1;  theSprings[1] = NULL; } else { fixedEnd[1] = 0; theSprings[1] = spring2.getCopy(); }  
  if ( &spring3 == NULL ) { fixedEnd[2] = 1;  theSprings[2] = NULL; } else { fixedEnd[2] = 0; theSprings[2] = spring3.getCopy(); }
  if ( &spring4 == NULL ) { fixedEnd[3] = 1;  theSprings[3] = NULL; } else { fixedEnd[3] = 0; theSprings[3] = spring4.getCopy(); }
  if ( &springC == NULL ) { opserr << "ERROR Joint2D::Joint2D(): The central node does not exist "; exit(-1); } else { fixedEnd[4] = 0; theSprings[4] = springC.getCopy(); }
  
  
  for ( i=0 ; i<5 ; i++ )
    {
      if ( fixedEnd[i] == 0  && theSprings[i] == NULL ) {
	opserr << "ERROR Joint2D::Joint2D(): Can not make copy of uniaxial materials, out of memory ";
	exit(-1);
      }
    }
  
  // Generate and add constraints to domain
  
  // get the constraint numbers
  int startMPtag = theDomain->getNumMPs();
  for ( i=0 ; i<4 ; i++ ) InternalConstraints(i) = startMPtag + i ;
  
  // create MP_Joint constraint node 1
  if ( addMP_Joint( TheDomain, InternalConstraints(0), ExternalNodes(4), ExternalNodes(0), 2, fixedEnd[0], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 1\n";
    return;
  }
  
  // create MP_Joint constraint node 2
  if ( addMP_Joint( TheDomain, InternalConstraints(1), ExternalNodes(4), ExternalNodes(1), 3, fixedEnd[1], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 2\n";
		return;
  }
  
  // create MP_Joint constraint node 3
  if ( addMP_Joint( TheDomain, InternalConstraints(2), ExternalNodes(4), ExternalNodes(2), 2, fixedEnd[2], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 3\n";
    return;
  }
  
  // create MP_Joint constraint node 4
  if ( addMP_Joint( TheDomain, InternalConstraints(3), ExternalNodes(4), ExternalNodes(3), 3, fixedEnd[3], LrgDisp ) != 0) {
    opserr << "WARNING Joint2D::Joint2D(): can not generate ForJoint MP at node 4\n";
    return;
  }
  	// Handle the damage models
	if ( &dmg1 == NULL ) { theDamages[0] = NULL; } else { theDamages[0] = dmg1.getCopy(); }
	if ( &dmg2 == NULL ) { theDamages[1] = NULL; } else { theDamages[1] = dmg2.getCopy(); }
	if ( &dmg3 == NULL ) { theDamages[2] = NULL; } else { theDamages[2] = dmg3.getCopy(); }
	if ( &dmg4 == NULL ) { theDamages[3] = NULL; } else { theDamages[3] = dmg4.getCopy(); }
	if ( &dmgC == NULL ) { theDamages[4] = NULL; } else { theDamages[4] = dmgC.getCopy(); }

	for ( i = 0 ; i < 5 ; i ++ ) if ( theDamages[i] != NULL ) theDamages[i]->revertToStart();
				
}


Joint2D::~Joint2D()
{

	if ( TheDomain != NULL)
	{
		MP_Constraint *Temp_MP;
		for ( int i=0 ; i < 4 ; i++ )
		{
			Temp_MP = TheDomain->getMP_Constraint( InternalConstraints(i) );
			
			if ( Temp_MP != NULL )
			{
				TheDomain->removeMP_Constraint( InternalConstraints(i) );
				delete Temp_MP;
			}
		}
		if ( theNodes[4] != NULL )
		{
			int intnodetag = theNodes[4]->getTag();
			TheDomain->removeNode( intnodetag );
			delete theNodes[4];
		}
	}

	for (int i=0 ; i<5 ; i++) {
		if ( theSprings[i] != NULL ) delete theSprings[i];
		if ( theDamages[i] != NULL ) delete theDamages[i];
	}
}


void Joint2D::setDomain(Domain *theDomain)
{
	//Ckeck domain not null - invoked when object removed from a domain
	if (theDomain == 0) {
		for(int i=0 ; i<4 ; i++) theNodes[i] = NULL;
	} else {
		
		TheDomain = theDomain;
		this->DomainComponent::setDomain(theDomain);
		
		for (int i=0 ; i<5 ; i++)
			if ( theNodes[i] ==0 )  theNodes[i] = TheDomain->getNode( ExternalNodes(i) );
	}
	
}//setDomain


int Joint2D::addMP_Joint(Domain *theDomain, int mpNum, 
				  int RnodeID, int CnodeID, 
				  int MainDOF, int FixedEnd, int LrgDispFlag )
{
	MP_Constraint *Temp_MP;

	// create MP_ForJoint constraint
	Temp_MP = new MP_Joint2D( theDomain, mpNum, RnodeID, CnodeID, MainDOF, FixedEnd, LrgDispFlag );
  
	if (Temp_MP == NULL)
	{
		opserr << "Joint2D::addMP_Joint - WARNING ran out of memory for ForJoint MP_Constraint ";
		return -1;
	}
	// Add the multi-point constraint to the domain
	if (theDomain->addMP_Constraint (Temp_MP) == false)
	{
		opserr << "Joint2D::addMP_Joint - WARNING could not add equalDOF MP_Constraint to domain ";
		delete Temp_MP;
		return -2;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////
// Public methods called, taken care of for 2D element subclasses
//////////////////////////////////////////////////////////////////////

int Joint2D::update(void)
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	const Vector &dispC = theNodes[4]->getTrialDisp();
	double Delta[5];
	Delta[0] = disp1(2) - dispC(3);
	Delta[1] = disp2(2) - dispC(2);
	Delta[2] = disp3(2) - dispC(3);
	Delta[3] = disp4(2) - dispC(2);
	Delta[4] = dispC(3) - dispC(2);
	int result = 0;

	for ( int i=0 ; i<5 ; i++ )
	{
		if ( theSprings[i] != NULL ) result = theSprings[i]->setTrialStrain(Delta[i]);
		if ( result != 0 ) break;
	}

	return result;
}

int Joint2D::commitState()
{
	int result = 0;


	// setting the trial state for the damage models

	Vector InforForDamage(3);


	for ( int i=0 ; i<5 ; i++ )
	{
		if ( theSprings[i] != NULL ) result = theSprings[i]->commitState();
		if ( result != 0 ) break;

		if ( theSprings[i] != NULL && theDamages[i] != NULL ) {
			InforForDamage(0) = theSprings[i]->getStrain();
			InforForDamage(1) = theSprings[i]->getStress();
			InforForDamage(2) = theSprings[i]->getInitialTangent();
				
			theDamages[i]->setTrial(InforForDamage);
			result = theDamages[i]->commitState();
			if ( result != 0 ) break;
		}

	}
	
	return result;
}

int Joint2D::revertToLastCommit()
{
	int result = 0;

	for ( int i=0 ; i<5 ; i++ )
	{
		if ( theSprings[i] != NULL ) result = theSprings[i]->revertToLastCommit();
		if ( result != 0 ) break;
		if ( theDamages[i] != NULL ) result = theDamages[i]->revertToLastCommit();
		if ( result != 0 ) break;
	}
	
	return result;
}

int Joint2D::revertToStart(void)
{
	int result = 0;

	for ( int i=0 ; i<5 ; i++ )
	{
		if ( theSprings[i] != NULL ) result = theSprings[i]->revertToStart();
		if ( result != 0 ) break;
		if ( theDamages[i] != NULL ) result = theDamages[i]->revertToStart();
		if ( result != 0 ) break;
	}
	
	return result;
}


int Joint2D::getNumExternalNodes(void) const
{
	return 5;
}

const ID &Joint2D::getExternalNodes(void)
{
	return ExternalNodes;
}

Node **Joint2D::getNodePtrs(void)
{
	return theNodes;
}

int Joint2D::getNumDOF(void)
{
  return numDof;
}

const Matrix &Joint2D::getTangentStiff(void)
{
	double Ktangent[5] ;
	for ( int i=0 ; i<5 ; i++ ) 
	{
		Ktangent[i] = 0;
		if ( theSprings[i] != NULL ) Ktangent[i] = theSprings[i]->getTangent();
	}

	K.Zero();

	K(2,2)  =  Ktangent[0];
	K(2,15) = -Ktangent[0];
	K(5,5)  =  Ktangent[1];
	K(5,14) = -Ktangent[1];
	K(8,8)  =  Ktangent[2];
	K(8,15) = -Ktangent[2];
	K(11,11)=  Ktangent[3];
	K(11,14)= -Ktangent[3];
	K(14,5) = -Ktangent[1];
	K(14,11)= -Ktangent[3];
	K(14,14)=  Ktangent[1] + Ktangent[3] + Ktangent[4];
	K(14,15)= -Ktangent[4];
	K(15,2) = -Ktangent[0];
	K(15,8) = -Ktangent[2];
	K(15,14)= -Ktangent[4];
	K(15,15)=  Ktangent[0] + Ktangent[2] + Ktangent[4];

	return K;
}


const Matrix &Joint2D::getInitialStiff(void)
{
	double Kintial[5] ;
	for ( int i=0 ; i<5 ; i++ ) 
	{
		Kintial[i] = 0;
		if ( theSprings[i] != NULL ) Kintial[i] = theSprings[i]->getTangent();
	}

	K.Zero();

	K(2,2)  =  Kintial[0];
	K(2,15) = -Kintial[0];
	K(5,5)  =  Kintial[1];
	K(5,14) = -Kintial[1];
	K(8,8)  =  Kintial[2];
	K(8,15) = -Kintial[2];
	K(11,11)=  Kintial[3];
	K(11,14)= -Kintial[3];
	K(14,5) = -Kintial[1];
	K(14,11)= -Kintial[3];
	K(14,14)=  Kintial[1] + Kintial[3] + Kintial[4];
	K(14,15)= -Kintial[4];
	K(15,2) = -Kintial[0];
	K(15,8) = -Kintial[2];
	K(15,14)= -Kintial[4];
	K(15,15)=  Kintial[0] + Kintial[2] + Kintial[4];

	return K;
}


const Matrix &Joint2D::getDamp(void)
{	
	K.Zero();
	return K;
}

const Matrix &Joint2D::getMass(void)
{
	K.Zero();
	return K;
}

void Joint2D::Print(OPS_Stream &s, int flag )
{
  s << "\nElement: " << getTag() << " type: Joint2D iNode: "
    << ExternalNodes(0) << " jNode: " << ExternalNodes(1) << "\n"
    << " kNode: " << ExternalNodes(2) << " lNode: " << ExternalNodes(3) << "\n"
	<< " Internal node: " << ExternalNodes(4) << "\n";
}

/////////////////////////////////////////////////////////////////////
// methods for applying and returning loads
//////////////////////////////////////////////////////////////////////

void Joint2D::zeroLoad(void)
{

}

int Joint2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}

int Joint2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}



const Vector &Joint2D::getResistingForce()
{
	double Force[5] ;
	for ( int i=0 ; i<5 ; i++ ) 
	{
		Force[i] = 0;
		if ( theSprings[i] != NULL ) Force[i] = theSprings[i]->getStress();
	}

	V.Zero();

	V(2) = Force[0];
	V(5) = Force[1];
	V(8) = Force[2];
	V(11)= Force[3];
	V(14)= -Force[4] - Force[1] - Force[3];
	V(15)= Force[4] - Force[0] - Force[2];

	return V;
}

const Vector &
Joint2D::getResistingForceIncInertia()
{
	return this->getResistingForce();
}



int Joint2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// first determine the four corner points of the element based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2

	const Vector &node1Crd = theNodes[0]->getCrds();
	const Vector &node2Crd = theNodes[1]->getCrds();	
	const Vector &node3Crd = theNodes[2]->getCrds();
	const Vector &node4Crd = theNodes[3]->getCrds();

	const Vector &node1Disp = theNodes[0]->getDisp();
	const Vector &node2Disp = theNodes[1]->getDisp();    
	const Vector &node3Disp = theNodes[2]->getDisp();
	const Vector &node4Disp = theNodes[3]->getDisp();  

	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	
	// calculate the current coordinates of four external nodes
	for (int i=0; i<2; i++) 
    {
		v1(i) = node1Crd(i)+node1Disp(i)*fact;
		v2(i) = node2Crd(i)+node2Disp(i)*fact;
		v3(i) = node3Crd(i)+node3Disp(i)*fact;
		v4(i) = node4Crd(i)+node4Disp(i)*fact;
	}

	// draw the center lines
	int dummy;
	dummy = theViewer.drawLine(v1, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v4, 1.0, 1.0);
	
	// calculate four corners of the element
	Vector vb(3);
	Vector vc(3);

	vb = v1 - v3;
	vc = v2 - v4;

	v1 = v3 - 0.5 * vc;
	v2 = v1 + vb;
	v3 = v2 + vc;
	v4 = v1 + vc;
	
	dummy = theViewer.drawLine(v1, v2, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v3, v4, 1.0, 1.0);
	dummy = theViewer.drawLine(v4, v1, 1.0, 1.0);

	return 0;

}


//most-probably requires to be overridden
Response* Joint2D::setResponse(const char **argv, int argc, Information &eleInformation, OPS_Stream &output)
{
//
// we compare argv[0] for known response types for the Truss
//
	if (strcmp(argv[0],"node") == 0 || strcmp(argv[0],"internalNode") == 0 )
    return new ElementResponse(this, 1, Vector(4));

	else if (strcmp(argv[0],"size") == 0 || strcmp(argv[0],"jointSize") == 0 )
    return new ElementResponse(this, 2, Vector(2));

	else if (strcmp(argv[0],"moment") == 0 || strcmp(argv[0],"-moment") == 0 
		|| strcmp(argv[0],"force") == 0 || strcmp(argv[0],"-force") == 0 )
    return new ElementResponse(this, 3, Vector(5));

	else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	   strcmp(argv[0],"deformation") == 0 )
    return new ElementResponse(this, 4, Vector(5));

	else if (strcmp(argv[0],"defoANDforce") == 0 || strcmp(argv[0],"deformationANDforce") == 0 ||
	   strcmp(argv[0],"deformationsANDforces") == 0 )
    return new ElementResponse(this, 5, Vector(10));

	else if ( strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0 )
    return new ElementResponse(this, 6, Matrix(16,16) );
	
	else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0)
		return new ElementResponse(this, 7, Vector(5));

	else if ( strcmp(argv[0],"damage") == 0 || strcmp(argv[0],"damages") == 0 ||
		strcmp(argv[0],"-damage") == 0 || strcmp(argv[0],"-damages") == 0)
		return new ElementResponse(this, 8, Vector(5));
	// material response
	else if ( (strcmp(argv[0],"spring")==0) || (strcmp(argv[0],"-spring") == 0) ||
		(strcmp(argv[0],"material")==0) || (strcmp(argv[0],"-material") == 0) ) {
		int materialNum = atoi(argv[1]) - 1;
		
		if (materialNum >= 0 && materialNum < 5)
			if (theSprings[materialNum] != 0)
			  return theSprings[materialNum]->setResponse(&argv[2], argc-2, eleInformation, output);
	}

	else 
		return 0;
  	
}

int Joint2D::getResponse(int responseID, Information &eleInformation)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		if(eleInformation.theVector!=0)
		{
			const Vector& disp = theNodes[4]->getTrialDisp();
			(*(eleInformation.theVector))(0) = disp(0);
			(*(eleInformation.theVector))(1) = disp(1);
			(*(eleInformation.theVector))(2) = disp(2);
			(*(eleInformation.theVector))(3) = disp(3);
		}
		return 0;

	case 2:
		if(eleInformation.theVector!=0)
		{
			const Vector &node1Crd = theNodes[0]->getCrds();
			const Vector &node2Crd = theNodes[1]->getCrds();	
			const Vector &node3Crd = theNodes[2]->getCrds();
			const Vector &node4Crd = theNodes[3]->getCrds();

			const Vector &node1Disp = theNodes[0]->getDisp();
			const Vector &node2Disp = theNodes[1]->getDisp();    
			const Vector &node3Disp = theNodes[2]->getDisp();
			const Vector &node4Disp = theNodes[3]->getDisp();  

			Vector v1(2);
			Vector v2(2);
			Vector v3(2);
			Vector v4(2);
	
			// calculate the current coordinates of four external nodes
			for (int i=0; i<2; i++) 
		    {
				v1(i) = node1Crd(i)+node1Disp(i);
				v2(i) = node2Crd(i)+node2Disp(i);
				v3(i) = node3Crd(i)+node3Disp(i);
				v4(i) = node4Crd(i)+node4Disp(i);
			}
			
			v3 = v3 - v1;
			v4 = v4 - v2;

			v1(0) = v3.Norm();
			v1(1) = v4.Norm();

			*(eleInformation.theVector) = v1;
		}
		return 0;

	case 3:
		if( eleInformation.theVector != 0 )
		{
			for ( int i =0 ; i<5 ; i++ )
			{
				(*(eleInformation.theVector))(i) = 0.0;
				if ( theSprings[i] != NULL ) 
					(*(eleInformation.theVector))(i) = theSprings[i]->getStress();
			}
		}
		return 0;

	case 4:
		if(eleInformation.theVector!=0)
		{
			for ( int i =0 ; i<5 ; i++ )
			{
				(*(eleInformation.theVector))(i) = 0.0;
				if ( theSprings[i] != NULL ) 
					(*(eleInformation.theVector))(i) = theSprings[i]->getStrain();
			}
		}
		return 0;

	case 5:
		if(eleInformation.theVector!=0)
		{
			for ( int i =0 ; i<5 ; i++ )
			{
				(*(eleInformation.theVector))(i) = 0.0;
				(*(eleInformation.theVector))(i+5) = 0.0;
				if ( theSprings[i] != NULL )
				{
					(*(eleInformation.theVector))(i) = theSprings[i]->getStrain();
					(*(eleInformation.theVector))(i+5) = theSprings[i]->getStress();
				}
			}
		}
		return 0;

	case 6:
		return eleInformation.setMatrix(this->getTangentStiff());
	
	case 7:
		if(eleInformation.theVector!=0)
		{
			for ( int i=0 ; i<5 ; i++ )
			{
				(*(eleInformation.theVector))(i) = 0.0;
				if ( theSprings[i] != NULL && theSprings[i]->getInitialTangent() != 0.0 )
				{
					(*(eleInformation.theVector))(i) = 
						theSprings[i]->getStrain() - theSprings[i]->getStress()/theSprings[i]->getInitialTangent();
				}
				
			}			
		}
		return 0;

	case 8:
		if(eleInformation.theVector!=0)
		{
			for ( int i=0 ; i<5 ; i++ )
			{
				(*(eleInformation.theVector))(i) = 0.0;
				if ( theDamages[i] != NULL ) {
					(*(eleInformation.theVector))(i) = theDamages[i]->getDamage();
				}
			}
		}
		return 0;
	
	default:
		return -1;
	}
}


int Joint2D::sendSelf(int commitTag, Channel &theChannel)
{
  int res;
  int i;
  int dataTag = this->getDbTag();
  
  static ID data(19);
  data(0) = this->getTag();
  data(1) = numDof;
  
  if (ExternalNodes.Size() != 0 && nodeDbTag == 0) nodeDbTag = theChannel.getDbTag();
  if (InternalConstraints.Size() != 0 && dofDbTag == 0) dofDbTag = theChannel.getDbTag();
  data(2) = nodeDbTag;
  data(3) = dofDbTag;
  
  // sending Sparing class and Db tags
  for (i=0 ; i<5 ; i++) {
    data( i+4 ) = fixedEnd[i];
    if ( theSprings[i] != NULL )
      {
	data( i+9 ) = theSprings[i]->getClassTag();
	int SpringDbTag = theSprings[i]->getDbTag();
	if (SpringDbTag == 0) {
	  SpringDbTag = theChannel.getDbTag();
	  if (SpringDbTag != 0)
	    theSprings[i]->setDbTag(SpringDbTag);
			}
	data(i+14) = SpringDbTag;
      } else
	{
	  data( i+9 ) = 0;
	  data( i+14 ) = 0;
	}
  }
  
  // send the ID vector
  
  res = theChannel.sendID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Joint2D::sendSelf() - " << this->getTag() << "failed to send ID\n";
    return -1;
  }
  
  // sends the tags of it's external nodes
  res = theChannel.sendID(nodeDbTag, commitTag, ExternalNodes);
  if (res < 0) {
    opserr << "WARNING Joint2D::sendSelf() - " << this->getTag()<< " failed to send Vector\n";
    return -2;
  }
  
  
  // sends the tags of it's internal constraints
  res = theChannel.sendID(dofDbTag, commitTag, InternalConstraints);
  if (res < 0) {
    opserr << "WARNING Joint2D::sendSelf() - %d failed to send Vector\n",this->getTag();
    return -2;
  }
  
  
  // finally send the materials one by one
  
  for ( i=0 ; i<5 ; i++ ) {
    if ( theSprings[i] != NULL ) {
      res = theSprings[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
	opserr << "WARNING Joint2D::sendSelf() - "<< this->getTag() << " failed to send its Spring " << (i+1) << " material\n";
	return -3;
      }
    }
  }
  
  return 0;
}

int Joint2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  
  int res;
  int dataTag = this->getDbTag();
  
  static ID data(19);
  res = theChannel.recvID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Joint2D::recvSelf() - failed to receive Vector\n";
    return -1;
  }
  
  this->setTag((int)data(0));
  numDof = data(1);
  nodeDbTag = data(2);
  dofDbTag = data(3);
  
  // Receving Springs
  for (int i=0 ; i<5 ; i++) {
    fixedEnd[i] = data( i+4 );
    int SpringClass = data( i+9 );
    int SpringDb = data( i+14 );
    
    if ( SpringClass != 0 && SpringDb != 0  && fixedEnd[i] == 0 )
      {
	// check if we have a material object already & if we do if of right type
	if ((theSprings[i] == 0) || (theSprings[i]->getClassTag() != SpringClass)) {
	  // if old one .. delete it
	  if (theSprings[i] != 0)
					delete theSprings[i];
	  // create a new material object
	  theSprings[i] = theBroker.getNewUniaxialMaterial(SpringClass);
	  if (theSprings[i] == 0) {
	    opserr << "WARNING Joint2D::recvSelf() - " << (i+1) << " failed to get a blank Material of type " << this->getTag() << " for Spring " << SpringClass << "\n";
	    return -3;
	  }
	}
	
	theSprings[i]->setDbTag(SpringDb); // note: we set the dbTag before we receive the material
	res = theSprings[i]->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "WARNING Joint2D::recvSelf() - " << this->getTag() << " failed to receive its Material for Spring " << (i+1) << "\n";
	  return -3;
	}
      }
    else theSprings[i] = NULL;
  }
  
  
  // receives the tags of it's external nodes
  res = theChannel.recvID(nodeDbTag, commitTag, ExternalNodes);
  if (res < 0) {
    opserr << "WARNING Joint2D::recvSelf() - " << this->getTag() << " failed to receive external nodes\n" ;
    return -2;
  }
  
  
  // receives the tags of it's constraint tags
  res = theChannel.recvID(dofDbTag, commitTag, InternalConstraints);
  if (res < 0) {
    opserr << "WARNING Joint2D::recvSelf() - " << this->getTag() << " failed to receive internal constraints\n";
    return -2;
  }
  
  return 0;
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int 
Joint2D::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
{

  if (theLoadSens == 0) {
    theLoadSens = new Vector(numDof);
  }
  else {
    theLoadSens->Zero();
  }

  return 0;
}


int
Joint2D::setParameter (const char **argv, int argc, Information &info)
{//
	// From the parameterID value it should be possible to extract
	// information about:
	//  1) Which parameter is in question. The parameter could
	//     be at element, section, or material level. 
	//  2) Which section and material number (tag) it belongs to. 
	//
	// To accomplish this the parameterID is given the following value:
	//     parameterID = type + 1000*matrTag + 100000*sectionTag
	// ...where 'type' is an integer in the range (1-99) and added 100
	// for each level (from material to section to element). 
	//
	// Example:
	//    If 'E0' (case 2) is random in material #3 of section #5
	//    the value of the parameterID at this (element) level would be:
	//    parameterID = 2 + 1000*3 + 100000*5 = 503002
	//    As seen, all given information can be extracted from this number. 
	//

	int ok = -1;
	if (argc < 3)
        return ok;

    // a material parameter
    if (strcmp(argv[0],"-material") == 0 || strcmp(argv[0],"material") == 0) {
		// Get material tag numbers from user input
		int paramMaterialTag = atoi(argv[1]);
		if ( paramMaterialTag<0 || paramMaterialTag>4 ) {
			opserr << "Joint2D::setParameter() - material number out of range, must be 0 to 4." << endln;
		return -1;
		}

		// Find the material and call its setParameter method
		if ( theSprings[paramMaterialTag] != NULL )
			ok = theSprings[paramMaterialTag]->setParameter(&argv[2], argc-2, info);
		
		// Check if the ok is valid
		if (ok < 0) {
			opserr << "Joint2D::setParameter() - could not set parameter. " << endln;
			return -1;
		}
		else {
			// Return the ok value (according to the above comments)
			return 1000*paramMaterialTag + ok;
		}
	}		
		 
    // otherwise parameter is unknown for the Joint2D class
      return -1;
}

int
Joint2D::updateParameter (int parameterID, Information &info)
{
	// Extract section and material tags from the passedParameterID
	int MaterialTag = (int)( floor( (double)parameterID / 1000.0 ) );
	int MaterialParameterID = parameterID-1000*MaterialTag;
	if ( MaterialTag<0 || MaterialTag>4 ) {
		opserr << "Joint2D::updateParameter() - material number out of range, must be 0 to 4." << endln;
		return -1;
	}
	
	// Go down to the material and set appropriate flags
	if ( theSprings[MaterialTag] != NULL )
		return theSprings[MaterialTag]->updateParameter(MaterialParameterID,info);
	
	return -1;
}


int
Joint2D::activateParameter(int passedParameterID)
{
	// Note that the parameteID that is stored here at the 
	// element level contains all information about section
	// and material tag number:
	parameterID = passedParameterID;

	if (passedParameterID == 0 ) {
		// "Zero out" all flags downwards through sections/materials 
		for (int i=0; i<5; i++) {
			if ( theSprings[i] !=0 )
				theSprings[i]->activateParameter(passedParameterID);
		}
	}
	else {
		
		// Extract section and material tags from the passedParameterID
		int activeMaterialTag = (int)( floor( (double)passedParameterID / 1000.0 ) );

		if ( activeMaterialTag<0 || activeMaterialTag>4 ) {
			opserr << "Joint2D::activateParameter() - material number out of range, must be 0 to 4." << endln;
			return -1;
		}
		
		// Go down to the material and set appropriate flags
		if ( theSprings[activeMaterialTag] != NULL ) {
			theSprings[activeMaterialTag]->activateParameter(passedParameterID-1000*activeMaterialTag);
		}
	}
	
	return 0;
}


const Matrix &
Joint2D::getKiSensitivity(int gradNumber)
{
	K.Zero();
    
	if (parameterID == 0) {
		// // Nothing here
	}
	
	else {
		double KtangentSensitivity[5] ;
		for ( int i=0 ; i<5 ; i++ ) 
		{
			KtangentSensitivity[i] = 0;
			if ( theSprings[i] != NULL ) 
				KtangentSensitivity[i] = theSprings[i]->getInitialTangentSensitivity(gradNumber);
		}
		
		K(2,2)  =  KtangentSensitivity[0];
		K(2,15) = -KtangentSensitivity[0];
		K(5,5)  =  KtangentSensitivity[1];
		K(5,14) = -KtangentSensitivity[1];
		K(8,8)  =  KtangentSensitivity[2];
		K(8,15) = -KtangentSensitivity[2];
		K(11,11)=  KtangentSensitivity[3];
		K(11,14)= -KtangentSensitivity[3];
		K(14,5) = -KtangentSensitivity[1];
		K(14,11)= -KtangentSensitivity[3];
		K(14,14)=  KtangentSensitivity[1] + KtangentSensitivity[3] + KtangentSensitivity[4];
		K(14,15)= -KtangentSensitivity[4];
		K(15,2) = -KtangentSensitivity[0];
		K(15,8) = -KtangentSensitivity[2];
		K(15,14)= -KtangentSensitivity[4];
		K(15,15)=  KtangentSensitivity[0] + KtangentSensitivity[2] + KtangentSensitivity[4];
	}
	
	return K;
}


const Matrix &
Joint2D::getMassSensitivity(int gradNumber)
{
  	K.Zero();
	return K;
}


const Vector &
Joint2D::getResistingForceSensitivity(int gradNumber)
{
	this->update();
	V.Zero();
	
	// Compute sensitivity depending on the material
	double ForceSensitivity[5] ;
	for ( int i=0 ; i<5 ; i++ ) 
	{
		ForceSensitivity[i] = 0;
		if ( theSprings[i] != NULL ) ForceSensitivity[i] = theSprings[i]->getStressSensitivity(gradNumber,true);
	}

	V(2) = ForceSensitivity[0];
	V(5) = ForceSensitivity[1];
	V(8) = ForceSensitivity[2];
	V(11)= ForceSensitivity[3];
	V(14)= -ForceSensitivity[4] - ForceSensitivity[1] - ForceSensitivity[3];
	V(15)= ForceSensitivity[4] - ForceSensitivity[0] - ForceSensitivity[2];

	return V;
}


int
Joint2D::commitSensitivity(int gradNumber, int numGrads)
{

	double strainSensitivity = 0.0;
	// Pass it down to the material
	for ( int i=0 ; i<5 ; i++ )
		if ( theSprings[i] != NULL )
			theSprings[i]->commitSensitivity(strainSensitivity, gradNumber, numGrads);

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
