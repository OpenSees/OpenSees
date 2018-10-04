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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-02-02 01:44:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint3d.cpp,v $
                                                                        
// Written: NM (nmitra@u.washington.edu)
// Created: February 2003
// Revised: September 2004
//
// Description: This file contains the class implementation for beam-column joint.
// This element is a 4 noded 24 dof (6 dof at each node) finite area super-element, being a slight
// variation of the 2d one. The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint. Though it has 6 dof per node 
// the out of the plane nodal dof are constrained or fixed and the inplane nodal dof are activated

#include <BeamColumnJoint3d.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>
#include <UniaxialMaterial.h>
#include <string.h>
#include <math.h>
#include <ElementResponse.h>
#include <elementAPI.h>

// class wide matrices 
Matrix BeamColumnJoint3d::Transf(12,24);
Matrix BeamColumnJoint3d::Tran(3,6);


void* OPS_BeamColumnJoint3d()
{
    if (OPS_GetNumRemainingInputArgs() < 18) {
	opserr << "WARNING insufficient arguments\n";
	  opserr << "Want: element beamColumnJoint eleTag? node1? node2? node3? node4? matTag1? matTag2? matTag3?\n";
	  opserr << "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? matTag10? matTag11? matTag12? matTag13?\n";
	  opserr << "<ElementHeightFactor? ElementWidthFactor?>\n";
	  return 0;
    }

    int idata[18];
    int numdata = 18;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    double data[2] = {1.0, 1.0};
    numdata = 2;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr<<"WARNING: invalid double inputs\n";
	    return 0;
	}
    }

    UniaxialMaterial* mats[13];
    for (int i = 0; i < 13; i++) {
	mats[i] = OPS_getUniaxialMaterial(idata[5+i]);
	if (mats[i] == 0) {
	    opserr<<"WARNING: material "<<idata[5+i]<<" is not defined\n";
	    return 0;
	}
    }

    return new BeamColumnJoint3d(idata[0],idata[1],idata[2],idata[3],idata[4],
				 *mats[0],*mats[1],*mats[2],*mats[3],*mats[4],
				 *mats[5],*mats[6],*mats[7],*mats[8],*mats[9],
				 *mats[10],*mats[11],*mats[12],data[0],data[1]);

}


// full constructors:
BeamColumnJoint3d::BeamColumnJoint3d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3,
				     UniaxialMaterial& theMat4,
				     UniaxialMaterial& theMat5,
				     UniaxialMaterial& theMat6,
				     UniaxialMaterial& theMat7,
				     UniaxialMaterial& theMat8,
				     UniaxialMaterial& theMat9,
				     UniaxialMaterial& theMat10,
				     UniaxialMaterial& theMat11,
				     UniaxialMaterial& theMat12,
				     UniaxialMaterial& theMat13):
  Element(tag,ELE_TAG_BeamColumnJoint3d),
  connectedExternalNodes(4), elemActHeight(0.0), elemActWidth(0.0), 
  elemWidth(0.0), elemHeight(0.0), HgtFac(1.0), WdtFac(1.0),
  Uecommit(24), UeIntcommit(4), UeprCommit(24), UeprIntCommit(4),
  BCJoint(13,16), dg_df(4,13), dDef_du(13,4), K(24,24), R(24), Node1(3), Node2(3), Node3(3), Node4(3)
{  
// ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 4)
         opserr << "ERROR : BeamColumnJoint::BeamColumnJoint - " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    connectedExternalNodes(2) = Nd3 ;
    connectedExternalNodes(3) = Nd4 ;

	MaterialPtr = new UniaxialMaterial*[13];
	for (int x = 0; x <13; x++)
	{	MaterialPtr[x] = 0; }

	Uecommit.Zero();
	UeIntcommit.Zero();
	UeprCommit.Zero();
	UeprIntCommit.Zero();

	BCJoint.Zero(); dg_df.Zero(); dDef_du.Zero();
	K.Zero();
	R.Zero();
	Node1.Zero(); Node2.Zero(); Node3.Zero(); Node4.Zero();

	nodePtr[0] = 0;
	nodePtr[1] = 0;

	// get a copy of the material and check we obtained a valid copy
    MaterialPtr[0] = theMat1.getCopy();
    if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
    MaterialPtr[1] = theMat2.getCopy();
    if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2" << endln;}
    MaterialPtr[2] = theMat3.getCopy();
    if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3" << endln;}
	MaterialPtr[3] = theMat4.getCopy();
	if (!MaterialPtr[3]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 4" << endln;}
	MaterialPtr[4] = theMat5.getCopy();
	if (!MaterialPtr[4]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 5" << endln;}
	MaterialPtr[5] = theMat6.getCopy();
	if (!MaterialPtr[5]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 6" << endln;}
	MaterialPtr[6] = theMat7.getCopy();
	if (!MaterialPtr[6]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 7" << endln;}
	MaterialPtr[7] = theMat8.getCopy();
	if (!MaterialPtr[7]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 8" << endln;}
	MaterialPtr[8] = theMat9.getCopy();
	if (!MaterialPtr[8]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 9" << endln;}
	MaterialPtr[9] = theMat10.getCopy();
	if (!MaterialPtr[9]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 10" << endln;}
	MaterialPtr[10] = theMat11.getCopy();
	if (!MaterialPtr[10]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 11" << endln;}
	MaterialPtr[11] = theMat12.getCopy();
	if (!MaterialPtr[11]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 12" << endln;}
	MaterialPtr[12] = theMat13.getCopy();
	if (!MaterialPtr[12]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 13" << endln;}	

}

// full constructors:
BeamColumnJoint3d::BeamColumnJoint3d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3,
				     UniaxialMaterial& theMat4,
				     UniaxialMaterial& theMat5,
				     UniaxialMaterial& theMat6,
				     UniaxialMaterial& theMat7,
				     UniaxialMaterial& theMat8,
				     UniaxialMaterial& theMat9,
				     UniaxialMaterial& theMat10,
				     UniaxialMaterial& theMat11,
				     UniaxialMaterial& theMat12,
				     UniaxialMaterial& theMat13,
					 double elHgtFac, double elWdtFac):
  Element(tag,ELE_TAG_BeamColumnJoint3d),
  connectedExternalNodes(4), elemActHeight(0.0), elemActWidth(0.0), 
  elemWidth(0), elemHeight(0), HgtFac(elHgtFac), WdtFac(elWdtFac),
  Uecommit(24), UeIntcommit(4), UeprCommit(24), UeprIntCommit(4),
  BCJoint(13,16), dg_df(4,13), dDef_du(13,4), K(24,24), R(24), Node1(3), Node2(3), Node3(3), Node4(3)
{  
// ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 4)
         opserr << "ERROR : BeamColumnJoint::BeamColumnJoint - " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    connectedExternalNodes(2) = Nd3 ;
    connectedExternalNodes(3) = Nd4 ;

	MaterialPtr = new UniaxialMaterial*[13];
	for (int x = 0; x <13; x++)
	{	MaterialPtr[x] = 0; }

	Uecommit.Zero();
	UeIntcommit.Zero();
	UeprCommit.Zero();
	UeprIntCommit.Zero();

	BCJoint.Zero(); dg_df.Zero(); dDef_du.Zero();
	K.Zero();
	R.Zero();
	Node1.Zero(); Node2.Zero(); Node3.Zero(); Node4.Zero();

	nodePtr[0] = 0;
	nodePtr[1] = 0;

	// get a copy of the material and check we obtained a valid copy
    MaterialPtr[0] = theMat1.getCopy();
    if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
    MaterialPtr[1] = theMat2.getCopy();
    if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2" << endln;}
    MaterialPtr[2] = theMat3.getCopy();
    if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3" << endln;}
	MaterialPtr[3] = theMat4.getCopy();
	if (!MaterialPtr[3]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 4" << endln;}
	MaterialPtr[4] = theMat5.getCopy();
	if (!MaterialPtr[4]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 5" << endln;}
	MaterialPtr[5] = theMat6.getCopy();
	if (!MaterialPtr[5]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 6" << endln;}
	MaterialPtr[6] = theMat7.getCopy();
	if (!MaterialPtr[6]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 7" << endln;}
	MaterialPtr[7] = theMat8.getCopy();
	if (!MaterialPtr[7]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 8" << endln;}
	MaterialPtr[8] = theMat9.getCopy();
	if (!MaterialPtr[8]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 9" << endln;}
	MaterialPtr[9] = theMat10.getCopy();
	if (!MaterialPtr[9]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 10" << endln;}
	MaterialPtr[10] = theMat11.getCopy();
	if (!MaterialPtr[10]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 11" << endln;}
	MaterialPtr[11] = theMat12.getCopy();
	if (!MaterialPtr[11]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 12" << endln;}
	MaterialPtr[12] = theMat13.getCopy();
	if (!MaterialPtr[12]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 13" << endln;}	

}

// default constructor:
BeamColumnJoint3d::BeamColumnJoint3d():
  Element(0,ELE_TAG_BeamColumnJoint3d),
  connectedExternalNodes(4), elemActHeight(0.0), elemActWidth(0.0), 
  elemWidth(0), elemHeight(0), HgtFac(1.0), WdtFac(1.0),
  Uecommit(24), UeIntcommit(4), UeprCommit(24), UeprIntCommit(4),
  BCJoint(13,16), dg_df(4,13), dDef_du(13,4), K(24,24), R(24), Node1(3), Node2(3), Node3(3), Node4(3)
{
    nodePtr[0] = 0;
	nodePtr[1] = 0;
	for (int x = 0; x <13; x++)
	{	MaterialPtr[x] = 0; }
   // does nothing (invoked by FEM_ObjectBroker)
}

//  destructor:
BeamColumnJoint3d::~BeamColumnJoint3d()
{
	for (int i =0; i<13; i++)
	{
	    if (MaterialPtr[i] != 0)
			delete MaterialPtr[i];
	}

	 if (MaterialPtr)
		 delete [] MaterialPtr;
}

// public methods

int
BeamColumnJoint3d::getNumExternalNodes(void) const
{
    return 4;
}

const ID &
BeamColumnJoint3d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **BeamColumnJoint3d::getNodePtrs(void)
{
	return nodePtr;
}

int
BeamColumnJoint3d::getNumDOF(void) 
{
    return 24;
}

void
BeamColumnJoint3d::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	opserr << "ERROR : BeamColumnJoint::setDomain -- Domain is null" << endln;

	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePtr[0] = 0;
	nodePtr[1] = 0;
    }

  //node pointers set
    for (int i = 0; i < 4; i++ ) {
		nodePtr[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR : BeamColumnJoint::setDomain -- node pointer is null" << endln;
			exit(-1); // donot go any further - otherwise segmentation fault
		}
	}

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

	// ensure connected nodes have correct dof's
	int dofNd1 = nodePtr[0]->getNumberDOF();
	int dofNd2 = nodePtr[1]->getNumberDOF();
	int dofNd3 = nodePtr[2]->getNumberDOF();
	int dofNd4 = nodePtr[3]->getNumberDOF();

	// check for proper nodal degrees of freedom
	if ((dofNd1 != 6) || (dofNd2 != 6) || (dofNd3 != 6) || (dofNd4 != 6))
	{
		opserr << "ERROR : BeamColumnJoint::setDomain -- number of DOF associated with the node incorrect" << endln;
		exit(-1); // donot go any further - otherwise segmentation fault
	}

    // obtain the nodal coordinates    
	const Vector &end1Crd = nodePtr[0]->getCrds();
	const Vector &end2Crd = nodePtr[1]->getCrds();
	const Vector &end3Crd = nodePtr[2]->getCrds();
	const Vector &end4Crd = nodePtr[3]->getCrds();

	// determine the element width and the element height. It should be noted that the element width is determined by nodal coordinates 2 and 4
	// whereas the element height is determined by the nodal coordinates 1 and 3
	// general order anticlockwise arrangement , but can be clockwise too. 
	Node1 = end1Crd;
	Node2 = end2Crd;
	Node3 = end3Crd;
	Node4 = end4Crd;
	
	Node3 = Node3 - Node1;
	Node2 = Node2 - Node4;

	double beam = 1.0; //HgtFac;
	double col = 1.0; //WdtFac;

	elemActHeight = fabs(Node3.Norm());
	elemActWidth = fabs(Node2.Norm());
	elemHeight = beam*elemActHeight;
	elemWidth = col*elemActWidth;

	if ((elemHeight <= 1e-12) || (elemWidth <= 1e-12))
	{
			opserr << "ERROR : BeamColumnJoint::setDomain -- length or width not correct, division by zero occurs"<< endln;
			exit(-1); // donot go any further - otherwise segmentation fault
	}

	getBCJoint();
	getdg_df();
	getdDef_du();
	formTransfMat();
}   

int
BeamColumnJoint3d::commitState()
{
	// store committed external nodal displacements
  Uecommit = UeprCommit;
  // store committed internal nodal displacements
  UeIntcommit = UeprIntCommit;

  // store material history data.
  int mcs = 0;
	for (int j=0; j<13; j++)
	{
			if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->commitState();
			if (mcs != 0) break;
	}
    
  return mcs;
}

int
BeamColumnJoint3d::revertToLastCommit()
{

	int mcs = 0;
	for (int j=0; j<13; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToLastCommit();
		if (mcs != 0) break;
	}
	UeprCommit = Uecommit;
	UeprIntCommit = UeIntcommit;   
	
	this->update();

  return mcs;
}

int
BeamColumnJoint3d::revertToStart()
{
	int mcs = 0;
	for (int j=0; j<13; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToStart();
		if (mcs != 0) break;
	}
	
	return mcs;
}

int
BeamColumnJoint3d::update()
{	

	Vector Ue(28);
	Ue.Zero();

	// determine committed displacements given trial displacements
	this->getGlobalDispls(Ue);
 
	// update displacements for the external nodes
	UeprCommit.Extract(Ue,0,1.0); 

	// update displacement for the internal nodes
	UeprIntCommit.Extract(Ue,24,1.0);  

	return 0;
}


const Matrix &
BeamColumnJoint3d::getTangentStiff(void)
{
	return K;
}

const Matrix &
BeamColumnJoint3d::getInitialStiff(void)
{
	return getTangentStiff();
}

const Vector &
BeamColumnJoint3d::getResistingForce()
{
  return R;
}

void BeamColumnJoint3d::getGlobalDispls(Vector &dg) 
{
	// local variables that will be used in this method
	int converge = 0;
	int linesearch = 0;
	int totalCount = 0;
	int dtConverge = 0;
	int incCount = 0;
	int count = 0;
	int maxTotalCount = 1000;
	int maxCount = 20;
	double loadStep = 0.0;
	double dLoadStep = 1.0;
	double stepSize;

	Vector uExtOld(24);       uExtOld.Zero();
	Vector uExt(12);          uExt.Zero();
	Vector duExt(12);         duExt.Zero();
	Vector uIntOld(4);        uIntOld.Zero();
	Vector uInt(4);           uInt.Zero();
	Vector duInt(4);          duInt.Zero();
	Vector duIntTemp(4);      duIntTemp.Zero();
	Vector intEq(4);          intEq.Zero();
	Vector intEqLast(4);      intEqLast.Zero();
	Vector Uepr(24);          Uepr.Zero();
	Vector UeprInt(4);        UeprInt.Zero();
	Vector Ut(24);            Ut.Zero();
	Vector duExtTemp(24);     duExtTemp.Zero();

    Vector disp1 = nodePtr[0]->getTrialDisp(); 
    Vector disp2 = nodePtr[1]->getTrialDisp();
    Vector disp3 = nodePtr[2]->getTrialDisp();
    Vector disp4 = nodePtr[3]->getTrialDisp();

	for (int i = 0; i < 6; i++)
    {
      Ut(i)     = disp1(i);
      Ut(i+6)   = disp2(i);
      Ut(i+12)   = disp3(i);
      Ut(i+18)   = disp4(i);
    }
	
	Uepr = Uecommit;   

	UeprInt = UeIntcommit;  

	uExtOld = Uepr;

	duExtTemp = Ut - Uepr;
	duExt.addMatrixVector(0.0,Transf,duExtTemp,1.0);
	uExt.addMatrixVector(0.0,Transf,uExtOld,1.0);  

	uIntOld = UeprInt;
	uInt = uIntOld;

	double tol = 1e-12;
	double tolIntEq = tol;
	double toluInt = (tol>tol*uInt.Norm())? tol:tol*uInt.Norm();
	double tolIntEqdU = tol;
	double ctolIntEqdU = tol;
	double ctolIntEq = tol;
	double normDuInt = toluInt;
	double normIntEq = tolIntEq;
	double normIntEqdU = tolIntEqdU;
	    
	Vector u(16);   u.Zero();

	double engrLast = 0.0;
	double engr = 0.0;

	Vector fSpring(13);          fSpring.Zero();
	Vector kSpring(13);          kSpring.Zero();
	Matrix dintEq_du(4,4);       dintEq_du.Zero();
	Matrix df_dDef(13,13);       df_dDef.Zero();
	Matrix tempintEq_du (4,13);  tempintEq_du.Zero();


	while ((loadStep < 1.0) && (totalCount < maxTotalCount))
	{
		count = 0;
		converge = 0;
		dtConverge = 0;
		while ((!converge) && (count < maxCount))
		{
			if (dLoadStep <= 1e-3)
			{
				dLoadStep = dLoadStep;
			}
			totalCount ++;
			count ++;
			
			for (int ic = 0; ic < 12; ic++ ) {
				u(ic) = uExt(ic) + duExt(ic);
			}
			u(12) = uInt(0);
			u(13) = uInt(1);
			u(14) = uInt(2);
			u(15) = uInt(3);

			getMatResponse(u,fSpring,kSpring);
		

		// performs internal equilibrium 

		intEq(0) = -fSpring(2)-fSpring(3)+fSpring(9)-fSpring(12)/elemHeight; 
		intEq(1) = fSpring(1)-fSpring(5)-fSpring(7)+fSpring(12)/elemWidth; 
		intEq(2) = -fSpring(4)-fSpring(8)+fSpring(10)+fSpring(12)/elemHeight; 
		intEq(3) = fSpring(0)-fSpring(6)-fSpring(11)-fSpring(12)/elemWidth; 

		matDiag(kSpring, df_dDef);

		//////////////////////// dintEq_du = dg_df*df_dDef*dDef_du
		tempintEq_du.addMatrixProduct(0.0,dg_df,df_dDef,1.0);
		dintEq_du.addMatrixProduct(0.0,tempintEq_du,dDef_du,1.0);
		normIntEq = intEq.Norm();
		normIntEqdU = 0.0;
		for (int jc = 0; jc<4 ; jc++)
		{
			normIntEqdU += intEq(jc)*duInt(jc);
		}
		normIntEqdU = fabs(normIntEqdU);

		if (totalCount == 1)
		{
			tolIntEq = (tol>tol*normIntEq) ? tol:tol*normIntEq;
			tolIntEqdU = tol;
		}
		else if (totalCount == 2)
		{
			tolIntEqdU = (tol>tol*normIntEqdU) ? tol:tol*normIntEqdU;
		}
		ctolIntEqdU = (tolIntEqdU*dLoadStep > tol) ? tolIntEqdU*dLoadStep:tol;
		ctolIntEq   = (tolIntEq*dLoadStep > tol) ? tolIntEq*dLoadStep:tol;

		// check for convergence starts  
		if ((normIntEq < tol) || ((normIntEqdU < tol) && (count >1)) || (normDuInt < toluInt) || (dLoadStep < 1e-3))
		{
		  	if ((normIntEq > ctolIntEq) || (normIntEqdU > tolIntEqdU) || (normDuInt > toluInt))
			{
				dtConverge = 1;
			}
			else
			{
				dtConverge = 0;
			}

			converge = 1;
			loadStep = loadStep + dLoadStep;
			if (fabs(1.0 - loadStep) < tol)
			{
				loadStep = 1.0;
			}
		}
		else
		{
			////////////// duInt = -dintEq_du/intEq
			dintEq_du.Solve(intEq,duInt);
			duInt *= -1;
			
			normDuInt = duInt.Norm();
			if (!linesearch)
			{
				uInt = uInt + duInt;
			}
			else
			{
				engrLast = 0.0;
				engr = 0.0;

				for (int jd = 0; jd<4 ; jd++)
				{
					engrLast += duInt(jd)*intEqLast(jd);
					engr += duInt(jd)*intEq(jd);
				}

				if (fabs(engr) > tol*engrLast)
				{
					duIntTemp = duInt;
					duIntTemp *= -1;
					// lineSearch algorithm requirement
					stepSize = getStepSize(engrLast,engr,uExt,duExt,uInt,duIntTemp,tol);
					
					if (fabs(stepSize) > 0.001)
					{
						uInt = uInt + stepSize*duInt;
					}
					else
					{
						uInt = uInt + duInt;
					}
				}
				else
				{
					uInt = uInt + duInt;
				}
				intEqLast = intEq;
			}
		}
	}

		if (!converge && loadStep < 1.0)
		{	
			incCount = 0;
			maxCount = 25;
			if (!linesearch)
			{
				linesearch = 1;
				uInt = uIntOld;
				duInt.Zero();
			}
			else
			{
				opserr << "WARNING : BeamColumnJoint::getGlobalDispls() - convergence problem in state determination" << endln;

				uInt = uIntOld;
				duInt.Zero();
				duExt = duExt*0.1;

				dLoadStep = dLoadStep*0.1;
			}
		}
		else if (loadStep < 1.0)
		{
			maxCount = 10;
			incCount ++;
			normDuInt = toluInt;
			if ((incCount < maxCount) || dtConverge)
			{
				uExt = uExt + duExt;
				if (loadStep + dLoadStep > 1.0)
				{
					duExt = duExt*(1.0 - loadStep)/dLoadStep;
					dLoadStep = 1.0 - loadStep;
					incCount = 9;
				}
			}
			else
			{
				incCount = 0;
				uExt = uExt + duExt;
				dLoadStep = dLoadStep*10;
				if (loadStep + dLoadStep > 1.0)
				{
					uExt = uExt + duExt*(1.0 - loadStep)/dLoadStep;
					dLoadStep = 1.0 - loadStep;
					incCount = 9;
				}
			}
		}
	}


	// determination of stiffness matrix and the residual force vector for the element

	formR(fSpring);

	formK(kSpring);

	for (int ig = 0; ig < 25; ig++ ) {
		if (ig<24)
		{
			dg(ig) = Ut(ig);
		}
	}


	dg(24) = uInt(0);
	dg(25) = uInt(1);
	dg(26) = uInt(2);
	dg(27) = uInt(3);

}

void BeamColumnJoint3d::getMatResponse(Vector U, Vector &fS, Vector &kS)
{
	double jh = HgtFac;            // factor for beams
	double jw = WdtFac;            // factor for column

	// obtains the material response from the material class
	Vector defSpring(13);
	defSpring.Zero();
	fS.Zero();
	kS.Zero();

	defSpring.addMatrixVector(0.0, BCJoint, U, 1.0);

	// slip @ bar = slip @ spring * tc couple

	defSpring(0) = defSpring(0)*jw;                   // new model applied on 
	defSpring(1) = defSpring(1)*jw;                   // 4th June 2004
	defSpring(6) = defSpring(6)*jw;                   // h, h_bar distinction was
	defSpring(7) = defSpring(7)*jw;                   // removed from the previous model

	defSpring(3)  = defSpring(3)*jh;                  // by making h = h_bar
	defSpring(4)  = defSpring(4)*jh;                  // similar case done for width
	defSpring(9) = defSpring(9)*jh;
	defSpring(10) = defSpring(10)*jh;

	for (int j=0; j<13; j++)
	{
		MaterialPtr[j]->setTrialStrain(defSpring[j]);
		kS[j] = MaterialPtr[j]->getTangent();
		fS[j] = MaterialPtr[j]->getStress();
	}

	// force @ spring = force @ bar * tc couple

	fS(0) = fS(0)*jw;
	fS(1) = fS(1)*jw;
	fS(6) = fS(6)*jw;
	fS(7) = fS(7)*jw;

	fS(3) = fS(3)*jh;
	fS(4) = fS(4)*jh;
	fS(9) = fS(9)*jh;
	fS(10) = fS(10)*jh;

	// stiffness @ spring = stiffness @ bar * (tc couple)^2

	kS(0) = kS(0)*jw*jw;
	kS(1) = kS(1)*jw*jw;
	kS(6) = kS(6)*jw*jw;
	kS(7) = kS(7)*jw*jw;

	kS(3) = kS(3)*jh*jh;
	kS(4) = kS(4)*jh*jh;
	kS(9) = kS(9)*jh*jh;
	kS(10) = kS(10)*jh*jh;
}


void BeamColumnJoint3d::getdDef_du()
{
	
	dDef_du.Zero();	
	for (int jb=0; jb<13; jb++)
	{
		dDef_du(jb,0) = BCJoint(jb,12);
		dDef_du(jb,1) = BCJoint(jb,13);
		dDef_du(jb,2) = BCJoint(jb,14);
		dDef_du(jb,3) = BCJoint(jb,15);
	}
}

void BeamColumnJoint3d::matDiag(Vector k,Matrix &dfd)
{
	dfd.Zero();
	// takes in a vector and converts it to a diagonal matrix (could have been placed as a method in matrix class)
	for (int ja=0; ja<13; ja++)
	{
		dfd(ja,ja) = k(ja);
	}

}

void BeamColumnJoint3d::formR(Vector f)
{
	
	// develops the element residual force vector
	Vector rForceTemp(16);
	Vector Rtempo(12);
	rForceTemp.Zero();
	
	rForceTemp.addMatrixTransposeVector(0.0,BCJoint,f,1.0);   // rForceTemp = BCJoint'*f
	Rtempo.Extract(rForceTemp,0,1.0);

	R.addMatrixTransposeVector(0.0,Transf,Rtempo,1.0);    // R = Transf'*Rtempo
}

void BeamColumnJoint3d::formK(Vector k)
{
    // develops the element stiffness matrix
	Matrix kSprDiag(13,13);
	kSprDiag.Zero();
    Matrix kRForce(16,16);
	kRForce.Zero();
	Matrix kRFT1(4,12);
	kRFT1.Zero();
	Matrix kRFT2(4,4);
	kRFT2.Zero();
	Matrix kRFT3(12,4);
	kRFT3.Zero();
	Matrix I(4,4);
	I.Zero();
	Matrix kRSTinv(4,4);
	kRSTinv.Zero();
    Matrix kRF(12,12);
	kRF.Zero();
	Matrix K2Temp(12,4);
	K2Temp.Zero();
	Matrix K2(12,12);
	K2.Zero();

	matDiag(k,kSprDiag);

	kRForce.addMatrixTripleProduct(0.0,BCJoint,kSprDiag,1.0);    // kRForce = BCJoint'*kSprDiag*BCJoint
	kRFT2.Extract(kRForce,12,12,1.0);
	kRFT1.Extract(kRForce,12,0,1.0);
	kRFT3.Extract(kRForce,0,12,1.0);
	kRF.Extract(kRForce,0,0,1.0);
	
    for (int ic=0; ic<4; ic++)
	{
		I(ic,ic) = 1.0;
	}
	kRFT2.Solve(I,kRSTinv);
	
	K2Temp.addMatrixProduct(0.0,kRFT3,kRSTinv,1.0);

	// done for some idiotic reason
	for(int i = 0; i < 12; ++i)
	{	
		for(int j = 0; j < 4; ++j)
		{
			if(fabs(K2Temp(i,j)) < 1e-15)
				K2Temp(i,j) = 0.;
		}
	}

	K2.addMatrixProduct(0.0,K2Temp,kRFT1,1.0);      // K2 = kRFT3*kRSTinv*kRFT1 

    // done for some idiotic reason
	for(int i1 = 0; i1 < 12; ++i1)
	{	
		for(int j1 = 0; j1 < 12; ++j1)
		{
			if(fabs(K2(i1,j1)) < 1e-15)
				K2(i1,j1) = 0.;
		}
	}

	kRF.addMatrix(1.0,K2,-1.0);                     // kRF = kRF - K2

	K.addMatrixTripleProduct(0.0,Transf,kRF,1.0);   // K = Transf'*kRF*Transf
}

void BeamColumnJoint3d::formTransfMat()
{
	Transf.Zero();
	Tran.Zero();
	
	double Crd24 = fabs(sqrt((Node2(0)*Node2(0)) + (Node2(1)*Node2(1)) + (Node2(2)*Node2(2))));
	double Crd31 = fabs(sqrt((Node3(0)*Node3(0)) + (Node3(1)*Node3(1)) + (Node3(2)*Node3(2))));

	double a24 = Node2(0)/Crd24; double b24 = Node2(1)/Crd24; double c24 = Node2(2)/Crd24; 
    double a31 = Node3(0)/Crd31; double b31 = Node3(1)/Crd31; double c31 = Node3(2)/Crd31; 

	Tran(0,0) = a24; Tran(0,1) =  b24; Tran(0,2) = c24;
	Tran(1,0) = a31; Tran(1,1) =  b31; Tran(1,2) = c31;
	Tran(2,3) = b24*c31 - b31*c24; Tran(2,4) = -a24*c31 + c24*a31; Tran(2,5) = a24*b31 - b24*a31;
  
	Transf.Assemble(Tran,0,0,1.0);
	Transf.Assemble(Tran,3,6,1.0);
	Transf.Assemble(Tran,6,12,1.0);
	Transf.Assemble(Tran,9,18,1.0);

}

void BeamColumnJoint3d::getdg_df()
{
		dg_df.Zero();        // basically derivative of intEq
		dg_df(0,2) = -1;
		dg_df(0,3) = -1;
		dg_df(0,9) = 1;
		dg_df(0,12) = -1/elemHeight; 
		dg_df(1,1) = 1;
		dg_df(1,5) = -1;
		dg_df(1,7) = -1;
		dg_df(1,12) = 1/elemWidth;
		dg_df(2,4) = -1;
		dg_df(2,8) = -1;
		dg_df(2,10) = 1;
		dg_df(2,12) = 1/elemHeight;
		dg_df(3,0) = 1;
		dg_df(3,6) = -1;
		dg_df(3,11) = -1;
		dg_df(3,12) = -1/elemWidth;

}

void BeamColumnJoint3d::getBCJoint() 
{
	BCJoint.Zero();                  // basically the transformation matrix for the element
	BCJoint(0,1) =  -1;
	BCJoint(0,2) =  elemWidth/2;
	BCJoint(0,15) = 1;
	BCJoint(1,1) =  -1;
	BCJoint(1,2) =  -elemWidth/2;
	BCJoint(1,13) = 1;
	BCJoint(2,0) =  1;
	BCJoint(2,12) = -1;
	BCJoint(3,3) =  1;
	BCJoint(3,5) =  elemHeight/2;
	BCJoint(3,12) = -1;
	BCJoint(4,3) = 1;
	BCJoint(4,5) = -elemHeight/2;
	BCJoint(4,14) = -1;
	BCJoint(5,4) = 1;
	BCJoint(5,13) = -1;
	BCJoint(6,7) =  1;
	BCJoint(6,8) =  -elemWidth/2;
	BCJoint(6,15) = -1;
	BCJoint(7,7) =  1;
	BCJoint(7,8) =  elemWidth/2;
	BCJoint(7,13) = -1;
	BCJoint(8,6) =  1;
	BCJoint(8,14) = -1;
	BCJoint(9,9) =  -1;
	BCJoint(9,11) =  -elemHeight/2;
	BCJoint(9,12) = 1;
	BCJoint(10,9) = -1;
	BCJoint(10,11) = elemHeight/2;
	BCJoint(10,14) = 1;
	BCJoint(11,10) = 1;
	BCJoint(11,15) = -1;
	BCJoint(12,12) = -1/elemHeight;
	BCJoint(12,13) = 1/elemWidth;
	BCJoint(12,14) = 1/elemHeight;
	BCJoint(12,15) = -1/elemWidth;

	// based upon the new formulation 
	BCJoint(2,2) = -(elemActHeight - elemHeight)/2;
	BCJoint(5,5) = -(elemActWidth - elemWidth)/2;
	BCJoint(8,8) = (elemActHeight - elemHeight)/2;;
	BCJoint(11,11) = (elemActWidth - elemWidth)/2;;

}

double BeamColumnJoint3d::getStepSize(double s0,double s1,Vector uExt,Vector duExt,Vector uInt,Vector duInt,double tol)
{
	// finds out the factor to be used for linesearch method
	Vector u(16);    u.Zero();
	Vector fSpr(13); fSpr.Zero();
	Vector kSpr(13); kSpr.Zero();
	Vector intEq(4); intEq.Zero();

	double r0 = 0.0;            // tolerance check for line-search
	double tolerance = 0.8;     // slack region tolerance set for line-search
    
	if (s0 != 0.0)
		r0 = fabs(s1/s0);

	if (r0 <= tolerance)
		return 1.0;   // Linsearch Not required residual decrease less than tolerance

	if (s1 == s0)
		return 1.0;   // Bisection would have divide by zero error if continued

	// set some variables
	double etaR;
	double eta = 1.0;
	double s = s1;
	double etaU = 1.0;
	double etaL = 0.0;
	double sU = s1;
	double sL = s0;
	double r = r0;
	double etaJ = 1.0;

	double minEta = 0.1;
	double maxEta = 10.0;
	int maxIter = 20;

	// bracket the region
	int count = 0;
	while ((sU*sL > 0.0) && (etaU < maxEta)) {
		count ++;
		etaU *= 2.0;
		etaR = etaU - etaJ;

		for (int i = 0; i < 12; i++) {
			u(i) = uExt(i) + duExt(i);
		}
		u(12) = uInt(0) - etaR*duInt(0);
		u(13) = uInt(1) - etaR*duInt(1);
		u(14) = uInt(2) - etaR*duInt(2);
		u(15) = uInt(3) - etaR*duInt(3);

		getMatResponse(u,fSpr,kSpr);

		intEq(0) = -fSpr(2) - fSpr(3) + fSpr(9) - fSpr(12)/elemHeight;
		intEq(1) = fSpr(1) - fSpr(5) - fSpr(7) + fSpr(12)/elemWidth;
		intEq(2) = -fSpr(4) - fSpr(8) + fSpr(10) + fSpr(12)/elemHeight;
		intEq(3) = fSpr(0) - fSpr(6) - fSpr(11) - fSpr(12)/elemWidth;

		sU = duInt^intEq;

		// check if we are happy with the solution
		r = fabs(sU/s0);
		if (r < tolerance)
			return etaU;

		etaJ = etaU;
	}

	if (sU*sL > 0.0)   // no bracketing could be done
		return 1.0;

	count = 0;
	while (r > tolerance && count < maxIter) {
		count ++;
		eta = (etaU + etaL)/2.0;

		if (r > r0) eta = 1.0;

		etaR = eta - etaJ;

		for (int i = 0; i < 12; i++) {
			u(i) = uExt(i) + duExt(i);
		}
		u(12) = uInt(0) - etaR*duInt(0);
		u(13) = uInt(1) - etaR*duInt(1);
		u(14) = uInt(2) - etaR*duInt(2);
		u(15) = uInt(3) - etaR*duInt(3);

		getMatResponse(u,fSpr,kSpr);

		intEq(0) = -fSpr(2) - fSpr(3) + fSpr(9) - fSpr(12)/elemHeight;
		intEq(1) = fSpr(1) - fSpr(5) - fSpr(7) + fSpr(12)/elemWidth;
		intEq(2) = -fSpr(4) - fSpr(8) + fSpr(10) + fSpr(12)/elemHeight;
		intEq(3) = fSpr(0) - fSpr(6) - fSpr(11) - fSpr(12)/elemWidth;

		s = duInt^intEq;

		// check if we are happy with the solution
		r = fabs(s/s0);

		// set variables for next iteration 
		etaJ = eta;
		
		if (s*sU < 0.0) {
			etaL = eta;
			sL = s;
		} else if (s*sU == 0.0) {
			count = maxIter;
		} else {
			etaU = eta;
			sU = s;
		}

		if (sL == sU)
			count = maxIter;

	}

	return eta;
}
    
const Matrix &
BeamColumnJoint3d::getDamp(void)
{
	// yet to do (stiffness being returned)
	return K;
}

const Matrix &
BeamColumnJoint3d::getMass(void)
{ 
	// yet to do (stiffness being returned)
	return K;
}

void 
BeamColumnJoint3d::zeroLoad(void)
{
  return;
}

int 
BeamColumnJoint3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return 0;
}

int 
BeamColumnJoint3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}


const Vector &
BeamColumnJoint3d::getResistingForceIncInertia()
{	
  //yet to do  (residual being returned)
	return R;
}

int
BeamColumnJoint3d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
BeamColumnJoint3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
BeamColumnJoint3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return -1; // yet to do
}

void
BeamColumnJoint3d::Print(OPS_Stream &s, int flag)
{
	s << "Element: " << this->getTag() << " Type: Beam Column Joint " << endln;
	for (int i = 0; i<4; i++)
	{
		s << "Node :" << connectedExternalNodes(i);
		s << "DOF :" << nodePtr[i]->getNumberDOF();
	}
	s << "\nResisting Forces :" <<this->getResistingForce();
	return;
}

Response*
BeamColumnJoint3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	// we will compare argv[0] to determine the type of response required
    
	if (strcmp(argv[0],"node1BarSlipL") == 0 || strcmp(argv[0],"node1BarslipL") ==0 || strcmp(argv[0],"Node1BarSlipL") == 0)
		return MaterialPtr[0]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node1BarSlipR") == 0 || strcmp(argv[0],"node1BarslipR") ==0 || strcmp(argv[0],"Node1BarSlipR") == 0)
		return MaterialPtr[1]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node1InterfaceShear") == 0 || strcmp(argv[0],"node1Interfaceshear") ==0 || strcmp(argv[0],"Node1InterfaceShear") ==0 )
		return MaterialPtr[2]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node2BarSlipB") == 0 || strcmp(argv[0],"node2BarslipB") ==0 || strcmp(argv[0],"Node2BarSlipB") == 0)
		return MaterialPtr[3]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node2BarSlipT") == 0 || strcmp(argv[0],"node2BarslipT") ==0 || strcmp(argv[0],"Node2BarSlipT") == 0)
		return MaterialPtr[4]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node2InterfaceShear") == 0 || strcmp(argv[0],"node2Interfaceshear") ==0 || strcmp(argv[0],"Node2InterfaceShear") ==0 )
		return MaterialPtr[5]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node3BarSlipL") == 0 || strcmp(argv[0],"node3BarslipL") ==0 || strcmp(argv[0],"Node3BarSlipL") == 0)
		return MaterialPtr[6]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node3BarSlipR") == 0 || strcmp(argv[0],"node3BarslipR") ==0 || strcmp(argv[0],"Node3BarSlipR") == 0)
		return MaterialPtr[7]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node3InterfaceShear") == 0 || strcmp(argv[0],"node3Interfaceshear") ==0 || strcmp(argv[0],"Node3InterfaceShear") ==0 )
		return MaterialPtr[8]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node4BarSlipB") == 0 || strcmp(argv[0],"node4BarslipB") ==0 || strcmp(argv[0],"Node4BarSlipB") == 0)
		return MaterialPtr[9]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node4BarSlipT") == 0 || strcmp(argv[0],"node4BarslipT") ==0 || strcmp(argv[0],"Node4BarSlipT") == 0)
		return MaterialPtr[10]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"node4InterfaceShear") == 0 || strcmp(argv[0],"node4Interfaceshear") ==0 || strcmp(argv[0],"Node4InterfaceShear") ==0 )
		return MaterialPtr[11]->setResponse(&argv[1], argc-1, output);

	else if (strcmp(argv[0],"shearpanel") == 0 || strcmp(argv[0],"shearPanel") ==0)
		return MaterialPtr[12]->setResponse(&argv[1], argc-1, output);
	

	else if (strcmp(argv[0],"externalDisplacement") == 0 || strcmp(argv[0],"externaldisplacement") == 0)
		return new ElementResponse(this,1,Vector(12));
    
	else if (strcmp(argv[0],"internalDisplacement") == 0 || strcmp(argv[0],"internaldisplacement") == 0)
		return new ElementResponse(this,2,Vector(4));

	else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"Deformation") == 0)
		return new ElementResponse(this,3,Vector(4));

	else
		return 0;
}

int 
BeamColumnJoint3d::getResponse(int responseID, Information &eleInfo)
{
	static Vector delta(13);
	static Vector def(4);
	static Vector U(16);
	static Vector Utemp(12);
	double bsFa, bsFb, bsFc, bsFd;
	double bsFac, bsFbd, isFac, isFbd; 

	switch (responseID) {
	case 1:       
		if(eleInfo.theVector!=0)
		{
			(*(eleInfo.theVector))(0) =  UeprCommit(0);
			(*(eleInfo.theVector))(1) =  UeprCommit(1);
			(*(eleInfo.theVector))(2) =  UeprCommit(2);
			(*(eleInfo.theVector))(3) =  UeprCommit(3);
			(*(eleInfo.theVector))(4) =  UeprCommit(4);
			(*(eleInfo.theVector))(5) =  UeprCommit(5);
			(*(eleInfo.theVector))(6) =  UeprCommit(6);
			(*(eleInfo.theVector))(7) =  UeprCommit(7);
			(*(eleInfo.theVector))(8) =  UeprCommit(8);
			(*(eleInfo.theVector))(9) =  UeprCommit(9);
			(*(eleInfo.theVector))(10) = UeprCommit(10);
			(*(eleInfo.theVector))(11) = UeprCommit(11);
			(*(eleInfo.theVector))(12) = UeprCommit(12);
			(*(eleInfo.theVector))(13) = UeprCommit(13);
			(*(eleInfo.theVector))(14) = UeprCommit(14);
			(*(eleInfo.theVector))(15) = UeprCommit(15);
			(*(eleInfo.theVector))(16) = UeprCommit(16);
			(*(eleInfo.theVector))(17) = UeprCommit(17);
			(*(eleInfo.theVector))(18) = UeprCommit(18);
			(*(eleInfo.theVector))(19) = UeprCommit(19);
			(*(eleInfo.theVector))(20) = UeprCommit(20);
			(*(eleInfo.theVector))(21) = UeprCommit(21);
			(*(eleInfo.theVector))(22) = UeprCommit(22);
			(*(eleInfo.theVector))(23) = UeprCommit(23);
		}
		return 0;

	case 2:
		if (eleInfo.theVector !=0) {
			(*(eleInfo.theVector))(0) = UeprIntCommit(0);
			(*(eleInfo.theVector))(1) = UeprIntCommit(1);
			(*(eleInfo.theVector))(2) = UeprIntCommit(2);
			(*(eleInfo.theVector))(3) = UeprIntCommit(3);
		}
		return 0;

	case 3:

		// modified 01.04.03  -------- determine which plane the joint lies
		Utemp.addMatrixVector(0.0,Transf,UeprCommit,1.0);
		U.Assemble(Utemp,0,1.0);
		U.Assemble(UeprIntCommit,12,1.0);

		delta.addMatrixVector(0.0,BCJoint,U,1.0);
		
		bsFa = fabs(delta(0) - delta(1))/elemWidth;
		bsFc = fabs(delta(7) - delta(6))/elemWidth;
		bsFac = bsFa + bsFc;
		bsFb = fabs(delta(4) - delta(3))/elemHeight;
		bsFd = fabs(delta(10) - delta(9))/elemHeight;
		bsFbd = bsFb + bsFd;

		def(0) = bsFac + bsFbd; // contribution of bar slip deformation
		
		
		isFac = (delta(2) + delta(8))/elemHeight;
		isFbd = (delta(5) + delta(11))/elemWidth;

		def(1) = isFac + isFbd;  // contribution due to interface shear spring

		def(2) = delta(12);  // contribution due to shear panel

		def(3) = def(0) + def(1) + def(2);  // total joint deformation


		return eleInfo.setVector(def);
	default:
		return -1;
	}
}

int
BeamColumnJoint3d::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
    
int
BeamColumnJoint3d::updateParameter (int parameterID, Information &info)
{
  return -1;
}
