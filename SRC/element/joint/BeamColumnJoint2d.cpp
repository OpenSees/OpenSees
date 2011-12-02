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
// $Date: 2004-06-07 23:19:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint2d.cpp,v $
                                                                        
// Written: NM (nmitra@u.washington.edu)
// Created: April 2002
//
// Description: This file contains the class implementation for beam-column joint.
// This element is a 4 noded 12 dof (3 dof at each node) finite area super-element introduced by
// Prof. Lowes and Arash. The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint.


#include <BeamColumnJoint2d.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>

#include <UniaxialMaterial.h>
#include <string.h>
#include <ElementResponse.h>



// full constructors:
BeamColumnJoint2d::BeamColumnJoint2d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
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
  Element(tag,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(4),
  nodeDbTag(0), dofDbTag(0), elemWidth(0), elemHeight(0),
  Uecommit(12), UeIntcommit(4), UeprCommit(12), UeprIntCommit(4), 
  BCJoint(13,16), dg_df(4,13), dDef_du(13,4), K(12,12), R(12)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 4)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJoint " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    connectedExternalNodes(2) = Nd3 ;
    connectedExternalNodes(3) = Nd4 ;

// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
  MaterialPtr[1] = theMat2.getCopy();
  if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2"<< endln;}
  MaterialPtr[2] = theMat3.getCopy();
  if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3"<< endln;}
  MaterialPtr[3] = theMat4.getCopy();
  if (!MaterialPtr[3]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 4"<< endln;}
  MaterialPtr[4] = theMat5.getCopy();
  if (!MaterialPtr[4]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 5"<< endln;}
  MaterialPtr[5] = theMat6.getCopy();
  if (!MaterialPtr[5]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 6"<< endln;}
  MaterialPtr[6] = theMat7.getCopy();
  if (!MaterialPtr[6]){
	    opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 7"<< endln;}
  MaterialPtr[7] = theMat8.getCopy();
  if (!MaterialPtr[7]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 8"<< endln;}
  MaterialPtr[8] = theMat9.getCopy();
  if (!MaterialPtr[8]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 9"<< endln;}
  MaterialPtr[9] = theMat10.getCopy();
  if (!MaterialPtr[9]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 10"<< endln;}
  MaterialPtr[10] = theMat11.getCopy();
  if (!MaterialPtr[10]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 11"<< endln;}
  MaterialPtr[11] = theMat12.getCopy();
  if (!MaterialPtr[11]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 12"<< endln;}
  MaterialPtr[12] = theMat13.getCopy();
  if (!MaterialPtr[12]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 13"<< endln;}	


}

/********************************************************************************************/

// default constructor:
BeamColumnJoint2d::BeamColumnJoint2d():
  Element(0,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(4),
  nodeDbTag(0), dofDbTag(0), elemWidth(0), elemHeight(0),
  Uecommit(12), UeIntcommit(4), UeprCommit(12), UeprIntCommit(4),
  BCJoint(13,16), dg_df(4,13), dDef_du(13,4), K(12,12), R(12)
{
   // does nothing (invoked by FEM_ObjectBroker)
}
/********************************************************************************************/

//  destructor:
BeamColumnJoint2d::~BeamColumnJoint2d()
{
	for (int i =0; i<13; i++)
	{
	    if (MaterialPtr[i] != 0)
		delete MaterialPtr[i];
	}

  for (int ia = 0 ;  ia < 4; ia++ )
  {

    nodePtr[ia] = 0 ;

  }

}
/********************************************************************************************/

// public methods

int
BeamColumnJoint2d::getNumExternalNodes(void) const
{
    return 4;
}
/********************************************************************************************/

const ID &
BeamColumnJoint2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}
/********************************************************************************************/
Node **BeamColumnJoint2d::getNodePtrs(void)
{
	return nodePtr;
}

/********************************************************************************************/

int
BeamColumnJoint2d::getNumDOF(void) 
{
    return 12;
}
/********************************************************************************************/

void
BeamColumnJoint2d::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	opserr << "ERROR : BeamColumnJoint::setDomain -- Domain is null" << endln;
	
  //node pointers set
    for (int i = 0; i < 4; i++ ) {
		nodePtr[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR : BeamColumnJoint::setDomain -- node pointer is null"<< endln;
			return; // donot go any further - otherwise segmentation fault
		}
	}

    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    this->revertToStart();

	// ensure connected nodes have correct dof's
	int dofNd1 = nodePtr[0]->getNumberDOF();
	int dofNd2 = nodePtr[1]->getNumberDOF();
	int dofNd3 = nodePtr[2]->getNumberDOF();
	int dofNd4 = nodePtr[3]->getNumberDOF();

	if ((dofNd1 != 3) || (dofNd2 != 3) || (dofNd3 != 3) || (dofNd4 != 3)) 
	{
			opserr << "ERROR : BeamColumnJoint::setDomain -- number of DOF associated with the node incorrect"<< endln;
			return; // donot go any further - otherwise segmentation fault
	}


    // obtain the nodal coordinates    
	const Vector &end1Crd = nodePtr[0]->getCrds();
	const Vector &end2Crd = nodePtr[1]->getCrds();
	const Vector &end3Crd = nodePtr[2]->getCrds();
	const Vector &end4Crd = nodePtr[3]->getCrds();

	nodeCrd[0][0] = end1Crd(0);
	nodeCrd[0][1] = end1Crd(1);
	nodeCrd[1][0] = end2Crd(0);
	nodeCrd[1][1] = end2Crd(1);
	nodeCrd[2][0] = end3Crd(0);
	nodeCrd[2][1] = end3Crd(1);
	nodeCrd[3][0] = end4Crd(0);
	nodeCrd[3][1] = end4Crd(1);
	

	// set the height and width of the element and perform check   
	elemHeight = fabs(nodeCrd[2][1] - nodeCrd[0][1]);
	elemWidth = fabs(nodeCrd[3][0] - nodeCrd[1][0]);

	if ((elemHeight == 0) || (elemWidth == 0))
	{
			opserr << "ERROR : BeamColumnJoint::setDomain -- length or width not correct, division by zero occurs"<< endln;
			return; // donot go any further - otherwise segmentation fault
	}


	getBCJoint();
	getdg_df();
	getdDef_du();
	
}   
/********************************************************************************************/

int
BeamColumnJoint2d::commitState()
{
	// store commited external nodal displacements
	Uecommit = UeprCommit;
	// store commited internal nodal displacements
	UeIntcommit = UeprIntCommit;

	// store material history data.
	int mcs;
		for (int j=0; j<13; j++)
		{
			mcs = MaterialPtr[j]->commitState();
		}
    
	return 0;
}
/********************************************************************************************/

int
BeamColumnJoint2d::revertToLastCommit()
{

	int mcs;
	for (int j=0; j<13; j++)
	{
		mcs = MaterialPtr[j]->revertToLastCommit();
	}
	UeprCommit = Uecommit;
	UeprIntCommit = UeIntcommit;   
	
	this->update();

  return 0;
}
/********************************************************************************************/

int
BeamColumnJoint2d::revertToStart()
{
	int mcs;
	for (int j=0; j<13; j++)
	{
		mcs = MaterialPtr[j]->revertToStart();
	}
	
	Uecommit.Zero();
	UeIntcommit.Zero();
	UeprCommit.Zero();
	UeprIntCommit.Zero();

	return 0;
}
/********************************************************************************************/

int
BeamColumnJoint2d::update()
{	

		Vector Ue(16);
		Ue.Zero();

	// determine commited displacements given trial displacements
		this->getGlobalDispls(Ue);
 
		// update displacements for the external nodes
		UeprCommit.Extract(Ue,0,1.0); 
		// update displacement for the internal nodes
		UeprIntCommit.Extract(Ue,12,1.0);

		return 0;
}
/********************************************************************************************/

const Matrix &
BeamColumnJoint2d::getTangentStiff(void)
{
	return K;
}

/********************************************************************************************/

const Matrix &
BeamColumnJoint2d::getInitialStiff(void)
{
	return getTangentStiff();
}

/********************************************************************************************/

const Vector &
BeamColumnJoint2d::getResistingForce()
{
  return R;
}

/*************************************************************************************************************************/
void BeamColumnJoint2d::getGlobalDispls(Vector &dg) 
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

	Vector uExtOld(12);
	Vector uExt(12);
	Vector duExt(12);
	Vector uIntOld(4);
	Vector uInt(4);
	Vector duInt(4);
	Vector duIntTemp(4);
	Vector intEq(4);
	Vector intEqLast(4);
	Vector Uepr(12);
	Vector UeprInt(4);
	Vector Ut(12);
	

    Vector disp1 = nodePtr[0]->getTrialDisp(); 
    Vector disp2 = nodePtr[1]->getTrialDisp();
    Vector disp3 = nodePtr[2]->getTrialDisp();
    Vector disp4 = nodePtr[3]->getTrialDisp();

	for (int i = 0; i < 3; i++)
    {
      Ut(i)     = disp1(i);
      Ut(i+3)   = disp2(i);
      Ut(i+6)   = disp3(i);
      Ut(i+9)   = disp4(i);

    }


	Uepr = Uecommit;   

	//UeprInt = UeIntcommit;
	UeprInt = UeprIntCommit;  

	uExtOld = Uepr;
	duExt = Ut -Uepr;

	uExt = uExtOld;

	uIntOld = UeprInt;
	uInt = uIntOld;
	duInt.Zero();


	double tol = 1e-12;
	double tolIntEq = tol;
	double toluInt = (tol>tol*uInt.Norm())? tol:tol*uInt.Norm();
	double tolIntEqdU = tol;
	double ctolIntEqdU = tol;
	double ctolIntEq = tol;
	double normDuInt = toluInt;
	double normIntEq = tolIntEq;
	double normIntEqdU = tolIntEqdU;
	
	intEq.Zero();
	intEqLast.Zero();
    
	Vector u(16);
	u.Zero();

	double engrLast;
	double engr;

	Vector fSpring(13);
	Vector kSpring(13);
    Matrix dintEq_du(4,4);
    Matrix df_dDef(13,13);
    Matrix tempintEq_du (4,13);

	
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
				u[ic] = uExt[ic] + duExt[ic];
			}
			u[12] = uInt[0];
			u[13] = uInt[1];
			u[14] = uInt[2];
			u[15] = uInt[3];


		getMatResponse(u,fSpring,kSpring);
		

		// performs internal equilibrium 

		intEq[0] = -fSpring[2]-fSpring[3]+fSpring[9]-fSpring[12]/elemHeight; 
		intEq[1] = fSpring[1]-fSpring[5]-fSpring[7]+fSpring[12]/elemWidth; 
		intEq[2] = -fSpring[4]-fSpring[8]+fSpring[10]+fSpring[12]/elemHeight; 
		intEq[3] = fSpring[0]-fSpring[6]-fSpring[11]-fSpring[12]/elemWidth; 


		df_dDef.Zero();               // the diagonal stiffness matrix for the element

		matDiag(kSpring, df_dDef);


		//////////////////////// dintEq_du = dg_df*df_dDef*dDef_du
		tempintEq_du.Zero();
		dintEq_du.Zero();

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
		ctolIntEq = (tolIntEqdU*dLoadStep > tol) ? tolIntEqdU*dLoadStep:tol;
		ctolIntEqdU = (tolIntEq*dLoadStep > tol) ? tolIntEq*dLoadStep:tol;

		// check for convergence starts  

      if ((normIntEq < ctolIntEq) || ((normIntEqdU < ctolIntEqdU) && (count >1)) || (normDuInt < toluInt) || (dLoadStep < 1e-3))
		{
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
				opserr << "WARNING : BeamColumnJoint2D::getGlobalDispls() - convergence problem in state determination" << endln;

				for (int je = 0; je<4 ; je++)
				{
					uInt[je] = uIntOld[je];
					duInt[je] = 0;
				}

				for (int jee = 0; jee<12; jee++)
				{
					duExt[jee] = 0.1*duExt[jee];
				}

				dLoadStep = dLoadStep*0.1;
			}
		}
		else if (loadStep < 1.0)
		{
			maxCount = 10;
			incCount ++;
			normDuInt = toluInt;
			if ((incCount < 10) || dtConverge)
			{
				for (int jf = 0; jf<12 ; jf++)
				{
					uExt[jf] = uExt[jf] + duExt[jf];
				}
				if (loadStep + dLoadStep > 1.0)
				{
					for (int jg = 0; jg<12 ; jg++)
					{
						duExt[jg] = duExt[jg]*(1.0 - loadStep)/dLoadStep;
					}
					dLoadStep = 1.0 - loadStep;
					incCount = 9;
				}
			}
			else
			{
				incCount = 0;
				for (int jh = 0; jh<12 ; jh++)
				{
					uExt[jh] = uExt[jh] + duExt[jh];
				}
				dLoadStep = dLoadStep*10;
				if (loadStep + dLoadStep > 1.0)
				{
					for (int ji = 0; ji<12 ; ji++)
					{
						uExt[ji] = uExt[ji] + duExt[ji]*(1.0 - loadStep)/dLoadStep;
					}
					dLoadStep = 1.0 - loadStep;
					incCount = 9;
				}
			}
		}
	}


	// determination of stiffness matrix and the residual force vector for the element

	formR(fSpring);

	formK(kSpring);


	
	// commmited external and internal displacement update
	for (int ig = 0; ig < 13; ig++ ) {
		if (ig<12)
		{
			dg[ig] = Ut[ig];
		}
	}

	dg[12] = uInt[0];
	dg[13] = uInt[1];
	dg[14] = uInt[2];
	dg[15] = uInt[3];

}

/*************************************************************************************************************************/
void BeamColumnJoint2d::getMatResponse(Vector U, Vector &fS, Vector &kS)
{
	// obtains the material response from the material class
	Vector defSpring(13);
	defSpring.Zero();
	fS.Zero();
	kS.Zero();


	defSpring.addMatrixVector(0.0, BCJoint, U, 1.0);

	for (int j=0; j<13; j++)
	{
		MaterialPtr[j]->setTrialStrain(defSpring[j]);
		kS[j] = MaterialPtr[j]->getTangent();
		fS[j] = MaterialPtr[j]->getStress();
	}

}

/*************************************************************************************************************************/
void BeamColumnJoint2d::getdDef_du()
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
/*************************************************************************************************************************/
void BeamColumnJoint2d::matDiag(Vector k,Matrix &df)
{
	// takes in a vector and converts it to a diagonal matrix (could have been placed as a method in matrix class)
	for (int ja=0; ja<13; ja++)
	{
		df(ja,ja) = k(ja);
	}

}
/*************************************************************************************************************************/
void BeamColumnJoint2d::formR(Vector f)
{
	// develops the element residual force vector
	Vector rForceTemp(16);
	rForceTemp.Zero();
	
	rForceTemp.addMatrixTransposeVector(0.0,BCJoint,f,1.0);   // rForceTemp = BCJoint'*f
	R.Extract(rForceTemp,0,1.0);

}
/*************************************************************************************************************************/
void BeamColumnJoint2d::formK(Vector k)
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
	K2.addMatrixProduct(0.0,K2Temp,kRFT1,1.0);


	kRF.addMatrix(1.0,K2,-1.0);

	
	K = kRF;    // K = kRF - kRFT3*kRSTinv*kRFT1


}
/*************************************************************************************************************************/

void BeamColumnJoint2d::getdg_df()
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
/*************************************************************************************************************************/

void BeamColumnJoint2d::getBCJoint() 
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

}

/*************************************************************************************************************************/
double BeamColumnJoint2d::getStepSize(double G0,double G,Vector uExt,Vector duExt,Vector uInt,Vector duInt,double tol)
{
	// finds out the factor to be used for linesearch method
	Vector u(16);
	Vector fSprng(13);
	Vector kSprng(13);
	Vector intEq(4);

	int count =0;
	int smax = 10;
	int sb = 0;
	int sa = 1;
	double Gb = G0;
	double Ga = G;


	while ((Ga*Gb > 0) && (sa<smax))
	{
		sb = sa;
		sa *= 2;
		Gb = Ga;

		for (int i = 0; i < 12; i++ ) {
			u[i] = uExt[i] + duExt[i];
			}
		u[12] = uInt[0] - sa*duInt[0];
		u[13] = uInt[1] - sa*duInt[1];
		u[14] = uInt[2] - sa*duInt[2];
		u[15] = uInt[3] - sa*duInt[3];

		getMatResponse(u,fSprng,kSprng);
		
		// performs internal equilibrium 

		intEq[0] = -fSprng[2]-fSprng[3]+fSprng[9]-fSprng[12]/elemHeight; 
		intEq[1] = fSprng[1]-fSprng[5]-fSprng[7]+fSprng[12]/elemWidth; 
		intEq[2] = -fSprng[4]-fSprng[8]+fSprng[10]+fSprng[12]/elemHeight; 
		intEq[3] = fSprng[0]-fSprng[6]-fSprng[11]-fSprng[12]/elemWidth; 

		Ga = 0.0;
		for (int ia=0; ia<4; ia++)
		{
			Ga += duInt[ia]*intEq[ia];
		}
	}

	double step = static_cast<double>(sa);
	G = Ga;

	while ((Ga*Gb < 0) && ((fabs(G) > tol*fabs(G0)) || (abs(sb-sa) > tol*0.5*(sa+sb))) && count<20)
	{
		count++;
		step = sa - Ga*(sa-sb)/(Ga-Gb);

		for (int ib = 0; ib < 12; ib++ ) {
			u[ib] = uExt[ib] + duExt[ib]; }

		u[12] = uInt[0] - sa*duInt[0];
		u[13] = uInt[1] - sa*duInt[1];
		u[14] = uInt[2] - sa*duInt[2];
		u[15] = uInt[3] - sa*duInt[3];
		
	
		getMatResponse(u,fSprng,kSprng);
		
		// performs internal equilibrium

		intEq[0] = -fSprng[2]-fSprng[3]+fSprng[9]-fSprng[12]/elemHeight; 
		intEq[1] = fSprng[1]-fSprng[5]-fSprng[7]+fSprng[12]/elemWidth; 
		intEq[2] = -fSprng[4]-fSprng[8]+fSprng[10]+fSprng[12]/elemHeight; 
		intEq[3] = fSprng[0]-fSprng[6]-fSprng[11]-fSprng[12]/elemWidth; 

		G = 0.0;
		for (int ic=0; ic<4; ic++)
		{
			G += duInt[ic]*intEq[ic];
		}
		if (G*Ga > 0)
		{
			Gb = 0.5*Gb;
		}
		else
		{
			sb = sa;
			Gb = Ga;
		}
		sa = static_cast<int>(step);
		Ga = G;
	}

	double stepSize = (step<1.0)? step:1.0;
	return stepSize;

}
/********************************************************************************************/
    
const Matrix &
BeamColumnJoint2d::getDamp(void)
{
	//not applicable (stiffness being returned)
	return K;
}
/********************************************************************************************/


const Matrix &
BeamColumnJoint2d::getMass(void)
{ 
	//not applicable  (stiffness being returned)
	return K;
}
/********************************************************************************************/

void 
BeamColumnJoint2d::zeroLoad(void)
{
	//not applicable  
	return;
}
/********************************************************************************************/

int 
BeamColumnJoint2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	//not applicable
	return 0;
}
/********************************************************************************************/

int 
BeamColumnJoint2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	//not applicable
	return 0;
}

/********************************************************************************************/

const Vector &
BeamColumnJoint2d::getResistingForceIncInertia()
{	
  //not applicable (residual being returned)
	return R;
}
/********************************************************************************************/

int
BeamColumnJoint2d::sendSelf(int commitTag, Channel &theChannel)
{
	// yet to do.
	return -1;
}
/********************************************************************************************/

int
BeamColumnJoint2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// yet to do.
	return -1;
}
/********************************************************************************************/

int
BeamColumnJoint2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	const Vector &node1Crd = nodePtr[0]->getCrds();
	const Vector &node2Crd = nodePtr[1]->getCrds();	
	const Vector &node3Crd = nodePtr[2]->getCrds();
	const Vector &node4Crd = nodePtr[3]->getCrds();

	const Vector &node1Disp = nodePtr[0]->getDisp();
	const Vector &node2Disp = nodePtr[1]->getDisp();    
	const Vector &node3Disp = nodePtr[2]->getDisp();
	const Vector &node4Disp = nodePtr[3]->getDisp();  

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
	
	// calculate four corners of the element
	Vector vb(3);
	Vector vc(3);

	vb = v1 - v3;
	vc = v2 - v4;

	v1 = v1 - 0.5 * vc;
	v2 = v1 + vc;
	v3 = v4 + 0.5 * vb;
	v4 = v4 + vb;
	
	int dummy;
	dummy = theViewer.drawLine(v1, v2, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v4, 1.0, 1.0);
	dummy = theViewer.drawLine(v4, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v3, v1, 1.0, 1.0);

	return 0;

}
/********************************************************************************************/

void
BeamColumnJoint2d::Print(OPS_Stream &s, int flag)
{
	s << "Element: " << this->getTag() << " Type: Beam Column Joint " << endln;
	for (int i = 0; i<4; i++)
	{
		s << "Node :" << connectedExternalNodes(i);
		s << "DOF :" << nodePtr[i]->getNumberDOF();
	}
	return;
}
/********************************************************************************************/

Response*
BeamColumnJoint2d::setResponse(const char **argv, int argc, Information &eleInfo)
{
	// we will compare argv[0] to determine the type of response required
    
	if (strcmp(argv[0],"columnbar1") == 0 || strcmp(argv[0],"columnBar1") ==0 )
		return MaterialPtr[0]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"columnbar2") == 0 || strcmp(argv[0],"columnBar2") ==0 )
		return MaterialPtr[1]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"columnbar3") == 0 || strcmp(argv[0],"columnBar3") ==0 )
		return MaterialPtr[2]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beambottombar1") == 0 || strcmp(argv[0],"beambottomBar1") ==0 || strcmp(argv[0],"beamBottomBar1") == 0)
		return MaterialPtr[3]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beamtopbar1") == 0 || strcmp(argv[0],"beamtopBar1") ==0 || strcmp(argv[0],"beamTopBar1") == 0)
		return MaterialPtr[4]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beammiddlebar1") == 0 || strcmp(argv[0],"beammiddleBar1") ==0 || strcmp(argv[0],"beamMiddleBar1") == 0)
		return MaterialPtr[5]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"columnbar4") == 0 || strcmp(argv[0],"columnBar4") ==0 )
		return MaterialPtr[6]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"columnbar5") == 0 || strcmp(argv[0],"columnBar5") ==0 )
		return MaterialPtr[7]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"columnbar6") == 0 || strcmp(argv[0],"columnBar6") ==0 )
		return MaterialPtr[8]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beambottombar2") == 0 || strcmp(argv[0],"beambottomBar2") ==0 || strcmp(argv[0],"beamBottomBar2") == 0)
		return MaterialPtr[9]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beamtopbar2") == 0 || strcmp(argv[0],"beamtopBar2") ==0 || strcmp(argv[0],"beamTopBar2") == 0)
		return MaterialPtr[10]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"beammiddlebar2") == 0 || strcmp(argv[0],"beammiddleBar2") ==0 || strcmp(argv[0],"beamMiddleBar2") == 0)
		return MaterialPtr[11]->setResponse(&argv[1], argc-1, eleInfo);

	else if (strcmp(argv[0],"shearpanel") == 0 || strcmp(argv[0],"shearPanel") ==0)
		return MaterialPtr[12]->setResponse(&argv[1], argc-1, eleInfo);
	

	else if (strcmp(argv[0],"externalDisplacement") == 0 || strcmp(argv[0],"externaldisplacement") == 0)
		return new ElementResponse(this,1,Vector(12));
    
	else if (strcmp(argv[0],"internalDisplacement") == 0 || strcmp(argv[0],"internaldisplacement") == 0)
		return new ElementResponse(this,2,Vector(4));

	else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"Deformation") == 0)
		return new ElementResponse(this,3,Vector(4));

	else
		return 0;
}
/********************************************************************************************/

int 
BeamColumnJoint2d::getResponse(int responseID, Information &eleInfo)
{
	static Vector delta(13);
	static Vector def(4);
	static Vector U(16);
	int tr,ty;
	double bsFa, bsFb, bsFc, bsFd;
	double bsFac, bsFbd, isFac, isFbd; 

	switch (responseID) {
	case 1:       
		if(eleInfo.theVector!=0)
		{
			(*(eleInfo.theVector))(0) = UeprCommit(0);
			(*(eleInfo.theVector))(1) = UeprCommit(1);
			(*(eleInfo.theVector))(2) = UeprCommit(2);
			(*(eleInfo.theVector))(3) = UeprCommit(3);
			(*(eleInfo.theVector))(4) = UeprCommit(4);
			(*(eleInfo.theVector))(5) = UeprCommit(5);
			(*(eleInfo.theVector))(6) = UeprCommit(6);
			(*(eleInfo.theVector))(7) = UeprCommit(7);
			(*(eleInfo.theVector))(8) = UeprCommit(8);
			(*(eleInfo.theVector))(9) = UeprCommit(9);
			(*(eleInfo.theVector))(10) = UeprCommit(10);
			(*(eleInfo.theVector))(11) = UeprCommit(11);
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

		for (ty =0; ty <12; ty++)
		{
			U(ty) = UeprCommit(ty);
		}
		for (tr = 0; tr <4 ; tr++)
		{
			U(12+tr) = UeprIntCommit(tr);
		}

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
/********************************************************************************************/

int
BeamColumnJoint2d::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
/********************************************************************************************/
    
int
BeamColumnJoint2d::updateParameter (int parameterID, Information &info)
{
  return -1;
}
