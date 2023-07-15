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
                                                                        
// $Revision: 1.6 $
// $Date: 2007/07/27 19:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/LehighJoint2d.cpp,v $
                                                                        
// Written: CY Seo
// Created August 2008
//
// // REFERENCES (added by CVMiculas)
// // C.Y. Seo, Y.C. Lin, R. Sause & J.M. Ricles (2009). Development of analytical models for 0.6 scale self-centering MRF with beam web friction devices. In: 6th International Conference for Steel Structures in Seismic Area (STESSA), Philadelphia. CRC Press, pp. 849-854.
// // The article from above is part of a conference proceedings book:
// // Mazzolani, F., Ricles, J.M., & Sause, R. (Eds.). (2009). Behaviour of Steel Structures in Seismic Areas: STESSA 2009 (1st ed.). CRC Press. https://doi.org/10.1201/9780203861592
// // Karavasilis, Theodore & Seo, Choungyeol & Ricles, James. (2008). HybridFEM: A PROGRAM FOR DYNAMIC TIME HISTORY ANALYSIS OF 2D INELASTIC FRAMED STRUCTURES AND REAL-TIME HYBRID SIMULATION HybridFEM Version 4.2.4 User's Manual.

#include <LehighJoint2d.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>
#include <math.h>

#include <UniaxialMaterial.h>
#include <string.h>
#include <ElementResponse.h>
#include <elementAPI.h>

void *OPS_LehighJoint2d()
{
  Domain* theDomain = OPS_GetDomain();
  if (theDomain == 0) return 0;

  // check no of arguments
  int numData = OPS_GetNumRemainingInputArgs();
  if (numData < 14) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element LehighJoint eleTag? node1? node2? node3? node4? matTag1? matTag2? matTag3? ";
    opserr << "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? \n";
    
    return 0;
  }

  // 1 ele tag, 4 node tags, 9 material tags
  int idata[14];
  int num = 14;
  if (OPS_GetIntInput(&num, idata) < 0) {
    opserr<<"WARNING: invalid integer data\n";
    return 0;
  }

  int eleTag = idata[0];
  int ndI = idata[1];
  int ndJ = idata[2];
  int ndK = idata[3];
  int ndL = idata[4];

  UniaxialMaterial *theMats[9];
  for (int i = 0; i < 9; i++) {
    theMats[i] = OPS_getUniaxialMaterial(idata[i+5]);
    if (theMats[i] == 0) {
      opserr << "WARNING: material not found\n";
      opserr << "Material: " << idata[i+5];
      opserr << "\nLehighJoint2d element: " << eleTag << endln;
      return 0;
    }
  }

  LehighJoint2d *theEle = 0;
  theEle = new LehighJoint2d(eleTag, ndI, ndJ, ndK, ndL,
			     *theMats[0], *theMats[1], *theMats[2],
			     *theMats[3], *theMats[4], *theMats[5],
			     *theMats[6], *theMats[7], *theMats[8]);
  return theEle;

}

// full constructors:
LehighJoint2d::LehighJoint2d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3,
				     UniaxialMaterial& theMat4,
				     UniaxialMaterial& theMat5,
				     UniaxialMaterial& theMat6,
				     UniaxialMaterial& theMat7,
				     UniaxialMaterial& theMat8,
				     UniaxialMaterial& theMat9):
  Element(tag,ELE_TAG_LehighJoint2d), connectedExternalNodes(4),
  nodeDbTag(0), dofDbTag(0), elemWidth(0.0), elemHeight(0.0), numDOF(12), numBasicDOF(9),
  vs(9), vt(9), avp(9,12), apq(12,12),  K(12,12), R(12)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 4)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJoint " << tag << "failed to create an ID of size 4" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    connectedExternalNodes(2) = Nd3 ;
    connectedExternalNodes(3) = Nd4 ;

	MaterialPtr = new UniaxialMaterial*[numBasicDOF];

	for (int x = 0; x <numBasicDOF; x++)
	{	MaterialPtr[x] = 0; }

	vs.Zero();
	vt.Zero();
	
	K.Zero();
	R.Zero();

	nodePtr[0] = 0;
	nodePtr[1] = 0;
	nodePtr[2] = 0;
	nodePtr[3] = 0;

	// transformation matrix
	avp.Zero();
	apq.Zero();

	
// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  MaterialPtr[1] = theMat2.getCopy();
  MaterialPtr[2] = theMat3.getCopy();
  MaterialPtr[3] = theMat4.getCopy();
  MaterialPtr[4] = theMat5.getCopy();
  MaterialPtr[5] = theMat6.getCopy();
  MaterialPtr[6] = theMat7.getCopy();
  MaterialPtr[7] = theMat8.getCopy();
  MaterialPtr[8] = theMat9.getCopy();

  for (int i=0; i<numBasicDOF; i++)
  {
      if (!MaterialPtr[i]){
		 opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material " 
			 << (i+1) << endln;
	  }
  }

 }



// default constructor:
LehighJoint2d::LehighJoint2d():
  Element(0,ELE_TAG_LehighJoint2d), connectedExternalNodes(4),
  nodeDbTag(0), dofDbTag(0), elemWidth(0.0), elemHeight(0.0), numDOF(12), numBasicDOF(9),
  vs(9), vt(9), avp(9,12), apq(12,12), K(12,12), R(12)
{
	nodePtr[0] = 0;
	nodePtr[1] = 0;
	nodePtr[2] = 0;
	nodePtr[3] = 0;	

	for (int x = 0; x <numBasicDOF; x++)
	{	MaterialPtr[x] = 0; }
   // does nothing (invoked by FEM_ObjectBroker)
}

//  destructor:
LehighJoint2d::~LehighJoint2d()
{
	for (int i =0; i<numBasicDOF; i++)
	{
	    if (MaterialPtr[i] != 0)
		delete MaterialPtr[i];
	}

	if (MaterialPtr)
		 delete [] MaterialPtr;
}

// public methods
int
LehighJoint2d::getNumExternalNodes(void) const
{
    return 4;
}

const ID &
LehighJoint2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **LehighJoint2d::getNodePtrs(void)
{
	return nodePtr;
}

int
LehighJoint2d::getNumDOF(void) 
{
    return numDOF;
}

void
LehighJoint2d::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	opserr << "ERROR : BeamColumnJoint::setDomain -- Domain is null" << endln;
	
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePtr[0] = 0;
	nodePtr[1] = 0;
	nodePtr[2] = 0;
	nodePtr[3] = 0;
    }

	//node pointers set
    for (int i = 0; i < 4; i++ ) {
		nodePtr[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR : BeamColumnJoint::setDomain -- node pointer is null"<< endln;
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

	if ((dofNd1 != 3) || (dofNd2 != 3) || (dofNd3 != 3) || (dofNd4 != 3)) 
	{
			opserr << "ERROR : BeamColumnJoint::setDomain -- number of DOF associated with the node incorrect"<< endln;
			exit(-1); // donot go any further - otherwise segmentation fault
	}


    // obtain the nodal coordinates    
	const Vector &end1Crd = nodePtr[0]->getCrds();
	const Vector &end2Crd = nodePtr[1]->getCrds();
	const Vector &end3Crd = nodePtr[2]->getCrds();
	const Vector &end4Crd = nodePtr[3]->getCrds();

	Vector Node1(end1Crd);
	Vector Node2(end2Crd);
	Vector Node3(end3Crd);
	Vector Node4(end4Crd);

	// set the height and width of the element and perform check   
	Node3 = Node3 - Node1;
	Node2 = Node2 - Node4;

	elemHeight = fabs(Node2.Norm());
	elemWidth  = fabs(Node3.Norm());
	
	if ((elemHeight <= 1e-12) || (elemWidth <= 1e-12))
	{
		opserr << "ERROR : BeamColumnJoint::setDomain -- length or width not correct, division by zero occurs"<< endln;
		exit(-1); // donot go any further - otherwise segmentation fault
	}

    // determine the length, cosines and fill the transformation
	double dx = end3Crd(0) - end1Crd(0);
	double dy = end3Crd(1) - end1Crd(1);

	double L = sqrt(dx*dx + dy*dy);
	
	// fill transformation coefficients, global coordinate to local coordinate
	apq.Zero();

    apq(0,0) =  dx/L;
    apq(0,1) =  dy/L;
    apq(1,0) = -apq(0,1);
    apq(1,1) =  apq(0,0);
    apq(2,2) =  1.0;
    apq(3,3) =  apq(0,0) ;
    apq(3,4) =  apq(0,1);
    apq(4,3) = -apq(0,1);
    apq(4,4) =  apq(0,0) ;
    apq(5,5) =  1.0;
    apq(6,6) =  apq(0,0) ;
    apq(6,7) =  apq(0,1);
    apq(7,6) = -apq(0,1);
    apq(7,7) =  apq(0,0) ;
    apq(8,8) =  1.0;
    apq(9,9) =   apq(0,0) ;
    apq(9,10) =  apq(0,1);
    apq(10,9) = -apq(0,1);
    apq(10,10)=  apq(0,0) ;
    apq(11,11)=  1.0;

	// fill transformation coefficients, local to basic system
	getAvp();

//	opserr<<" Apq" <<endln << apq<< endln;
    
//	opserr<<" Avp" <<endln << avp<< endln;

}   

int
LehighJoint2d::commitState(void)
{
	// following code is not necessary for linear elastic element
	
	// store committed external nodal displacements
	vs = vt;

	// store material history data.
	int mcs = 0;
		for (int j=0; j<numBasicDOF; j++)
		{
			if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->commitState();
			if (mcs != 0) break;
		}
    
	return mcs;
}

int
LehighJoint2d::revertToLastCommit(void)
{
	int mcs = 0;
	for (int j=0; j<numBasicDOF; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToLastCommit();
		if (mcs != 0) break;
	}
	vt = vs;
	
//	this->update();

    return mcs;
}

int
LehighJoint2d::revertToStart(void)
{
	int mcs = 0;
	for (int j=0; j<numBasicDOF; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToStart();
		if (mcs != 0) break;
	}
	
	return mcs;
}

int
LehighJoint2d::update(void)
{	
	// update material state
	getBasicTrialDisp();

	for (int j=0; j<numBasicDOF; j++)
	{
		MaterialPtr[j]->setTrialStrain(vt(j));
	}

	return 0;
}

const Matrix &
LehighJoint2d::getTangentStiff(void)
{
	// determine the stiffness coefficient
	// assuming diagonal stiffness
	static Matrix kb(numBasicDOF,numBasicDOF);
	kb.Zero();
	
	for (int j=0; j<numBasicDOF; j++)
	{
		kb(j,j) = MaterialPtr[j]->getTangent();
	}
//	opserr << kb << endln;

    static Matrix Kp(numDOF,numDOF);
	Kp.Zero();
	Kp.addMatrixTripleProduct(0.0,avp,kb,1.0);    // kp = avp'*kb*avp
	
//	opserr << Kp << endln;

	// transform stiffness in local coordinate to global coordinate
	K.addMatrixTripleProduct(0.0,apq,Kp,1.0);    // K = apq'*kp*apq

//	opserr << K << endln;

	return K;
}

const Matrix &
LehighJoint2d::getInitialStiff(void)
{
	// assuming linear elastic element
	// assuming diagonal stiffness
	static Matrix kb(numBasicDOF,numBasicDOF);
	kb.Zero();
	
	for (int j=0; j<numBasicDOF; j++)
	{
		kb(j,j) = MaterialPtr[j]->getInitialTangent();
	}

    static Matrix Kp(numDOF,numDOF);
	Kp.Zero();
	Kp.addMatrixTripleProduct(0.0,avp,kb,1.0);    // kp = avp'*kb*avp
	
	// transform stiffness in local coordinate to global coordinate
	K.addMatrixTripleProduct(0.0,apq,Kp,1.0);    // K = apq'*kp*apq

	return K;
}

const Vector &
LehighJoint2d::getResistingForce(void)
{
	// basic forces
	static Vector fs(numBasicDOF);
  
	for (int j=0; j<numBasicDOF; j++)
	{
		fs(j) = MaterialPtr[j]->getStress();
	}

//	opserr << "fs" << endln << fs << endln;
	// transform the basic force to force in local coordinates
    static Vector P(numDOF);
    P.Zero();

	// determine deformation in the basic system eliminating rigid body modes 
    P.addMatrixTransposeVector(0.0,avp,fs,1.0);  // P = avp'*fs 

//	opserr << "P" << endln << P << endln;

	R.addMatrixTransposeVector(0.0,apq,P,1.0);  // P = apq'*P 
	
//	opserr << "R" << endln << R << endln;
	
	return R;
}

void LehighJoint2d::getBasicTrialDisp(void) 
{

	static Vector q(numDOF);
	q.Zero();
	
	// determine the trial displacement
    const Vector &disp1 = nodePtr[0]->getTrialDisp(); 
    const Vector &disp2 = nodePtr[1]->getTrialDisp();
    const Vector &disp3 = nodePtr[2]->getTrialDisp();
    const Vector &disp4 = nodePtr[3]->getTrialDisp();

	// assuming no initial displacement
	for (int i = 0; i < 3; i++)
    {
      q(i)     = disp1(i);
      q(i+3)   = disp2(i);
      q(i+6)   = disp3(i);
      q(i+9)   = disp4(i);
	  
    }

    // transform global end displacements to local coordinates
    static Vector p(numDOF);
	p.Zero();
	
    // determine deformation in the basic system eliminating rigid body modes 
    p.addMatrixVector(0.0,apq,q,1.0);   // p = apq*q

//	opserr<< "p :" << endln << p<< endln;

    // determine deformation in the basic system eliminating rigid body modes 
    vt.addMatrixVector(0.0,avp,p,1.0);   // v = avp*p
//	opserr<< "v :" << endln << vt << endln;
}


void LehighJoint2d::getAvp() 
{
	avp.Zero();                  // basically the transformation matrix for the element
                                 // from local to basic system

	avp(0,0) =  -1;
	avp(0,6) =  -avp(0,0);
	avp(1,4) =  -1;
	avp(1,10)=  -avp(1,4); 
	avp(2,1) =  -elemHeight/elemWidth;
	avp(2,3) =  -1;
	avp(2,7) =   -avp(2,1);
	avp(2,9) =  -avp(2,3);
	avp(3,2) =  -1;
	avp(3,8) =  -avp(3,2);
	avp(4,5) =  -1;
	avp(4,11) = -avp(4,5);
	avp(5,2) =   1;
	avp(5,3) =  -2/elemHeight;
	avp(5,8) =   1;
	avp(5,9) =  -avp(5,3);
	avp(6,1) =   2/elemWidth;
	avp(6,5) =   1;
	avp(6,7) =   -avp(6,1);
	avp(6,11) =  1;
	avp(7,0) =   1;
	avp(7,3) =  -avp(7,0);
	avp(7,6) =   1;
	avp(7,9) =  -avp(7,6);
	avp(8,1) =  -1;
	avp(8,4) =  -avp(8,1);
	avp(8,7) =  -1;
	avp(8,10) = -avp(8,7);

//	opserr<< avp << endln;
}

    
const Matrix &
LehighJoint2d::getDamp(void)
{
	//not applicable (stiffness being returned)
	K.Zero();
	return K;
}

const Matrix &
LehighJoint2d::getMass(void)
{ 
	//not applicable  (stiffness being returned)
	K.Zero();
	
	return K;
}

void 
LehighJoint2d::zeroLoad(void)
{
	//not applicable  
	return;
}

int 
LehighJoint2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	//not applicable
	return 0;
}

int 
LehighJoint2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	//not applicable
	return 0;
}


const Vector &
LehighJoint2d::getResistingForceIncInertia()
{	
  //not applicable (residual being returned)
	return this->getResistingForce();
	//return R;
}

int
LehighJoint2d::sendSelf(int commitTag, Channel &theChannel)
{
	// yet to do.

  int res;
  int i;
  int dataTag = this->getDbTag();
  
  static ID data(22);
  data(0) = this->getTag();
  data(1) = numDOF;
  
  if (nodeDbTag == 0) nodeDbTag = theChannel.getDbTag();
  if (dofDbTag == 0) dofDbTag = theChannel.getDbTag();
  data(2) = nodeDbTag;
  data(3) = dofDbTag;
  
   // sending Spring class and Db tags
  for (i=0 ; i<numBasicDOF ; i++) {
  //  data( i+4 ) = fixedEnd[i];
    if ( MaterialPtr[i] != NULL )
    {
		data(i+4) = MaterialPtr[i]->getClassTag();
		int SpringDbTag = MaterialPtr[i]->getDbTag();
		if (SpringDbTag == 0) {
			SpringDbTag = theChannel.getDbTag();
			if (SpringDbTag != 0)
				MaterialPtr[i]->setDbTag(SpringDbTag);
		}
		data(i+13) = SpringDbTag;
    } 
	else
	{
	  data( i+4 ) = 0;
	  data( i+13 ) = 0;
	}
  }
  
  // send the ID vector
  
  res = theChannel.sendID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LehighJoint::sendSelf() - " << this->getTag() << "failed to send ID\n";
    return -1;
  }
  
  // sends the tags of it's external nodes
  res = theChannel.sendID(nodeDbTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING LehighJoint::sendSelf() - " << this->getTag()<< " failed to send Vector\n";
    return -2;
  }
  
   
  // finally send the materials one by one
  
  for ( i=0 ; i<numBasicDOF ; i++ ) {
    if ( MaterialPtr[i] != NULL ) {
      res = MaterialPtr[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
	opserr << "WARNING LehighJoint::sendSelf() - "<< this->getTag() << " failed to send its Spring " << (i+1) << " material\n";
	return -3;
      }
    }
  }
  
  return 0;
}

int
LehighJoint2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// yet to do.
  int res;
  int dataTag = this->getDbTag();
  
  static ID data(22);
  res = theChannel.recvID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Joint2D::recvSelf() - failed to receive Vector\n";
    return -1;
  }
  
  this->setTag((int)data(0));
  numDOF = data(1);
  nodeDbTag = data(2);
  dofDbTag = data(3);
  
  // Receiving Springs
  for (int i=0 ; i<numBasicDOF ; i++) {
    int SpringClass = data( i+4 );
    int SpringDb = data( i+13 );
    
    if ( SpringClass != 0 && SpringDb != 0  )
      {
	// check if we have a material object already & if we do if of right type
	if ((MaterialPtr[i] == 0) || (MaterialPtr[i]->getClassTag() != SpringClass)) {
	  // if old one .. delete it
	  if (MaterialPtr[i] != 0)
					delete MaterialPtr[i];
	  // create a new material object
	  MaterialPtr[i] = theBroker.getNewUniaxialMaterial(SpringClass);
	  if (MaterialPtr[i] == 0) {
	    opserr << "WARNING Joint2D::recvSelf() - " << (i+1) << " failed to get a blank Material of type " << this->getTag() << " for Spring " << SpringClass << "\n";
	    return -3;
	  }
	}
	
	MaterialPtr[i]->setDbTag(SpringDb); // note: we set the dbTag before we receive the material
	res = MaterialPtr[i]->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "WARNING Joint2D::recvSelf() - " << this->getTag() << " failed to receive its Material for Spring " << (i+1) << "\n";
	  return -3;
	}
      }
    else MaterialPtr[i] = NULL;
  }
  
  
  // receives the tags of it's external nodes
  res = theChannel.recvID(nodeDbTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING Joint2D::recvSelf() - " << this->getTag() << " failed to receive external nodes\n" ;
    return -2;
  }
   
  
  return 0;
}

int
LehighJoint2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// get coordinates for connecting nodes
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);

	nodePtr[0]->getDisplayCrds(v1, fact, displayMode);
	nodePtr[1]->getDisplayCrds(v2, fact, displayMode);
	nodePtr[2]->getDisplayCrds(v3, fact, displayMode);
	nodePtr[3]->getDisplayCrds(v4, fact, displayMode);

	// calculate four corners of the element
	Vector w(3); // width vector
	Vector c1(3);
	Vector c2(3);
	Vector c3(3);
	Vector c4(3);

	w = v2 - v4;
	c1 = v1 - 0.5 * w;
	c2 = v1 + 0.5 * w;
	c3 = v3 + 0.5 * w;
	c4 = v3 - 0.5 * w;

	int res = 0;
	res += theViewer.drawLine(c1, c2, 1.0, 1.0, this->getTag());
	res += theViewer.drawLine(c2, c3, 1.0, 1.0, this->getTag());
	res += theViewer.drawLine(c3, c4, 1.0, 1.0, this->getTag());
	res += theViewer.drawLine(c4, c1, 1.0, 1.0, this->getTag());

	return res;
}

void
LehighJoint2d::Print(OPS_Stream &s, int flag)
{
	s << "Element: " << this->getTag() << " Type: Beam Column Joint " << endln;
	for (int i = 0; i<4; i++)
	{
		s << "Node :" << connectedExternalNodes(i);
		s << "DOF :" << nodePtr[i]->getNumberDOF();
	}
	return;
}

Response*
LehighJoint2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  // we will compare argv[0] to determine the type of response required
  
  if (strcmp(argv[0],"globalForces") == 0 || strcmp(argv[0],"globalForce") == 0)
    return new ElementResponse(this,1,Vector(12));

  else if (strcmp(argv[0],"localForces") == 0 || strcmp(argv[0],"localForce") == 0)
    return new ElementResponse(this,2,Vector(12));  
  
  else if (strcmp(argv[0],"elementForces") == 0 || strcmp(argv[0],"basicForces") == 0)
  	return new ElementResponse(this,3,Vector(9));

  else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"Deformation") == 0)
  	return new ElementResponse(this,4,Vector(9));
  
  else
    return 0;
}

int 
LehighJoint2d::getResponse(int responseID, Information &eleInfo)
{

	switch (responseID) {
	case 1:       
		if(eleInfo.theVector!=0)
		{
			static Vector R(numDOF);
            R = this->getResistingForce();

			(*(eleInfo.theVector))(0) = R(0);
			(*(eleInfo.theVector))(1) = R(1);
			(*(eleInfo.theVector))(2) = R(2);
			(*(eleInfo.theVector))(3) = R(3);
			(*(eleInfo.theVector))(4) = R(4);
			(*(eleInfo.theVector))(5) = R(5);
			(*(eleInfo.theVector))(6) = R(6);
			(*(eleInfo.theVector))(7) = R(7);
			(*(eleInfo.theVector))(8) = R(8);
			(*(eleInfo.theVector))(9) = R(9);
			(*(eleInfo.theVector))(10) = R(10);
			(*(eleInfo.theVector))(11) = R(11);

		}
		return 0;

	case 2:
		if (eleInfo.theVector !=0) {
			
			static Vector fs(numBasicDOF);
			static Vector P(numDOF);			
			
			for (int j=0; j<numBasicDOF; j++)
			{
				fs(j) = MaterialPtr[j]->getStress();
			}

			// transform the basic force to force in local coordinates
			P.Zero();

			// determine element force in the local coord 
			P.addMatrixTransposeVector(0.0,avp,fs,1.0);  // P = avp'*fs 

			(*(eleInfo.theVector))(0) = P(0);
			(*(eleInfo.theVector))(1) = P(1);
			(*(eleInfo.theVector))(2) = P(2);
			(*(eleInfo.theVector))(3) = P(3);
			(*(eleInfo.theVector))(4) = P(4);
			(*(eleInfo.theVector))(5) = P(5);
			(*(eleInfo.theVector))(6) = P(6);
			(*(eleInfo.theVector))(7) = P(7);
			(*(eleInfo.theVector))(8) = P(8);
			(*(eleInfo.theVector))(9) = P(9);
			(*(eleInfo.theVector))(10) = P(10);
			(*(eleInfo.theVector))(11) = P(11);
		}
		return 0;

    case 3:
		if(eleInfo.theVector!=0) {

			for ( int i =0 ; i<numBasicDOF ; i++ )	{
				
				(*(eleInfo.theVector))(i) = 0.0;

				if ( MaterialPtr[i] != NULL ) 
					(*(eleInfo.theVector))(i) = MaterialPtr[i]->getStress();
			}
		}
		return 0;

	case 4:
		if(eleInfo.theVector!=0) {

			for ( int i =0 ; i<numBasicDOF ; i++ )	{
				
				(*(eleInfo.theVector))(i) = 0.0;

				if ( MaterialPtr[i] != NULL ) 
					(*(eleInfo.theVector))(i) = MaterialPtr[i]->getStrain();
			}
		}
		return 0;
	
	

	default:
		return -1;
	}
}

int
LehighJoint2d::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
    
int
LehighJoint2d::updateParameter (int parameterID, Information &info)
{
  return -1;
}
