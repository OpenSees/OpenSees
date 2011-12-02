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
// $Date: 2009-04-17 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elastomericBearing/ElastomericBearing3d.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the implementation of the
// ElastomericBearing3d class.
//
// What: "@(#) ElastomericBearing3d.cpp, revA"

#include "ElastomericBearing3d.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


// initialize the class wide variables
Matrix ElastomericBearing3d::theMatrix(12,12);
Vector ElastomericBearing3d::theVector(12);
Vector ElastomericBearing3d::theLoad(12);


ElastomericBearing3d::ElastomericBearing3d(int tag, int Nd1, int Nd2,
    double ke, double fy, double alpha, UniaxialMaterial **materials,
    const Vector _y, const Vector _x, double m)
    : Element(tag, ELE_TAG_ElastomericBearing3d),
    connectedExternalNodes(2),
    k0(0.0), qYield(0.0), k2(0.0), x(_x), y(_y), mass(m),
    L(0.0), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElastomericBearing3d::setUp() - element: "
            << this->getTag() << " failed to create an ID of size 2\n";
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // initialize parameters
    k0 = (1.0-alpha)*ke;
    qYield = (1.0-alpha)*fy;
    k2 = alpha*ke;
    
    // check material input
    if (materials == 0)  {
        opserr << "ElastomericBearing3d::ElastomericBearing3d() - "
            << "null material array passed.\n";
        exit(-1);
    }
    
    // get copies of the uniaxial materials
    for (int i=0; i<4; i++)  {
        if (materials[i] == 0) {
            opserr << "ElastomericBearing3d::ElastomericBearing3d() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0) {
            opserr << "ElastomericBearing3d::ElastomericBearing3d() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = ke;
    kbInit(2,2) = ke;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
    
    // initialize variables
    this->revertToStart();
}


ElastomericBearing3d::ElastomericBearing3d()
    : Element(0, ELE_TAG_ElastomericBearing3d),
    connectedExternalNodes(2),
    k0(0.0), qYield(0.0), k2(0.0), x(0), y(0), mass(0.0),
    L(0.0), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(6),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6)
{	
    // ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2)  {
		opserr << "ElastomericBearing3d::ElastomericBearing3d() - "
			<<  "failed to create an ID of size 2\n";
		exit(-1);
    }
    
    // set node pointers to NULL
	for (int i=0; i<2; i++)
		theNodes[i] = 0;
    
    // set material pointers to NULL
	for (int i=0; i<4; i++)
		theMaterials[i] = 0;
}


ElastomericBearing3d::~ElastomericBearing3d()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int ElastomericBearing3d::getNumExternalNodes() const
{
    return 2;
}


const ID& ElastomericBearing3d::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** ElastomericBearing3d::getNodePtrs() 
{
	return theNodes;
}


int ElastomericBearing3d::getNumDOF() 
{
    return 12;
}


void ElastomericBearing3d::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain)  {
		theNodes[0] = 0;
		theNodes[1] = 0;
        
		return;
    }

    // first set the node pointers
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));	
	
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
		if (!theNodes[0])  {
			opserr << "WARNING ElastomericBearing3d::setDomain() - Nd1: " 
				<< connectedExternalNodes(0) << " does not exist in the model for ";
		} else  {
			opserr << "WARNING ElastomericBearing3d::setDomain() - Nd2: " 
				<< connectedExternalNodes(1) << " does not exist in the model for ";
		}
		opserr << "ElastomericBearing3d ele: " << this->getTag() << endln;
		
		return;
    }
	
	// now determine the number of dof and the dimension    
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();	
	
	// if differing dof at the ends - print a warning message
    if (dofNd1 != 6)  {
		opserr << "ElastomericBearing3d::setDomain() - node 1: "
			<< connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
		return;
    }
    if (dofNd2 != 6)  {
		opserr << "ElastomericBearing3d::setDomain() - node 2: "
			<< connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
		return;
    }
	
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}   	 


int ElastomericBearing3d::commitState()
{
	int errCode = 0;
    
    // commit trial history variables
    ubPlasticC = ubPlastic;
	
    // commit material models
    for (int i=0; i<4; i++)
	    errCode += theMaterials[i]->commitState();
    
	return errCode;
}


int ElastomericBearing3d::revertToLastCommit()
{
    int errCode = 0;
    
    // revert material models
    for (int i=0; i<4; i++)
	    errCode += theMaterials[i]->revertToLastCommit();
    
    return errCode;
}


int ElastomericBearing3d::revertToStart()
{   
    int errCode=0;
    
    // reset trial history variables
    ub.Zero();
    ubPlastic.Zero();
    qb.Zero();
    
    // reset committed history variables
    ubPlasticC.Zero();
    
    // reset stiffness matrix in basic system
    kb = kbInit;
    
    // revert material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->revertToStart();
    
    return errCode;
}


int ElastomericBearing3d::update()
{
    // get global trial displacements and velocities
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    
    static Vector ug(12), ugdot(12), uldot(12), ubdot(6);
    for (int i=0; i<6; i++)  {
        ug(i)   = dsp1(i);  ugdot(i)   = vel1(i);
        ug(i+6) = dsp2(i);  ugdot(i+6) = vel2(i);
    }
    
    // transform response from the global to the local system
    ul = Tgl*ug;
    uldot = Tgl*ugdot;
    
    // transform response from the local to the basic system
    ub = Tlb*ul;
    ubdot = Tlb*uldot;
    
    // 1) get axial force and stiffness in basic x-direction
    theMaterials[0]->setTrialStrain(ub(0),ubdot(0));
    qb(0) = theMaterials[0]->getStress();
    kb(0,0) = theMaterials[0]->getTangent();
    
    // 2) calculate shear forces and stiffnesses in basic y- and z-direction
    // get trial shear forces of hysteretic component
    Vector qTrial(2);
    qTrial(0) = k0*(ub(1) - ubPlasticC(0));
    qTrial(1) = k0*(ub(2) - ubPlasticC(1));
    
    // compute yield criterion of hysteretic component
    double qTrialNorm = qTrial.Norm();
    double Y = qTrialNorm - qYield;
    
    // elastic step -> no updates required
    if (Y <= 0.0)  {
        // set shear forces
        qb(1) = k2*ub(1) + qTrial(0);
        qb(2) = k2*ub(2) + qTrial(1);
        // set tangent stiffnesses
        kb(1,1) = kb(2,2) = k2 + k0;
        kb(1,2) = kb(2,1) = 0.0;
    }
    // plastic step -> return mapping
    else  {
        // compute consistency parameters
        double dGamma = Y/k0;
        // update plastic displacements
        ubPlastic(0) = ubPlasticC(0) + dGamma*qTrial(0)/qTrialNorm;
        ubPlastic(1) = ubPlasticC(1) + dGamma*qTrial(1)/qTrialNorm;
        // set shear forces
        qb(1) = k2*ub(1) + qYield*qTrial(0)/qTrialNorm;
        qb(2) = k2*ub(2) + qYield*qTrial(1)/qTrialNorm;
        // set tangent stiffnesses
        kb(1,1) = k2 + qYield*k0*pow(qTrial(1),2)/pow(qTrialNorm,3);
        kb(1,2) = -qYield*k0*qTrial(0)*qTrial(1)/pow(qTrialNorm,3);
        kb(2,1) = -qYield*k0*qTrial(0)*qTrial(1)/pow(qTrialNorm,3);
        kb(2,2) = k2 + qYield*k0*pow(qTrial(0),2)/pow(qTrialNorm,3);
    }
    
    // 3) get moment and stiffness in basic x-direction
    theMaterials[1]->setTrialStrain(ub(3),ubdot(3));
    qb(3) = theMaterials[1]->getStress();
    kb(3,3) = theMaterials[1]->getTangent();
    
    // 4) get moment and stiffness in basic y-direction
    theMaterials[2]->setTrialStrain(ub(4),ubdot(4));
    qb(4) = theMaterials[2]->getStress();
    kb(4,4) = theMaterials[2]->getTangent();
    
    // 5) get moment and stiffness in basic z-direction
    theMaterials[3]->setTrialStrain(ub(5),ubdot(5));
    qb(5) = theMaterials[3]->getStress();
    kb(5,5) = theMaterials[3]->getTangent();
    
    return 0;
}


const Matrix& ElastomericBearing3d::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    kl(5,1) -= 0.5*qb(0);
    kl(5,7) += 0.5*qb(0);
    kl(11,1) -= 0.5*qb(0);
    kl(11,7) += 0.5*qb(0);
    kl(4,2) += 0.5*qb(0);
    kl(4,8) -= 0.5*qb(0);
    kl(10,2) += 0.5*qb(0);
    kl(10,8) -= 0.5*qb(0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearing3d::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearing3d::getMass()
{
	// zero the matrix
    theMatrix.Zero();
    
	// check for quick return
	if (mass == 0.0)  {
		return theMatrix;
	}    
    
	double m = 0.5*mass;
	for (int i = 0; i < 3; i++)  {
		theMatrix(i,i)     = m;
		theMatrix(i+3,i+3) = m;
	}
	
    return theMatrix; 
}


void ElastomericBearing3d::zeroLoad()
{
    theLoad.Zero();
}


int ElastomericBearing3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
	opserr <<"ElastomericBearing3d::addLoad() - "
		<< "load type unknown for element: "
		<< this->getTag() << endln;
    
	return -1;
}


int ElastomericBearing3d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// check for quick return
	if (mass == 0.0)  {
		return 0;
	}
    
	// get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
	if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
		opserr << "ElastomericBearing3d::addInertiaLoadToUnbalance() - "
			<< "matrix and vector sizes are incompatible\n";
		return -1;
	}
    
	// want to add ( - fact * M R * accel ) to unbalance
	// take advantage of lumped mass matrix
	double m = 0.5*mass;
    for (int i = 0; i < 3; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+3) -= m * Raccel2(i);
    }
    
	return 0;
}


const Vector& ElastomericBearing3d::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(12);
    ql = Tlb^qb;
    
    // add P-Delta moments to local forces
    double MpDelta;
    MpDelta = qb(0)*(ul(7)-ul(1));
    ql(5)  += 0.5*MpDelta;
    ql(11) += 0.5*MpDelta;
    MpDelta = qb(0)*(ul(8)-ul(2));
    ql(4)  -= 0.5*MpDelta;
    ql(10) -= 0.5*MpDelta;
    
    // determine resisting forces in global system
    theVector = Tgl^ql;
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    return theVector;
}


const Vector& ElastomericBearing3d::getResistingForceIncInertia()
{	
	theVector = this->getResistingForce();
	
	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		theVector += this->getRayleighDampingForces();
    
	// now include the mass portion
	if (mass != 0.0)  {
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();    
		
		double m = 0.5*mass;
		for (int i = 0; i < 3; i++)  {
			theVector(i)   += m * accel1(i);
			theVector(i+3) += m * accel2(i);
		}
	}
	
	return theVector;
}


int ElastomericBearing3d::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(7);
    data(0) = this->getTag();
    data(1) = k0;
    data(2) = qYield;
    data(3) = k2;
    data(4) = mass;
    data(5) = x.Size();
    data(6) = y.Size();
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send the material class tags
    ID matClassTags(4);
    for (int i=0; i<4; i++)
        matClassTags(i) = theMaterials[i]->getClassTag();
    sChannel.sendID(0, commitTag, matClassTags);
    
    // send the material models
    for (int i=0; i<4; i++)
        theMaterials[i]->sendSelf(commitTag, sChannel);
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    
    return 0;
}


int ElastomericBearing3d::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete material memory
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
    
    // receive element parameters
    static Vector data(7);
    rChannel.recvVector(0, commitTag, data);    
    this->setTag((int)data(0));
    k0 = data(1);
    qYield = data(2);
    k2 = data(3);
    mass = data(4);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive the material class tags
    ID matClassTags(4);
    rChannel.recvID(0, commitTag, matClassTags);
    
    // receive the material models
    for (int i=0; i<4; i++)  {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "ElastomericBearing2d::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
        theMaterials[i]->recvSelf(commitTag, rChannel, theBroker);
    }
    
    // receive remaining data
    if ((int)data(5) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(6) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = k2 + k0;
    kbInit(2,2) = k2 + k0;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
    
    // initialize variables
    this->revertToStart();
    
    return 0;
}


int ElastomericBearing3d::displaySelf(Renderer &theViewer,
    int displayMode, float fact)
{
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();
    
    static Vector v1(3);
    static Vector v2(3);
    
    for (int i = 0; i < 3; i++)  {
        v1(i) = end1Crd(i) + end1Disp(i)*fact;
        v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}


void ElastomericBearing3d::Print(OPS_Stream &s, int flag)
{
    if (flag == 0)  {
        // print everything
		s << "Element: " << this->getTag(); 
		s << "  type: ElastomericBearing3d  iNode: " << connectedExternalNodes(0);
		s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  k0: " << k0 << "  qYield: " << qYield << "  k2: " << k2 << endln;
        s << "  Material ux: " << theMaterials[0]->getTag() << endln;
        s << "  Material rx: " << theMaterials[1]->getTag() << endln;
        s << "  Material ry: " << theMaterials[2]->getTag() << endln;
        s << "  Material rz: " << theMaterials[3]->getTag() << endln;
        s << "  mass: " << mass << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    } else if (flag == 1)  {
		// does nothing
    }
}


Response* ElastomericBearing3d::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ElastomericBearing3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Pz_1");
        output.tag("ResponseType","Mx_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","Mx_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
        output.tag("ResponseType","N_ 1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Tz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 2, theVector);
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
        output.tag("ResponseType","qb1");
        output.tag("ResponseType","qb2");
        output.tag("ResponseType","qb3");
        output.tag("ResponseType","qb4");
        output.tag("ResponseType","qb5");
        output.tag("ResponseType","qb6");
        
        theResponse = new ElementResponse(this, 3, Vector(6));
    }
	// local displacements
    else if (strcmp(argv[0],"localDisplacement") == 0 ||
        strcmp(argv[0],"localDisplacements") == 0)
    {
        output.tag("ResponseType","ux_1");
        output.tag("ResponseType","uy_1");
        output.tag("ResponseType","uz_1");
        output.tag("ResponseType","rx_1");
        output.tag("ResponseType","ry_1");
        output.tag("ResponseType","rz_1");
        output.tag("ResponseType","ux_2");
        output.tag("ResponseType","uy_2");
        output.tag("ResponseType","uz_2");
        output.tag("ResponseType","rx_2");
        output.tag("ResponseType","ry_2");
        output.tag("ResponseType","rz_2");
        
        theResponse = new ElementResponse(this, 4, theVector);
    }
	// basic displacements
    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
        output.tag("ResponseType","ub1");
        output.tag("ResponseType","ub2");
        output.tag("ResponseType","ub3");
        output.tag("ResponseType","ub4");
        output.tag("ResponseType","ub5");
        output.tag("ResponseType","ub6");
        
        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // material output
    else if (strcmp(argv[0],"material") == 0)  {
        if (argc > 2)  {
            int matNum = atoi(argv[1]);
            if (matNum >= 1 && matNum <= 4)
                theResponse =  theMaterials[matNum-1]->setResponse(&argv[2], argc-2, output);
        }
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int ElastomericBearing3d::getResponse(int responseID, Information &eleInfo)
{
    double MpDelta;
    
    switch (responseID)  {
	case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
	case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector = Tlb^qb;
        // add P-Delta moments
        MpDelta = qb(0)*(ul(7)-ul(1));
        theVector(5)  += 0.5*MpDelta;
        theVector(11) += 0.5*MpDelta;
        MpDelta = qb(0)*(ul(8)-ul(2));
        theVector(4)  -= 0.5*MpDelta;
        theVector(10) -= 0.5*MpDelta;
        
        return eleInfo.setVector(theVector);
        
	case 3:  // basic forces
        return eleInfo.setVector(qb);
        
	case 4:  // local displacements
        return eleInfo.setVector(ul);
        
	case 5:  // basic displacements
        return eleInfo.setVector(ub);
        
    default:
		return -1;
	}
}


// set up the transformation matrix for orientation
void ElastomericBearing3d::setUp()
{ 
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    if (L > DBL_EPSILON)  {
		if (x.Size() == 0)  {
		    x.resize(3);
		    x = xp;
        } else  {
            opserr << "WARNING ElastomericBearing3d::setUp() - " 
                << "element: " << this->getTag() << endln
                << "ignoring nodes and using specified "
                << "local x vector to determine orientation\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "ElastomericBearing3d::setUp() - "
            << "element: " << this->getTag() << endln
            << "incorrect dimension of orientation vectors\n";
        exit(-1);
    }
    
    // establish orientation of element for the tranformation matrix
    // z = x cross y
    Vector z(3);
    z(0) = x(1)*y(2) - x(2)*y(1);
    z(1) = x(2)*y(0) - x(0)*y(2);
    z(2) = x(0)*y(1) - x(1)*y(0);
    
    // y = z cross x
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);
    
    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();
    
    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0)  {
        opserr << "ElastomericBearing3d::setUp() - "
            << "element: " << this->getTag() << endln
            << "invalid orientation vectors\n";
        exit(-1);
    }
    
    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = x(0)/xn;
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = x(1)/xn;
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = x(2)/xn;
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = y(0)/yn;
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = y(1)/yn;
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = y(2)/yn;
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = z(0)/zn;
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = z(1)/zn;
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = z(2)/zn;
    
    // create transformation matrix from local to basic system (linear)
    Tlb.Zero();
    Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
    Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;
    Tlb(1,5) = Tlb(1,11) = -0.5*L;
    Tlb(2,4) = Tlb(2,10) = 0.5*L;
}


double ElastomericBearing3d::sgn(double x)
{ 
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
