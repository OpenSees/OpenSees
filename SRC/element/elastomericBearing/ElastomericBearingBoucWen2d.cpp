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

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the implementation of the
// ElastomericBearingBoucWen2d class.

#include "ElastomericBearingBoucWen2d.h"

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

#include <elementAPI.h>

// initialize the class wide variables
Matrix ElastomericBearingBoucWen2d::theMatrix(6,6);
Vector ElastomericBearingBoucWen2d::theVector(6);

void * OPS_ADD_RUNTIME_VPV(OPS_ElastomericBearingBoucWen2d)
{
    int ndf = OPS_GetNDF();
    if (ndf != 3)  {
	opserr << "WARNING invalid ndf: " << ndf;
	opserr << ", for plane problem need 3 - elastomericBearingBoucWen\n";
	return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 15) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: elastomericBearingBoucWen eleTag iNode jNode kInit qd alpha1 alpha2 mu eta beta gamma -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
	return 0;
    }

    // tags
    int idata[3];
    int num = 3;
    if (OPS_GetIntInput(&num, idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    // data
    double data[8];
    num = 8;
    if (OPS_GetDoubleInput(&num, data) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    // materials
    UniaxialMaterial* mats[2] = {0,0};
    const char* type = OPS_GetString();
    if (strcmp(type, "-P") != 0) {
	opserr<<"WARNING: want -P\n";
	return 0;
    }
    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num, &matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }
    mats[0] = OPS_getUniaxialMaterial(matTag);
    if (mats[0] == 0) {
	opserr<<"WARNING: material not found\n";
	return 0;
    }

    type = OPS_GetString();
    if (strcmp(type, "-Mz") != 0) {
	opserr<<"WARNING: want -Mz\n";
	return 0;
    }
    num = 1;
    if (OPS_GetIntInput(&num, &matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }
    mats[1] = OPS_getUniaxialMaterial(matTag);
    if (mats[1] == 0) {
	opserr<<"WARNING: material not found\n";
	return 0;
    }
    
    // options
    Vector x,y;
    double sDistI = 0.5;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1e-12;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	return new ElastomericBearingBoucWen2d(idata[0],idata[1], idata[2],data[0],data[1],data[2],
					       mats,y,x,data[3],data[4],data[5],data[6],data[7],
					       sDistI,doRayleigh,mass,maxIter,tol);
    }
    while (OPS_GetNumRemainingInputArgs() > 0) {
	type = OPS_GetString();
	if (strcmp(type, "-orient") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 6) {
		opserr<<"WARNING: insufficient arguments after -orient\n";
		return 0;
	    }
	    num = 3;
	    x.resize(3);
	    if (OPS_GetDoubleInput(&num, &x(0)) < 0) {
		opserr<<"WARNING: invalid orient value\n";
		return 0;
	    }
	    y.resize(3);
	    if (OPS_GetDoubleInput(&num, &y(0)) < 0) {
		opserr<<"WARNING: invalid orient value\n";
		return 0;
	    }
	} else if (strcmp(type, "-shearDist") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return 0;
	    }
	    num = 1;
	    if (OPS_GetDoubleInput(&num, &sDistI) < 0) {
		opserr<<"WARNING: invalid shearDist\n";
		return 0;
	    }
	} else if (strcmp(type, "-doRayleigh") == 0) {
	    doRayleigh = 1;
	} else if (strcmp(type, "-mass") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return 0;
	    }
	    num = 1;
	    if (OPS_GetDoubleInput(&num, &mass) < 0) {
		opserr<<"WARNING: invalid mass\n";
		return 0;
	    }
	} else if (strcmp(type, "-iter") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr<<"WARNING: insufficient args\n";
		return 0;
	    }
	    num = 1;
	    if (OPS_GetIntInput(&num,&maxIter) < 0) {
		opserr<<"WARNING: invalid maxIter\n";
		return 0;
	    }
	    if (OPS_GetDoubleInput(&num,&tol) < 0) {
		opserr<<"WARNING: invalid tol\n";
		return 0;
	    }
	}
    }

    return new ElastomericBearingBoucWen2d(idata[0],idata[1], idata[2],data[0],data[1],data[2],
					   mats,y,x,data[3],data[4],data[5],data[6],data[7],
					   sDistI,doRayleigh,mass,maxIter,tol);
}


ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d(int tag,
    int Nd1, int Nd2, double kInit, double qd, double alpha1,
    UniaxialMaterial **materials, const Vector _y, const Vector _x,
    double alpha2, double _mu, double _eta, double _beta, double _gamma,
    double sdI, int addRay, double m, int maxiter, double _tol)
    : Element(tag, ELE_TAG_ElastomericBearingBoucWen2d),
    connectedExternalNodes(2), k0(0.0), qYield(qd), k2(0.0), k3(0.0),
    mu(_mu), eta(_eta), beta(_beta), gamma(_gamma), A(1.0), x(_x), y(_y),
    shearDistI(sdI), addRayleigh(addRay), mass(m), maxIter(maxiter), tol(_tol),
    L(0.0), onP0(true), ub(3), z(0.0), dzdu(0.0), qb(3), kb(3,3), ul(6),
    Tgl(6,6), Tlb(3,6), ubC(3), zC(0.0), kbInit(3,3), theLoad(6)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // initialize stiffnesses
    k0 = (1.0-alpha1)*kInit;
    k2 = alpha1*kInit;
    k3 = alpha2*kInit;
    
    // check material input
    if (materials == 0)  {
        opserr << "ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d() - "
            << "null material array passed.\n";
        exit(-1);
    }
    
    // get copies of the uniaxial materials
    for (int i=0; i<2; i++)  {
        if (materials[i] == 0) {
            opserr << "ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0) {
            opserr << "ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = A*k0 + k2;
    kbInit(2,2) = theMaterials[1]->getInitialTangent();
    
    // initialize other variables
    this->revertToStart();
}


ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d()
    : Element(0, ELE_TAG_ElastomericBearingBoucWen2d),
    connectedExternalNodes(2), k0(0.0), qYield(0.0), k2(0.0), k3(0.0),
    mu(2.0), eta(1.0), beta(0.5), gamma(0.5), A(1.0), x(0), y(0),
    shearDistI(0.5), addRayleigh(0), mass(0.0), maxIter(25), tol(1E-12),
    L(0.0), onP0(false), ub(3), z(0.0), dzdu(0.0), qb(3), kb(3,3), ul(6),
    Tgl(6,6), Tlb(3,6), ubC(3), zC(0.0), kbInit(3,3), theLoad(6)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElastomericBearingBoucWen2d::ElastomericBearingBoucWen2d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;    
    
    // set material pointers to NULL
    for (int i=0; i<2; i++)
        theMaterials[i] = 0;
}


ElastomericBearingBoucWen2d::~ElastomericBearingBoucWen2d()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    for (int i=0; i<2; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int ElastomericBearingBoucWen2d::getNumExternalNodes() const
{
    return 2;
}


const ID& ElastomericBearingBoucWen2d::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** ElastomericBearingBoucWen2d::getNodePtrs()
{
    return theNodes;
}


int ElastomericBearingBoucWen2d::getNumDOF()
{
    return 6;
}


void ElastomericBearingBoucWen2d::setDomain(Domain *theDomain)
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
            opserr << "WARNING ElastomericBearingBoucWen2d::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING ElastomericBearingBoucWen2d::setDomain() - Nd2: "
                << connectedExternalNodes(1)
                << " does not exist in the model for";
        }
        opserr << " element: " << this->getTag() << ".\n";
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != 3)  {
        opserr << "ElastomericBearingBoucWen2d::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 3).\n";
        return;
    }
    if (dofNd2 != 3)  {
        opserr << "ElastomericBearingBoucWen2d::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 3).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int ElastomericBearingBoucWen2d::commitState()
{
    int errCode = 0;
    
    // commit trial history variables
    ubC = ub;
    zC  = z;
    
    // commit material models
    for (int i=0; i<2; i++)
        errCode += theMaterials[i]->commitState();
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int ElastomericBearingBoucWen2d::revertToLastCommit()
{
    int errCode = 0;
    
    // revert material models
    for (int i=0; i<2; i++)
        errCode += theMaterials[i]->revertToLastCommit();
    
    return errCode;
}


int ElastomericBearingBoucWen2d::revertToStart()
{
    int errCode = 0;
    
    // reset trial history variables
    ub.Zero();
    z = 0.0;
    qb.Zero();
    
    // reset committed history variables
    ubC.Zero();
    zC = 0.0;
    
    // reset derivative of hysteretic evolution parameter * uy
    dzdu = A;
    
    // reset stiffness matrix in basic system
    kb = kbInit;
    
    // revert material models
    for (int i=0; i<2; i++)
        errCode += theMaterials[i]->revertToStart();
    
    return errCode;
}


int ElastomericBearingBoucWen2d::update()
{
    // get global trial displacements and velocities
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    
    static Vector ug(6), ugdot(6), uldot(6), ubdot(3);
    for (int i=0; i<3; i++)  {
        ug(i)   = dsp1(i);  ugdot(i)   = vel1(i);
        ug(i+3) = dsp2(i);  ugdot(i+3) = vel2(i);
    }
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
    
    // 1) get axial force and stiffness in basic x-direction
    theMaterials[0]->setTrialStrain(ub(0), ubdot(0));
    qb(0) = theMaterials[0]->getStress();
    kb(0,0) = theMaterials[0]->getTangent();
    
    // 2) calculate shear force and stiffness in basic y-direction
    // get displacement increment (trial - committed)
    double delta_ub = ub(1) - ubC(1);
    if (fabs(delta_ub) > 0.0)  {
        
        // get yield displacement
        double uy = qYield/k0;
        
        // calculate hysteretic evolution parameter z using Newton-Raphson
        int iter = 0;
        double zAbs, tmp1, f, Df, delta_z;
        do  {
            zAbs = fabs(z);
            if (zAbs == 0.0)    // check because of negative exponents
                zAbs = DBL_EPSILON;
            tmp1 = gamma + beta*sgn(z*delta_ub);
            
            // function and derivative
            f  = z - zC - delta_ub/uy*(A - pow(zAbs,eta)*tmp1);
            Df = 1.0 + delta_ub/uy*eta*pow(zAbs,eta-1.0)*sgn(z)*tmp1;
            
            // issue warning if derivative Df is zero
            if (fabs(Df) <= DBL_EPSILON)  {
                opserr << "WARNING: ElastomericBearingBoucWen2d::update() - "
                    << "zero derivative in Newton-Raphson scheme for "
                    << "hysteretic evolution parameter z.\n";
                return -1;
            }
            
            // advance one step
            delta_z = f/Df;
            z -= delta_z;
            iter++;
        } while ((fabs(delta_z) >= tol) && (iter < maxIter));
        
        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter)   {
            opserr << "WARNING: ElastomericBearingBoucWen2d::update() - "
                << "did not find the hysteretic evolution parameter z after "
                << iter << " iterations and norm: " << fabs(delta_z) << endln;
            return -2;
        }
        
        // get derivative of hysteretic evolution parameter * uy
        dzdu = A - pow(fabs(z),eta)*(gamma + beta*sgn(z*delta_ub));
        // set shear force
        qb(1) = qYield*z + k2*ub(1) + k3*sgn(ub(1))*pow(fabs(ub(1)),mu);
        // set tangent stiffness
        kb(1,1) = k0*dzdu + k2 + k3*mu*pow(fabs(ub(1)),mu-1.0);
    }
    
    // 3) get moment and stiffness about basic z-direction
    theMaterials[1]->setTrialStrain(ub(2), ubdot(2));
    qb(2) = theMaterials[1]->getStress();
    kb(2,2) = theMaterials[1]->getTangent();
    
    return 0;
}


const Matrix& ElastomericBearingBoucWen2d::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(6,6);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    double kGeo1 = 0.5*qb(0);
    kl(2,1) -= kGeo1;
    kl(2,4) += kGeo1;
    kl(5,1) -= kGeo1;
    kl(5,4) += kGeo1;
    double kGeo2 = kGeo1*shearDistI*L;
    kl(2,2) += kGeo2;
    kl(5,2) -= kGeo2;
    double kGeo3 = kGeo1*(1.0 - shearDistI)*L;
    kl(2,5) -= kGeo3;
    kl(5,5) += kGeo3;
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWen2d::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix klInit(6,6);
    klInit.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWen2d::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    // call base class to setup Rayleigh damping
    double factThis = 0.0;
    if (addRayleigh == 1)  {
        theMatrix = this->Element::getDamp();
        factThis = 1.0;
    }
    
    // now add damping tangent from materials
    static Matrix cb(3,3);
    cb.Zero();
    cb(0,0) = theMaterials[0]->getDampTangent();
    cb(2,2) = theMaterials[1]->getDampTangent();
    
    // transform from basic to local system
    static Matrix cl(6,6);
    cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);
    
    // transform from local to global system and add to cg
    theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWen2d::getMass()
{
    // zero the matrix
    theMatrix.Zero();
    
    // check for quick return
    if (mass == 0.0)  {
        return theMatrix;
    }
    
    double m = 0.5*mass;
    for (int i=0; i<2; i++)  {
        theMatrix(i,i)     = m;
        theMatrix(i+3,i+3) = m;
    }
    
    return theMatrix;
}


void ElastomericBearingBoucWen2d::zeroLoad()
{
    theLoad.Zero();
}


int ElastomericBearingBoucWen2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"ElastomericBearingBoucWen2d::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << ".\n";
    
    return -1;
}


int ElastomericBearingBoucWen2d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)
        return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (3 != Raccel1.Size() || 3 != Raccel2.Size())  {
        opserr << "ElastomericBearingBoucWen2d::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible.\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5*mass;
    for (int i=0; i<2; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+3) -= m * Raccel2(i);
    }
    
    return 0;
}


const Vector& ElastomericBearingBoucWen2d::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(6);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
    
    // add P-Delta moments to local forces
    double kGeo1 = 0.5*qb(0);
    double MpDelta1 = kGeo1*(ul(4)-ul(1));
    ql(2) += MpDelta1;
    ql(5) += MpDelta1;
    double MpDelta2 = kGeo1*shearDistI*L*ul(2);
    ql(2) += MpDelta2;
    ql(5) -= MpDelta2;
    double MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(5);
    ql(2) -= MpDelta3;
    ql(5) += MpDelta3;
    
    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    return theVector;
}


const Vector& ElastomericBearingBoucWen2d::getResistingForceIncInertia()
{
    // this already includes damping forces from materials
    theVector = this->getResistingForce();
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    // add the damping forces from rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add inertia forces from element mass
    if (mass != 0.0)  {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();
        
        double m = 0.5*mass;
        for (int i=0; i<2; i++)  {
            theVector(i)   += m * accel1(i);
            theVector(i+3) += m * accel2(i);
        }
    }
    
    return theVector;
}


int ElastomericBearingBoucWen2d::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(21);
    data(0) = this->getTag();
    data(1) = k0;
    data(2) = qYield;
    data(3) = k2;
    data(4) = k3;
    data(5) = mu;
    data(6) = eta;
    data(7) = beta;
    data(8) = gamma;
    data(9) = A;
    data(10) = shearDistI;
    data(11) = addRayleigh;
    data(12) = mass;
    data(13) = maxIter;
    data(14) = tol;
    data(15) = x.Size();
    data(16) = y.Size();
    data(17) = alphaM;
    data(18) = betaK;
    data(19) = betaK0;
    data(20) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send the material class tags
    ID matClassTags(2);
    for (int i=0; i<2; i++)
        matClassTags(i) = theMaterials[i]->getClassTag();
    sChannel.sendID(0, commitTag, matClassTags);
    
    // send the material models
    for (int i=0; i<2; i++)
        theMaterials[i]->sendSelf(commitTag, sChannel);
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    
    return 0;
}


int ElastomericBearingBoucWen2d::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete material memory
    for (int i=0; i<2; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
    
    // receive element parameters
    static Vector data(21);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    k0 = data(1);
    qYield = data(2);
    k2 = data(3);
    k3 = data(4);
    mu = data(5);
    eta = data(6);
    beta = data(7);
    gamma = data(8);
    A = data(9);
    shearDistI = data(10);
    addRayleigh = (int)data(11);
    mass = data(12);
    maxIter = (int)data(13);
    tol = data(14);
    alphaM = data(17);
    betaK = data(18);
    betaK0 = data(19);
    betaKc = data(20);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive the material class tags
    ID matClassTags(2);
    rChannel.recvID(0, commitTag, matClassTags);
    
    // receive the material models
    for (int i=0; i<2; i++)  {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "ElastomericBearingBoucWen2d::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
        theMaterials[i]->recvSelf(commitTag, rChannel, theBroker);
    }
    
    // receive remaining data
    if ((int)data(15) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(16) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    onP0 = false;
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = A*k0 + k2;
    kbInit(2,2) = theMaterials[1]->getInitialTangent();
    
    // initialize other variables
    this->revertToStart();
    
    return 0;
}


int ElastomericBearingBoucWen2d::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)					     
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void ElastomericBearingBoucWen2d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln; 
        s << "  type: ElastomericBearingBoucWen2d\n";
        s << "  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  k0: " << k0 << "  qYield: " << qYield << "  k2: " << k2 << endln;
        s << "  k3: " << k3 << "  mu: " << mu << endln;
        s << "  eta: " << eta << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  Material ux: " << theMaterials[0]->getTag();
        s << "  Material rz: " << theMaterials[1]->getTag() << endln;
        s << "  shearDistI: " << shearDistI << "  addRayleigh: "
            << addRayleigh << "  mass: " << mass << endln;
        s << "  maxIter: " << maxIter << "  tol: " << tol << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ElastomericBearingBoucWen2d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"k0\": " << k0 << ", ";
        s << "\"qYield\": " << qYield << ", ";
        s << "\"k2\": " << k2 << ", ";
        s << "\"k3\": " << k3 << ", ";
        s << "\"mu\": " << mu << ", ";
        s << "\"eta\": " << eta << ", ";
        s << "\"beta\": " << beta << ", ";
        s << "\"gamma\": " << gamma << ", ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\"], ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"addRayleigh\": " << addRayleigh << ", ";
        s << "\"mass\": " << mass << "}";
    }
}


Response* ElastomericBearingBoucWen2d::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ElastomericBearingBoucWen2d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
        output.tag("ResponseType","N_1");
        output.tag("ResponseType","V_1");
        output.tag("ResponseType","M_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","V_2");
        output.tag("ResponseType","M_2");
        
        theResponse = new ElementResponse(this, 2, theVector);
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
        output.tag("ResponseType","qb1");
        output.tag("ResponseType","qb2");
        output.tag("ResponseType","qb3");
        
        theResponse = new ElementResponse(this, 3, Vector(3));
    }
    // local displacements
    else if (strcmp(argv[0],"localDisplacement") == 0 ||
        strcmp(argv[0],"localDisplacements") == 0)
    {
        output.tag("ResponseType","ux_1");
        output.tag("ResponseType","uy_1");
        output.tag("ResponseType","rz_1");
        output.tag("ResponseType","ux_2");
        output.tag("ResponseType","uy_2");
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
        
        theResponse = new ElementResponse(this, 5, Vector(3));
    }
    // hysteretic evolution parameter
    else if (strcmp(argv[0],"hystereticParameter") == 0 || strcmp(argv[0],"hystParameter") == 0 || 
        strcmp(argv[0],"hystereticParam") == 0 || strcmp(argv[0],"hystParam") == 0 ||
        strcmp(argv[0],"z") == 0)
    {
        output.tag("ResponseType","z");
        
        theResponse = new ElementResponse(this, 6, z);
    }
    // dzdu
    else if (strcmp(argv[0],"dzdu") == 0)
    {
        output.tag("ResponseType","dzdu");
        
        theResponse = new ElementResponse(this, 7, dzdu);
    }
    // basic stiffness
    else if (strcmp(argv[0],"kb") == 0 || strcmp(argv[0],"basicStiff") == 0 ||
        strcmp(argv[0],"basicStiffness") == 0)
    {
        output.tag("ResponseType","kb22");
        
        theResponse = new ElementResponse(this, 8, k0);
    }
    // material output
    else if (strcmp(argv[0],"material") == 0)  {
        if (argc > 2)  {
            int matNum = atoi(argv[1]);
            if (matNum >= 1 && matNum <= 2)
                theResponse =  theMaterials[matNum-1]->setResponse(&argv[2], argc-2, output);
        }
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int ElastomericBearingBoucWen2d::getResponse(int responseID, Information &eleInfo)
{
    double kGeo1, MpDelta1, MpDelta2, MpDelta3;
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta moments
        kGeo1 = 0.5*qb(0);
        MpDelta1 = kGeo1*(ul(4)-ul(1));
        theVector(2) += MpDelta1;
        theVector(5) += MpDelta1;
        MpDelta2 = kGeo1*shearDistI*L*ul(2);
        theVector(2) += MpDelta2;
        theVector(5) -= MpDelta2;
        MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(5);
        theVector(2) -= MpDelta3;
        theVector(5) += MpDelta3;
        
        return eleInfo.setVector(theVector);
        
    case 3:  // basic forces
        return eleInfo.setVector(qb);
        
    case 4:  // local displacements
        return eleInfo.setVector(ul);
        
    case 5:  // basic displacements
        return eleInfo.setVector(ub);
        
    case 6:  // hysteretic evolution parameter
        return eleInfo.setDouble(z);
        
    case 7:  // dzdu
        return eleInfo.setDouble(dzdu);
        
    case 8:  // basic stiffness
        return eleInfo.setDouble(kb(1,1));
        
    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void ElastomericBearingBoucWen2d::setUp()
{
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    if (L > DBL_EPSILON)  {
        if (x.Size() == 0)  {
            x.resize(3);
            x(0) = xp(0);  x(1) = xp(1);  x(2) = 0.0;
            y.resize(3);
            y(0) = -x(1);  y(1) = x(0);  y(2) = 0.0;
        } else if (onP0)  {
            opserr << "WARNING ElastomericBearingBoucWen2d::setUp() - " 
                << "element: " << this->getTag()
                << " - ignoring nodes and using specified "
                << "local x vector to determine orientation.\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "ElastomericBearingBoucWen2d::setUp() - "
            << "element: " << this->getTag()
            << " - incorrect dimension of orientation vectors.\n";
        exit(-1);
    }
    
    // establish orientation of element for the transformation matrix
    // z = x cross y
    static Vector z(3);
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
        opserr << "ElastomericBearingBoucWen2d::setUp() - "
            << "element: " << this->getTag()
            << " - invalid orientation vectors.\n";
        exit(-1);
    }
    
    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = x(0)/xn;
    Tgl(0,1) = Tgl(3,4) = x(1)/xn;
    Tgl(1,0) = Tgl(4,3) = y(0)/yn;
    Tgl(1,1) = Tgl(4,4) = y(1)/yn;
    Tgl(2,2) = Tgl(5,5) = z(2)/zn;
    
    // create transformation matrix from local to basic system (linear)
    Tlb.Zero();
    Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = -1.0;
    Tlb(0,3) = Tlb(1,4) = Tlb(2,5) = 1.0;
    Tlb(1,2) = -shearDistI*L;
    Tlb(1,5) = -(1.0 - shearDistI)*L;
}


double ElastomericBearingBoucWen2d::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
