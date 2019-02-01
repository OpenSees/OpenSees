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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the implementation of the
// FlatSliderSimple3d class.

#include "FlatSliderSimple3d.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <FrictionModel.h>
#include <UniaxialMaterial.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

void* OPS_FlatSliderSimple3d()
{
    int ndf = OPS_GetNDF();
    if (ndf != 6)  {
	opserr << "WARNING invalid ndf: " << ndf;
	opserr << ", for space problem need 6 - flatSliderBearing \n";    
	return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 13) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: flatSliderBearing eleTag iNode jNode frnMdlTag kInit -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
	return 0;
    }

    // tags
    int idata[4];
    int num = 4;
    if (OPS_GetIntInput(&num, idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    FrictionModel* theFrnMdl = OPS_getFrictionModel(idata[3]);
    if (theFrnMdl == 0) {
	opserr << "WARNING friction model not found\n";
	opserr << "frictionModel: " << idata[3] << endln;
	return 0;
    }

    // data
    double kInit;
    num = 1;
    if (OPS_GetDoubleInput(&num, &kInit) < 0) {
	opserr<<"WARNING: invalid double kInit\n";
	return 0;
    }

    // materials
    UniaxialMaterial* mats[4] = {0,0,0,0};
    const char* type = OPS_GetString();
    if (strcmp(type,"-P") != 0) {
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
    if (strcmp(type,"-T") != 0) {
	opserr<<"WARNING: want -T\n";
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

    type = OPS_GetString();
    if (strcmp(type,"-My") != 0) {
	opserr<<"WARNING: want -My\n";
	return 0;
    }
    num = 1;
    if (OPS_GetIntInput(&num, &matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }
    mats[2] = OPS_getUniaxialMaterial(matTag);
    if (mats[2] == 0) {
	opserr<<"WARNING: material not found\n";
	return 0;
    }

    type = OPS_GetString();
    if (strcmp(type,"-Mz") != 0) {
	opserr<<"WARNING: want -Mz\n";
	return 0;
    }
    num = 1;
    if (OPS_GetIntInput(&num, &matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }
    mats[3] = OPS_getUniaxialMaterial(matTag);
    if (mats[3] == 0) {
	opserr<<"WARNING: material not found\n";
	return 0;
    }

    // options
    Vector x,y(3);
    y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
    double sDistI = 0.0;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1e-12;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	type = OPS_GetString();
	if (strcmp(type,"-orient") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 3) {
		opserr<<"WARNING: insufficient arguments after -orient\n";
		return 0;
	    }
	    num = 3;
	    x.resize(3);
	    if (OPS_GetDoubleInput(&num, &x(0)) < 0) {
		opserr<<"WARNING: invalid orient value\n";
		return 0;
	    }
	    if (OPS_GetNumRemainingInputArgs() < 3) {
		y = x;
		x = Vector();
		continue;
	    }
	    y.resize(3);
	    if (OPS_GetDoubleInput(&num, &y(0)) < 0) {
		y = x;
		x = Vector();
		continue;
	    }
	} else if (strcmp(type,"-shearDist") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return 0;
	    }
	    num = 1;
	    if (OPS_GetDoubleInput(&num, &sDistI) < 0) {
		opserr<<"WARNING: invalid shearDist\n";
		return 0;
	    }
	} else if (strcmp(type,"-doRayleigh") == 0) {
	    doRayleigh = 1;
	} else if (strcmp(type,"-mass") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"WARNING: insufficient args\n";
		return 0;
	    }
	    num = 1;
	    if (OPS_GetDoubleInput(&num, &mass) < 0) {
		opserr<<"WARNING: invalid mass\n";
		return 0;
	    }
	} else if (strcmp(type,"-iter") == 0) {
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

    return new FlatSliderSimple3d(idata[0],idata[1],idata[2],*theFrnMdl,
				  kInit,mats,y,x,sDistI,doRayleigh,mass,
				  maxIter,tol);
}


// initialize the class wide variables
Matrix FlatSliderSimple3d::theMatrix(12,12);
Vector FlatSliderSimple3d::theVector(12);


FlatSliderSimple3d::FlatSliderSimple3d(int tag, int Nd1, int Nd2,
    FrictionModel &thefrnmdl, double kInit, UniaxialMaterial **materials,
    const Vector _y, const Vector _x, double sdI, int addRay, double m,
    int maxiter, double _tol, double kfactuplift)
    : Element(tag, ELE_TAG_FlatSliderSimple3d),
    connectedExternalNodes(2), theFrnMdl(0), k0(kInit), x(_x), y(_y),
    shearDistI(sdI), addRayleigh(addRay), mass(m), maxIter(maxiter), tol(_tol),
    kFactUplift(kfactuplift),
    L(0.0), onP0(true), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // get a copy of the friction model
    theFrnMdl = thefrnmdl.getCopy();
    if (!theFrnMdl)  {
        opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - element: "
            << this->getTag() << " - failed to get copy of the "
            << "friction model.\n";
        exit(-1);
    }
    
    // check material input
    if (materials == 0)  {
        opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - "
            << "null material array passed.\n";
        exit(-1);
    }
    
    // get copies of the uniaxial materials
    for (int i=0; i<4; i++)  {
        if (materials[i] == 0)  {
            opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0)  {
            opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = kbInit(2,2) = k0;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
    
    // initialize other variables
    this->revertToStart();
}


FlatSliderSimple3d::FlatSliderSimple3d()
    : Element(0, ELE_TAG_FlatSliderSimple3d),
    connectedExternalNodes(2), theFrnMdl(0), k0(0.0), x(0), y(0),
    shearDistI(0.0), addRayleigh(0), mass(0.0), maxIter(25), tol(1E-12),
    kFactUplift(1E-12),
    L(0.0), onP0(false), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "FlatSliderSimple3d::FlatSliderSimple3d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // set material pointers to NULL
    for (int i=0; i<4; i++)
        theMaterials[i] = 0;
}


FlatSliderSimple3d::~FlatSliderSimple3d()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theFrnMdl)
        delete theFrnMdl;
    
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int FlatSliderSimple3d::getNumExternalNodes() const
{
    return 2;
}


const ID& FlatSliderSimple3d::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** FlatSliderSimple3d::getNodePtrs()
{
    return theNodes;
}


int FlatSliderSimple3d::getNumDOF()
{
    return 12;
}


void FlatSliderSimple3d::setDomain(Domain *theDomain)
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
            opserr << "WARNING FlatSliderSimple3d::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING FlatSliderSimple3d::setDomain() - Nd2: "
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
    if (dofNd1 != 6)  {
        opserr << "FlatSliderSimple3d::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "FlatSliderSimple3d::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int FlatSliderSimple3d::commitState()
{
    int errCode = 0;
    
    // commit trial history variables
    ubPlasticC = ubPlastic;
    
    // commit friction model
    errCode += theFrnMdl->commitState();
    
    // commit material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->commitState();
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int FlatSliderSimple3d::revertToLastCommit()
{
    int errCode = 0;
    
    // revert friction model
    errCode += theFrnMdl->revertToLastCommit();
    
    // revert material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->revertToLastCommit();
    
    return errCode;
}


int FlatSliderSimple3d::revertToStart()
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
    
    // revert friction model
    errCode += theFrnMdl->revertToStart();
    
    // revert material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->revertToStart();
    
    return errCode;
}


int FlatSliderSimple3d::update()
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
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
    
    // get absolute velocity
    double ubdotAbs = sqrt(pow(ubdot(1),2) + pow(ubdot(2),2));
    
    // 1) get axial force and stiffness in basic x-direction
    double ub0Old = theMaterials[0]->getStrain();
    theMaterials[0]->setTrialStrain(ub(0), ubdot(0));
    qb(0) = theMaterials[0]->getStress();
    kb(0,0) = theMaterials[0]->getTangent();
    
    // check for uplift
    if (qb(0) >= 0.0)  {
        // update plastic displacements
        ubPlastic(0) = ub(1);
        ubPlastic(1) = ub(2);
        // set basic forces
        qb.Zero();
        // set tangent stiffnesses
        kb = kbInit;
        if (qb(0) > 0.0)  {
            theMaterials[0]->setTrialStrain(ub0Old, 0.0);
            kb = kFactUplift*kbInit;  // kb = DBL_EPSILON*kbInit;
        }
        return 0;
    }
    
    // 2) calculate shear forces and stiffnesses in basic y- and z-direction
    int iter = 0;
    Vector qbOld(2);
    do  {
        // save old shear forces
        iter++;
        qbOld(0) = qb(1);
        qbOld(1) = qb(2);
        
        // get normal and friction (yield) forces
        double N = -qb(0) - qb(1)*ul(5) + qb(2)*ul(4);
        N = N > 0.0 ? N : 0.0;  // can not be negative
        theFrnMdl->setTrial(N, ubdotAbs);
        double qYield = (theFrnMdl->getFrictionForce());
        
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
            qb(1) = qTrial(0) - N*ul(5);
            qb(2) = qTrial(1) + N*ul(4);
            // set tangent stiffnesses
            kb(1,1) = kb(2,2) = k0;
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
            qb(1) = qYield*qTrial(0)/qTrialNorm - N*ul(5);
            qb(2) = qYield*qTrial(1)/qTrialNorm + N*ul(4);
            // set tangent stiffnesses
            double D = pow(qTrialNorm,3);
            kb(1,1) = qYield*k0*qTrial(1)*qTrial(1)/D;
            kb(1,2) = kb(2,1) = -qYield*k0*qTrial(1)*qTrial(0)/D;
            kb(2,2) = qYield*k0*qTrial(0)*qTrial(0)/D;
        }
        
    } while ((sqrt(pow(qb(1)-qbOld(0),2)+pow(qb(2)-qbOld(1),2)) >= tol) && (iter <= maxIter));
    
    // issue warning if iteration did not converge
    if (iter >= maxIter)   {
        opserr << "WARNING: FlatSliderSimple3d::update() - element: "
            << this->getTag() << " - did not find the shear force after "
            << iter << " iterations and norm: "
            << sqrt(pow(qb(1)-qbOld(0),2)+pow(qb(2)-qbOld(1),2)) << ".\n";
        return -1;
    }
    
    // 3) get moment and stiffness in basic x-direction
    theMaterials[1]->setTrialStrain(ub(3), ubdot(3));
    qb(3) = theMaterials[1]->getStress();
    kb(3,3) = theMaterials[1]->getTangent();
    
    // 4) get moment and stiffness in basic y-direction
    theMaterials[2]->setTrialStrain(ub(4), ubdot(4));
    qb(4) = theMaterials[2]->getStress();
    kb(4,4) = theMaterials[2]->getTangent();
    
    // 5) get moment and stiffness in basic z-direction
    theMaterials[3]->setTrialStrain(ub(5), ubdot(5));
    qb(5) = theMaterials[3]->getStress();
    kb(5,5) = theMaterials[3]->getTangent();
    
    return 0;
}


const Matrix& FlatSliderSimple3d::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    double Ls = (1.0 - shearDistI)*L;
    // add P-Delta moment stiffness terms
    kl(5,1)   -= qb(0);
    kl(5,7)   += qb(0);
    kl(5,11)  -= qb(0)*Ls;
    kl(11,11) += qb(0)*Ls;
    kl(4,2)   += qb(0);
    kl(4,8)   -= qb(0);
    kl(4,10)  -= qb(0)*Ls;
    kl(10,10) += qb(0)*Ls;
    // add V-Delta torsion stiffness terms
    kl(3,1)   += qb(2);
    kl(3,2)   -= qb(1);
    kl(3,7)   -= qb(2);
    kl(3,8)   += qb(1);
    kl(3,10)  += qb(1)*Ls;
    kl(3,11)  += qb(2)*Ls;
    kl(9,10)  -= qb(1)*Ls;
    kl(9,11)  -= qb(2)*Ls;
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& FlatSliderSimple3d::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix klInit(12,12);
    klInit.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);
    
    return theMatrix;
}


const Matrix& FlatSliderSimple3d::getDamp()
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
    static Matrix cb(6,6);
    cb.Zero();
    cb(0,0) = theMaterials[0]->getDampTangent();
    cb(3,3) = theMaterials[1]->getDampTangent();
    cb(4,4) = theMaterials[2]->getDampTangent();
    cb(5,5) = theMaterials[3]->getDampTangent();
    
    // transform from basic to local system
    static Matrix cl(12,12);
    cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);
    
    // transform from local to global system and add to cg
    theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    
    return theMatrix;
}


const Matrix& FlatSliderSimple3d::getMass()
{
    // zero the matrix
    theMatrix.Zero();
    
    // check for quick return
    if (mass == 0.0)  {
        return theMatrix;
    }    
    
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theMatrix(i,i)     = m;
        theMatrix(i+6,i+6) = m;
    }
    
    return theMatrix; 
}


void FlatSliderSimple3d::zeroLoad()
{
    theLoad.Zero();
}


int FlatSliderSimple3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"FlatSliderSimple3d::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << ".\n";
    
    return -1;
}


int FlatSliderSimple3d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
        opserr << "FlatSliderSimple3d::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible.\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+6) -= m * Raccel2(i);
    }
    
    return 0;
}


const Vector& FlatSliderSimple3d::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(12);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
    
    // add P-Delta moments to local forces
    double MpDelta1 = qb(0)*(ul(7)-ul(1));
    ql(5)  += MpDelta1;
    double MpDelta2 = qb(0)*(1.0 - shearDistI)*L*ul(11);
    ql(5)  -= MpDelta2;
    ql(11) += MpDelta2;
    double MpDelta3 = qb(0)*(ul(8)-ul(2));
    ql(4)  -= MpDelta3;
    double MpDelta4 = qb(0)*(1.0 - shearDistI)*L*ul(10);
    ql(4)  -= MpDelta4;
    ql(10) += MpDelta4;
    
    // add V-Delta torsion to local forces
    double Vdelta1 = qb(1)*(ul(8)-ul(2)) - qb(2)*(ul(7)-ul(1));
    ql(3) += Vdelta1;
    double Vdelta2 = (1.0 - shearDistI)*L*(qb(1)*ul(10) + qb(2)*ul(11));
    ql(3) += Vdelta2;
    ql(9) -= Vdelta2;
    
    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    return theVector;
}


const Vector& FlatSliderSimple3d::getResistingForceIncInertia()
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
        for (int i=0; i<3; i++)  {
            theVector(i)   += m * accel1(i);
            theVector(i+6) += m * accel2(i);
        }
    }
    
    return theVector;
}


int FlatSliderSimple3d::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(14);
    data(0)  = this->getTag();
    data(1)  = k0;
    data(2)  = shearDistI;
    data(3)  = addRayleigh;
    data(4)  = mass;
    data(5)  = maxIter;
    data(6)  = tol;
    data(7)  = kFactUplift;
    data(8)  = x.Size();
    data(9)  = y.Size();
    data(10) = alphaM;
    data(11) = betaK;
    data(12) = betaK0;
    data(13) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send the friction model class tag
    ID frnClassTag(1);
    frnClassTag(0) = theFrnMdl->getClassTag();
    sChannel.sendID(0, commitTag, frnClassTag);
    
    // send the friction model
    theFrnMdl->sendSelf(commitTag, sChannel);
    
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


int FlatSliderSimple3d::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete material memory
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
    
    // receive element parameters
    static Vector data(14);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    k0 = data(1);
    shearDistI = data(2);
    addRayleigh = (int)data(3);
    mass = data(4);
    maxIter = (int)data(5);
    tol = data(6);
    kFactUplift = data(7);
    alphaM = data(10);
    betaK = data(11);
    betaK0 = data(12);
    betaKc = data(13);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive the friction model class tag
    ID frnClassTag(1);
    rChannel.recvID(0, commitTag, frnClassTag);
    
    // receive the friction model
    theFrnMdl = theBroker.getNewFrictionModel(frnClassTag(0));
    if (theFrnMdl == 0) {
        opserr << "FlatSliderSimple3d::recvSelf() - "
            << "failed to get blank friction model.\n";
        return -1;
    }
    theFrnMdl->recvSelf(commitTag, rChannel, theBroker);
    
    // receive the material class tags
    ID matClassTags(4);
    rChannel.recvID(0, commitTag, matClassTags);
    
    // receive the material models
    for (int i=0; i<4; i++)  {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "FlatSliderSimple3d::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
        theMaterials[i]->recvSelf(commitTag, rChannel, theBroker);
    }
    
    // receive remaining data
    if ((int)data(8) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(9) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    onP0 = false;
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = kbInit(2,2) = k0;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
    
    // initialize variables
    this->revertToStart();
    
    return 0;
}


int FlatSliderSimple3d::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    int errCode = 0;
    
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;
    
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    
    if (displayMode >= 0)  {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i=0; i<3; i++)  {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v3(i) = end2Crd(i) + end2Disp(i)*fact;
        }
        v2(0) = end1Crd(0) + (end2Disp(0) + xp(1)*end2Disp(5) - xp(2)*end2Disp(4))*fact;
        v2(1) = end1Crd(1) + (end2Disp(1) - xp(0)*end2Disp(5) + xp(2)*end2Disp(3))*fact;
        v2(2) = end1Crd(2) + (end2Disp(2) + xp(0)*end2Disp(4) - xp(1)*end2Disp(3))*fact;
    } else  {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        
        if (eigen1.noCols() >= mode)  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v3(i) = end2Crd(i) + eigen2(i,mode-1)*fact;
            }
            v2(0) = end1Crd(0) + (eigen2(0,mode-1) + xp(1)*eigen2(5,mode-1) - xp(2)*eigen2(4,mode-1))*fact;
            v2(1) = end1Crd(1) + (eigen2(1,mode-1) - xp(0)*eigen2(5,mode-1) + xp(2)*eigen2(3,mode-1))*fact;
            v2(2) = end1Crd(2) + (eigen2(2,mode-1) + xp(0)*eigen2(4,mode-1) - xp(1)*eigen2(3,mode-1))*fact;
        } else  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i);
                v2(i) = end1Crd(i);
                v3(i) = end2Crd(i);
            }
        }
    }
    
    errCode += theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag(), 0);
    errCode += theViewer.drawLine (v2, v3, 1.0, 1.0, this->getTag(), 0);
    
    return errCode;
}


void FlatSliderSimple3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag(); 
        s << "  type: FlatSliderSimple3d  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  FrictionModel: " << theFrnMdl->getTag() << endln;
        s << "  kInit: " << k0 << endln;
        s << "  Material ux: " << theMaterials[0]->getTag() << endln;
        s << "  Material rx: " << theMaterials[1]->getTag() << endln;
        s << "  Material ry: " << theMaterials[2]->getTag() << endln;
        s << "  Material rz: " << theMaterials[3]->getTag() << endln;
        s << "  shearDistI: " << shearDistI << "  addRayleigh: "
            << addRayleigh << "  mass: " << mass << endln;
        s << "  maxIter: " << maxIter << "  tol: " << tol << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FlatSliderSimple3d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"frictionModel\": \"" << theFrnMdl->getTag() << "\", ";
        s << "\"kInit\": " << k0 << ", ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\", \"";
        s << theMaterials[2]->getTag() << "\", \"";
        s << theMaterials[3]->getTag() << "\"], ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"addRayleigh\": " << addRayleigh << ", ";
        s << "\"mass\": " << mass << ", ";
        s << "\"maxIter\": " << maxIter << ", ";
        s << "\"tol\": " << tol << "}";
    }
}


Response* FlatSliderSimple3d::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","FlatSliderSimple3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 ||
        strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0)
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
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
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
    else if (strcmp(argv[0],"basicForce") == 0 ||
        strcmp(argv[0],"basicForces") == 0)
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
    else if (strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || 
        strcmp(argv[0],"basicDeformation") == 0 ||
        strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 ||
        strcmp(argv[0],"basicDisplacements") == 0)
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
                theResponse = theMaterials[matNum-1]->setResponse(&argv[2], argc-2, output);
        }
    }
    // friction model output
    else if (strcmp(argv[0],"frictionModel") == 0 || strcmp(argv[0],"frnMdl") == 0 ||
        strcmp(argv[0],"frictionMdl") == 0 || strcmp(argv[0],"frnModel") == 0)  {
            if (argc > 1)
                theResponse = theFrnMdl->setResponse(&argv[1], argc-1, output);
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int FlatSliderSimple3d::getResponse(int responseID, Information &eleInfo)
{
    double MpDelta1, MpDelta2, MpDelta3, MpDelta4, Vdelta1, Vdelta2;
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta moments
        MpDelta1 = qb(0)*(ul(7)-ul(1));
        theVector(5)  += MpDelta1;
        MpDelta2 = qb(0)*(1.0 - shearDistI)*L*ul(11);
        theVector(5)  -= MpDelta2;
        theVector(11) += MpDelta2;
        MpDelta3 = qb(0)*(ul(8)-ul(2));
        theVector(4)  -= MpDelta3;
        MpDelta4 = qb(0)*(1.0 - shearDistI)*L*ul(10);
        theVector(4)  -= MpDelta4;
        theVector(10) += MpDelta4;
        // add V-Delta torsion
        Vdelta1 = qb(1)*(ul(8)-ul(2)) - qb(2)*(ul(7)-ul(1));
        theVector(3)  += Vdelta1;
        Vdelta2 = (1.0 - shearDistI)*L*(qb(1)*ul(10) + qb(2)*ul(11));
        theVector(3)  += Vdelta2;
        theVector(9)  -= Vdelta2;
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


// Establish the external nodes and set up the transformation matrix for orientation
void FlatSliderSimple3d::setUp()
{
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    if (L > DBL_EPSILON)  {
        if (x.Size() == 0)  {
            x.resize(3);
            x = xp;
        } else if (onP0)  {
            opserr << "WARNING FlatSliderSimple3d::setUp() - "
                << "element: " << this->getTag()
                << " - ignoring nodes and using specified "
                << "local x vector to determine orientation.\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "FlatSliderSimple3d::setUp() - "
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
        opserr << "FlatSliderSimple3d::setUp() - "
            << "element: " << this->getTag()
            << " - invalid orientation vectors.\n";
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
    Tlb(1,5) = -shearDistI*L;
    Tlb(1,11) = -(1.0 - shearDistI)*L;
    Tlb(2,4) = -Tlb(1,5);
    Tlb(2,10) = -Tlb(1,11);
}
