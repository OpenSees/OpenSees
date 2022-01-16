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
// Created: 03/19
// Revision: A
//
// Description: This file contains the implementation of the LinearElasticSpring class.

#include <LinearElasticSpring.h>
#include <Information.h>

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
#include <vector>

// initialize the class wide variables
Matrix LinearElasticSpring::LinearElasticSpringM2(2,2);
Matrix LinearElasticSpring::LinearElasticSpringM4(4,4);
Matrix LinearElasticSpring::LinearElasticSpringM6(6,6);
Matrix LinearElasticSpring::LinearElasticSpringM12(12,12);
Vector LinearElasticSpring::LinearElasticSpringV2(2);
Vector LinearElasticSpring::LinearElasticSpringV4(4);
Vector LinearElasticSpring::LinearElasticSpringV6(6);
Vector LinearElasticSpring::LinearElasticSpringV12(12);

void * OPS_ADD_RUNTIME_VPV(OPS_LinearElasticSpring)
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (OPS_GetNumRemainingInputArgs() < 7) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: linearElasticSpring eleTag iNode jNode -dir dirs -stif kb <-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-doRayleigh> <-damp cb>\n";
        return 0;
    }
    
    // tags
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "WARNING: invalid integer data\n";
        return 0;
    }
    
    // dirs
    const char* type = OPS_GetString();
    if (strcmp(type, "-dir") != 0 && strcmp(type, "-dof") != 0) {
        opserr << "WARNING expecting -dir dirs\n";
        return 0;
    }
    ID dirs(ndf);
    int numDIR = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        int dir;
        numdata = 1;
        int numArgs = OPS_GetNumRemainingInputArgs();
        if (OPS_GetIntInput(&numdata, &dir) < 0) {
            if (numArgs > OPS_GetNumRemainingInputArgs()) {
                // move current arg back by one
                OPS_ResetCurrentInputArg(-1);
            }
            break;
        }
        if (dir < 1 || ndf < dir) {
            opserr << "WARNING invalid direction ID\n";
            return 0;
        }
        dirs(numDIR++) = dir - 1;
    }
    dirs.resize(numDIR);
    
    // stiffness matrix terms
    type = OPS_GetString();
    if (strcmp(type, "-stif") != 0 && strcmp(type, "-stiff") != 0) {
        opserr << "WARNING expecting -stif kij\n";
        return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < numDIR*numDIR) {
        opserr << "WARNING wrong number of kij specified\n";
        return 0;
    }
    numdata = 1;
    Matrix kb(numDIR, numDIR);
    for (int i = 0; i < numDIR; i++) {
        for (int j = 0; j < numDIR; j++) {
            if (OPS_GetDoubleInput(&numdata, &kb(i, j)) < 0) {
                opserr << "WARNING invalid stiffness value\n";
                    return 0;
            }
        }
    }
    
    // options
    Vector x, y, Mratio;
    int doRayleigh = 0;
    Matrix *cb = 0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        return new LinearElasticSpring(idata[0], ndm, idata[1], idata[2],
            dirs, kb);
    }
    
    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        if (strcmp(type, "-orient") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING: insufficient arguments after -orient\n";
                return 0;
            }
            numdata = 3;
            x.resize(3);
            if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
                opserr << "WARNING: invalid -orient values\n";
                return 0;
            }
            if (OPS_GetNumRemainingInputArgs() < 3) {
                y = x;
                x = Vector();
                continue;
            }
            y.resize(3);
            if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
                y = x;
                x = Vector();
                continue;
            }
        }
        else if (strcmp(type, "-pDelta") == 0) {
            Mratio.resize(4);
            Mratio.Zero();
            numdata = 4;
            double* ptr = &Mratio(0);
            if (ndm == 2) {
                numdata = 2;
                ptr += 2;
            }
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING: insufficient data for -pDelta\n";
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, ptr) < 0) {
                opserr << "WARNING: invalid -pDelta value\n";
                return 0;
            }
        }
        else if (strcmp(type, "-doRayleigh") == 0) {
            doRayleigh = 1;
        }
        else if (strcmp(type, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < numDIR*numDIR) {
                opserr << "WARNING wrong number of cij specified\n";
                return 0;
            }
            numdata = 1;
            double cij;
            cb = new Matrix(numDIR, numDIR);
            for (int i = 0; i < numDIR; i++) {
                for (int j = 0; j < numDIR; j++) {
                    if (OPS_GetDoubleInput(&numdata, &cij) < 0) {
                        opserr << "WARNING invalid damping value\n";
                        delete cb;
                        return 0;
                    }
                    (*cb)(i, j) = cij;
                }
            }
        }
    }
    
    // create object
    Element *theEle = new LinearElasticSpring(idata[0], ndm, idata[1], idata[2],
        dirs, kb, y, x, Mratio,
        doRayleigh, cb);
    
    // clean up memory
    delete cb;
    
    return theEle;
}


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
LinearElasticSpring::LinearElasticSpring(int tag, int dim,
    int Nd1, int Nd2, const ID &direction, const Matrix &stif,
    const Vector _y, const Vector _x, const Vector Mr,
    int addRay, const Matrix *damp)
    : Element(tag, ELE_TAG_LinearElasticSpring),
    numDIM(dim), numDOF(0), connectedExternalNodes(2),
    numDIR(direction.Size()), dir(direction), kb(stif), cb(0),
    x(_x), y(_y), Mratio(Mr), addRayleigh(addRay), L(0.0), onP0(true),
    trans(3, 3), ub(0), ubdot(0), qb(0), ul(0),
    Tgl(0,0), Tlb(0,0), theMatrix(0), theVector(0), theLoad(0)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "LinearElasticSpring::LinearElasticSpring() - element: "
            << this->getTag() << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // check the number of directions
    if (numDIR < 1 || numDIR > 6)  {
        opserr << "LinearElasticSpring::LinearElasticSpring() - element: "
            << this->getTag() << " wrong number of directions\n";
        exit(-1);
    }
    
    // initialize directions and check for valid values
    for (int i=0; i< numDIR; i++)  {
        if (dir(i) < 0 ||
            (numDIM == 1 && dir(i) > 0) ||
            (numDIM == 2 && dir(i) > 2) ||
            (numDIM == 3 && dir(i) > 5))  {
            opserr << "LinearElasticSpring::LinearElasticSpring() - "
                << "incorrect direction " << dir(i)
                << " is set to 0\n";
            dir(i) = 0;
        }
    }
    
    // check p-delta moment distribution ratios
    if (Mratio.Size() == 4) {
        if (Mratio(0) < 0.0 || Mratio(1) < 0.0 ||
            Mratio(2) < 0.0 || Mratio(3) < 0.0) {
            opserr << "LinearElasticSpring::LinearElasticSpring() - "
                << "p-delta moment ratios can not be negative\n";
            exit(-1);
        }
        if (Mratio(0) + Mratio(1) > 1.0) {
            opserr << "LinearElasticSpring::LinearElasticSpring() - "
                << "incorrect p-delta moment ratios:\nrMy1 + rMy2 = "
                << Mratio(0) + Mratio(1) << " > 1.0\n";
            exit(-1);
        }
        if (Mratio(2) + Mratio(3) > 1.0) {
            opserr << "LinearElasticSpring::LinearElasticSpring() - "
                << "incorrect p-delta moment ratios:\nrMz1 + rMz2 = "
                << Mratio(2) + Mratio(3) << " > 1.0\n";
            exit(-1);
        }
    }
    
    // initialize optional damping matrix
    if (damp != 0)
        cb = new Matrix(*damp);
    
    // initialize response vectors in basic system
    ub.resize(numDIR);
    ubdot.resize(numDIR);
    qb.resize(numDIR);
    this->revertToStart();
}


LinearElasticSpring::LinearElasticSpring()
    : Element(0, ELE_TAG_LinearElasticSpring),
    numDIM(0), numDOF(0), connectedExternalNodes(2),
    numDIR(0), dir(0), kb(1,1), cb(0), x(0), y(0), Mratio(0),
    addRayleigh(0), L(0.0), onP0(false), trans(3, 3),
    ub(0), ubdot(0), qb(0), ul(0), Tgl(0,0), Tlb(0,0),
    theMatrix(0), theVector(0), theLoad(0)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "LinearElasticSpring::LinearElasticSpring() - "
            << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}


// delete must be invoked on any objects created by the object.
LinearElasticSpring::~LinearElasticSpring()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (cb != 0)
        delete cb;
    if (theLoad != 0)
        delete theLoad;
}


int LinearElasticSpring::getNumExternalNodes() const
{
    return 2;
}


const ID& LinearElasticSpring::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** LinearElasticSpring::getNodePtrs() 
{
    return theNodes;
}


int LinearElasticSpring::getNumDOF() 
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void LinearElasticSpring::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        
        return;
    }
    
    // set default values for error conditions
    numDOF = 2;
    theMatrix = &LinearElasticSpringM2;
    theVector = &LinearElasticSpringV2;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "LinearElasticSpring::setDomain() - Nd1: "
                << Nd1 << " does not exist in the model for ";
        } else  {
            opserr << "LinearElasticSpring::setDomain() - Nd2: " 
                << Nd2 << " does not exist in the model for ";
        }
        opserr << "LinearElasticSpring ele: " << this->getTag() << endln;
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2)  {
        opserr << "LinearElasticSpring::setDomain(): nodes " << Nd1 << " and " << Nd2
            << "have differing dof at ends for element: " << this->getTag() << endln;
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 1 && dofNd1 == 1)  {
        numDOF = 2;
        theMatrix = &LinearElasticSpringM2;
        theVector = &LinearElasticSpringV2;
        elemType  = D1N2;
    }
    else if (numDIM == 2 && dofNd1 == 2)  {
        numDOF = 4;
        theMatrix = &LinearElasticSpringM4;
        theVector = &LinearElasticSpringV4;
        elemType  = D2N4;
    }
    else if (numDIM == 2 && dofNd1 == 3)  {
        numDOF = 6;
        theMatrix = &LinearElasticSpringM6;
        theVector = &LinearElasticSpringV6;
        elemType  = D2N6;
    }
    else if (numDIM == 3 && dofNd1 == 3)  {
        numDOF = 6;
        theMatrix = &LinearElasticSpringM6;
        theVector = &LinearElasticSpringV6;
        elemType  = D3N6;
    }
    else if (numDIM == 3 && dofNd1 == 6)  {
        numDOF = 12;
        theMatrix = &LinearElasticSpringM12;
        theVector = &LinearElasticSpringV12;
        elemType  = D3N12;
    }
    else  {
        opserr << "LinearElasticSpring::setDomain() can not handle "
            << numDIM << "dofs at nodes in " << dofNd1 << " d problem\n";
        return;
    }
    
    // set the local displacement vector size
    ul.resize(numDOF);
    ul.Zero();
    
    // allocate memory for the load vector
    if (theLoad == 0)
        theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF)  {
        delete theLoad;
        theLoad = new Vector(numDOF);
    }
    if (theLoad == 0)  {
        opserr << "LinearElasticSpring::setDomain() - element: " << this->getTag()
            << " out of memory creating vector of size: " << numDOF << endln;
        return;
    }
    
    // setup the transformation matrix for orientation
    this->setUp();
    
    // set transformation matrix from global to local system
    this->setTranGlobalLocal();
    
    // set transformation matrix from local to basic system
    this->setTranLocalBasic();
}


int LinearElasticSpring::commitState()
{
    int errCode = 0;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int LinearElasticSpring::revertToLastCommit()
{
    return 0;
}


int LinearElasticSpring::revertToStart()
{   
    // reset trial history variables
    ub.Zero();
    ubdot.Zero();
    qb.Zero();
    
    return 0;
}


int LinearElasticSpring::update()
{
    int errCode = 0;
    
    // get global trial response
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    
    int numDOF2 = numDOF/2;
    Vector ug(numDOF), ugdot(numDOF), uldot(numDOF);
    for (int i=0; i<numDOF2; i++)  {
        ug(i)         = dsp1(i);  ugdot(i)         = vel1(i);
        ug(i+numDOF2) = dsp2(i);  ugdot(i+numDOF2) = vel2(i);
    }
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
    
    return errCode;
}


const Matrix& LinearElasticSpring::getTangentStiff()
{
    // zero the matrix
    theMatrix->Zero();
    
    // get resisting force and stiffness
    qb.addMatrixVector(0.0, kb, ub, 1.0);
    
    // transform stiffness from basic to local system
    Matrix kl(numDOF,numDOF);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add P-Delta effects to local stiffness
    if (Mratio.Size() == 4)
        this->addPDeltaStiff(kl, qb);
    
    // transform stiffness from local to global system
    theMatrix->addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return *theMatrix;
}


const Matrix& LinearElasticSpring::getInitialStiff()
{
    // zero the matrix
    theMatrix->Zero();
    
    // transform stiffness from basic to local system
    Matrix klInit(numDOF,numDOF);
    klInit.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // transform stiffness from local to global system
    theMatrix->addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);
    
    return *theMatrix;
}


const Matrix& LinearElasticSpring::getDamp()
{
    // zero the matrix
    theMatrix->Zero();
    
    // call base class to setup Rayleigh damping
    double factThis = 0.0;
    if (addRayleigh == 1)  {
        (*theMatrix) = this->Element::getDamp();
        factThis = 1.0;
    }
    
    // add damping from optional user input
    if (cb != 0) {
        // transform damping from basic to local system
        Matrix cl(numDOF,numDOF);
        cl.addMatrixTripleProduct(0.0, Tlb, *cb, 1.0);

        // transform damping from local to global system
        theMatrix->addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    }
    
    return *theMatrix;
}


void LinearElasticSpring::zeroLoad()
{
    theLoad->Zero();
}


int LinearElasticSpring::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"LinearElasticSpring::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int LinearElasticSpring::addInertiaLoadToUnbalance(const Vector &accel)
{
    // this element has no mass
    return 0;
}


const Vector& LinearElasticSpring::getResistingForce()
{
    // zero the residual
    theVector->Zero();
    
    // get resisting force
    qb.addMatrixVector(0.0, kb, ub, 1.0);
    
    // determine resisting force in local system
    Vector ql(numDOF);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
    
    // add P-Delta effects to local force
    if (Mratio.Size() == 4)
        this->addPDeltaForces(ql, qb);
    
    // determine resisting force in global system
    theVector->addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    return *theVector;
}


const Vector& LinearElasticSpring::getResistingForceIncInertia()
{
    this->getResistingForce();
    
    // subtract external load
    theVector->addVector(1.0, *theLoad, -1.0);
    
    // add the damping force from Rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add damping force from optional user input
    if (cb != 0) {
        // get damping force
        Vector qdb(numDIR);
        qdb.addMatrixVector(0.0, *cb, ubdot, 1.0);
        
        // transform damping force from basic to local system
        Vector qdl(numDOF);
        qdl.addMatrixTransposeVector(0.0, Tlb, qdb, 1.0);
        
        // add P-Delta effects to local force
        if (Mratio.Size() == 4)
            this->addPDeltaForces(qdl, qdb);

        // transform damping force from local to global system
        theVector->addMatrixTransposeVector(1.0, Tgl, qdl, 1.0);
    }
    
    return *theVector;
}


int LinearElasticSpring::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(13);
    data(0) = this->getTag();
    data(1) = numDIM;
    data(2) = numDOF;
    data(3) = numDIR;
    data(4) = x.Size();
    data(5) = y.Size();
    data(6) = Mratio.Size();
    data(7) = addRayleigh;
    data(8) = (cb != 0);
    data(9) = alphaM;
    data(10) = betaK;
    data(11) = betaK0;
    data(12) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send the direction array
    sChannel.sendID(0, commitTag, dir);
    
    // send the stiffness matrix
    sChannel.sendMatrix(0, commitTag, kb);
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    if (Mratio.Size() == 4)
        sChannel.sendVector(0, commitTag, Mratio);
    if (cb != 0)
        sChannel.sendMatrix(0, commitTag, *cb);
    
    return 0;
}


int LinearElasticSpring::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete dynamic memory
    if (cb != 0)
        delete cb;
    
    // receive element parameters
    static Vector data(13);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    numDIM = (int)data(1);
    numDOF = (int)data(2);
    numDIR = (int)data(3);
    addRayleigh = (int)data(7);
    alphaM = data(9);
    betaK = data(10);
    betaK0 = data(11);
    betaKc = data(12);
   
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive the direction array
    rChannel.recvID(0, commitTag, dir);
    
    // receive the stiffness matrix
    rChannel.recvMatrix(0, commitTag, kb);
    
    // receive remaining data
    if ((int)data(4) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(5) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    if ((int)data(6) == 4)  {
        Mratio.resize(4);
        rChannel.recvVector(0, commitTag, Mratio);
        // check p-delta moment distribution ratios
        if (Mratio(0) < 0.0 || Mratio(1) < 0.0 ||
            Mratio(2) < 0.0 || Mratio(3) < 0.0) {
            opserr << "LinearElasticSpring::recvSelf() - "
                << "p-delta moment ratios can not be negative\n";
            return -1;
        }
        if (Mratio(0)+Mratio(1) > 1.0)  {
            opserr << "LinearElasticSpring::recvSelf() - "
                << "incorrect p-delta moment ratios:\nrMy1 + rMy2 = "
                << Mratio(0)+Mratio(1) << " > 1.0\n";
            return -1;
        }
        if (Mratio(2)+Mratio(3) > 1.0)  {
            opserr << "LinearElasticSpring::recvSelf() - "
                << "incorrect p-delta moment ratios:\nrMz1 + rMz2 = "
                << Mratio(2)+Mratio(3) << " > 1.0\n";
            return -1;
        }
    }
    // allocate memory for the damping matrix and receive it
    if ((bool)data(8) == true) {
        cb = new Matrix(numDIR, numDIR);
        if (cb == 0) {
            opserr << "LinearElasticSpring::recvSelf() - "
                << "failed to create damping matrix\n";
            return -2;
        }
        rChannel.recvMatrix(0, commitTag, *cb);
    }
    onP0 = false;
    
    // initialize response vectors in basic system
    ub.resize(numDIR);
    ubdot.resize(numDIR);
    qb.resize(numDIR);
    this->revertToStart();
    
    return 0;
}


int LinearElasticSpring::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void LinearElasticSpring::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: LinearElasticSpring" << endln;
        s << "  iNode: " << connectedExternalNodes(0)
            << ", jNode: " << connectedExternalNodes(1) << endln;
        s << "  kb: " << kb << endln;
        s << "  Mratio: " << Mratio << endln;
        s << "  addRayleigh: " << addRayleigh << endln;
        if (cb != 0)
            s << "  cb: " << *cb << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"LinearElasticSpring\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"dof\": [";
        for (int i = 0; i < numDIR - 1; i++) {
            if (dir(i) == 0)
                s << "\"P\", ";
            else if (dir(i) == 1)
                s << "\"Vy\", ";
            else if (dir(i) == 2)
                s << "\"Vz\", ";
            else if (dir(i) == 3)
                s << "\"T\", ";
            else if (dir(i) == 4)
                s << "\"My\", ";
            else if (dir(i) == 5)
                s << "\"Mz\", ";
        }
        if (dir(numDIR - 1) == 0)
            s << "\"P\"], ";
        else if (dir(numDIR - 1) == 1)
            s << "\"Vy\"], ";
        else if (dir(numDIR - 1) == 2)
            s << "\"Vz\"], ";
        else if (dir(numDIR - 1) == 3)
            s << "\"T\"], ";
        else if (dir(numDIR - 1) == 4)
            s << "\"My\"], ";
        else if (dir(numDIR - 1) == 5)
            s << "\"Mz\"], ";
        s << "\"transMatrix\": [[";
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (j < 2)
                    s << trans(i, j) << ", ";
                else if (j == 2 && i < 2)
                    s << trans(i, j) << "], [";
                else if (j == 2 && i == 2)
                    s << trans(i, j) << "]],";
            }
        }
        s << "\"addRayleigh\": " << addRayleigh << "}";
    }
}


Response* LinearElasticSpring::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","LinearElasticSpring");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    char outputData[10];
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
        for (int i=0; i<numDOF; i++)  {
            sprintf(outputData,"P%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 1, *theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
        for (int i=0; i<numDOF; i++)  {
            sprintf(outputData,"p%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, *theVector);
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
        for (int i=0; i<numDIR; i++)  {
            sprintf(outputData,"q%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 3, Vector(numDIR));
    }
    // local displacements
    else if (strcmp(argv[0],"localDisplacement") == 0 ||
        strcmp(argv[0],"localDisplacements") == 0)
    {
        for (int i=0; i<numDOF; i++)  {
            sprintf(outputData,"dl%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 4, Vector(numDOF));
    }
    // basic displacements
    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
        for (int i=0; i<numDIR; i++)  {
            sprintf(outputData,"db%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 5, Vector(numDIR));
    }
    // basic deformations and basic forces
    else if (strcmp(argv[0],"defoANDforce") == 0 || strcmp(argv[0],"deformationANDforce") == 0 ||
        strcmp(argv[0],"deformationsANDforces") == 0)
    {
        int i;
        for (i=0; i<numDIR; i++)  {
            sprintf(outputData,"db%d",i+1);
            output.tag("ResponseType",outputData);
        }
        for (i=0; i<numDIR; i++)  {
            sprintf(outputData,"q%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 6, Vector(numDIR*2));
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int LinearElasticSpring::getResponse(int responseID, Information &eleInfo)
{
    Vector defoAndForce(numDIR*2);
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector->Zero();
        // determine resisting forces in local system
        theVector->addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta effects to local forces
        if (Mratio.Size() == 4)
            this->addPDeltaForces(*theVector, qb);
        
        return eleInfo.setVector(*theVector);
        
    case 3:  // basic forces
        return eleInfo.setVector(qb);
        
    case 4:  // local displacements
        return eleInfo.setVector(ul);
        
    case 5:  // basic displacements
        return eleInfo.setVector(ub);
        
    case 6:  // basic deformations and basic forces
        defoAndForce.Zero();
        defoAndForce.Assemble(ub,0);
        defoAndForce.Assemble(qb,numDIR);
        
        return eleInfo.setVector(defoAndForce);
        
    default:
        return 0;
    }
}


// set up the transformation matrix for orientation
void LinearElasticSpring::setUp()
{
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    // setup x and y orientation vectors
    if (L > DBL_EPSILON)  {
        if (x.Size() == 0)  {
            x.resize(3);
            x.Zero();
            x(0) = xp(0);
            if (xp.Size() > 1)
                x(1) = xp(1);
            if (xp.Size() > 2)
                x(2) = xp(2);
        } else if (onP0)  {
            opserr << "WARNING LinearElasticSpring::setUp() - " 
                << "element: " << this->getTag() << endln
                << "ignoring nodes and using specified "
                << "local x vector to determine orientation\n";
        }
        if (y.Size() == 0)  {
            y.resize(3);
            y.Zero();
            y(0) = -xp(1);
            if (xp.Size() > 1)
                y(1) = xp(0);
            if (xp.Size() > 2)
                opserr << "WARNING LinearElasticSpring::setUp() - " 
                    << "element: " << this->getTag() << endln
                    << "no local y vector specified\n";
        }
    } else  {
        if (x.Size() == 0)  {
            x.resize(3);
            x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
        }
        if (y.Size() == 0)  {
            y.resize(3);
            y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
        }
    }
    
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "LinearElasticSpring::setUp() - "
            << "element: " << this->getTag() << endln
            << "incorrect dimension of orientation vectors\n";
        exit(-1);
    }
    
    // establish orientation of element for the transformation matrix
    // z = x cross yp
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
        opserr << "LinearElasticSpring::setUp() - "
            << "element: " << this->getTag() << endln
            << "invalid orientation vectors\n";
        exit(-1);
    }
    
    // create transformation matrix of direction cosines
    for (int i=0; i<3; i++)  {
        trans(0,i) = x(i)/xn;
        trans(1,i) = y(i)/yn;
        trans(2,i) = z(i)/zn;
    }
}


// set transformation matrix from global to local system
void LinearElasticSpring::setTranGlobalLocal()
{
    // resize transformation matrix and zero it
    Tgl.resize(numDOF,numDOF);
    Tgl.Zero();
    
    // switch on dimensionality of element
    switch (elemType)  {
    case D1N2:
        Tgl(0,0) = Tgl(1,1) = trans(0,0);
        break;
    case D2N4:
        Tgl(0,0) = Tgl(2,2) = trans(0,0);
        Tgl(0,1) = Tgl(2,3) = trans(0,1);
        Tgl(1,0) = Tgl(3,2) = trans(1,0);
        Tgl(1,1) = Tgl(3,3) = trans(1,1);
        break;
    case D2N6:
        Tgl(0,0) = Tgl(3,3) = trans(0,0);
        Tgl(0,1) = Tgl(3,4) = trans(0,1);
        Tgl(1,0) = Tgl(4,3) = trans(1,0);
        Tgl(1,1) = Tgl(4,4) = trans(1,1);
        Tgl(2,2) = Tgl(5,5) = trans(2,2);
        break;
    case D3N6:
        Tgl(0,0) = Tgl(3,3) = trans(0,0);
        Tgl(0,1) = Tgl(3,4) = trans(0,1);
        Tgl(0,2) = Tgl(3,5) = trans(0,2);
        Tgl(1,0) = Tgl(4,3) = trans(1,0);
        Tgl(1,1) = Tgl(4,4) = trans(1,1);
        Tgl(1,2) = Tgl(4,5) = trans(1,2);
        Tgl(2,0) = Tgl(5,3) = trans(2,0);
        Tgl(2,1) = Tgl(5,4) = trans(2,1);
        Tgl(2,2) = Tgl(5,5) = trans(2,2);
        break;
    case D3N12:
        Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = trans(0,0);
        Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = trans(0,1);
        Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = trans(0,2);
        Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = trans(1,0);
        Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = trans(1,1);
        Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = trans(1,2);
        Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = trans(2,0);
        Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = trans(2,1);
        Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = trans(2,2);
        break;
    }
}


// set transformation matrix from local to basic system
void LinearElasticSpring::setTranLocalBasic()
{
    // resize transformation matrix and zero it
    Tlb.resize(numDIR,numDOF);
    Tlb.Zero();
    
    for (int i=0; i<numDIR; i++)  {
        
        int dirID = dir(i);     // direction 0 to 5;
        Tlb(i,dirID) = -1.0;
        Tlb(i,dirID+numDOF/2) = 1.0;
    }
}


void LinearElasticSpring::addPDeltaForces(Vector &pLocal, const Vector& qBasic)
{
    int dirID;
    double N = 0.0;
    double deltal1 = 0.0;
    double deltal2 = 0.0;
    
    for (int i=0; i<numDIR; i++)  {
        dirID = dir(i);  // direction 0 to 5;
        
        // get axial force and local disp differences
        if (dirID == 0)
            N = qBasic(i);
        else if (dirID == 1 && numDIM > 1)
            deltal1 = ul(1+numDOF/2) - ul(1);
        else if (dirID == 2 && numDIM > 2)
            deltal2 = ul(2+numDOF/2) - ul(2);
    }
    
    if (N != 0.0 && (deltal1 != 0.0 || deltal2 != 0.0))  {
        for (int i=0; i<numDIR; i++)  {
            dirID = dir(i);  // direction 0 to 5;
            
            // switch on dimensionality of element
            switch (elemType)  {
            case D2N4:
                if (dirID == 1)  {
                    double VpDelta = N*deltal1/L;
                    VpDelta *= 1.0-Mratio(2)-Mratio(3);
                    pLocal(1) -= VpDelta;
                    pLocal(3) += VpDelta;
                }
                break;
            case D2N6: 
                if (dirID == 1)  {
                    double VpDelta = N*deltal1/L;
                    VpDelta *= 1.0-Mratio(2)-Mratio(3);
                    pLocal(1) -= VpDelta;
                    pLocal(4) += VpDelta;
                }
                else if (dirID == 2)  {
                    double MpDelta = N*deltal1;
                    pLocal(2) += Mratio(2)*MpDelta;
                    pLocal(5) += Mratio(3)*MpDelta;
                }
                break;
            case D3N6:
                if (dirID == 1)  {
                    double VpDelta = N*deltal1/L;
                    VpDelta *= 1.0-Mratio(2)-Mratio(3);
                    pLocal(1) -= VpDelta;
                    pLocal(4) += VpDelta;
                }
                else if (dirID == 2)  {
                    double VpDelta = N*deltal2/L;
                    VpDelta *= 1.0-Mratio(0)-Mratio(1);
                    pLocal(2) -= VpDelta;
                    pLocal(5) += VpDelta;
                }
                break;
            case D3N12:
                if (dirID == 1)  {
                    double VpDelta = N*deltal1/L;
                    VpDelta *= 1.0-Mratio(2)-Mratio(3);
                    pLocal(1) -= VpDelta;
                    pLocal(7) += VpDelta;
                }
                else if (dirID == 2)  {
                    double VpDelta = N*deltal2/L;
                    VpDelta *= 1.0-Mratio(0)-Mratio(1);
                    pLocal(2) -= VpDelta;
                    pLocal(8) += VpDelta;
                }
                else if (dirID == 4)  {
                    double MpDelta = N*deltal2;
                    pLocal(4) -= Mratio(0)*MpDelta;
                    pLocal(10) -= Mratio(1)*MpDelta;
                }
                else if (dirID == 5)  {
                    double MpDelta = N*deltal1;
                    pLocal(5) += Mratio(2)*MpDelta;
                    pLocal(11) += Mratio(3)*MpDelta;
                }
                break;
            default :
                // do nothing
                break;
            }
        }
    }
}


void LinearElasticSpring::addPDeltaStiff(Matrix &kLocal, const Vector& qBasic)
{
    int dirID;
    double N = 0.0;
    
    // get axial force
    for (int i=0; i<numDIR; i++)  {
        if (dir(i) == 0)
            N = qBasic(i);
    }
    
    if (N != 0.0)  {
        for (int i=0; i<numDIR; i++)  {
            dirID = dir(i);  // direction 0 to 5;
            
            // switch on dimensionality of element
            switch (elemType)  {
            case D2N4:
                if (dirID == 1)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(2)-Mratio(3);
                    kLocal(1,1) += NoverL;
                    kLocal(1,3) -= NoverL;
                    kLocal(3,1) -= NoverL;
                    kLocal(3,3) += NoverL;
                }
                break;
            case D2N6: 
                if (dirID == 1)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(2)-Mratio(3);
                    kLocal(1,1) += NoverL;
                    kLocal(1,4) -= NoverL;
                    kLocal(4,1) -= NoverL;
                    kLocal(4,4) += NoverL;
                }
                else if (dirID == 2)  {
                    kLocal(2,1) -= Mratio(2)*N;
                    kLocal(2,4) += Mratio(2)*N;
                    kLocal(5,1) -= Mratio(3)*N;
                    kLocal(5,4) += Mratio(3)*N;
                }
                break;
            case D3N6:
                if (dirID == 1)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(2)-Mratio(3);
                    kLocal(1,1) += NoverL;
                    kLocal(1,4) -= NoverL;
                    kLocal(4,1) -= NoverL;
                    kLocal(4,4) += NoverL;
                }
                else if (dirID == 2)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(0)-Mratio(1);
                    kLocal(2,2) += NoverL;
                    kLocal(2,5) -= NoverL;
                    kLocal(5,2) -= NoverL;
                    kLocal(5,5) += NoverL;
                }
                break;
            case D3N12:
                if (dirID == 1)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(2)-Mratio(3);
                    kLocal(1,1) += NoverL;
                    kLocal(1,7) -= NoverL;
                    kLocal(7,1) -= NoverL;
                    kLocal(7,7) += NoverL;
                }
                else if (dirID == 2)  {
                    double NoverL = N/L;
                    NoverL *= 1.0-Mratio(0)-Mratio(1);
                    kLocal(2,2) += NoverL;
                    kLocal(2,8) -= NoverL;
                    kLocal(8,2) -= NoverL;
                    kLocal(8,8) += NoverL;
                }
                else if (dirID == 4)  {
                    kLocal(4,2) += Mratio(0)*N;
                    kLocal(4,8) -= Mratio(0)*N;
                    kLocal(10,2) += Mratio(1)*N;
                    kLocal(10,8) -= Mratio(1)*N;
                }
                else if (dirID == 5)  {
                    kLocal(5,1) -= Mratio(2)*N;
                    kLocal(5,7) += Mratio(2)*N;
                    kLocal(11,1) -= Mratio(3)*N;
                    kLocal(11,7) += Mratio(3)*N;
                }
                break;
            default :
                // do nothing
                break;
            }
        }
    }
}
