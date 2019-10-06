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
// Created: 09/07
// Revision: A
//
// Description: This file contains the implementation of the Adapter class.

#include "Adapter.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <TCP_Socket.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

void* OPS_Adapter()
{
    int ndf = OPS_GetNDF();
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element adapter eleTag -node Ndi Ndj ... -dof dofNdi -dof dofNdj ... -stif Kij ipPort <-doRayleigh> <-mass Mij>\n";
        return 0;
    }
    
    // tags
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING: invalid tag\n";
        return 0;
    }
    
    // nodes
    const char* type = OPS_GetString();
    if (strcmp(type, "-node") != 0) {
        opserr << "WARNING expecting -node Ndi Ndj ...\n";
        return 0;
    }
    ID nodes(32);
    int numNodes = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        int node;
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &node) < 0) {
            break;
        }
        nodes(numNodes++) = node;
    }
    nodes.resize(numNodes);
    
    // dofs
    int numDOF = 0;
    ID *dofs = new ID[numNodes];
    for (int i = 0; i < numNodes; i++) {
        type = OPS_GetString();
        if (strcmp(type, "-dof") != 0) {
            opserr << "WARNING expecting -dof dofNdi\n";
            return 0;
        }
        ID dofsi(ndf);
        int numDOFi = 0;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            int dof;
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &dof) < 0) {
                break;
            }
            if (dof < 1 || ndf < dof) {
                opserr << "WARNING invalid dof ID\n";
                return 0;
            }
            dofsi(numDOFi++) = dof - 1;
            numDOF++;
        }
        dofsi.resize(numDOFi);
        dofs[i] = dofsi;
    }
    
    // stiffness matrix terms
    type = OPS_GetString();
    if (strcmp(type, "-stif") != 0 && strcmp(type, "-stiff") != 0) {
        opserr << "WARNING expecting -stif kij\n";
        return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < numDOF*numDOF) {
        opserr << "WARNING wrong number of kij specified\n";
        return 0;
    }
    Matrix kb(numDOF, numDOF);
    numdata = 1;
    for (int i = 0; i < numDOF; i++) {
        for (int j = 0; j < numDOF; j++) {
            if (OPS_GetDoubleInput(&numdata, &kb(i, j)) < 0) {
                opserr << "WARNING invalid stiffness value\n";
                return 0;
            }
        }
    }
    // ipPort
    int ipPort;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &ipPort) < 0) {
        opserr << "WARNING: invalid ipPort\n";
        return 0;
    }
    
    // options
    int doRayleigh = 0;
    Matrix *mb = 0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        return new Adapter(tag, nodes, dofs, kb, ipPort);
    }
    
    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        if (strcmp(type, "-doRayleigh") == 0) {
            doRayleigh = 1;
        }
        else if (strcmp(type, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < numDOF*numDOF) {
                opserr << "WARNING wrong number of mij specified\n";
                return 0;
            }
            double mij;
            numdata = 1;
            mb = new Matrix(numDOF, numDOF);
            for (int i = 0; i < numDOF; i++) {
                for (int j = 0; j < numDOF; j++) {
                    if (OPS_GetDoubleInput(&numdata, &mij) < 0) {
                        opserr << "WARNING invalid damping value\n";
                        delete mb;
                        return 0;
                    }
                    (*mb)(i, j) = mij;
                }
            }
        }
    }
    
    // create object
    Element *theEle = new Adapter(tag, nodes, dofs, kb, ipPort,
        doRayleigh, mb);
    
    // clean up memory
    delete mb;
    
    return theEle;
}


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
Adapter::Adapter(int tag, ID nodes, ID *dof,
    const Matrix &_kb, int ipport, int addRay, const Matrix *_mb)
    : Element(tag, ELE_TAG_Adapter),
    connectedExternalNodes(nodes), basicDOF(1), numExternalNodes(0),
    numDOF(0), numBasicDOF(0), kb(_kb), ipPort(ipport), addRayleigh(addRay),
    mb(0), tPast(0.0), theMatrix(1,1), theVector(1), theLoad(1), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlVel(0), ctrlAccel(0), ctrlForce(0), ctrlTime(0),
    daqDisp(0), daqVel(0), daqAccel(0), daqForce(0), daqTime(0)
{
    // initialize nodes
    numExternalNodes = connectedExternalNodes.Size();
    theNodes = new Node* [numExternalNodes];
    if (!theNodes)  {
        opserr << "Adapter::Adapter() "
            << "- failed to create node array\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    int i;
    for (i=0; i<numExternalNodes; i++)
        theNodes[i] = 0;
    
    // initialize dof
    theDOF = new ID [numExternalNodes];
    if (!theDOF)  {
        opserr << "Adapter::Adapter() "
            << "- failed to create dof array\n";
        exit(-1);
    }
    numBasicDOF = 0;
    for (i=0; i<numExternalNodes; i++)  {
        theDOF[i] = dof[i];
        numBasicDOF += theDOF[i].Size();
    }
    
    // initialize mass matrix
    if (_mb != 0)
        mb = new Matrix(*_mb);
    
    // set the vector sizes and zero them
    basicDOF.resize(numBasicDOF);
    basicDOF.Zero();
    db.resize(numBasicDOF);
    db.Zero();
    q.resize(numBasicDOF);
    q.Zero();
}


// invoked by a FEM_ObjectBroker - blank object that recvSelf
// needs to be invoked upon
Adapter::Adapter()
    : Element(0, ELE_TAG_Adapter),
    connectedExternalNodes(1), basicDOF(1), numExternalNodes(0),
    numDOF(0), numBasicDOF(0), kb(1,1), ipPort(0), addRayleigh(0), mb(0),
    tPast(0.0), theMatrix(1,1), theVector(1), theLoad(1), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlVel(0), ctrlAccel(0), ctrlForce(0), ctrlTime(0),
    daqDisp(0), daqVel(0), daqAccel(0), daqForce(0), daqTime(0)
{
    // initialize variables
    theNodes = 0;
    theDOF = 0;
}


// delete must be invoked on any objects created by the object.
Adapter::~Adapter()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theNodes != 0)
        delete [] theNodes;
    if (theDOF != 0)
        delete [] theDOF;
    if (mb != 0)
        delete mb;
    
    if (daqDisp != 0)
        delete daqDisp;
    if (daqVel != 0)
        delete daqVel;
    if (daqAccel != 0)
        delete daqAccel;
    if (daqForce != 0)
        delete daqForce;
    if (daqTime != 0)
        delete daqTime;
    
    if (ctrlDisp != 0)
        delete ctrlDisp;
    if (ctrlVel != 0)
        delete ctrlVel;
    if (ctrlAccel != 0)
        delete ctrlAccel;
    if (ctrlForce != 0)
        delete ctrlForce;
    if (ctrlTime != 0)
        delete ctrlTime;
    
    if (sendData != 0)
        delete sendData;
    if (sData != 0)
        delete [] sData;
    if (recvData != 0)
        delete recvData;
    if (rData != 0)
        delete [] rData;
    if (theChannel != 0)
        delete theChannel;
}


int Adapter::getNumExternalNodes() const
{
    return numExternalNodes;
}


const ID& Adapter::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** Adapter::getNodePtrs()
{
    return theNodes;
}


int Adapter::getNumDOF()
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void Adapter::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    int i;
    if (!theDomain)  {
        for (i=0; i<numExternalNodes; i++)
            theNodes[i] = 0;
        return;
    }
    
    // first set the node pointers
    for (i=0; i<numExternalNodes; i++)
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    
    // if can't find all - send a warning message
    for (i=0; i<numExternalNodes; i++)  {
        if (!theNodes[i])  {
            opserr << "Adapter::setDomain() - Nd" << i << ": " 
                << connectedExternalNodes(i) << " does not exist in the "
                << "model for Adapter ele: " << this->getTag() << endln;
            return;
        }
    }
    
    // now determine the number of dof
    numDOF = 0;
    for (i=0; i<numExternalNodes; i++)  {
        numDOF += theNodes[i]->getNumberDOF();
    }
    
    // set the basicDOF ID
    int j, k = 0, ndf = 0;
    for (i=0; i<numExternalNodes; i++)  {
        for (j=0; j<theDOF[i].Size(); j++)  {
            basicDOF(k) = ndf + theDOF[i](j);
            k++;
        }
        ndf += theNodes[i]->getNumberDOF();
    }
    
    // set the matrix and vector sizes and zero them
    theMatrix.resize(numDOF,numDOF);
    theMatrix.Zero();
    theVector.resize(numDOF);
    theVector.Zero();
    theLoad.resize(numDOF);
    theLoad.Zero();
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
}


int Adapter::commitState()
{
    int errCode = 0;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int Adapter::revertToLastCommit()
{
    opserr << "Adapter::revertToLastCommit() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to last commit. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int Adapter::revertToStart()
{
    opserr << "Adapter::revertToStart() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to start. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int Adapter::update()
{
    if (theChannel == 0)  {
        if (this->setupConnection() != 0)  {
            opserr << "Adapter::update() - "
                << "failed to setup connection\n";
            return -1;
        }
    }
    
    // assemble dsp in basic system
    int ndim = 0;
    db.Zero();
    for (int i=0; i<numExternalNodes; i++)  {
        Vector disp = theNodes[i]->getTrialDisp();
        db.Assemble(disp(theDOF[i]), ndim);
        ndim += theDOF[i].Size();
    }
    
    return 0;
}


const Matrix& Adapter::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // assemble stiffness matrix
    theMatrix.Assemble(kb, basicDOF, basicDOF);
    
    return theMatrix;
}


const Matrix& Adapter::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // assemble stiffness matrix
    theMatrix.Assemble(kb, basicDOF, basicDOF);
    
    return theMatrix;
}


const Matrix& Adapter::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    // call base class to setup Rayleigh damping
    if (addRayleigh == 1)
        theMatrix = this->Element::getDamp();
    
    return theMatrix;
}


const Matrix& Adapter::getMass()
{
    // zero the matrix
    theMatrix.Zero();
    
    // assemble mass matrix
    if (mb != 0)
        theMatrix.Assemble(*mb, basicDOF, basicDOF);
    
    return theMatrix;
}


void Adapter::zeroLoad()
{
    theLoad.Zero();
}


int Adapter::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"Adapter::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int Adapter::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mb == 0)
        return 0;
    
    int ndim = 0, i;
    Vector Raccel(numDOF);
    
    // get mass matrix
    Matrix M = this->getMass();
    // assemble Raccel vector
    for (i=0; i<numExternalNodes; i++ )  {
        Raccel.Assemble(theNodes[i]->getRV(accel), ndim);
        ndim += theNodes[i]->getNumberDOF();
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    theLoad.addMatrixVector(1.0, M, Raccel, -1.0);
    
    return 0;
}


const Vector& Adapter::getResistingForce()
{
    // get current time
    Domain *theDomain = this->getDomain();
    double t = theDomain->getCurrentTime();
    
    // update response if time has advanced
    if (t > tPast)  {
        // receive data
        theChannel->recvVector(0, 0, *recvData, 0);
        
        // check if force request was received
        if (rData[0] == RemoteTest_getForce)  {
            // send daq displacements and forces
            theChannel->sendVector(0, 0, *sendData, 0);
            
            // receive new trial response
            theChannel->recvVector(0, 0, *recvData, 0);
        }
        
        if (rData[0] != RemoteTest_setTrialResponse)  {
            if (rData[0] == RemoteTest_DIE)  {
                opserr << "\nThe Simulation has successfully completed.\n";
            } else  {
                opserr << "Adapter::getResistingForce() - "
                    << "wrong action received: expecting 3 but got "
                    << rData[0] << endln;
            }
            exit(-1);
        }
        
        // set velocities at nodes
        if (ctrlVel != 0)  {
            int i, j, ndim = 0;
            for (i=0; i<numExternalNodes; i++ )  {
                Vector vel = theNodes[i]->getTrialVel();
                for (j=0; j<theDOF[i].Size(); j++)  {
                    vel(theDOF[i](j)) = (*ctrlVel)(ndim);
                    ndim++;
                }
                theNodes[i]->setTrialVel(vel);
            }
        }
        
        // set accelerations at nodes
        if (ctrlAccel != 0)  {
            int i, j, ndim = 0;
            for (i=0; i<numExternalNodes; i++ )  {
                Vector accel = theNodes[i]->getTrialAccel();
                for (j=0; j<theDOF[i].Size(); j++)  {
                    accel(theDOF[i](j)) = (*ctrlAccel)(ndim);
                    ndim++;
                }
                theNodes[i]->setTrialAccel(accel);
            }
        }
        
        // save current time
        tPast = t;
    }
    
    // get resisting force in basic system q = k*db + q0 = k*(db - db0)
    q.addMatrixVector(0.0, kb, (db - *ctrlDisp), 1.0);
    //q = kb*(db - *ctrlDisp);
    
    // assign daq values for feedback
    *daqDisp  = db;
    *daqForce = -1.0*q;
    
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in global system
    theVector.Assemble(q, basicDOF);
    
    return theVector;
}


const Vector& Adapter::getResistingForceIncInertia()
{
    theVector = this->getResistingForce();
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    // add the damping forces from rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add inertia forces from element mass
    if (mb != 0)  {
        int ndim = 0, i;
        Vector accel(numDOF);
        
        // get mass matrix
        Matrix M = this->getMass();
        // assemble accel vector
        for (i=0; i<numExternalNodes; i++ )  {
            accel.Assemble(theNodes[i]->getTrialAccel(), ndim);
            ndim += theNodes[i]->getNumberDOF();
        }
        
        theVector.addMatrixVector(1.0, M, accel, 1.0);
    }
    
    return theVector;
}


int Adapter::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(9);
    data(0) = this->getTag();
    data(1) = numExternalNodes;
    data(2) = ipPort;
    data(3) = addRayleigh;
    data(4) = (mb==0) ? 0 : 1;
    data(5) = alphaM;
    data(6) = betaK;
    data(7) = betaK0;
    data(8) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the end nodes and dofs
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    for (int i=0; i<numExternalNodes; i++)
        sChannel.sendID(0, commitTag, theDOF[i]);
    
    // send the stiffness and mass matrices
    sChannel.sendMatrix(0, commitTag, kb);
    if ((int)data(4))
        sChannel.sendMatrix(0, commitTag, *mb);
    
    return 0;
}


int Adapter::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete dynamic memory
    if (theNodes != 0)
        delete [] theNodes;
    if (theDOF != 0)
        delete [] theDOF;
    if (mb != 0)
        delete mb;
    
    // receive element parameters
    static Vector data(9);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    numExternalNodes = (int)data(1);
    ipPort = (int)data(2);
    addRayleigh = (int)data(3);
    alphaM = data(5);
    betaK = data(6);
    betaK0 = data(7);
    betaKc = data(8);
    
    // initialize nodes and receive them
    connectedExternalNodes.resize(numExternalNodes);
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    theNodes = new Node* [numExternalNodes];
    if (theNodes == 0)  {
        opserr << "GenericClient::recvSelf() "
            << "- failed to create node array\n";
        return -1;
    }
    
    // set node pointers to NULL
    int i;
    for (i=0; i<numExternalNodes; i++)
        theNodes[i] = 0;
    
    // initialize dof
    theDOF = new ID [numExternalNodes];
    if (theDOF == 0)  {
        opserr << "GenericClient::recvSelf() "
            << "- failed to create dof array\n";
        return -2;
    }
    
    // initialize number of basic dof
    numBasicDOF = 0;
    for (i=0; i<numExternalNodes; i++)  {
        rChannel.recvID(0, commitTag, theDOF[i]);
        numBasicDOF += theDOF[i].Size();
    }
    
    // receive the stiffness and mass matrices
    kb.resize(numBasicDOF,numBasicDOF);
    rChannel.recvMatrix(0, commitTag, kb);
    if ((int)data(4))  {
        mb = new Matrix(numBasicDOF,numBasicDOF);
        rChannel.recvMatrix(0, commitTag, *mb);
    }
    
    // set the vector sizes and zero them
    basicDOF.resize(numBasicDOF);
    basicDOF.Zero();
    db.resize(numBasicDOF);
    db.Zero();
    q.resize(numBasicDOF);
    q.Zero();
    
    return 0;
}


int Adapter::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    int rValue = 0, i, j;
    
    if (numExternalNodes > 1)  {
        if (displayMode >= 0)  {
            for (i=0; i<numExternalNodes-1; i++)  {
                const Vector &end1Crd = theNodes[i]->getCrds();
                const Vector &end2Crd = theNodes[i+1]->getCrds();
                
                const Vector &end1Disp = theNodes[i]->getDisp();
                const Vector &end2Disp = theNodes[i+1]->getDisp();
                
                int end1NumCrds = end1Crd.Size();
                int end2NumCrds = end2Crd.Size();
                
                static Vector v1(3), v2(3);
                
                for (j=0; j<end1NumCrds; j++)
                    v1(j) = end1Crd(j) + end1Disp(j)*fact;
                for (j=0; j<end2NumCrds; j++)
                    v2(j) = end2Crd(j) + end2Disp(j)*fact;
                
                rValue += theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(), 0);
            }
        } else  {
            int mode = displayMode * -1;
            for (i=0; i<numExternalNodes-1; i++)  {
                const Vector &end1Crd = theNodes[i]->getCrds();
                const Vector &end2Crd = theNodes[i+1]->getCrds();
                
                const Matrix &eigen1 = theNodes[i]->getEigenvectors();
                const Matrix &eigen2 = theNodes[i+1]->getEigenvectors();
                
                int end1NumCrds = end1Crd.Size();
                int end2NumCrds = end2Crd.Size();
                
                static Vector v1(3), v2(3);
                
                if (eigen1.noCols() >= mode)  {
                    for (j=0; j<end1NumCrds; j++)
                        v1(j) = end1Crd(j) + eigen1(j,mode-1)*fact;
                    for (j=0; j<end2NumCrds; j++)
                        v2(j) = end2Crd(j) + eigen2(j,mode-1)*fact;
                } else  {
                    for (j=0; j<end1NumCrds; j++)
                        v1(j) = end1Crd(j);
                    for (j=0; j<end2NumCrds; j++)
                        v2(j) = end2Crd(j);
                }
                
                rValue += theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(), 0);
            }
        }
    }
    
    return rValue;
}


void Adapter::Print(OPS_Stream &s, int flag)
{
    int i;
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: Adapter";
        for (i=0; i<numExternalNodes; i++ )
            s << ", Node" << i+1 << ": " << connectedExternalNodes(i);
        s << endln;
        s << "  kb: " << kb << endln;
        s << "  ipPort: " << ipPort << endln;
        s << "  addRayleigh: " << addRayleigh << endln;
        if (mb != 0)
            s << "  mb: " << *mb << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"Adapter\", ";
        s << "\"nodes\": [";
        for (i = 0; i < numExternalNodes - 1; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(numExternalNodes) << "], ";
        s << "\"kb\": [" << kb << "], ";
        s << "\"ipPort\": " << ipPort << ", ";
        s << "\"addRayleigh\": " << addRayleigh;
        if (mb != 0)
            s << ", \"mb\": [" << *mb << "]}";
        else
            s << "}";
    }
}


Response* Adapter::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    int i;
    char outputData[10];
    
    output.tag("ElementOutput");
    output.attr("eleType","Adapter");
    output.attr("eleTag",this->getTag());
    for (i=0; i<numExternalNodes; i++ )  {
        sprintf(outputData,"node%d",i+1);
        output.attr(outputData,connectedExternalNodes[i]);
    }
    
    // global forces
    if (strcmp(argv[0],"force") == 0 ||
        strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0)
    {
         for (i=0; i<numDOF; i++)  {
            sprintf(outputData,"P%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, theVector);
    }
    
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
    {
        for (i=0; i<numDOF; i++)  {
            sprintf(outputData,"p%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 3, theVector);
    }
    
    // forces in basic system
    else if (strcmp(argv[0],"basicForce") == 0 ||
        strcmp(argv[0],"basicForces") == 0 ||
        strcmp(argv[0],"daqForce") == 0 ||
        strcmp(argv[0],"daqForces") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"q%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 4, Vector(numBasicDOF));
    }
    
    // ctrl basic displacements
    else if (strcmp(argv[0],"defo") == 0 ||
        strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 ||
        strcmp(argv[0],"basicDefo") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 ||
        strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"ctrlDisp") == 0 ||
        strcmp(argv[0],"ctrlDisplacement") == 0 ||
        strcmp(argv[0],"ctrlDisplacements") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"db%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 5, Vector(numBasicDOF));
    }
    
    // ctrl basic velocities
    else if (strcmp(argv[0],"basicVel") == 0 ||
        strcmp(argv[0],"basicVelocity") == 0 ||
        strcmp(argv[0],"basicVelocities") == 0 ||
        strcmp(argv[0],"ctrlVel") == 0 ||
        strcmp(argv[0],"ctrlVelocity") == 0 ||
        strcmp(argv[0],"ctrlVelocities") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"vb%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 6, Vector(numBasicDOF));
    }
    
    // ctrl basic accelerations
    else if (strcmp(argv[0],"basicAccel") == 0 ||
        strcmp(argv[0],"basicAcceleration") == 0 ||
        strcmp(argv[0],"basicAccelerations") == 0 ||
        strcmp(argv[0],"ctrlAccel") == 0 ||
        strcmp(argv[0],"ctrlAcceleration") == 0 ||
        strcmp(argv[0],"ctrlAccelerations") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"ab%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 7, Vector(numBasicDOF));
    }
    
    // daq basic displacements
    else if (strcmp(argv[0],"daqDisp") == 0 ||
        strcmp(argv[0],"daqDisplacement") == 0 ||
        strcmp(argv[0],"daqDisplacements") == 0)
    {
        for (int i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"dbm%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 8, Vector(numBasicDOF));
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int Adapter::getResponse(int responseID, Information &eleInformation)
{
    switch (responseID)  {
    case -1:
        return -1;
        
    case 1:  // stiffness
        if (eleInformation.theMatrix != 0)  {
            *(eleInformation.theMatrix) = this->getTangentStiff();
        }
        return 0;
        
    case 2:  // global forces
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = this->getResistingForce();
        }
        return 0;
        
    case 3:  // local forces
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = this->getResistingForce();
        }
        return 0;
        
    case 4:  // basic forces
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = q;
        }
        return 0;
        
    case 5:  // ctrl basic displacements
        if (eleInformation.theVector != 0  &&  ctrlDisp != 0)  {
            *(eleInformation.theVector) = *ctrlDisp;
        }
        return 0;
        
    case 6:  // ctrl basic velocities
        if (eleInformation.theVector != 0  &&  ctrlVel != 0)  {
            *(eleInformation.theVector) = *ctrlVel;
        }
        return 0;
        
    case 7:  // ctrl basic accelerations
        if (eleInformation.theVector != 0  &&  ctrlAccel != 0)  {
            *(eleInformation.theVector) = *ctrlAccel;
        }
        return 0;
        
    case 8:  // daq basic displacements
        if (eleInformation.theVector != 0  &&  daqDisp != 0)  {
            *(eleInformation.theVector) = *daqDisp;
        }
        return 0;
        
    default:
        return -1;
    }
}


int Adapter::setupConnection()
{
    // setup the connection
    theChannel = new TCP_Socket(ipPort);
    if (theChannel != 0) {
        opserr << "\nChannel successfully created: "
            << "Waiting for ECSimAdapter experimental control...\n";
    } else {
        opserr << "Adapter::setupConnection() - "
            << "could not create channel\n";
        return -1;
    }
    if (theChannel->setUpConnection() != 0)  {
        opserr << "Adapter::setupConnection() - "
            << "failed to setup connection\n";
        return -2;
    }
    
    // get the data sizes
    // sizes = {ctrlDisp, ctrlVel, ctrlAccel, ctrlForce, ctrlTime,
    //          daqDisp,  daqVel,  daqAccel,  daqForce,  daqTime,  dataSize}
    ID sizes(11);
    theChannel->recvID(0, 0, sizes, 0);
    for (int i=0; i<10; i++)  {
        if (sizes(i) != 0 && sizes(i) != numBasicDOF)  {
            opserr << "Adapter::Adapter() - wrong data sizes != "
                << numBasicDOF << " received\n";
            return -3;
        }
    }
    
    // allocate memory for the receive vectors
    int id = 1;
    rData = new double [sizes(10)];
    recvData = new Vector(rData, sizes(10));
    if (sizes(0) != 0)  {
        ctrlDisp = new Vector(&rData[id], sizes(0));
        id += sizes(0);
    }
    if (sizes(1) != 0)  {
        ctrlVel = new Vector(&rData[id], sizes(1));
        id += sizes(1);
    }
    if (sizes(2) != 0)  {
        ctrlAccel = new Vector(&rData[id], sizes(2));
        id += sizes(2);
    }
    if (sizes(3) != 0)  {
        ctrlForce = new Vector(&rData[id], sizes(3));
        id += sizes(3);
    }
    if (sizes(4) != 0)  {
        ctrlTime = new Vector(&rData[id], sizes(4));
        id += sizes(4);
    }
    recvData->Zero();
    
    // allocate memory for the send vectors
    id = 0;
    sData = new double [sizes(10)];
    sendData = new Vector(sData, sizes(10));
    if (sizes(5) != 0)  {
        daqDisp = new Vector(&sData[id], sizes(5));
        id += sizes(5);
    }
    if (sizes(6) != 0)  {
        daqVel = new Vector(&sData[id], sizes(6));
        id += sizes(6);
    }
    if (sizes(7) != 0)  {
        daqAccel = new Vector(&sData[id], sizes(7));
        id += sizes(7);
    }
    if (sizes(8) != 0)  {
        daqForce = new Vector(&sData[id], sizes(8));
        id += sizes(8);
    }
    if (sizes(9) != 0)  {
        daqTime = new Vector(&sData[id], sizes(9));
        id += sizes(9);
    }
    sendData->Zero();
    
    opserr << "\nAdapter element " << this->getTag()
        << " now running...\n";
    
    return 0;
}
