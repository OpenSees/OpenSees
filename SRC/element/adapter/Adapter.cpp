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

// $Revision: 1.4 $
// $Date: 2009-06-02 21:10:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/adapter/Adapter.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
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


// initialize the class wide variables
Matrix Adapter::theMatrix(1,1);
Vector Adapter::theVector(1);
Vector Adapter::theLoad(1);


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
Adapter::Adapter(int tag, ID nodes, ID *dof,
    const Matrix &_kb, int ipport, const Matrix *_mb)
    : Element(tag, ELE_TAG_Adapter),
    connectedExternalNodes(nodes), basicDOF(1),
    numExternalNodes(0), numDOF(0), numBasicDOF(0),
    kb(_kb), ipPort(ipport), mb(0), tPast(0.0), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
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
    connectedExternalNodes(1), basicDOF(1),
    numExternalNodes(0), numDOF(0), numBasicDOF(0),
    kb(1,1), ipPort(0), mb(0), tPast(0.0), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
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
    if (daqForce != 0)
        delete daqForce;
    if (ctrlDisp != 0)
        delete ctrlDisp;
    if (ctrlForce != 0)
        delete ctrlForce;
    
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
    return 0;
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
    theMatrix.Assemble(kb,basicDOF,basicDOF);
    
    return theMatrix;
}


const Matrix& Adapter::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // assemble stiffness matrix
    theMatrix.Assemble(kb,basicDOF,basicDOF);
    
    return theMatrix;
}


const Matrix& Adapter::getMass()
{
    // zero the matrix
    theMatrix.Zero();

    // assemble mass matrix
    if (mb != 0)
        theMatrix.Assemble(*mb,basicDOF,basicDOF);
    
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
    static Vector Raccel(numDOF);
    Raccel.Zero();

    // get mass matrix
    Matrix M = this->getMass();
    // assemble Raccel vector
    for (i=0; i<numExternalNodes; i++ )  {
        Raccel.Assemble(theNodes[i]->getRV(accel), ndim);
        ndim += theNodes[i]->getNumberDOF();
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    theLoad -= M * Raccel;
    
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
        
        // save current time
        tPast = t;
    }
    
    // get resisting force in basic system q = k*db + q0 = k*(db - db0)
    q = kb*(db - *ctrlDisp);
    
    // assign daq values for feedback
    *daqDisp  = db;
    *daqForce = -1.0*q;
    
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in global system
    theVector.Assemble(q, basicDOF);
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    return theVector;
}


const Vector& Adapter::getResistingForceIncInertia()
{	
    theVector = this->getResistingForce();
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        theVector += this->getRayleighDampingForces();
    
    // now include the mass portion
    int ndim = 0, i;
    static Vector accel(numDOF);
    accel.Zero();

    // get mass matrix
    Matrix M = this->getMass();
    // assemble accel vector
    for (i=0; i<numExternalNodes; i++ )  {
        accel.Assemble(theNodes[i]->getTrialAccel(), ndim);
        ndim += theNodes[i]->getNumberDOF();
    }
    
    theVector += M * accel;
    
    return theVector;
}


int Adapter::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static ID idData(4);
    idData(0) = this->getTag();
    idData(1) = numExternalNodes;
    idData(2) = ipPort;
    idData(3) = (mb==0) ? 0 : 1;
    sChannel.sendID(0, commitTag, idData);
    
    // send the end nodes and dofs
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    for (int i=0; i<numExternalNodes; i++)
        sChannel.sendID(0, commitTag, theDOF[i]);
    
    // send the stiffness and mass matrices
    sChannel.sendMatrix(0, commitTag, kb);
    if (idData(3))
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
    static ID idData(4);
    rChannel.recvID(0, commitTag, idData);
    this->setTag(idData(0));
    numExternalNodes = idData(1);
    ipPort = idData(2);
    
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
    if (idData(3))  {
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
    int displayMode, float fact)
{
    int rValue = 0, i, j;

    if (numExternalNodes > 1)  {
        Vector *v = new Vector [numExternalNodes];

        // first determine the end points of the element based on
        // the display factor (a measure of the distorted image)
        for (i=0; i<numExternalNodes; i++)  {
            const Vector &endCrd = theNodes[i]->getCrds();
            const Vector &endDisp = theNodes[i]->getDisp();
            int numCrds = endCrd.Size();
            for (j=0; j<numCrds; i++)
                v[i](j) = endCrd(j) + endDisp(j)*fact;
        }

        for (i=0; i<numExternalNodes-1; i++)
            rValue += theViewer.drawLine (v[i], v[i+1], 1.0, 1.0);
        //rValue += theViewer.drawLine (v[i+1], v[0], 1.0, 1.0);
    }
    
    return rValue;
}


void Adapter::Print(OPS_Stream &s, int flag)
{
    int i;
    if (flag == 0)  {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: Adapter";
        for (i=0; i<numExternalNodes; i++ )
            s << ", Node" << i+1 << ": " << connectedExternalNodes(i);
        s << endln;
        s << "  kb: " << kb << endln;
        s << "  ipPort: " << ipPort << endln;
        if (mb != 0)
            s << "  mb: " << *mb << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    } else if (flag == 1)  {
        // does nothing
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
    
    // daq basic displacements
    else if (strcmp(argv[0],"daqDisp") == 0 ||
        strcmp(argv[0],"daqDisplacement") == 0 ||
        strcmp(argv[0],"daqDisplacements") == 0)
    {
        for (int i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"dbm%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 6, Vector(numBasicDOF));
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
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = *ctrlDisp;
        }
        return 0;      
        
    case 6:  // daq basic displacements
        if (eleInformation.theVector != 0)  {
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
    if ((sizes(0) != 0 && sizes(0) != numBasicDOF) ||
        (sizes(3) != 0 && sizes(3) != numBasicDOF) ||
        (sizes(5) != 0 && sizes(5) != numBasicDOF) ||
        (sizes(8) != 0 && sizes(8) != numBasicDOF))  {
        opserr << "Adapter::Adapter() - wrong data sizes != "
            << numBasicDOF << " received\n";
        return -3;
    }
    
    // allocate memory for the receive vectors
    int id = 1;
    rData = new double [sizes(10)];
    recvData = new Vector(rData, sizes(10));
    if (sizes(0) != 0)  {
        ctrlDisp = new Vector(&rData[id], sizes(0));
        id += sizes(0);
    }
    if (sizes(3) != 0)  {
        ctrlForce = new Vector(&rData[id], sizes(3));
        id += sizes(3);
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
    if (sizes(8) != 0)  {
        daqForce = new Vector(&sData[id], sizes(8));
        id += sizes(8);
    }
    sendData->Zero();
    
    opserr << "\nAdapter element " << this->getTag()
        << " now running...\n";
    
    return 0;
}
