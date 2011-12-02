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
// $Date: 2009-06-02 21:09:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/generic/GenericClient.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/06
// Revision: A
//
// Description: This file contains the implementation of the GenericClient class.

#include "GenericClient.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <TCP_Socket.h>
#ifdef SSL
    #include <TCP_SocketSSL.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>


// initialize the class wide variables
Matrix GenericClient::theMatrix(1,1);
Matrix GenericClient::theInitStiff(1,1);
Matrix GenericClient::theMass(1,1);
Vector GenericClient::theVector(1);
Vector GenericClient::theLoad(1);


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
GenericClient::GenericClient(int tag, ID nodes, ID *dof,
    int _port, char *machineinetaddr, int _ssl, int datasize)
    : Element(tag, ELE_TAG_GenericClient),
    connectedExternalNodes(nodes), basicDOF(1),
    numExternalNodes(0), numDOF(0), numBasicDOF(0),
    port(_port), machineInetAddr(0), ssl(_ssl), dataSize(datasize),
    theChannel(0), sData(0), sendData(0), rData(0), recvData(0),
    db(0), vb(0), ab(0), t(0), qDaq(0), rMatrix(0),
    dbCtrl(1), vbCtrl(1), abCtrl(1), dbPast(1), tPast(0.0),
    initStiffFlag(false), massFlag(false)
{    
    // initialize nodes
    numExternalNodes = connectedExternalNodes.Size();
    theNodes = new Node* [numExternalNodes];
    if (!theNodes)  {
        opserr << "GenericClient::GenericClient() "
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
        opserr << "GenericClient::GenericClient() "
            << "- failed to create dof array\n";
        exit(-1);
    }
    numBasicDOF = 0;
    for (i=0; i<numExternalNodes; i++)  {
        theDOF[i] = dof[i];
        numBasicDOF += theDOF[i].Size();
    }
    
    // save ipAddress
    machineInetAddr = machineinetaddr;
    
    // set the vector sizes and zero them
    basicDOF.resize(numBasicDOF);
    basicDOF.Zero();
    dbCtrl.resize(numBasicDOF);
    dbCtrl.Zero();
    vbCtrl.resize(numBasicDOF);
    vbCtrl.Zero();
    abCtrl.resize(numBasicDOF);
    abCtrl.Zero();
    dbPast.resize(numBasicDOF);
    dbPast.Zero();
}


// invoked by a FEM_ObjectBroker - blank object that recvSelf
// needs to be invoked upon
GenericClient::GenericClient()
    : Element(0, ELE_TAG_GenericClient),
    connectedExternalNodes(1), basicDOF(1),
    numExternalNodes(0), numDOF(0), numBasicDOF(0),
    port(0), machineInetAddr(0), ssl(0), dataSize(0),
    theChannel(0), sData(0), sendData(0), rData(0), recvData(0),
    db(0), vb(0), ab(0), t(0), qDaq(0), rMatrix(0),
    dbCtrl(1), vbCtrl(1), abCtrl(1), dbPast(1), tPast(0.0),
    initStiffFlag(false), massFlag(false)
{    
    // initialize variables
    theNodes = 0;
    theDOF = 0;
}


// delete must be invoked on any objects created by the object.
GenericClient::~GenericClient()
{
    // terminate remote process
    if (theChannel != 0)  {
        sData[0] = RemoteTest_DIE;
        theChannel->sendVector(0, 0, *sendData, 0);
    }
    
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theNodes != 0)
        delete [] theNodes;
    if (theDOF != 0)
        delete [] theDOF;
    if (machineInetAddr != 0)
        delete [] machineInetAddr;
    
    if (db != 0)
        delete db;
    if (vb != 0)
        delete vb;
    if (ab != 0)
        delete ab;
    if (t != 0)
        delete t;
    
    if (qDaq != 0)
        delete qDaq;
    if (rMatrix != 0)
        delete rMatrix;
    
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


int GenericClient::getNumExternalNodes() const
{
    return numExternalNodes;
}


const ID& GenericClient::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** GenericClient::getNodePtrs() 
{
    return theNodes;
}


int GenericClient::getNumDOF() 
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void GenericClient::setDomain(Domain *theDomain)
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
            opserr << "GenericClient::setDomain() - Nd" << i << ": " 
                << connectedExternalNodes(i) << " does not exist in the "
                << "model for GenericClient ele: " << this->getTag() << endln;
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
    theInitStiff.resize(numDOF,numDOF);
    theInitStiff.Zero();
    theMass.resize(numDOF,numDOF);
    theMass.Zero();
    theVector.resize(numDOF);
    theVector.Zero();
    theLoad.resize(numDOF);
    theLoad.Zero();
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
}   	 


int GenericClient::commitState()
{
    int rValue = 0;
    
    sData[0] = RemoteTest_commitState;
    rValue += theChannel->sendVector(0, 0, *sendData, 0);
    
    return rValue;
}


int GenericClient::revertToLastCommit()
{
    opserr << "GenericClient::revertToLastCommit() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to last commit. This element "
        << "shadows an experimental element." 
        << endln;
    
    return -1;
}


int GenericClient::revertToStart()
{
    opserr << "GenericClient::revertToStart() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to start. This element "
        << "shadows an experimental element." 
        << endln;
    
    return -1;
}


int GenericClient::update()
{
    int rValue = 0;

    if (theChannel == 0)  {
        if (this->setupConnection() != 0)  {
            opserr << "GenericClient::update() - "
                << "failed to setup connection\n";
            return -1;
        }
    }

    // get current time
    Domain *theDomain = this->getDomain();
    (*t)(0) = theDomain->getCurrentTime();

    // assemble response vectors
    int ndim = 0, i;
    db->Zero(); vb->Zero(); ab->Zero();

    for (i=0; i<numExternalNodes; i++)  {
        Vector disp = theNodes[i]->getTrialDisp();
        Vector vel = theNodes[i]->getTrialVel();
        Vector accel = theNodes[i]->getTrialAccel();
        db->Assemble(disp(theDOF[i]), ndim);
        vb->Assemble(vel(theDOF[i]), ndim);
        ab->Assemble(accel(theDOF[i]), ndim);
        ndim += theDOF[i].Size();
    }
 
    if ((*db) != dbPast || (*t)(0) != tPast)  {
        // save the displacements and the time
        dbPast = (*db);
        tPast = (*t)(0);
        // set the trial response at the element
        sData[0] = RemoteTest_setTrialResponse;
        rValue += theChannel->sendVector(0, 0, *sendData, 0);
    }
    
    return rValue;
}


const Matrix& GenericClient::getTangentStiff()
{
    // zero the matrices
    theMatrix.Zero();
    rMatrix->Zero();

    sData[0] = RemoteTest_getTangentStiff;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    theMatrix.Assemble(*rMatrix,basicDOF,basicDOF);
    
    return theMatrix;
}


const Matrix& GenericClient::getInitialStiff()
{
    if (initStiffFlag == false)  {
        // zero the matrices
        theInitStiff.Zero();
        rMatrix->Zero();

        sData[0] = RemoteTest_getInitialStiff;
        theChannel->sendVector(0, 0, *sendData, 0);
        theChannel->recvVector(0, 0, *recvData, 0);
        
        theInitStiff.Assemble(*rMatrix,basicDOF,basicDOF);        
        initStiffFlag = true;
    }
    
    return theInitStiff;
}


/*const Matrix& GenericClient::getDamp()
{
    // zero the matrices
    theMatrix.Zero();
    rMatrix->Zero();

    sData[0] = RemoteTest_getDamp;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    theMatrix.Assemble(*rMatrix,basicDOF,basicDOF);
    
    return theMatrix;
}*/


const Matrix& GenericClient::getMass()
{
    if (massFlag == false)  {
        // zero the matrices
        theMass.Zero();
        rMatrix->Zero();

        sData[0] = RemoteTest_getMass;
        theChannel->sendVector(0, 0, *sendData, 0);
        theChannel->recvVector(0, 0, *recvData, 0);
        
        theMass.Assemble(*rMatrix,basicDOF,basicDOF);
        massFlag = true;
    }
    
    return theMass;
}


void GenericClient::zeroLoad()
{
    theLoad.Zero();
}


int GenericClient::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
    opserr <<"GenericClient::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int GenericClient::addInertiaLoadToUnbalance(const Vector &accel)
{    
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


const Vector& GenericClient::getResistingForce()
{    
    // zero the residual
    theVector.Zero();
    
    // get daq resisting forces
    sData[0] = RemoteTest_getForce;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    // save corresponding ctrl response for recorder
    dbCtrl = (*db);
    vbCtrl = (*vb);
    abCtrl = (*ab);
    
    // determine resisting forces in global system
    theVector.Assemble(*qDaq, basicDOF);
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    return theVector;
}


const Vector& GenericClient::getResistingForceIncInertia()
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


/*const Vector& GenericClient::getTime()
{	
    sData[0] = RemoteTest_getTime;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    return *tDaq;
}


const Vector& GenericClient::getBasicDisp()
{	
    sData[0] = RemoteTest_getDisp;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    return *dbDaq;
}


const Vector& GenericClient::getBasicVel()
{	
    sData[0] = RemoteTest_getVel;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    return *vbDaq;
}


const Vector& GenericClient::getBasicAccel()
{	
    sData[0] = RemoteTest_getAccel;
    theChannel->sendVector(0, 0, *sendData, 0);
    theChannel->recvVector(0, 0, *recvData, 0);
    
    return *abDaq;
}*/


int GenericClient::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static ID idData(6);
    idData(0) = this->getTag();
    idData(1) = numExternalNodes;
    idData(2) = port;
    idData(3) = strlen(machineInetAddr);
    idData(4) = ssl;
    idData(5) = dataSize;
    sChannel.sendID(0, commitTag, idData);

    // send the end nodes and dofs
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    for (int i=0; i<numExternalNodes; i++)
        sChannel.sendID(0, commitTag, theDOF[i]);
    
    // send the ip-address
    Message theMessage(machineInetAddr, strlen(machineInetAddr));
    sChannel.sendMsg(0, commitTag, theMessage);
    
    return 0;
}


int GenericClient::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete dynamic memory
    if (theNodes != 0)
        delete [] theNodes;
    if (theDOF != 0)
        delete [] theDOF;
    if (machineInetAddr != 0)
        delete [] machineInetAddr;
    
    // receive element parameters
    static ID idData(6);
    rChannel.recvID(0, commitTag, idData);    
    this->setTag(idData(0));
    numExternalNodes = idData(1);
    port = idData(2);
    machineInetAddr = new char [idData(3) + 1];
    ssl = idData(4);
    dataSize = idData(5);
    
    // initialize nodes and receive them
    connectedExternalNodes.resize(numExternalNodes);
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    theNodes = new Node* [numExternalNodes];
    if (!theNodes)  {
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
    if (!theDOF)  {
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
    
    // receive the ip-address
    Message theMessage(machineInetAddr, strlen(machineInetAddr));  
    rChannel.recvMsg(0, commitTag, theMessage);
    
    // set the vector sizes and zero them
    basicDOF.resize(numBasicDOF);
    basicDOF.Zero();
    dbCtrl.resize(numBasicDOF);
    dbCtrl.Zero();
    vbCtrl.resize(numBasicDOF);
    vbCtrl.Zero();
    abCtrl.resize(numBasicDOF);
    abCtrl.Zero();
    dbPast.resize(numBasicDOF);
    dbPast.Zero();
    
    return 0;
}


int GenericClient::displaySelf(Renderer &theViewer,
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


void GenericClient::Print(OPS_Stream &s, int flag)
{
    int i;
    if (flag == 0)  {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: GenericClient" << endln;
        for (i=0; i<numExternalNodes; i++ )
            s << "  Node" << i+1 << ": " << connectedExternalNodes(i);
        s << endln;
        s << "  ipAddress: " << machineInetAddr
            << ", ipPort: " << port << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    } else if (flag == 1)  {
        // does nothing
    }
}


Response* GenericClient::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    int i;
    char outputData[10];
    
    output.tag("ElementOutput");
    output.attr("eleType","GenericClient");
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
        theResponse = new ElementResponse(this, 1, theVector);
    }
    
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
    {
        for (i=0; i<numDOF; i++)  {
            sprintf(outputData,"p%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, theVector);
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
        theResponse = new ElementResponse(this, 3, Vector(numBasicDOF));
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
        theResponse = new ElementResponse(this, 4, Vector(numBasicDOF));
    }
    
    // ctrl basic velocities
    else if (strcmp(argv[0],"ctrlVel") == 0 ||
        strcmp(argv[0],"ctrlVelocity") == 0 ||
        strcmp(argv[0],"ctrlVelocities") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"vb%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 5, Vector(numBasicDOF));
    }
    
    // ctrl basic accelerations
    else if (strcmp(argv[0],"ctrlAccel") == 0 ||
        strcmp(argv[0],"ctrlAcceleration") == 0 ||
        strcmp(argv[0],"ctrlAccelerations") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"ab%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 6, Vector(numBasicDOF));
    }
    
    /* daq basic displacements
    else if (strcmp(argv[0],"daqDisp") == 0 ||
        strcmp(argv[0],"daqDisplacement") == 0 ||
        strcmp(argv[0],"daqDisplacements") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"dbm%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 7, Vector(numBasicDOF));
    }
    
    // daq basic velocities
    else if (strcmp(argv[0],"daqVel") == 0 ||
        strcmp(argv[0],"daqVelocity") == 0 ||
        strcmp(argv[0],"daqVelocities") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"vbm%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 8, Vector(numBasicDOF));
    }
    
    // daq basic accelerations
    else if (strcmp(argv[0],"daqAccel") == 0 ||
        strcmp(argv[0],"daqAcceleration") == 0 ||
        strcmp(argv[0],"daqAccelerations") == 0)
    {
        for (i=0; i<numBasicDOF; i++)  {
            sprintf(outputData,"abm%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 9, Vector(numBasicDOF));
    }*/
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int GenericClient::getResponse(int responseID, Information &eleInfo)
{    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 3:  // basic forces
        return eleInfo.setVector(*qDaq);
        
    case 4:  // ctrl basic displacements
        return eleInfo.setVector(dbCtrl);
        
    case 5:  // ctrl basic velocities
        return eleInfo.setVector(vbCtrl);
        
    case 6:  // ctrl basic accelerations
        return eleInfo.setVector(abCtrl);
        
    /*case 7:  // daq basic displacements
        return eleInfo.setVector(this->getBasicDisp());
        
    case 8:  // daq basic velocities
        return eleInfo.setVector(this->getBasicVel());
        
    case 9:  // daq basic accelerations
        return eleInfo.setVector(this->getBasicAccel());*/
        
    default:
        return -1;
    }
}


int GenericClient::setupConnection()
{
    // setup the connection
    if (!ssl)  {
        if (machineInetAddr == 0)
            theChannel = new TCP_Socket(port, "127.0.0.1");
        else
            theChannel = new TCP_Socket(port, machineInetAddr);
    }
#ifdef SSL
    else  {
        if (machineInetAddr == 0)
            theChannel = new TCP_SocketSSL(port, "127.0.0.1");
        else
            theChannel = new TCP_SocketSSL(port, machineInetAddr);
    }
#endif
    if (!theChannel)  {
        opserr << "GenericClient::GenericClient() "
            << "- failed to create channel\n";
        return -1;
    }
    if (theChannel->setUpConnection() != 0)  {
        opserr << "GenericClient::GenericClient() "
            << "- failed to setup connection\n";
        return -2;
    }
    
    // set the data size for the experimental element
    ID idData(2*5+1);
    idData.Zero();
    
    idData(0) = numBasicDOF;  // sizeCtrl->disp
    idData(1) = numBasicDOF;  // sizeCtrl->vel
    idData(2) = numBasicDOF;  // sizeCtrl->accel
    idData(4) = 1;            // sizeCtrl->time
    
    idData(8) = numBasicDOF;  // sizeDaq->force
    
    if (dataSize < 1+3*numBasicDOF+1) dataSize = 1+3*numBasicDOF+1;
    if (dataSize < numBasicDOF*numBasicDOF) dataSize = numBasicDOF*numBasicDOF;
    idData(10) = dataSize;
    
    theChannel->sendID(0, 0, idData, 0);
    
    // allocate memory for the send vectors
    int id = 1;
    sData = new double [dataSize];
    sendData = new Vector(sData, dataSize);
    db = new Vector(&sData[id], numBasicDOF);
    id += numBasicDOF;
    vb = new Vector(&sData[id], numBasicDOF);
    id += numBasicDOF;
    ab = new Vector(&sData[id], numBasicDOF);
    id += numBasicDOF;
    t = new Vector(&sData[id], 1);
    sendData->Zero();
    
    // allocate memory for the receive vectors
    id = 0;
    rData = new double [dataSize];
    recvData = new Vector(rData, dataSize);
    qDaq = new Vector(&rData[id], numBasicDOF);
    recvData->Zero();
    
    // allocate memory for the receive matrix
    rMatrix = new Matrix(rData, numBasicDOF, numBasicDOF);
    
    return 0;
}
