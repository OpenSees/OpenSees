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

// $Date: 2007-12-06 20:33:18 $

// $URL: $



// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)

// Created: 11/06

// Revision: A

//

// Description: This file contains the implementation of the GenericClient class.



#include "GenericClient.h"



#include <Domain.h>

#include <Node.h>

#include <Channel.h>

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

    int port, char *machineInetAddr, int ssl, int dataSize)

    : Element(tag, ELE_TAG_GenericClient),

    connectedExternalNodes(nodes), basicDOF(1),

    numExternalNodes(0), numDOF(0), numBasicDOF(0),

    theChannel(0), sData(0), sendData(0), rData(0), recvData(0),

    db(0), vb(0), ab(0), t(0), qMeas(0), rMatrix(0),

    dbTarg(1), dbPast(1), initStiffFlag(false), massFlag(false)

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

        numBasicDOF += dof[i].Size();

        theDOF[i] = dof[i];

    }

    

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

    theChannel->setUpConnection();



    // set the data size for the experimental element

    int intData[2*5+1];

    ID idData(intData, 2*5+1);

    ID *sizeCtrl = new ID(intData, 5);

    ID *sizeDaq = new ID(&intData[5], 5);

    idData.Zero();

        

    (*sizeCtrl)[0]  = numBasicDOF;

    (*sizeCtrl)[1]   = numBasicDOF;

    (*sizeCtrl)[2] = numBasicDOF;

    (*sizeCtrl)[4]  = 1;

    

    (*sizeDaq)[3]  = numBasicDOF;

    

    if (dataSize < 1+3*numBasicDOF+1) dataSize = 1+3*numBasicDOF+1;

    if (dataSize < numBasicDOF*numBasicDOF) dataSize = numBasicDOF*numBasicDOF;

    intData[2*5] = dataSize;



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

    qMeas = new Vector(&rData[id], numBasicDOF);

    recvData->Zero();



    // allocate memory for the receive matrix

    rMatrix = new Matrix(rData, numBasicDOF, numBasicDOF);



    // set the vector and matrix sizes and zero them

    basicDOF.resize(numBasicDOF);

    basicDOF.Zero();

    dbTarg.resize(numBasicDOF);

    dbTarg.Zero();

    dbPast.resize(numBasicDOF);

    dbPast.Zero();

}





// delete must be invoked on any objects created by the object.

GenericClient::~GenericClient()

{

    sData[0] = 99;

    theChannel->sendVector(0, 0, *sendData, 0);



    // invoke the destructor on any objects created by the object

    // that the object still holds a pointer to

    if (theDOF != 0)

        delete [] theDOF;

    if (theNodes != 0)

        delete [] theNodes;



    if (db != 0)

        delete db;

    if (vb != 0)

        delete vb;

    if (ab != 0)

        delete ab;

    if (t != 0)

        delete t;



    if (qMeas != 0)

        delete qMeas;

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

    

    sData[0] = 6;

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

 

    if ((*db) != dbPast)  {

        // save the displacements

        dbPast = (*db);

        // set the trial response at the element

        sData[0] = 3;

        rValue += theChannel->sendVector(0, 0, *sendData, 0);

    }

    

    return rValue;

}





const Matrix& GenericClient::getTangentStiff()

{

    // zero the matrices

    theMatrix.Zero();

    rMatrix->Zero();



    sData[0] = 12;

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



        sData[0] = 13;

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



    sData[0] = 14;

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



        sData[0] = 15;

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

    

    // get measured resisting forces

    sData[0] = 10;

    theChannel->sendVector(0, 0, *sendData, 0);

    theChannel->recvVector(0, 0, *recvData, 0);

   

    // save corresponding target displacements for recorder

    dbTarg = (*db);



    // determine resisting forces in global system

    theVector.Assemble(*qMeas, basicDOF);

    

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





int GenericClient::sendSelf(int commitTag, Channel &theChannel)

{

    // has not been implemented yet.....

    return 0;

}





int GenericClient::recvSelf(int commitTag, Channel &theChannel,

    FEM_ObjectBroker &theBroker)

{

    // has not been implemented yet.....

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

        s << "  type: GenericClient";

        for (i=0; i<numExternalNodes; i++ )

            s << "  Node" << i+1 << ": " << connectedExternalNodes(i);

        s << endln;

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

    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||

        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)

    {

         for (i=0; i<numDOF; i++)  {

            sprintf(outputData,"P%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 2, theVector);

    }

    // local forces

    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)

    {

        for (i=0; i<numDOF; i++)  {

            sprintf(outputData,"p%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 3, theVector);

    }

    // forces in basic system

    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)

    {

        for (i=0; i<numBasicDOF; i++)  {

            sprintf(outputData,"q%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 4, Vector(numBasicDOF));

    }

    // target basic displacements

    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 

        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||

        strcmp(argv[0],"targetDisplacement") == 0 || strcmp(argv[0],"targetDisplacements") == 0)

    {

        for (i=0; i<numBasicDOF; i++)  {

            sprintf(outputData,"db%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 5, Vector(numBasicDOF));

    }

    // measured basic displacements

    /*else if (strcmp(argv[0],"measuredDisplacement") == 0 || 

        strcmp(argv[0],"measuredDisplacements") == 0)

    {

        for (int i=0; i<numBasicDOF; i++)  {

            sprintf(outputData,"dbm%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 6, Vector(numBasicDOF));

    }

    // measured basic velocities

    else if (strcmp(argv[0],"measuredVelocity") == 0 || 

        strcmp(argv[0],"measuredVelocities") == 0)

    {

        for (int i=0; i<numBasicDOF; i++)  {

            sprintf(outputData,"vbm%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 7, Vector(numBasicDOF));

    }

    // measured basic accelerations

    else if (strcmp(argv[0],"measuredAcceleration") == 0 || 

        strcmp(argv[0],"measuredAccelerations") == 0)

    {

        for (int i=0; i<numBasicDOF; i++)  {

            sprintf(outputData,"abm%d",i+1);

            output.tag("ResponseType",outputData);

        }

        theResponse = new ElementResponse(this, 8, Vector(numBasicDOF));

    }*/



    output.endTag(); // ElementOutput



    return theResponse;

}





int GenericClient::getResponse(int responseID, Information &eleInformation)

{    

    switch (responseID)  {

    case -1:

        return -1;

        

    case 1:  // initial stiffness

        if (eleInformation.theMatrix != 0)  {

            *(eleInformation.theMatrix) = this->getInitialStiff();

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

        

    case 4:  // forces in basic system

        if (eleInformation.theVector != 0)  {

            *(eleInformation.theVector) = (*qMeas);

        }

        return 0;      

        

    case 5:  // target basic displacements

        if (eleInformation.theVector != 0)  {

            *(eleInformation.theVector) = dbTarg;

        }

        return 0;      

        

    /*case 6:  // measured basic displacements

        if (eleInformation.theVector != 0)  {

            sData[0] = OF_RemoteTest_getDisp;

            theChannel->sendVector(0, 0, *sendData, 0);

            theChannel->recvVector(0, 0, *recvData, 0);

            *(eleInformation.theVector) = (*dbMeas);

        }

        return 0;



    case 7:  // measured basic velocities

        if (eleInformation.theVector != 0)  {

            sData[0] = OF_RemoteTest_getVel;

            theChannel->sendVector(0, 0, *sendData, 0);

            theChannel->recvVector(0, 0, *recvData, 0);

            *(eleInformation.theVector) = (*vbMeas);

        }

        return 0;



    case 8:  // measured basic accelerations

        if (eleInformation.theVector != 0)  {

            sData[0] = OF_RemoteTest_getAccel;

            theChannel->sendVector(0, 0, *sendData, 0);

            theChannel->recvVector(0, 0, *recvData, 0);

            *(eleInformation.theVector) = (*abMeas);

        }

        return 0;*/



    default:

        return -1;

    }

}

