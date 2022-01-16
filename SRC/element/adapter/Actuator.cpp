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
// Description: This file contains the implementation of the Actuator class.

#include "Actuator.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <TCP_Socket.h>
#include <UDP_Socket.h>
#ifdef SSL
    #include <TCP_SocketSSL.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>


// initialize the class wide variables
Matrix Actuator::ActuatorM2(2,2);
Matrix Actuator::ActuatorM4(4,4);
Matrix Actuator::ActuatorM6(6,6);
Matrix Actuator::ActuatorM12(12,12);
Vector Actuator::ActuatorV2(2);
Vector Actuator::ActuatorV4(4);
Vector Actuator::ActuatorV6(6);
Vector Actuator::ActuatorV12(12);

void * OPS_ADD_RUNTIME_VPV(OPS_Actuator)
{
    // check the number of arguments is correct
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element actuator eleTag iNode jNode EA ipPort <-ssl> <-udp> <-doRayleigh> <-rho rho>\n";
        return 0;
    }
    
    int ndm = OPS_GetNDM();
    
    // get the id and end nodes
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid actuator int inputs" << endln;
	return 0;
    }
    int tag = idata[0];
    int iNode = idata[1];
    int jNode = idata[2];
    
    // stiffness
    double EA;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &EA) < 0) {
	opserr << "WARNING invalid actuator EA" << endln;
	return 0;
    }
    
    // ipPort
    int ipPort;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &ipPort) < 0) {
	opserr << "WARNING invalid actuator ipPort" << endln;
	return 0;
    }
    
    // options
    int ssl = 0, udp = 0;
    int doRayleigh = 0;
    double rho = 0.0;
    
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* flag = OPS_GetString();
		if (strcmp(flag, "-ssl") == 0) {
			ssl = 1; udp = 0;
		}
		else if (strcmp(flag, "-udp") == 0) {
			udp = 1; ssl = 0;
		}
		else if (strcmp(flag, "-doRayleigh") == 0) {
			doRayleigh = 1;
		}
		else if (strcmp(flag, "-rho") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				numdata = 1;
				if (OPS_GetDoubleInput(&numdata, &rho) < 0) {
					opserr << "WARNING invalid rho\n";
					opserr << "actuator element: " << tag << endln;
					return 0;
				}
			}
		}
	}
    
    // now create the actuator and add it to the Domain
    return new Actuator(tag, ndm, iNode, jNode, EA, ipPort,
			ssl, udp, doRayleigh, rho);
}


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
Actuator::Actuator(int tag, int dim, int Nd1, int Nd2,
    double ea, int ipport, int _ssl, int _udp, int addRay, double r)
    : Element(tag, ELE_TAG_Actuator), numDIM(dim), numDOF(0),
    connectedExternalNodes(2), EA(ea), ipPort(ipport), ssl(_ssl),
    udp(_udp), addRayleigh(addRay), rho(r), L(0.0),
    tPast(0.0), theMatrix(0), theVector(0), theLoad(0), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
{    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "Actuator::Actuator() - element: "
            <<  tag << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    // zero direction cosines
    cosX[0] = 0.0;
    cosX[1] = 0.0;
    cosX[2] = 0.0;
}

// invoked by a FEM_ObjectBroker - blank object that recvSelf
// needs to be invoked upon
Actuator::Actuator()
    : Element(0, ELE_TAG_Actuator), numDIM(0), numDOF(0),
    connectedExternalNodes(2), EA(0.0), ipPort(0), ssl(0),
    udp(0), addRayleigh(0), rho(0.0), L(0.0), tPast(0.0),
    theMatrix(0), theVector(0), theLoad(0), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "Actuator::Actuator() - "
            <<  "failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    // zero direction cosines
    cosX[0] = 0.0;
    cosX[1] = 0.0;
    cosX[2] = 0.0;
}


// delete must be invoked on any objects created by the object.
Actuator::~Actuator()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theLoad != 0)
        delete theLoad;
    
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


int Actuator::getNumExternalNodes() const
{
    return 2;
}


const ID& Actuator::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** Actuator::getNodePtrs() 
{
    return theNodes;
}


int Actuator::getNumDOF() 
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void Actuator::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        L = 0.0;
        return;
    }
    
    // set default values for error conditions
    numDOF = 2;
    theMatrix = &ActuatorM2;
    theVector = &ActuatorV2;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "Actuator::setDomain() - Nd1: "
                << Nd1 << "does not exist in the model for ";
        } else  {
            opserr << "Actuator::setDomain() - Nd2: "
                << Nd2 << "does not exist in the model for ";
        }
        opserr << "Actuator ele: " << this->getTag() << endln;
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2)  {
        opserr <<"Actuator::setDomain(): nodes " << Nd1 << " and " << Nd2
            << "have differing dof at ends for element: " << this->getTag() << endln;
        
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 1 && dofNd1 == 1)  {
        numDOF = 2;    
        theMatrix = &ActuatorM2;
        theVector = &ActuatorV2;
    }
    else if (numDIM == 2 && dofNd1 == 2)  {
        numDOF = 4;
        theMatrix = &ActuatorM4;
        theVector = &ActuatorV4;
    }
    else if (numDIM == 2 && dofNd1 == 3)  {
        numDOF = 6;	
        theMatrix = &ActuatorM6;
        theVector = &ActuatorV6;
    }
    else if (numDIM == 3 && dofNd1 == 3)  {
        numDOF = 6;	
        theMatrix = &ActuatorM6;
        theVector = &ActuatorV6;
    }
    else if (numDIM == 3 && dofNd1 == 6)  {
        numDOF = 12;	    
        theMatrix = &ActuatorM12;
        theVector = &ActuatorV12;
    }
    else  {
        opserr <<"Actuator::setDomain() - can not handle "
            << numDIM << " dofs at nodes in " << dofNd1  << " d problem\n";
        
        return;
    }
    
    if (!theLoad)
        theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF)  {
        delete theLoad;
        theLoad = new Vector(numDOF);
    }
    
    if (!theLoad)  {
        opserr << "Actuator::setDomain() - element: " << this->getTag()
            << " out of memory creating vector of size: " << numDOF << endln;
        
        return;
    }
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    // initialize the cosines
    cosX[0] = cosX[1] = cosX[2] = 0.0;
    for (int i=0; i<numDIM; i++)
        cosX[i] = end2Crd(i)-end1Crd(i);
    
    // get initial length
    L = sqrt(cosX[0]*cosX[0] + cosX[1]*cosX[1] + cosX[2]*cosX[2]);
    if (L == 0.0)  {
        opserr <<"Actuator::setDomain() - element: "
            << this->getTag() << " has zero length\n";
        return;
    }
    
    // set global orientations
    cosX[0] /= L;
    cosX[1] /= L;
    cosX[2] /= L;
}


int Actuator::commitState()
{
    int errCode = 0;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int Actuator::revertToLastCommit()
{
    opserr << "Actuator::revertToLastCommit() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to last commit. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int Actuator::revertToStart()
{
    opserr << "Actuator::revertToStart() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to start. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int Actuator::update()
{
    if (theChannel == 0)  {
        if (this->setupConnection() != 0)  {
            opserr << "Actuator::update() - "
                << "failed to setup connection\n";
            return -1;
        }
    }
    
    // determine dsp in basic system
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    db(0) = 0.0;
    for (int i=0; i<numDIM; i++)
        db(0) += (dsp2(i)-dsp1(i))*cosX[i];
    
    return 0;
}


const Matrix& Actuator::getTangentStiff()
{
    // zero the matrix
    theMatrix->Zero();
    
    // transform the stiffness from the basic to the global system
    int numDOF2 = numDOF/2;
    double temp;
    for (int i=0; i<numDIM; i++)  {
        for (int j=0; j<numDIM; j++)  {
            temp = cosX[i]*cosX[j]*EA/L;
            (*theMatrix)(i,j) = temp;
            (*theMatrix)(i+numDOF2,j) = -temp;
            (*theMatrix)(i,j+numDOF2) = -temp;
            (*theMatrix)(i+numDOF2,j+numDOF2) = temp;
        }
    }
    
    return *theMatrix;
}


const Matrix& Actuator::getInitialStiff()
{
    // zero the matrix
    theMatrix->Zero();
    
    // transform the stiffness from the basic to the global system
    int numDOF2 = numDOF/2;
    double temp;
    for (int i=0; i<numDIM; i++)  {
        for (int j=0; j<numDIM; j++)  {
            temp = cosX[i]*cosX[j]*EA/L;
            (*theMatrix)(i,j) = temp;
            (*theMatrix)(i+numDOF2,j) = -temp;
            (*theMatrix)(i,j+numDOF2) = -temp;
            (*theMatrix)(i+numDOF2,j+numDOF2) = temp;
        }
    }
    
    return *theMatrix;
}


const Matrix& Actuator::getDamp()
{
    // zero the matrix
    theMatrix->Zero();
    
    // call base class to setup Rayleigh damping
    if (addRayleigh == 1)
        (*theMatrix) = this->Element::getDamp();
    
    return *theMatrix;
}


const Matrix& Actuator::getMass()
{   
    // zero the matrix
    theMatrix->Zero();
    
    // form mass matrix
    if (L != 0.0 && rho != 0.0)  {
        double m = 0.5*rho*L;
        int numDOF2 = numDOF/2;
        for (int i=0; i<numDIM; i++)  {
            (*theMatrix)(i,i) = m;
            (*theMatrix)(i+numDOF2,i+numDOF2) = m;
        }
    }
    
    return *theMatrix;
}


void Actuator::zeroLoad()
{
    theLoad->Zero();
}


int Actuator::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
    opserr <<"Actuator::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int Actuator::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for a quick return
    if (L == 0.0 || rho == 0.0)
        return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    int nodalDOF = numDOF/2;
    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size())  {
        opserr <<"Actuator::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    double m = 0.5*rho*L;
    for (int i=0; i<numDIM; i++)  {
        double val1 = Raccel1(i);
        double val2 = Raccel2(i);
        
        // perform - fact * M*(R * accel) // remember M a diagonal matrix
        val1 *= -m;
        val2 *= -m;
        
        (*theLoad)(i) += val1;
        (*theLoad)(i+nodalDOF) += val2;
    }
    
    return 0;
}


const Vector& Actuator::getResistingForce()
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
                opserr << "Actuator::getResistingForce() - "
                    << "wrong action received: expecting 3 but got "
                    << rData[0] << endln;
            }
            exit(-1);
        }
        
        // save current time
        tPast = t;
    }
    
    // get resisting force in basic system q = k*db + q0 = k*(db - db0)
    q(0) = EA/L*(db(0) - (*ctrlDisp)(0));
    
    // assign daq values for feedback
    (*daqDisp)(0)  = db(0);
    (*daqForce)(0) = -q(0);
    
    // zero the residual
    theVector->Zero();
    
    // determine resisting forces in global system
    int numDOF2 = numDOF/2;
    for (int i=0; i<numDIM; i++)  {
        (*theVector)(i) = -cosX[i]*q(0);
        (*theVector)(i+numDOF2) = cosX[i]*q(0);
    }
    
    return *theVector;
}


const Vector& Actuator::getResistingForceIncInertia()
{
    this->getResistingForce();
    
    // subtract external load
    (*theVector) -= *theLoad;
    
    // add the damping forces from rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add inertia forces from element mass
    if (L != 0.0 && rho != 0.0)  {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();
        
        int numDOF2 = numDOF/2;
        double m = 0.5*rho*L;
        for (int i=0; i<numDIM; i++)  {
            (*theVector)(i) += m * accel1(i);
            (*theVector)(i+numDOF2) += m * accel2(i);
        }
    }
    
    return *theVector;
}


int Actuator::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(13);
    data(0) = this->getTag();
    data(1) = numDIM;
    data(2) = numDOF;
    data(3) = EA;
    data(4) = ipPort;
    data(5) = ssl;
    data(6) = udp;
    data(7) = addRayleigh;
    data(8) = rho;
    data(9) = alphaM;
    data(10) = betaK;
    data(11) = betaK0;
    data(12) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    return 0;
}


int Actuator::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // receive element parameters
    static Vector data(13);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    numDIM = (int)data(1);
    numDOF = (int)data(2);
    EA = data(3);
    ipPort = (int)data(4);
    ssl = (int)data(5);
    udp = (int)data(6);
    addRayleigh = (int)data(7);
    rho = data(8);
    alphaM = data(9);
    betaK = data(10);
    betaK0 = data(11);
    betaKc = data(12);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    return 0;
}


int Actuator::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void Actuator::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE)  {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: Actuator, iNode: " << connectedExternalNodes(0)
            << ", jNode: " << connectedExternalNodes(1) << endln;
        s << "  EA: " << EA << ", L: " << L << endln;
        s << "  ipPort: " << ipPort << endln;
        s << "  addRayleigh: " << addRayleigh;
        s << "  mass per unit length: " << rho << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"Actuator\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"EA\": " << EA << ", ";
        s << "\"L\": " << L << ", ";
        s << "\"ipPort\": " << ipPort << ", ";
        s << "\"addRayleigh\": " << addRayleigh << ", ";
        s << "\"massperlength\": " << rho << "}";
    }
}


Response* Actuator::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","Actuator");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    char outputData[10];
    
    // global forces
    if (strcmp(argv[0],"force") == 0 ||
        strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0)
    {
        for (int i=0; i<numDOF; i++)  {
            sprintf(outputData,"P%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, *theVector);
    }
    
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
    {
        for (int i=0; i<numDOF; i++)  {
            sprintf(outputData,"p%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 3, *theVector);
    }
    
    // basic force
    else if (strcmp(argv[0],"basicForce") == 0 ||
        strcmp(argv[0],"basicForces") == 0 ||
        strcmp(argv[0],"daqForce") == 0 ||
        strcmp(argv[0],"daqForces") == 0)
    {
        output.tag("ResponseType","q1");
        
        theResponse = new ElementResponse(this, 4, Vector(1));
    }
    
    // ctrl basic displacement
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
        output.tag("ResponseType","db1");
        
        theResponse = new ElementResponse(this, 5, Vector(1));
    }
    
    // daq basic displacement
    else if (strcmp(argv[0],"daqDisp") == 0 ||
        strcmp(argv[0],"daqDisplacement") == 0 ||
        strcmp(argv[0],"daqDisplacements") == 0)
    {
        output.tag("ResponseType","dbm1");
        
        theResponse = new ElementResponse(this, 6, Vector(1));
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int Actuator::getResponse(int responseID, Information &eleInformation)
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
            theVector->Zero();
            // Axial
            (*theVector)(0)        = -q(0);
            (*theVector)(numDOF/2) =  q(0);
            
            *(eleInformation.theVector) = *theVector;
        }
        return 0;
        
    case 4:  // basic force
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = q;
        }
        return 0;
        
    case 5:  // ctrl basic displacement
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = *ctrlDisp;
        }
        return 0;
        
    case 6:  // daq basic displacement
        if (eleInformation.theVector != 0)  {
            *(eleInformation.theVector) = *daqDisp;
        }
        return 0;
        
    default:
        return 0;
    }
}


int Actuator::setupConnection()
{
    // setup the connection
    if (udp)
        theChannel = new UDP_Socket(ipPort);
#ifdef SSL
    else if (ssl)
        theChannel = new TCP_SocketSSL(ipPort);
#endif
    else
        theChannel = new TCP_Socket(ipPort);
    
    if (theChannel != 0)  {
        opserr << "\nChannel successfully created: "
            << "Waiting for ECSimAdapter experimental control...\n";
    } else {
        opserr << "Actuator::setupConnection() - "
            << "could not create channel\n";
        return -1;
    }
    if (theChannel->setUpConnection() != 0)  {
        opserr << "Actuator::setupConnection() - "
            << "failed to setup connection\n";
        return -2;
    }
    
    // get the data sizes and check values
    // sizes = {ctrlDisp, ctrlVel, ctrlAccel, ctrlForce, ctrlTime,
    //          daqDisp,  daqVel,  daqAccel,  daqForce,  daqTime,  dataSize}
    ID sizes(11);
    theChannel->recvID(0, 0, sizes, 0);
    if (sizes(0) > 1 || sizes(3) > 1 || sizes(5) > 1 || sizes(8) > 1)  {
        opserr << "Actuator::setupConnection() - "
            << "wrong data sizes > 1 received\n";
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
    
    opserr << "\nActuator element " << this->getTag()
        << " now running...\n";
    
    return 0;
}
