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
// $Source: /usr/local/cvs/OpenSees/SRC/element/adapter/ActuatorCorot.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/07
// Revision: A
//
// Description: This file contains the implementation of the ActuatorCorot class.

#include "ActuatorCorot.h"

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
Matrix ActuatorCorot::ActuatorCorotM2(2,2);
Matrix ActuatorCorot::ActuatorCorotM4(4,4);
Matrix ActuatorCorot::ActuatorCorotM6(6,6);
Matrix ActuatorCorot::ActuatorCorotM12(12,12);
Vector ActuatorCorot::ActuatorCorotV2(2);
Vector ActuatorCorot::ActuatorCorotV4(4);
Vector ActuatorCorot::ActuatorCorotV6(6);
Vector ActuatorCorot::ActuatorCorotV12(12);


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
ActuatorCorot::ActuatorCorot(int tag, int dim, int Nd1, int Nd2,
    double ea, int ipport, double r)
    : Element(tag, ELE_TAG_ActuatorCorot), numDIM(dim), numDOF(0),
    connectedExternalNodes(2), EA(ea), ipPort(ipport), rho(r), L(0.0), Ln(0.0),
    tPast(0.0), theMatrix(0), theVector(0), theLoad(0), R(3,3), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
{    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ActuatorCorot::ActuatorCorot() - element: "
            <<  tag << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
}


// invoked by a FEM_ObjectBroker - blank object that recvSelf
// needs to be invoked upon
ActuatorCorot::ActuatorCorot()
    : Element(0, ELE_TAG_ActuatorCorot), numDIM(0), numDOF(0),
    connectedExternalNodes(2), EA(0.0), ipPort(0), rho(0.0), L(0.0), Ln(0.0),
    tPast(0.0), theMatrix(0), theVector(0), theLoad(0), R(3,3), db(1), q(1),
    theChannel(0), rData(0), recvData(0), sData(0), sendData(0),
    ctrlDisp(0), ctrlForce(0), daqDisp(0), daqForce(0)
{    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ActuatorCorot::ActuatorCorot() - "
            <<  "failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
}


// delete must be invoked on any objects created by the object.
ActuatorCorot::~ActuatorCorot()
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


int ActuatorCorot::getNumExternalNodes() const
{
    return 2;
}


const ID& ActuatorCorot::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** ActuatorCorot::getNodePtrs() 
{
    return theNodes;
}


int ActuatorCorot::getNumDOF() 
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void ActuatorCorot::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        L  = 0.0;
        Ln = 0.0;
        return;
    }
    
    // set default values for error conditions
    numDOF = 2;
    theMatrix = &ActuatorCorotM2;
    theVector = &ActuatorCorotV2;	
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "ActuatorCorot::setDomain() - Nd1: "
                << Nd1 << "does not exist in the model for ";
        } else  {
            opserr << "ActuatorCorot::setDomain() - Nd2: "
                << Nd2 << "does not exist in the model for ";
        }
        opserr << "ActuatorCorot ele: " << this->getTag() << endln;
        
        return;
    }
    
    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2)  {
        opserr <<"ActuatorCorot::setDomain(): nodes " << Nd1 << " and " << Nd2
            << "have differing dof at ends for element: " << this->getTag() << endln;
        
        return;
    }	
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 1 && dofNd1 == 1)  {
        numDOF = 2;    
        theMatrix = &ActuatorCorotM2;
        theVector = &ActuatorCorotV2;
    }
    else if (numDIM == 2 && dofNd1 == 2)  {
        numDOF = 4;
        theMatrix = &ActuatorCorotM4;
        theVector = &ActuatorCorotV4;	
    }
    else if (numDIM == 2 && dofNd1 == 3)  {
        numDOF = 6;	
        theMatrix = &ActuatorCorotM6;
        theVector = &ActuatorCorotV6;		
    }
    else if (numDIM == 3 && dofNd1 == 3)  {
        numDOF = 6;	
        theMatrix = &ActuatorCorotM6;
        theVector = &ActuatorCorotV6;			
    }
    else if (numDIM == 3 && dofNd1 == 6)  {
        numDOF = 12;	    
        theMatrix = &ActuatorCorotM12;
        theVector = &ActuatorCorotV12;			
    }
    else  {
        opserr <<"ActuatorCorot::setDomain() - can not handle "
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
        opserr << "ActuatorCorot::setDomain() - element: " << this->getTag()
            << " out of memory creating vector of size: " << numDOF << endln;
        
        return;
    }          
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    // initalize the cosines
    double cosX[3];
    cosX[0] = cosX[1] = cosX[2] = 0.0;
    for (int i=0; i<numDIM; i++)
        cosX[i] = end2Crd(i)-end1Crd(i);

    // get initial length
    L = sqrt(cosX[0]*cosX[0] + cosX[1]*cosX[1] + cosX[2]*cosX[2]);
    if (L == 0.0)  {
        opserr <<"ActuatorCorot::setDomain() - element: "
            << this->getTag() << " has zero length\n";
        return;
    }
    Ln = L;
    
    // set global orientations
	cosX[0] /= L;
	cosX[1] /= L;
	cosX[2] /= L;

    // initialize rotation matrix
    R(0,0) = cosX[0];
    R(0,1) = cosX[1];
    R(0,2) = cosX[2];
    if (fabs(cosX[0]) > 0.0)  {  // element outside the YZ plane
        R(1,0) = -cosX[1];
        R(1,1) =  cosX[0];
        R(1,2) =  0.0;

        R(2,0) = -cosX[0]*cosX[2];
        R(2,1) = -cosX[1]*cosX[2];
        R(2,2) =  cosX[0]*cosX[0] + cosX[1]*cosX[1];
    }
    else  {  // element in the YZ plane
        R(1,0) =  0.0;
        R(1,1) = -cosX[2];
        R(1,2) =  cosX[1];

        R(2,0) =  1.0;
        R(2,1) =  0.0;
        R(2,2) =  0.0;
    }

    // orthonormalize last two rows of R
    double norm;
    for (int i=1; i<3; i++)  {
        norm = sqrt(R(i,0)*R(i,0) + R(i,1)*R(i,1) + R(i,2)*R(i,2));
        R(i,0) /= norm;
        R(i,1) /= norm;
        R(i,2) /= norm;
    }

    // set initial offsets
   	d21[0] = L;
	d21[1] = 0.0;
	d21[2] = 0.0;
}   	 


int ActuatorCorot::commitState()
{
    return 0;
}


int ActuatorCorot::revertToLastCommit()
{
    opserr << "ActuatorCorot::revertToLastCommit() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to last commit. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int ActuatorCorot::revertToStart()
{
    opserr << "ActuatorCorot::revertToStart() - "
        << "Element: " << this->getTag() << endln
        << "Can't revert to start. This element "
        << "is connected to an external process." 
        << endln;
    
    return -1;
}


int ActuatorCorot::update()
{
    if (theChannel == 0)  {
        if (this->setupConnection() != 0)  {
            opserr << "ActuatorCorot::update() - "
                << "failed to setup connection\n";
            return -1;
        }
    }
    
    // determine dsp in basic system
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    
    // initial offsets
    d21[0] = L;
    d21[1] = 0.0;
    d21[2] = 0.0;
    
    // update offsets in basic system due to nodal displacements
    for (int i=0; i<numDIM; i++) {
        d21[0] += (dsp2(i)-dsp1(i))*R(0,i);
        d21[1] += (dsp2(i)-dsp1(i))*R(1,i);
        d21[2] += (dsp2(i)-dsp1(i))*R(2,i);
    }
    
    // compute new length and deformation
    Ln = sqrt(d21[0]*d21[0] + d21[1]*d21[1] + d21[2]*d21[2]);
    db(0) = Ln - L;

    return 0;
}


const Matrix& ActuatorCorot::getTangentStiff(void)
{
    // zero the matrix
    theMatrix->Zero();
    
    // local stiffness matrix
    static Matrix kl(3,3);
    
    // material stiffness portion
    int i,j;
    double EAoverL3 = EA/(Ln*Ln*L);
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            kl(i,j) = EAoverL3*d21[i]*d21[j];
    
    // geometric stiffness portion
    q(0) = EA/L*(db(0) - (*ctrlDisp)(0));
    double SL = q(0)/Ln;
    double SA = q(0)/(Ln*Ln*Ln);
    for (i=0; i<3; i++)  {
        kl(i,i) += SL;
        for (j=0; j<3; j++)
            kl(i,j) -= SA*d21[i]*d21[j];
    }
    
    // compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);
    
    // copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (i=0; i<numDIM; i++)  {
        for (j=0; j<numDIM; j++)  {
            (*theMatrix)(i,j) = kg(i,j);
            (*theMatrix)(i,j+numDOF2) = -kg(i,j);
            (*theMatrix)(i+numDOF2,j) = -kg(i,j);
            (*theMatrix)(i+numDOF2,j+numDOF2) = kg(i,j);
        }
    }
    
    return *theMatrix;
}


const Matrix& ActuatorCorot::getInitialStiff(void)
{
    // zero the matrix
    theMatrix->Zero();
    
    // local stiffness matrix
    static Matrix kl(3,3);
    
    // material stiffness portion
    kl.Zero();
    kl(0,0) = EA/L;
    
    // compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);
    
    // copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (int i=0; i<numDIM; i++)  {
        for (int j=0; j<numDIM; j++)  {
            (*theMatrix)(i,j) = kg(i,j);
            (*theMatrix)(i,j+numDOF2) = -kg(i,j);
            (*theMatrix)(i+numDOF2,j) = -kg(i,j);
            (*theMatrix)(i+numDOF2,j+numDOF2) = kg(i,j);
        }
    }
    
    return *theMatrix;
}


const Matrix& ActuatorCorot::getMass()
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


void ActuatorCorot::zeroLoad()
{
    theLoad->Zero();
}


int ActuatorCorot::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
    opserr <<"ActuatorCorot::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int ActuatorCorot::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for a quick return
    if (L == 0.0 || rho == 0.0)
        return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size())  {
        opserr <<"ActuatorCorot::addInertiaLoadToUnbalance() - "
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


const Vector& ActuatorCorot::getResistingForce()
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
                opserr << "ActuatorCorot::getResistingForce() - "
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
    
    // local resisting force
    static Vector ql(3);
    
	ql(0) = d21[0]/Ln*q(0);
	ql(1) = d21[1]/Ln*q(0);
	ql(2) = d21[2]/Ln*q(0);
    
    static Vector qg(3);
    qg.addMatrixTransposeVector(0.0, R, ql, 1.0);
    
    // zero the residual
    theVector->Zero();
    
    // copy forces into appropriate places
    int numDOF2 = numDOF/2;
    for (int i=0; i<numDIM; i++)  {
        (*theVector)(i) = -qg(i);
        (*theVector)(i+numDOF2) = qg(i);
    }
    
    // subtract external load
    (*theVector) -= *theLoad;
    
    return *theVector;
}


const Vector& ActuatorCorot::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        (*theVector) += this->getRayleighDampingForces();
    
    // now include the mass portion
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


int ActuatorCorot::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(6);
    data(0) = this->getTag();
    data(1) = numDIM;
    data(2) = numDOF;
    data(3) = EA;
    data(4) = ipPort;
    data(5) = rho;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    return 0;
}


int ActuatorCorot::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // receive element parameters
    static Vector data(6);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    numDIM = (int)data(1);
    numDOF = (int)data(2);
    EA     = data(3);
    ipPort = (int)data(4);
    rho    = data(5);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    return 0;
}


int ActuatorCorot::displaySelf(Renderer &theViewer,
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
    
    for (int i=0; i<numDIM; i++)  {
        v1(i) = end1Crd(i) + end1Disp(i)*fact;
        v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}


void ActuatorCorot::Print(OPS_Stream &s, int flag)
{    
    if (flag == 0)  {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: ActuatorCorot, iNode: " << connectedExternalNodes(0)
            << ", jNode: " << connectedExternalNodes(1) << endln;
        s << "  EA: " << EA << ", L: " << L << ", Ln: " << Ln << endln;
        s << "  ipPort: " << ipPort << endln;
        s << "  mass per unit length: " << rho << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    } else if (flag == 1)  {
        // does nothing
    }
}


Response* ActuatorCorot::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ActuatorCorot");
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


int ActuatorCorot::getResponse(int responseID, Information &eleInformation)
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


int ActuatorCorot::setupConnection()
{
    // setup the connection
    theChannel = new TCP_Socket(ipPort);
    if (theChannel != 0)  {
        opserr << "\nChannel successfully created: "
            << "Waiting for ECSimAdapter experimental control...\n";
    } else {
        opserr << "ActuatorCorot::setupConnection() - "
            << "could not create channel\n";
        return -1;
    }
    if (theChannel->setUpConnection() != 0)  {
        opserr << "ActuatorCorot::setupConnection() - "
            << "failed to setup connection\n";
        return -2;
    }

    // get the data sizes and check values
    // sizes = {ctrlDisp, ctrlVel, ctrlAccel, ctrlForce, ctrlTime,
    //          daqDisp,  daqVel,  daqAccel,  daqForce,  daqTime,  dataSize}
    ID sizes(11);
    theChannel->recvID(0, 0, sizes, 0);
    if (sizes(0) > 1 || sizes(3) > 1 || sizes(5) > 1 || sizes(8) > 1)  {
        opserr << "ActuatorCorot::setupConnection() - "
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

    opserr << "\nActuatorCorot element " << this->getTag()
        << " now running...\n";

    return 0;
}
