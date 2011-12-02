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

// Description: This file contains the implementation of the GenericCopy class.



#include "GenericCopy.h"



#include <Domain.h>

#include <Node.h>

#include <Channel.h>

#include <FEM_ObjectBroker.h>

#include <Renderer.h>

#include <Information.h>

#include <ElementResponse.h>



#include <math.h>

#include <stdlib.h>

#include <string.h>





// initialize the class wide variables

Matrix GenericCopy::theMatrix(1,1);

Matrix GenericCopy::theInitStiff(1,1);

Matrix GenericCopy::theMass(1,1);

Vector GenericCopy::theVector(1);

Vector GenericCopy::theLoad(1);



// responsible for allocating the necessary space needed

// by each object and storing the tags of the end nodes.

GenericCopy::GenericCopy(int tag, ID nodes, int srctag)

    : Element(tag, ELE_TAG_GenericCopy),

    connectedExternalNodes(nodes),

    numExternalNodes(0), numDOF(0),

    srcTag(srctag), theSource(0),

    initStiffFlag(false), massFlag(false)

{   

    // initialize nodes

    numExternalNodes = connectedExternalNodes.Size();

    theNodes = new Node* [numExternalNodes];

    if (!theNodes)  {

        opserr << "GenericCopy::GenericCopy() "

            << "- failed to create node array\n";

        exit(-1);

    }



    // set node pointers to NULL

    int i;

    for (i=0; i<numExternalNodes; i++)

        theNodes[i] = 0;

}





// delete must be invoked on any objects created by the object.

GenericCopy::~GenericCopy()

{

    // invoke the destructor on any objects created by the object

    // that the object still holds a pointer to

    if (theNodes != 0)

        delete [] theNodes;

}





int GenericCopy::getNumExternalNodes() const

{

    return numExternalNodes;

}





const ID& GenericCopy::getExternalNodes() 

{

    return connectedExternalNodes;

}





Node** GenericCopy::getNodePtrs() 

{

    return theNodes;

}





int GenericCopy::getNumDOF() 

{

    return numDOF;

}





// to set a link to the enclosing Domain and to set the node pointers.

void GenericCopy::setDomain(Domain *theDomain)

{

    // check Domain is not null - invoked when object removed from a domain

    int i;

    if (!theDomain)  {

        for (i=0; i<numExternalNodes; i++)

            theNodes[i] = 0;

        return;

    }



    // get a pointer to the source element

    theSource = theDomain->getElement(srcTag);

    if (theSource == 0) {

        opserr << "GenericCopy::setDomain() "

            << "- failed to get a pointer to the source "

            << "element with tag " << srcTag << endln;

        return;

    }

    

    // check we got correct number of nodes

    if (numExternalNodes != theSource->getNumExternalNodes()) {

        opserr << "GenericCopy::setDomain() "

            << "- number of external nodes of copy do not "

            << "agree with source\n";

        return;

    }



    // now set the node pointers

    for (i=0; i<numExternalNodes; i++)

        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));

    

    // if can't find all - send a warning message

    for (i=0; i<numExternalNodes; i++)  {

        if (!theNodes[i])  {

            opserr << "GenericCopy::setDomain() - Nd" << i << ": " 

                << connectedExternalNodes(i) << " does not exist in the "

                << "model for GenericCopy ele: " << this->getTag() << endln;

            return;

        }

    }

    

    // now determine the number of dof

    numDOF = 0;

    for (i=0; i<numExternalNodes; i++)  {

        numDOF += theNodes[i]->getNumberDOF();

    }

    if (numDOF != theSource->getNumDOF()) {

        opserr << "GenericCopy::setDomain() "

            << "- number of DOFs of copy do not "

            << "agree with source\n";

        return;

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





int GenericCopy::commitState()

{   

    // does nothing

    return 0;

}





int GenericCopy::revertToLastCommit()

{ 

    // does nothing

    return 0;

}





int GenericCopy::revertToStart()

{

    // does nothing

    return 0;

}





int GenericCopy::update()

{

    // does nothing

    return 0;

}





const Matrix& GenericCopy::getTangentStiff()

{

    // zero the matrix

    theMatrix.Zero();

    

    // get tangent stiffness matrix from source element

    theMatrix = theSource->getTangentStiff();

    

    return theMatrix;

}





const Matrix& GenericCopy::getInitialStiff()

{

    if (initStiffFlag == false)  {

        // zero the matrix

        theInitStiff.Zero();

    

        // get initial stiffness matrix from source element

        theInitStiff = theSource->getInitialStiff();

        initStiffFlag = true;

    }

    

    return theInitStiff;

}





/*const Matrix& GenericCopy::getDamp()

{

    // zero the matrix

    theMatrix.Zero();

    

    // get damping matrix from source element

    theMatrix = theSource->getDamp();

    

    return theMatrix;

}*/





const Matrix& GenericCopy::getMass()

{

    if (massFlag == false)  {

        // zero the matrix

        theMass.Zero();

        

        // get mass matrix from source element

        theMatrix = theSource->getMass();

        massFlag = true;

    }

    

    return theMass;

}





void GenericCopy::zeroLoad()

{

    theLoad.Zero();

}





int GenericCopy::addLoad(ElementalLoad *theLoad, double loadFactor)

{  

    opserr <<"GenericCopy::addLoad() - "

        << "load type unknown for element: "

        << this->getTag() << endln;

    

    return -1;

}





int GenericCopy::addInertiaLoadToUnbalance(const Vector &accel)

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





const Vector& GenericCopy::getResistingForce()

{

    // zero the residual

    theVector.Zero();

    

    // determine resisting forces in global system

    theVector = theSource->getResistingForce();

    

    // subtract external load

    theVector.addVector(1.0, theLoad, -1.0);

    

    return theVector;

}





const Vector& GenericCopy::getResistingForceIncInertia()

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





int GenericCopy::sendSelf(int commitTag, Channel &theChannel)

{

    // has not been implemented yet.....

    return 0;

}





int GenericCopy::recvSelf(int commitTag, Channel &theChannel,

    FEM_ObjectBroker &theBroker)

{

    // has not been implemented yet.....

    return 0;

}





int GenericCopy::displaySelf(Renderer &theViewer,

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





void GenericCopy::Print(OPS_Stream &s, int flag)

{

    int i;

    if (flag == 0)  {

        // print everything

        s << "Element: " << this->getTag() << endln;

        s << "  type: GenericCopy";

        for (i=0; i<numExternalNodes; i++ )

            s << "  Node" << i+1 << ": " << connectedExternalNodes(i);

        s << "\n  source element: " << srcTag << endln;

        // determine resisting forces in global system

        s << "  resisting force: " << this->getResistingForce() << endln;

    } else if (flag == 1)  {

        // does nothing

    }

}





Response* GenericCopy::setResponse(const char **argv, int argc,

    OPS_Stream &output)

{

    Response *theResponse = 0;



    int i;

    char outputData[10];



    output.tag("ElementOutput");

    output.attr("eleType","GenericCopy");

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



    output.endTag(); // ElementOutput



    return theResponse;

}





int GenericCopy::getResponse(int responseID, Information &eleInformation)

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

        

    default:

        return -1;

    }

}

