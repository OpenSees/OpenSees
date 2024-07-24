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

#include <TwoNodeLinkSection.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <SectionForceDeformation.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <elementAPI.h>
#include <vector>

// initialize the class wide variables
Matrix TwoNodeLinkSection::TwoNodeLinkSectionM2(2,2);
Matrix TwoNodeLinkSection::TwoNodeLinkSectionM4(4,4);
Matrix TwoNodeLinkSection::TwoNodeLinkSectionM6(6,6);
Matrix TwoNodeLinkSection::TwoNodeLinkSectionM12(12,12);
Vector TwoNodeLinkSection::TwoNodeLinkSectionV2(2);
Vector TwoNodeLinkSection::TwoNodeLinkSectionV4(4);
Vector TwoNodeLinkSection::TwoNodeLinkSectionV6(6);
Vector TwoNodeLinkSection::TwoNodeLinkSectionV12(12);

void* OPS_TwoNodeLinkSection()
{
    int ndm = OPS_GetNDM();
    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: twoNodeLink eleTag iNode jNode secTag <-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-shearDist sDratios> <-doRayleigh> <-mass m>\n";
        return 0;
    }
    
    // tags
    int idata[4];
    int numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "WARNING: invalid integer data\n";
        return 0;
    }

    // get section
    SectionForceDeformation* theSection = OPS_getSectionForceDeformation(idata[3]);
    if(theSection == 0) {
      opserr << "twoNodeLinkSection -- no section with tag " << idata[3] << " exists in Domain" << endln;
      return 0;
    }    

    // options
    Vector x, y, Mratio, sDistI;
    int doRayleigh = 0;
    double mass = 0.0;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        return new TwoNodeLinkSection(idata[0], ndm, idata[1], idata[2], *theSection);
    }
    
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *type = OPS_GetString();
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
        else if (strcmp(type, "-shearDist") == 0) {
            sDistI.resize(2);
            numdata = 2;
            if (ndm == 2) {
                numdata = 1;
                sDistI(1) = 0.5;
            }
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING: insufficient data for -shearDist\n";
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, &sDistI(0)) < 0) {
                opserr << "WARNING: invalid -shearDist value\n";
                return 0;
            }
        }
        else if (strcmp(type, "-doRayleigh") == 0) {
            doRayleigh = 1;
        }
        else if (strcmp(type, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WANRING: insufficient mass value\n";
                return 0;
            }
            numdata = 1;
            if (OPS_GetDoubleInput(&numdata, &mass) < 0) {
                opserr << "WANRING: invalid -mass value\n";
                return 0;
            }
        }
    }
    
    // create object
    return new TwoNodeLinkSection(idata[0], ndm, idata[1], idata[2],
        *theSection, y, x, Mratio, sDistI, doRayleigh, mass);
}


// responsible for allocating the necessary space needed
// by each object and storing the tags of the end nodes.
TwoNodeLinkSection::TwoNodeLinkSection(int tag, int dim, int Nd1, int Nd2, 
				       SectionForceDeformation &theSec,
				       const Vector _y, const Vector _x, const Vector Mr,
				       const Vector sdI, int addRay, double m)
  : Element(tag, ELE_TAG_TwoNodeLinkSection),
    numDIM(dim), numDOF(0), connectedExternalNodes(2),
    theSection(0), trans(3,3),
    x(_x), y(_y), Mratio(Mr), shearDistI(sdI), addRayleigh(addRay),
    mass(m), L(0.0), onP0(true), ub(0), ubdot(0), qb(0), ul(0),
    Tgl(0,0), Tlb(0,0), theMatrix(0), theVector(0), theLoad(0)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - element: "
            << this->getTag() << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;

    theSection = theSec.getCopy();
    if (theSection == 0)  {
      opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
	     << "failed to copy of section force-deformation.\n";
      exit(-1);
    }

    int numDIR = theSection->getOrder();
    
    // check the number of directions
    if (numDIR < 1) {
        opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - element: "
            << this->getTag() << " wrong number of directions\n";
        exit(-1);
    }
    
    // check p-delta moment distribution ratios
    if (Mratio.Size() == 4)  {
      for (int i = 0; i < 4; i++) {
	if (Mratio(i) < 0.0) {
	  opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
		 << "P-Delta moment ratio " << i+1 << " is negative\n";
	  Mratio(i) = -Mratio(i);
	  opserr << "Making the value positive " << Mratio(i) << endln;
	}
      }
      double sumRatio = Mratio(0)+Mratio(1);
      if (sumRatio > 1.0)  {
	opserr << "TwoNodeLink::TwoNodeLink() - "
	       << "incorrect P-Delta moment ratios:\nrMy1 + rMy2 = "
	       << sumRatio << " > 1.0\n";
	Mratio(0) = Mratio(0)/sumRatio;
	Mratio(1) = Mratio(1)/sumRatio;	    
	opserr << "Scaling ratios down to " << Mratio(0) << " and " << Mratio(1) << endln;
      }
      sumRatio = Mratio(2) + Mratio(3);
      if (sumRatio > 1.0)  {
	opserr << "TwoNodeLink::TwoNodeLink() - "
	       << "incorrect P-Delta moment ratios:\nrMz1 + rMz2 = "
	       << sumRatio << " > 1.0\n";
	Mratio(2) = Mratio(2)/sumRatio;
	Mratio(3) = Mratio(3)/sumRatio;	    
	opserr << "Scaling ratios down to " << Mratio(2) << " and " << Mratio(3) << endln;	    
      }      
    }
    
    // check or initialize shear distance ratios
    if (shearDistI.Size() == 2)  {
      if (shearDistI(0) < 0.0) {
	opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
	       << "Shear distance ratio shearDistIy is negative\n";
	shearDistI(0) = -shearDistI(0);
	opserr << "Making the value positive " << shearDistI(0) << endln;
      }
      if (shearDistI(0) > 1.0) {
	opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
	       << "Shear distance ratio shearDistIy is greater than 1\n";
	shearDistI(0) = 1.0;
	opserr << "Making the value 1.0" << endln;
      }
      if (shearDistI(1) < 0.0) {
	opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
	       << "Shear distance ratio shearDistIz is negative\n";
	shearDistI(1) = -shearDistI(1);
	opserr << "Making the value positive " << shearDistI(1) << endln;
      }
      if (shearDistI(1) > 1.0) {
	opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
	       << "Shear distance ratio shearDistIz is greater than 1\n";
	shearDistI(1) = 1.0;
	opserr << "Making the value 1.0" << endln;
      }		

    } else  {
        shearDistI.resize(2);
        shearDistI(0) = 0.5;
        shearDistI(1) = 0.5;
    }
    
    // initialize response vectors in basic system
    ub.resize(numDIR);
    ubdot.resize(numDIR);
    qb.resize(numDIR);
    this->revertToStart();
}


TwoNodeLinkSection::TwoNodeLinkSection()
    : Element(0, ELE_TAG_TwoNodeLinkSection),
    numDIM(0), numDOF(0), connectedExternalNodes(2),
    theSection(0), trans(3,3), x(0), y(0),
    Mratio(0), shearDistI(0), addRayleigh(0), mass(0.0), L(0.0),
    onP0(false), ub(0), ubdot(0), qb(0), ul(0), Tgl(0,0), Tlb(0,0),
    theMatrix(0), theVector(0), theLoad(0)
{
    // ensure the connectedExternalNode ID is of correct size
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "TwoNodeLinkSection::TwoNodeLinkSection() - "
            << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}


// delete must be invoked on any objects created by the object.
TwoNodeLinkSection::~TwoNodeLinkSection()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theLoad != 0)
        delete theLoad;
    
    // delete the section
    if (theSection != 0)
      delete theSection;
}


int TwoNodeLinkSection::getNumExternalNodes() const
{
    return 2;
}


const ID& TwoNodeLinkSection::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** TwoNodeLinkSection::getNodePtrs() 
{
    return theNodes;
}


int TwoNodeLinkSection::getNumDOF() 
{
    return numDOF;
}


// to set a link to the enclosing Domain and to set the node pointers.
void TwoNodeLinkSection::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        
        return;
    }
    
    // set default values for error conditions
    numDOF = 2;
    theMatrix = &TwoNodeLinkSectionM2;
    theVector = &TwoNodeLinkSectionV2;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "TwoNodeLinkSection::setDomain() - Nd1: "
                << Nd1 << " does not exist in the model for ";
        } else  {
            opserr << "TwoNodeLinkSection::setDomain() - Nd2: " 
                << Nd2 << " does not exist in the model for ";
        }
        opserr << "TwoNodeLinkSection ele: " << this->getTag() << endln;
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2)  {
        opserr << "TwoNodeLinkSection::setDomain(): nodes " << Nd1 << " and " << Nd2
            << "have differing dof at ends for element: " << this->getTag() << endln;
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 1 && dofNd1 == 1)  {
        numDOF = 2;
        theMatrix = &TwoNodeLinkSectionM2;
        theVector = &TwoNodeLinkSectionV2;
        elemType  = D1N2;
    }
    else if (numDIM == 2 && dofNd1 == 2)  {
        numDOF = 4;
        theMatrix = &TwoNodeLinkSectionM4;
        theVector = &TwoNodeLinkSectionV4;
        elemType  = D2N4;
    }
    else if (numDIM == 2 && dofNd1 == 3)  {
        numDOF = 6;
        theMatrix = &TwoNodeLinkSectionM6;
        theVector = &TwoNodeLinkSectionV6;
        elemType  = D2N6;
    }
    else if (numDIM == 3 && dofNd1 == 3)  {
        numDOF = 6;
        theMatrix = &TwoNodeLinkSectionM6;
        theVector = &TwoNodeLinkSectionV6;
        elemType  = D3N6;
    }
    else if (numDIM == 3 && dofNd1 == 6)  {
        numDOF = 12;
        theMatrix = &TwoNodeLinkSectionM12;
        theVector = &TwoNodeLinkSectionV12;
        elemType  = D3N12;
    }
    else  {
        opserr << "TwoNodeLinkSection::setDomain() can not handle "
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
        opserr << "TwoNodeLinkSection::setDomain() - element: " << this->getTag()
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


int TwoNodeLinkSection::commitState()
{
    int errCode = 0;

    errCode += theSection->commitState();
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int TwoNodeLinkSection::revertToLastCommit()
{
    int errCode = 0;

    errCode += theSection->revertToLastCommit();
    
    return errCode;
}


int TwoNodeLinkSection::revertToStart()
{   
    int errCode = 0;
    
    // reset trial history variables
    ub.Zero();
    ubdot.Zero();
    qb.Zero();
    
    errCode += theSection->revertToStart();
    
    return errCode;
}


int TwoNodeLinkSection::update()
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
    //ub = (Tlb*Tgl)*ug;
    //ubdot = (Tlb*Tgl)*ugdot;

    // MHS hack
    ub(0) /= L;
    ub(1) /= L;
    
    errCode += theSection->setTrialSectionDeformation(ub);
    
    return errCode;
}


const Matrix& TwoNodeLinkSection::getTangentStiff()
{
    // zero the matrix
    theMatrix->Zero();

    int numDIR = theSection->getOrder();
    
    // get resisting force and stiffness
    //const Matrix &kb = theSection->getSectionTangent();
    qb = theSection->getStressResultant();

    // MHS hack
    Matrix kb = theSection->getSectionTangent();
    kb(0,0) /= L;
    kb(1,1) /= L;
    //kb(2,2) *= L;

    //qb(2) *= L;
    
    // transform stiffness from basic to local system
    Matrix kl(numDOF,numDOF);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add P-Delta effects to local stiffness
    if (Mratio.Size() == 4)
        this->addPDeltaStiff(kl, qb);
    
    // transform stiffness from local to global system
    theMatrix->addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    //Matrix kg(numDOF,numDOF);
    //kg.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    //theMatrix->addMatrixTranspose(0.5, kg, 0.5);
    
    return *theMatrix;
}


const Matrix& TwoNodeLinkSection::getInitialStiff()
{
    // zero the matrix
    theMatrix->Zero();

    int numDIR = theSection->getOrder();
    
    // get initial stiffness
    const Matrix &kbInit = theSection->getInitialTangent();
    
    // transform stiffness from basic to local system
    Matrix klInit(numDOF,numDOF);
    klInit.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform stiffness from local to global system
    theMatrix->addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);
    //Matrix kgInit(numDOF,numDOF);
    //kgInit.addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);
    //theMatrix->addMatrixTranspose(0.5, kgInit, 0.5);
    
    return *theMatrix;
}


const Matrix& TwoNodeLinkSection::getDamp()
{
    // zero the matrix
    theMatrix->Zero();

    // call base class to set up Rayleigh damping
    double factThis = 0.0;
    if (addRayleigh == 1)  {
        (*theMatrix) = this->Element::getDamp();
        factThis = 1.0;
    }

    int numDIR = theSection->getOrder();
        
    // add damping tangent from materials
    Matrix cb(numDIR,numDIR);
    //for (int i=0; i<numDIR; i++)  {
    //    cb(i,i) = theMaterials[i]->getDampTangent();
    //}
    
    // transform damping from basic to local system
    Matrix cl(numDOF,numDOF);
    cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);
    
    // transform damping from local to global system
    theMatrix->addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    //Matrix cg(numDOF,numDOF);
    //cg.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    //theMatrix->addMatrixTranspose(0.5, cg, 0.5);
    
    return *theMatrix;
}


const Matrix& TwoNodeLinkSection::getMass()
{
    // zero the matrix
    theMatrix->Zero();
    
    // form mass matrix
    if (mass != 0.0)  {
        double m = 0.5*mass;
        int numDOF2 = numDOF/2;
        for (int i=0; i<numDIM; i++)  {
            (*theMatrix)(i,i)                 = m;
            (*theMatrix)(i+numDOF2,i+numDOF2) = m;
        }
    }
    
    return *theMatrix; 
}


void TwoNodeLinkSection::zeroLoad()
{
    theLoad->Zero();
}


int TwoNodeLinkSection::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"TwoNodeLinkSection::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int TwoNodeLinkSection::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    int numDOF2 = numDOF/2;
    if (numDOF2 != Raccel1.Size() || numDOF2 != Raccel2.Size())  {
        opserr << "TwoNodeLinkSection::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5*mass;
    for (int i=0; i<numDIM; i++)  {
        (*theLoad)(i)         -= m * Raccel1(i);
        (*theLoad)(i+numDOF2) -= m * Raccel2(i);
    }
    
    return 0;
}


const Vector& TwoNodeLinkSection::getResistingForce()
{
    // zero the residual
    theVector->Zero();
    
    // get resisting force
    qb = theSection->getStressResultant();

    // MHS hack
    //qb(2) *= L;
    
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


const Vector& TwoNodeLinkSection::getResistingForceIncInertia()
{
    // this already includes damping force from materials
    this->getResistingForce();
    
    // subtract external load
    theVector->addVector(1.0, *theLoad, -1.0);
    
    // add the damping force from Rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add inertia force from element mass
    if (mass != 0.0)  {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();
        
        int numDOF2 = numDOF/2;
        double m = 0.5*mass;
        for (int i=0; i<numDIM; i++)  {
            (*theVector)(i)         += m * accel1(i);
            (*theVector)(i+numDOF2) += m * accel2(i);
        }
    }
    
    return *theVector;
}


int TwoNodeLinkSection::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(14);
    data(0) = this->getTag();
    data(1) = numDIM;
    data(2) = numDOF;
    //data(3) = numDIR;
    data(4) = x.Size();
    data(5) = y.Size();
    data(6) = Mratio.Size();
    data(7) = shearDistI.Size();
    data(8) = addRayleigh;
    data(9) = mass;
    data(10) = alphaM;
    data(11) = betaK;
    data(12) = betaK0;
    data(13) = betaKc;
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send the material class tags
    ID matClassTags(1);
    matClassTags(0) = theSection->getClassTag();
    sChannel.sendID(0, commitTag, matClassTags);
    
    theSection->sendSelf(commitTag, sChannel);
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    if (Mratio.Size() == 4)
        sChannel.sendVector(0, commitTag, Mratio);
    if (shearDistI.Size() == 2)
        sChannel.sendVector(0, commitTag, shearDistI);
    
    return 0;
}


int TwoNodeLinkSection::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete dynamic memory
    if (theSection != 0) 
      delete theSection;
    
    // receive element parameters
    static Vector data(14);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    numDIM = (int)data(1);
    numDOF = (int)data(2);
    //numDIR = (int)data(3);
    addRayleigh = (int)data(8);
    mass = data(9);
    alphaM = data(10);
    betaK = data(11);
    betaK0 = data(12);
    betaKc = data(13);
   
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive the material class tags
    ID matClassTags(1);
    rChannel.recvID(0, commitTag, matClassTags);
    
    theSection = theBroker.getNewSection(matClassTags(0));
    if (theSection == 0) {
      opserr << "TwoNodeLinkSection::recvSelf() - "
	     << "failed to get blank section.\n";
      return -3;
    }
    theSection->recvSelf(commitTag, rChannel, theBroker);

    
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
            opserr << "TwoNodeLinkSection::recvSelf() - "
                << "p-delta moment ratios can not be negative\n";
            return -4;
        }
        if (Mratio(0)+Mratio(1) > 1.0)  {
            opserr << "TwoNodeLinkSection::recvSelf() - "
                << "incorrect p-delta moment ratios:\nrMy1 + rMy2 = "
                << Mratio(0)+Mratio(1) << " > 1.0\n";
            return -4;
        }
        if (Mratio(2)+Mratio(3) > 1.0)  {
            opserr << "TwoNodeLinkSection::recvSelf() - "
                << "incorrect p-delta moment ratios:\nrMz1 + rMz2 = "
                << Mratio(2)+Mratio(3) << " > 1.0\n";
            return -4;
        }
    }
    if ((int)data(7) == 2)  {
        shearDistI.resize(2);
        rChannel.recvVector(0, commitTag, shearDistI);
        // check shear distance ratios
        if (shearDistI(0) < 0.0 || shearDistI(0) > 1.0)  {
            opserr << "TwoNodeLinkSection::recvSelf() - "
                << "incorrect shear distance ratio:\n shearDistIy = "
                << shearDistI(0) << " < 0.0 or > 1.0\n";
            return -5;
        }
        if (shearDistI(1) < 0.0 || shearDistI(1) > 1.0)  {
            opserr << "TwoNodeLinkSection::recvSelf() - "
                << "incorrect shear distance ratio:\n shearDistIz = "
                << shearDistI(1) << " < 0.0 or > 1.0\n";
            return -5;
        }
    } else  {
        // initialize shear distance ratios
        shearDistI.resize(2);
        shearDistI(0) = 0.5;
        shearDistI(1) = 0.5;
    }
    onP0 = false;
    
    // initialize response vectors in basic system
    int numDIR = theSection->getOrder();
    ub.resize(numDIR);
    ubdot.resize(numDIR);
    qb.resize(numDIR);
    this->revertToStart();
    
    return 0;
}


int TwoNodeLinkSection::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void TwoNodeLinkSection::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: TwoNodeLinkSection" << endln;
        s << "  iNode: " << connectedExternalNodes(0)
	  << ", jNode: " << connectedExternalNodes(1) << endln;
	s << theSection->getTag() << endln;
        s << "  Mratio: " << Mratio << "  shearDistI: " << shearDistI << endln;
        s << "  addRayleigh: " << addRayleigh << "  mass: " << mass << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
}


Response* TwoNodeLinkSection::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","TwoNodeLinkSection");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    char outputData[80];
    
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
      int numDIR = theSection->getOrder();
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
      int numDIR = theSection->getOrder();
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
      int numDIR = theSection->getOrder();
      for (int i=0; i<numDIR; i++)  {
	sprintf(outputData,"db%d",i+1);
	output.tag("ResponseType",outputData);
      }
      for (int i=0; i<numDIR; i++)  {
	sprintf(outputData,"q%d",i+1);
	output.tag("ResponseType",outputData);
        }
      theResponse = new ElementResponse(this, 6, Vector(numDIR*2));
    }
    // material output
    else if (strcmp(argv[0],"section") == 0)  {
        if (argc > 1)  {
	  theResponse =  theSection->setResponse(&argv[1], argc-1, output);
        }
    }

    if (strcmp(argv[0],"xaxis") == 0) {
      theResponse = new ElementResponse(this, 20, Vector(3));
    }
    if (strcmp(argv[0],"yaxis") == 0) {
      theResponse = new ElementResponse(this, 21, Vector(3));
    }
    if (strcmp(argv[0],"zaxis") == 0) {
      theResponse = new ElementResponse(this, 22, Vector(3));
    }

    output.endTag(); // ElementOutput
    
    return theResponse;
}


int TwoNodeLinkSection::getResponse(int responseID, Information &eleInfo)
{
  int numDIR = theSection->getOrder();
  
    Vector defoAndForce(numDIR*2);
    Vector &theVec = *(eleInfo.theVector);
    ID &theID = *(eleInfo.theID);    
    
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

    case 20:
      theVec(0) = trans(0,0);
      theVec(1) = trans(0,1);
      theVec(2) = trans(0,2);
      return 0;
    case 21:
      theVec(0) = trans(1,0);
      theVec(1) = trans(1,1);
      theVec(2) = trans(1,2);
      return 0;
    case 22:
      theVec(0) = trans(2,0);
      theVec(1) = trans(2,1);
      theVec(2) = trans(2,2);
      return 0;

    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void TwoNodeLinkSection::setUp()
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
            opserr << "WARNING TwoNodeLinkSection::setUp() - " 
                << "element: " << this->getTag() << endln
                << "ignoring nodes and using specified "
                << "local x vector to determine orientation\n";
        }
        if (y.Size() == 0)  {
            y.resize(3);
            y.Zero();
	    if (xp.Size() == 1) // Make a 1D model work
	      y(1) = 1.0;
	    else
	      y(0) = -xp(1);
            if (xp.Size() > 1)
                y(1) = xp(0);
            if (xp.Size() > 2)
                opserr << "WARNING TwoNodeLinkSection::setUp() - " 
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
        opserr << "TwoNodeLinkSection::setUp() - "
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
        opserr << "TwoNodeLinkSection::setUp() - "
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
void TwoNodeLinkSection::setTranGlobalLocal()
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
void TwoNodeLinkSection::setTranLocalBasic()
{
  int numDIR = theSection->getOrder();
  const ID &type = theSection->getType();
  
  // resize transformation matrix and zero it
  Tlb.resize(numDIR,numDOF);
  Tlb.Zero();

  for (int i=0; i<numDIR; i++)  {
    int dirID = this->getDirID(type(i));

    Tlb(i,dirID) = -1.0;
    Tlb(i,dirID+numDOF/2) = 1.0;
    
    // switch on dimensionality of element
    switch (elemType)  {
    case D2N6:
      if (dirID == 1)  {
	Tlb(i,2) = -shearDistI(0)*L;
	Tlb(i,5) = -(1.0 - shearDistI(0))*L;
      }
      break;
    case D3N12:
      if (dirID == 1)  {
	Tlb(i,5)  = -shearDistI(0)*L;
	Tlb(i,11) = -(1.0-shearDistI(0))*L;
      }
      else if (dirID == 2)  {
	Tlb(i,4)  = shearDistI(1)*L;
	Tlb(i,10) = (1.0-shearDistI(1))*L;
      }
      break;
    default :
      // do nothing
      break;
    }
  }
}

int TwoNodeLinkSection::getDirID(int sectionCode)
{
  int dirID = -1;
  
  switch (sectionCode) {
  case SECTION_RESPONSE_P:
    dirID = 0;
    break;
  case SECTION_RESPONSE_VY:
    dirID = 1;
    break;
  case SECTION_RESPONSE_VZ:
    dirID = 2;
    break;
  case SECTION_RESPONSE_T:
    dirID = 3;
    break;
  case SECTION_RESPONSE_MY:
    dirID = 4;
    break;
  case SECTION_RESPONSE_MZ:
    dirID = numDIM == 2 ? 2 : 5;
    break;
  default:
    break;
  }
  return dirID;
}

void TwoNodeLinkSection::addPDeltaForces(Vector &pLocal, const Vector& qBasic)
{
    int dirID;
    double N = 0.0;
    double deltal1 = 0.0;
    double deltal2 = 0.0;

    int numDIR = theSection->getOrder();
    const ID &type = theSection->getType();
    
    for (int i=0; i<numDIR; i++)  {
      int dirID = this->getDirID(type(i));
      
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
	  int dirID = this->getDirID(type(i));
	  
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


void TwoNodeLinkSection::addPDeltaStiff(Matrix &kLocal, const Vector& qBasic)
{
    int dirID;
    double N = 0.0;

    int numDIR = theSection->getOrder();
    const ID &type = theSection->getType();
    
    // get axial force
    for (int i=0; i<numDIR; i++)  {
      int dirID = this->getDirID(type(i));      
      if (dirID == 0)
	N = qBasic(i);
    }
    
    if (N != 0.0)  {
        for (int i=0; i<numDIR; i++)  {
	  int dirID = this->getDirID(type(i));      
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


int TwoNodeLinkSection::setParameter(const char **argv, int argc, Parameter &param)
{
    int result = -1;
    
    if (argc < 1)
        return -1;
    
    if (strcmp(argv[0], "section") == 0) {
      if (argc > 1) {
	return theSection->setParameter(&argv[1], argc - 1, param);
      }
      else {
	return -1;
      }
    }
    
    int res = theSection->setParameter(argv, argc, param);
    if (res != -1) {
      result = res;
    }
    
    return result;
}
