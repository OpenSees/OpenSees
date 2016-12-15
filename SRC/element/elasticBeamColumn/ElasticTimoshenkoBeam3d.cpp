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
// Created: 03/13
// Revision: A
//
// Purpose: This file contains the class definition for ElasticTimoshenkoBeam3d.
// ElasticTimoshenkoBeam3d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <ElasticTimoshenkoBeam3d.h>

#include <Domain.h>
#include <Node.h>
#include <CrdTransf.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>

#include <math.h>
#include <stdlib.h>


// initialize the class wide variables
Matrix ElasticTimoshenkoBeam3d::theMatrix(12,12);
Vector ElasticTimoshenkoBeam3d::theVector(12);


void *OPS_ElasticTimoshenkoBeam3d()
{
    Element *theElement = 0;
    
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    
    if (numRemainingArgs == 0)  { // parallel processing
        theElement = new ElasticTimoshenkoBeam3d();
        return theElement;
    }
    
    if (numRemainingArgs < 11)  {
        opserr << "ERROR not enough args provided, want: element ElasticTimoshenkoBeam3d $tag $iNode $jNode $E $G $A $Jx $Iy $Iz $Avy $Avz $transTag <-mass $m> <-cMass> \n";
        return 0;
    }
    
    int numData;
    int iData[5];     // tag, iNode, jNode, transTag, cMass
    double dData[9];  // E, G, A, Jx, Iy, Iz, Avy, Avz, mass
    
    iData[4] = 0;     // cMass
    dData[8] = 0.0;   // mass per unit length
    
    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0)  {
        opserr << "WARNING invalid element data (tag, iNode, jNode) element ElasticTimoshenkoBeam3d.\n";
        return 0;
    }
    
    numData = 8;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element data (E, G, A, Jx, Iy, Iz, Avy, Avz) element ElasticTimoshenkoBeam3d " << iData[0] << endln;
        return 0;
    }
    
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData[3]) != 0)  {
        opserr << "WARNING invalid element data (transTag) element ElasticTimoshenkoBeam3d " << iData[0] << endln;
        return 0;
    }

    CrdTransf *theTrans = OPS_GetCrdTransf(iData[3]);
    if (theTrans == 0)  {
	    opserr << "WARNING transformation object not found for ElasticTimoshenkoBeam3d " << iData[0] << endln;
	    return 0;
    }
    
    numRemainingArgs = OPS_GetNumRemainingInputArgs();
    while (numRemainingArgs > 0)  {
      const char *argvLoc = OPS_GetString();
        numData = 1;
        
        if ((strcmp(argvLoc, "-mass") == 0) || (strcmp(argvLoc, "mass") == 0))  {
            if (OPS_GetDoubleInput(&numData, &dData[8]) != 0)  {
                opserr << "WARNING error reading element data (mass) element ElasticTimoshenkoBeam3d " << iData[0] << endln;
                return 0;
            }
        }
        if ((strcmp(argvLoc, "-lMass") == 0) || (strcmp(argvLoc, "lMass") == 0))  {
            iData[4] = 0;  // lumped mass matrix (default)
        }
        if ((strcmp(argvLoc, "-cMass") == 0) || (strcmp(argvLoc, "cMass") == 0))  {
            iData[4] = 1;  // consistent mass matrix
        }
        numRemainingArgs = OPS_GetNumRemainingInputArgs();      
    }
    
    theElement = new ElasticTimoshenkoBeam3d(iData[0], iData[1], iData[2],
        dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
        dData[7], *theTrans, dData[8], iData[4]);
    
    return theElement;
}


ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d(int tag, int Nd1, int Nd2, 
    double e, double g, double a, double jx, double iy, double iz, double avy,
    double avz, CrdTransf &coordTransf, double r, int cm)
    : Element(tag, ELE_TAG_ElasticTimoshenkoBeam3d),
    connectedExternalNodes(2), theCoordTransf(0), E(e), G(g), A(a), Jx(jx),
    Iy(iy), Iz(iz), Avy(avy), Avz(avz), rho(r), cMass(cm), nlGeo(0), phiY(0.0),
    phiZ(0.0), L(0.0), ul(12), ql(12), ql0(12), kl(12,12), klgeo(12,12),
    Tgl(12,12), Ki(12,12), M(12,12), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // get a copy of the coordinate transformation
    theCoordTransf = coordTransf.getCopy3d();
    if (!theCoordTransf)  {
        opserr << "ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d() - "
            << "failed to get copy of coordinate transformation.\n";
        exit(-1);
    }
    
    // get coordinate transformation type and save flag
    if (strncmp(theCoordTransf->getClassType(),"Linear",6) == 0)  {
        nlGeo = 0;
    } else if (strncmp(theCoordTransf->getClassType(),"PDelta",6) == 0)  {
        nlGeo = 1;
    } else if (strncmp(theCoordTransf->getClassType(),"Corot",5) == 0)  {
        nlGeo = 1;
        opserr << "\nWARNING ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d()"
            << " - Element: " << this->getTag() << endln
            << "Unsupported Corotational transformation assigned.\n"
            << "Using PDelta transformation instead.\n";
    }
    
    // zero fixed end forces vector
    ql0.Zero();
}


ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d()
    : Element(0, ELE_TAG_ElasticTimoshenkoBeam3d),
    connectedExternalNodes(2), theCoordTransf(0), E(0.0), G(0.0), A(0.0),
    Jx(0.0), Iy(0.0), Iz(0.0), Avy(0.0), Avz(0.0), rho(0.0), cMass(0),
    nlGeo(0), phiY(0.0), phiZ(0.0), L(0.0), ul(12), ql(12), ql0(12),
    kl(12,12), klgeo(12,12), Tgl(12,12), Ki(12,12), M(12,12), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElasticTimoshenkoBeam3d::ElasticTimoshenkoBeam3d() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // zero fixed end forces vector
    ql0.Zero();
}


ElasticTimoshenkoBeam3d::~ElasticTimoshenkoBeam3d()
{
    if (theCoordTransf)
	    delete theCoordTransf;
}


int ElasticTimoshenkoBeam3d::getNumExternalNodes() const
{
    return 2;
}


const ID& ElasticTimoshenkoBeam3d::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** ElasticTimoshenkoBeam3d::getNodePtrs() 
{
    return theNodes;
}


int ElasticTimoshenkoBeam3d::getNumDOF()
{
    return 12;
}


void ElasticTimoshenkoBeam3d::setDomain(Domain *theDomain)
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
            opserr << "WARNING ElasticTimoshenkoBeam3d::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING ElasticTimoshenkoBeam3d::setDomain() - Nd2: "
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
        opserr << "ElasticTimoshenkoBeam3d::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "ElasticTimoshenkoBeam3d::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // initialize the coordinate transformation
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0)  {
	    opserr << "ElasticTimoshenkoBeam3d::setDomain() - "
            << "error initializing coordinate transformation.\n";
	    return;
    }
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int ElasticTimoshenkoBeam3d::commitState()
{
    int errCode = 0;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    // no need to commit coordinate transformation
    // since it is only used to get type of transf
    // errCode += theCoordTransf->commitState();
    
    return errCode;
}


int ElasticTimoshenkoBeam3d::revertToLastCommit()
{
    return 0;
}


int ElasticTimoshenkoBeam3d::revertToStart()
{
    return 0;
}


int ElasticTimoshenkoBeam3d::update()
{
    return 0;
}


const Matrix& ElasticTimoshenkoBeam3d::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    if (nlGeo == 0)  {
        // transform from local to global system
        theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
        
    } else  {
        // initialize local stiffness matrix
        static Matrix klTot(12,12);
        klTot.addMatrix(0.0, kl, 1.0);
        
        // get global trial displacements
        const Vector &dsp1 = theNodes[0]->getTrialDisp();
        const Vector &dsp2 = theNodes[1]->getTrialDisp();
        static Vector ug(12);
        for (int i=0; i<6; i++)  {
            ug(i)   = dsp1(i);
            ug(i+6) = dsp2(i);
        }
        
        // transform response from the global to the local system
        ul.addMatrixVector(0.0, Tgl, ug, 1.0);
        
        // get the resisting forces in local system
        ql.addMatrixVector(0.0, kl, ul, 1.0);
        
        // add geometric stiffness to local stiffness
        if (ql(0) != 0.0)
            klTot.addMatrix(1.0, klgeo, ql(0));
        
        // transform from local to global system
        theMatrix.addMatrixTripleProduct(0.0, Tgl, klTot, 1.0);
    }
    
    return theMatrix;
}


const Matrix& ElasticTimoshenkoBeam3d::getInitialStiff()
{
    return Ki;
}


const Matrix& ElasticTimoshenkoBeam3d::getMass()
{ 
    return M;
}


void ElasticTimoshenkoBeam3d::zeroLoad()
{
    theLoad.Zero();
    ql0.Zero();
    
    return;
}


int ElasticTimoshenkoBeam3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    
    if (type == LOAD_TAG_Beam3dUniformLoad) {
        double wy = data(0)*loadFactor;  // Transverse
        double wz = data(1)*loadFactor;  // Transverse
        double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)
        
        double Vy = 0.5*wy*L;
        double Mz = Vy*L/6.0; // wy*L*L/12
        double Vz = 0.5*wz*L;
        double My = Vz*L/6.0; // wz*L*L/12
        double P  = 0.5*wx*L;
        
        // fixed end forces in local system
        ql0(0)  -= P;
        ql0(1)  -= Vy;
        ql0(2)  -= Vz;
        ql0(4)  += My;
        ql0(5)  -= Mz;
        ql0(6)  -= P;
        ql0(7)  -= Vy;
        ql0(8)  -= Vz;
        ql0(10) -= My;
        ql0(11) += Mz;
    }
    
    else {
        opserr << "ElasticTimoshenkoBeam3d::addLoad() - "
            << "load type unknown for element: " 
            << this->getTag() << ".\n";
        return -1;
    }
    
    return 0;
}


int ElasticTimoshenkoBeam3d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (rho == 0.0)
        return 0;
    
    // assemble Raccel vector
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    static Vector Raccel(12);
    for (int i=0; i<6; i++)  {
        Raccel(i)   = Raccel1(i);
        Raccel(i+6) = Raccel2(i);
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    theLoad.addMatrixVector(1.0, M, Raccel, -1.0);
    
    return 0;
}


const Vector& ElasticTimoshenkoBeam3d::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // get global trial displacements
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    static Vector ug(12);
    for (int i=0; i<6; i++)  {
        ug(i)   = dsp1(i);
        ug(i+6) = dsp2(i);
    }
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    
    // get the resisting forces in local system
    ql.addMatrixVector(0.0, kl, ul, 1.0);
    
    // add P-Delta moments to local forces
    if ((ql(0) != 0.0) && (nlGeo == 1))
        ql.addMatrixVector(1.0, klgeo, ul, ql(0));
    
    // add effects of element loads, ql = ql(ul) + ql0
    ql.addVector(1.0, ql0, 1.0);
    
    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    return theVector;
}


const Vector& ElasticTimoshenkoBeam3d::getResistingForceIncInertia()
{	
    // first get the resisting forces
    theVector = this->getResistingForce();
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    // add the damping forces from rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
    // check for quick return
    if (rho == 0.0)
        return theVector;
    
    // add inertia forces from element mass
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();    
    static Vector accel(12);
    for (int i=0; i<6; i++)  {
        accel(i)   = accel1(i);
        accel(i+6) = accel2(i);
    }
    theVector.addMatrixVector(1.0, M, accel, 1.0);
    
    return theVector;
}


int ElasticTimoshenkoBeam3d::sendSelf(int commitTag, Channel &sChannel)
{
    int res = 0;
    
    static Vector data(19);
    data(0) = this->getTag();
    data(1) = connectedExternalNodes(0);
    data(2) = connectedExternalNodes(1);
    data(3) = E;
    data(4) = G;
    data(5) = A;
    data(6) = Jx;
    data(7) = Iy;
    data(8) = Iz;
    data(9) = Avy;
    data(10) = Avz;
    data(11) = rho;
    data(12) = cMass;
    data(13) = alphaM;
    data(14) = betaK;
    data(15) = betaK0;
    data(16) = betaKc;
    data(17) = theCoordTransf->getClassTag();
    
    int dbTag = theCoordTransf->getDbTag();
    if (dbTag == 0)  {
        dbTag = sChannel.getDbTag();
        if (dbTag != 0)
            theCoordTransf->setDbTag(dbTag);
    }
    data(18) = dbTag;
    
    // send the data vector
    res += sChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTimoshenkoBeam3d::sendSelf() - could not send data Vector.\n";
        return res;
    }
    
    // ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(commitTag, sChannel);
    if (res < 0) {
        opserr << "ElasticTimoshenkoBeam3d::sendSelf() - could not send CoordTransf.\n";
        return res;
    }
    
    return res;
}


int ElasticTimoshenkoBeam3d::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
    static Vector data(19);
    res += rChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "ElasticTimoshenkoBeam3d::recvSelf() - could not receive data Vector.\n";
        return res;
    }
    
    this->setTag((int)data(0));
    connectedExternalNodes(0) = (int)data(1);
    connectedExternalNodes(1) = (int)data(2);
    E = data(3);
    G = data(4);
    A = data(5);
    Jx = data(6);
    Iy = data(7);
    Iz = data(8);
    Avy = data(9);
    Avz = data(10);
    rho = data(11);
    cMass = (int)data(12);
    alphaM = data(13);
    betaK  = data(14);
    betaK0 = data(15);
    betaKc = data(16);
    
    // check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(17);
    if (theCoordTransf == 0)  {
        theCoordTransf = theBroker.getNewCrdTransf(crdTag);
        if (theCoordTransf == 0)  {
            opserr << "ElasticTimoshenkoBeam3d::recvSelf() - could not get a CrdTransf3d.\n";
            return -1;
        }
    }
    
    // check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag)  {
        delete theCoordTransf;
        theCoordTransf = theBroker.getNewCrdTransf(crdTag);
        if (theCoordTransf == 0)  {
            opserr << "ElasticTimoshenkoBeam3d::recvSelf() - could not get a CrdTransf3d.\n";
            return -1;
        }
    }
    
    // receive the CoordTransf
    theCoordTransf->setDbTag((int)data(18));
    res += theCoordTransf->recvSelf(commitTag, rChannel, theBroker);
    if (res < 0) {
        opserr << "ElasticTimoshenkoBeam3d::recvSelf() - could not receive CoordTransf.\n";
        return res;
    }
    
    // get coordinate transformation type and save flag
    if (strncmp(theCoordTransf->getClassType(),"Linear",6) == 0)  {
        nlGeo = 0;
    } else if (strncmp(theCoordTransf->getClassType(),"PDelta",6) == 0)  {
        nlGeo = 1;
    } else if (strncmp(theCoordTransf->getClassType(),"Corot",5) == 0)  {
        nlGeo = 1;
        opserr << "\nWARNING ElasticTimoshenkoBeam3d::recvSelf()"
            << " - Element: " << this->getTag() << endln
            << "Unsupported Corotational transformation assigned.\n"
            << "Using PDelta transformation instead.\n";
    }
    
    // revert the CoordTransf to its last committed state
    theCoordTransf->revertToLastCommit();
    
    return res;
}


int ElasticTimoshenkoBeam3d::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)				 
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    static Vector v1(3);
    static Vector v2(3);
    
    if (displayMode >= 0) {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i = 0; i < 3; i++) {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v2(i) = end2Crd(i) + end2Disp(i)*fact;    
        }
    } else {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
            }    
        } else {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i);
                v2(i) = end2Crd(i);
            }    
        }
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag(), 0);
}


void ElasticTimoshenkoBeam3d::Print(OPS_Stream &s, int flag)
{
    if (flag == 0)  {
        // print everything
        s << "Element: " << this->getTag(); 
        s << "  type: ElasticTimoshenkoBeam3d";
        s << "  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  E: " << E << "  G: " << G << endln;
        s << "  A: " << A << "  Jx: " << Jx << "  Iy: " << Iy;
        s << "  Iz: " << Iz << "  Avy: " << Avy << "  Avz: " << Avz << endln;
        s << "  coordTransf: " << theCoordTransf->getClassType() << endln;
        s << "  rho: " << rho << "  cMass: " << cMass << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    } else if (flag == 1)  {
        // does nothing
    }
}


Response* ElasticTimoshenkoBeam3d::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ElasticTimoshenkoBeam3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
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
        
        theResponse =  new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
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
        
        theResponse =  new ElementResponse(this, 2, theVector);
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int ElasticTimoshenkoBeam3d::getResponse (int responseID, Information &eleInfo)
{
    switch (responseID) {
    case 1: // global forces
        return eleInfo.setVector(this->getResistingForce());

    case 2: // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector = ql;
        
        return eleInfo.setVector(theVector);

    default:
        return -1;
    }
}


int ElasticTimoshenkoBeam3d::setParameter(const char **argv,
    int argc, Parameter &param)
{
    if (argc < 1)
        return -1;
    
    // E of the beam
    if (strcmp(argv[0],"E") == 0)
        return param.addObject(1, this);
    
    // G of the beam
    if (strcmp(argv[0],"G") == 0)
        return param.addObject(2, this);
    
    // A of the beam
    if (strcmp(argv[0],"A") == 0)
        return param.addObject(3, this);
    
    // J of the beam
    if (strcmp(argv[0],"J") == 0)
        return param.addObject(4, this);
    
    // Iy of the beam
    if (strcmp(argv[0],"Iy") == 0)
        return param.addObject(5, this);
    
    // Iz of the beam
    if (strcmp(argv[0],"Iz") == 0)
        return param.addObject(6, this);
    
    // Avy of the beam
    if (strcmp(argv[0],"Avy") == 0)
        return param.addObject(7, this);
    
    // Avz of the beam
    if (strcmp(argv[0],"Avz") == 0)
        return param.addObject(8, this);
    
    return -1;
}


int ElasticTimoshenkoBeam3d::updateParameter (int parameterID,
    Information &info)
{
    switch (parameterID) {
    case -1:
        return -1;
    case 1:
        E = info.theDouble;
        return 0;
    case 2:
        G = info.theDouble;
        return 0;
    case 3:
        A = info.theDouble;
        return 0;
    case 4:
        Jx = info.theDouble;
        return 0;
    case 5:
        Iy = info.theDouble;
        return 0;
    case 6:
        Iz = info.theDouble;
        return 0;
    case 7:
        Avy = info.theDouble;
        return 0;
    case 8:
        Avz = info.theDouble;
        return 0;
    default:
        return -1;
    }
}


void ElasticTimoshenkoBeam3d::setUp()
{
    // element projection
    static Vector dx(3);
    
    const Vector &ndICoords = theNodes[0]->getCrds();
    const Vector &ndJCoords = theNodes[1]->getCrds();
    
    dx = ndJCoords - ndICoords;
    
    /*if (nodeJOffset != 0) {
        dx(0) += nodeJOffset[0];
        dx(1) += nodeJOffset[1];
        dx(2) += nodeJOffset[2];
    }
    
    if (nodeIOffset != 0) {
        dx(0) -= nodeIOffset[0];
        dx(1) -= nodeIOffset[1];
        dx(2) -= nodeIOffset[2];
    }
    
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }*/
    
    // determine the element length
    L = dx.Norm();
    if (L == 0.0)  {
        opserr << "WARNING ElasticTimoshenkoBeam3d::setUp()  - "
            << "element: " << this->getTag() << " has zero length.\n";
        return;
    }
    
    // get local axes vectors (these are already normalized)
    static Vector xAxis(3);
    static Vector yAxis(3);
    static Vector zAxis(3);
    theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);
    
    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = xAxis(0);
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = xAxis(1);
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = xAxis(2);
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = yAxis(0);
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = yAxis(1);
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = yAxis(2);
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = zAxis(0);
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = zAxis(1);
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = zAxis(2);
    
    // determine ratios of bending to shear stiffness
    phiY = 12.0*E*Iy/(L*L*G*Avz);
    phiZ = 12.0*E*Iz/(L*L*G*Avy);
    
    // compute initial stiffness matrix in local system
    kl.Zero();
    kl(0,0) = kl(6,6) = E*A/L;
    kl(0,6) = kl(6,0) = -kl(0,0);
    kl(3,3) = kl(9,9) = G*Jx/L;
    kl(3,9) = kl(9,3) = -kl(3,3);
    double a1y = E*Iy/(1.0 + phiY);
    kl(2,2) = kl(8,8) = 12.0*a1y/(L*L*L);
    kl(2,8) = kl(8,2) = -kl(2,2);
    kl(4,4) = kl(10,10) = (4.0 + phiY)*a1y/L;
    kl(4,10) = kl(10,4) = (2.0 - phiY)*a1y/L;
    kl(2,4) = kl(4,2) = kl(2,10) = kl(10,2) = -6.0*a1y/(L*L);
    kl(4,8) = kl(8,4) = kl(8,10) = kl(10,8) = -kl(2,4);
    double a1z = E*Iz/(1.0 + phiZ);
    kl(1,1) = kl(7,7) = 12.0*a1z/(L*L*L);
    kl(1,7) = kl(7,1) = -kl(1,1);
    kl(5,5) = kl(11,11) = (4.0 + phiZ)*a1z/L;
    kl(5,11) = kl(11,5) = (2.0 - phiZ)*a1z/L;
    kl(1,5) = kl(5,1) = kl(1,11) = kl(11,1) = 6.0*a1z/(L*L);
    kl(5,7) = kl(7,5) = kl(7,11) = kl(11,7) = -kl(1,5);
    
    // compute geometric stiffness matrix in local system
    klgeo.Zero();
    if (nlGeo == 1)  {
        double b1y = 1.0/(30.0*L*pow(1.0 + phiY,2));
        klgeo(2,2) = klgeo(8,8) = b1y*(30.0*phiY*phiY + 60.0*phiY + 36.0);
        klgeo(2,8) = klgeo(8,2) = -klgeo(2,2);
        klgeo(4,4) = klgeo(10,10) = b1y*L*L*(2.5*phiY*phiY + 5.0*phiY + 4.0);
        klgeo(4,10) = klgeo(10,4) = -b1y*L*L*(2.5*phiY*phiY + 5.0*phiY + 1.0);
        klgeo(2,4) = klgeo(4,2) = klgeo(2,10) = klgeo(10,2) = -3.0*L;
        klgeo(4,8) = klgeo(8,4) = klgeo(8,10) = klgeo(10,8) = -klgeo(2,4);
        double b1z = 1.0/(30.0*L*pow(1.0 + phiZ,2));
        klgeo(1,1) = klgeo(7,7) = b1z*(30.0*phiZ*phiZ + 60.0*phiZ + 36.0);
        klgeo(1,7) = klgeo(7,1) = -klgeo(1,1);
        klgeo(5,5) = klgeo(11,11) = b1z*L*L*(2.5*phiZ*phiZ + 5.0*phiZ + 4.0);
        klgeo(5,11) = klgeo(11,5) = -b1z*L*L*(2.5*phiZ*phiZ + 5.0*phiZ + 1.0);
        klgeo(1,5) = klgeo(5,1) = klgeo(1,11) = klgeo(11,1) = 3.0*L;
        klgeo(5,7) = klgeo(7,5) = klgeo(7,11) = klgeo(11,7) = -klgeo(1,5);
    }
    
    // compute initial stiffness matrix in global system
    Ki.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    // compute mass matrix in global system
    M.Zero();
    if (rho > 0.0)  {
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*L;
            for (int i=0; i<3; i++)  {
                M(i,i)     = m;
                M(i+6,i+6) = m;
            }
        } else  {
            // consistent mass matrix
            Matrix mlTrn(12,12), mlRot(12,12), ml(12,12);
            mlTrn.Zero(); mlRot.Zero(); ml.Zero();
            double c1x = rho*L/210.0;
            mlTrn(0,0) = mlTrn(6,6) = c1x*70.0;
            mlTrn(0,6) = mlTrn(6,0) = c1x*35.0;
            double c2x = rho/A*Jx*L/210.0;
            mlTrn(3,3) = mlTrn(9,9) = c2x*70.0;
            mlTrn(3,9) = mlTrn(9,3) = c2x*35.0;
            
            double c1y = c1x/pow(1.0 + phiY,2);
            mlTrn(2,2) = mlTrn(8,8) = c1y*(70.0*phiY*phiY + 147.0*phiY + 78.0);
            mlTrn(2,8) = mlTrn(8,2) = c1y*(35.0*phiY*phiY + 63.0*phiY + 27.0);
            mlTrn(4,4) = mlTrn(10,10) = c1y*L*L/4.0*(7.0*phiY*phiY + 14.0*phiY + 8.0);
            mlTrn(4,10) = mlTrn(10,4) = -c1y*L*L/4.0*(7.0*phiY*phiY + 14.0*phiY + 6.0);
            mlTrn(2,4) = mlTrn(4,2) = -c1y*L/4.0*(35.0*phiY*phiY + 77.0*phiY + 44.0);
            mlTrn(8,10) = mlTrn(10,8) = -mlTrn(2,4);
            mlTrn(2,10) = mlTrn(10,2) = c1y*L/4.0*(35.0*phiY*phiY + 63.0*phiY + 26.0);
            mlTrn(4,8) = mlTrn(8,4) = -mlTrn(2,10);
            double c2y = rho/A*Iy/(30.0*L*pow(1.0 + phiY,2));
            mlRot(2,2) = mlRot(8,8) = c2y*36.0;
            mlRot(2,8) = mlRot(8,2) = -mlRot(2,2);
            mlRot(4,4) = mlRot(10,10) = c2y*L*L*(10.0*phiY*phiY + 5.0*phiY + 4.0);
            mlRot(4,10) = mlRot(10,4) = c2y*L*L*(5.0*phiY*phiY - 5.0*phiY - 1.0);
            mlRot(2,4) = mlRot(4,2) = mlRot(2,10) = mlRot(10,2) = c2y*L*(15.0*phiY - 3.0);
            mlRot(4,8) = mlRot(8,4) = mlRot(8,10) = mlRot(10,8) = -mlRot(2,4);
            
            double c1z = c1x/pow(1.0 + phiZ,2);
            mlTrn(1,1) = mlTrn(7,7) = c1z*(70.0*phiZ*phiZ + 147.0*phiZ + 78.0);
            mlTrn(1,7) = mlTrn(7,1) = c1z*(35.0*phiZ*phiZ + 63.0*phiZ + 27.0);
            mlTrn(5,5) = mlTrn(11,11) = c1z*L*L/4.0*(7.0*phiZ*phiZ + 14.0*phiZ + 8.0);
            mlTrn(5,11) = mlTrn(11,5) = -c1z*L*L/4.0*(7.0*phiZ*phiZ + 14.0*phiZ + 6.0);
            mlTrn(1,5) = mlTrn(5,1) = c1z*L/4.0*(35.0*phiZ*phiZ + 77.0*phiZ + 44.0);
            mlTrn(7,11) = mlTrn(11,7) = -mlTrn(1,5);
            mlTrn(1,11) = mlTrn(11,1) = -c1z*L/4.0*(35.0*phiZ*phiZ + 63.0*phiZ + 26.0);
            mlTrn(5,7) = mlTrn(7,5) = -mlTrn(1,11);
            double c2z = rho/A*Iz/(30.0*L*pow(1.0 + phiZ,2));
            mlRot(1,1) = mlRot(7,7) = c2z*36.0;
            mlRot(1,7) = mlRot(7,1) = -mlRot(1,1);
            mlRot(5,5) = mlRot(11,11) = c2z*L*L*(10.0*phiZ*phiZ + 5.0*phiZ + 4.0);
            mlRot(5,11) = mlRot(11,5) = c2z*L*L*(5.0*phiZ*phiZ - 5.0*phiZ - 1.0);
            mlRot(1,5) = mlRot(5,1) = mlRot(1,11) = mlRot(11,1) = -c2z*L*(15.0*phiZ - 3.0);
            mlRot(5,7) = mlRot(7,5) = mlRot(7,11) = mlRot(11,7) = -mlRot(1,5);
            // add translational and rotational parts
            ml = mlTrn + mlRot;
            // transform from local to global system
            M.addMatrixTripleProduct(0.0, Tgl, ml, 1.0);
        }
    }
}
