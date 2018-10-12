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

// $Revision $
// $Date$
// $URL$
                                                                        
// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: This MINI element with quasi-incompressible form

#include "PFEMElement2DQuasi.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Renderer.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <cmath>
#include <ElementIter.h>
#include <NodeIter.h>
#include <Parameter.h>

Matrix PFEMElement2DQuasi::K;
Vector PFEMElement2DQuasi::P;

void* OPS_PFEMElement2DQuasi()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 8) {
	opserr<<"WARNING: insufficient number of arguments\n";
	return 0;
    }

    // tag, nd1, nd2, nd3
    numdata = 4;
    int idata[4];
    if(OPS_GetIntInput(&numdata,idata)<0) {
	opserr << "WARNING: failed to read integers -- PFEMElement2DQuasi\n";
	return 0;
    }

    // rho, mu, b1, b2, (thickness, kappa)
    numdata = OPS_GetNumRemainingInputArgs();
    if(numdata > 6) numdata = 6;
    double data[6] = {0,0,0,0,1.0,2.15e9};
    if(OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr << "WARNING: failed to read doubles -- PFEMElement2DQuasi\n";
	return 0;
    }

    return new PFEMElement2DQuasi(idata[0],idata[1],idata[2],idata[3],
				  data[0],data[1],data[2],data[3],data[4],
				  data[5]);
}

int OPS_PFEMElement2DQuasi(Domain& theDomain, const ID& elenodes, ID& eletags)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 4) {
	opserr<<"WARNING: insufficient number of arguments\n";
	return 0;
    }

    // rho, mu, b1, b2, (thickness, kappa)
    numdata = OPS_GetNumRemainingInputArgs();
    if(numdata > 6) numdata = 6;
    double data[6] = {0,0,0,0,1.0,2.15e9};
    if(OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr << "WARNING: failed to read doubles -- PFEMElement2DQuasi\n";
	return 0;
    }

    // create elements
    ElementIter& theEles = theDomain.getElements();
    Element* theEle = theEles();
    int currTag = 0;
    if (theEle != 0) {
	currTag = theEle->getTag();
    }

    eletags.resize(elenodes.Size()/3);

    for (int i=0; i<eletags.Size(); i++) {
	theEle = new PFEMElement2DQuasi(--currTag,elenodes(3*i),elenodes(3*i+1),elenodes(3*i+2),data[0],data[1],data[2],data[3],data[4],data[5]);
	if (theEle == 0) {
	    opserr<<"WARING: run out of memory for creating element\n";
	    return -1;
	}
	if (theDomain.addElement(theEle) == false) {
	    opserr<<"WARNING: failed to add element to domain\n";
	    delete theEle;
	    return -1;
	}
	eletags(i) = currTag;
    }

    return 0;
}

// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2DQuasi::PFEMElement2DQuasi()
    :Element(0, ELE_TAG_PFEMElement2DQuasi), ntags(7),
     rho(0), mu(0), b1(0), b2(0), thickness(1.0), kappa(2.15e9),
     ndf(0), J(0.0), tangent(true)
{
    for(int i=0;i<4;i++)
    {
	nodes[2*i] = 0;
	ntags(2*i) = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	if (i < 3) {
	    nodes[2*i+1] = 0;
	    ntags(2*i+1) = 0;
	    thePCs[i] = 0;
	    cc[i] = 0;
	    dd[i] = 0;
	    pdof[i] = 0;
	}

    }
}

// for object
PFEMElement2DQuasi::PFEMElement2DQuasi(
    int tag, int nd1, int nd2, int nd3,
    double r, double m, double bx, double by, double thk,
    double ka)
    :Element(tag, ELE_TAG_PFEMElement2DQuasi), ntags(7),
     rho(r), mu(m), b1(bx), b2(by), thickness(thk), kappa(ka),
     ndf(0), J(0.0), tangent(true)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3; ntags(6) = nd3;
    for(int i=0;i<4;i++)
    {
	nodes[2*i] = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	if (i < 3) {
	    nodes[2*i+1] = 0;
	    ntags(2*i+1) = ntags(2*i);
	    thePCs[i] = 0;
	    cc[i] = 0;
	    dd[i] = 0;
	    pdof[i] = 0;
	}
    }
    if (kappa <= 0.0) {
	kappa = 2.15e9;
    }
}


PFEMElement2DQuasi::~PFEMElement2DQuasi()
{
    for(int i=0; i<3; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
    Domain* domain = this->getDomain();
    if (domain != 0) {
	if (nodes[6] != 0) {
	    domain->removeNode(nodes[6]->getTag());
	    delete nodes[6];
	    nodes[6] = 0;
	}
    }
}

int
PFEMElement2DQuasi::commitState()
{
    Vector bdisp(2);
    for(int a=0; a<3; a++) {
        const Vector& disp = nodes[2*a]->getTrialDisp();
	bdisp(0) += disp(0);
	bdisp(1) += disp(1);
    }
    bdisp /= 3.0;
    nodes[6]->setTrialDisp(bdisp);
    nodes[6]->commitState();
    return Element::commitState();
}

int
PFEMElement2DQuasi::update()
{
    // get nodal coordinates
    double x[3], y[3];
    for(int a=0; a<3; a++) {
        const Vector& coord = nodes[2*a]->getCrds();
        const Vector& disp = nodes[2*a]->getTrialDisp();
        x[a] = coord(0) + disp(0);
        y[a] = coord(1) + disp(1);
    }

    // get c and d
    cc[0] = y[1]-y[2];
    dd[0] = x[2]-x[1];
    cc[1] = y[2]-y[0];
    dd[1] = x[0]-x[2];
    cc[2] = y[0]-y[1];
    dd[2] = x[1]-x[0];

    // get Jacobi
    J = cc[0]*dd[1]-dd[0]*cc[1];

    // check Jacobi
    if(fabs(J)<1e-15) {
	//if(J < 0) {
        opserr<<"WARING: element area is negative";
        opserr<<" -- PFEMElement2DQuasi::update\n";
        return -1;
    }

    return 0;
}

const Matrix&
PFEMElement2DQuasi::getMass()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M and Mp
    double M = rho*J*thickness/24.0;
    double M2 = 2.0*M;
    double M3 = 3.0*rho*J*thickness/40.0;
    double M4 = 81.0*rho*J*thickness/560.0;
    double Mp = J*thickness/24.0/kappa;
    double Mp2 = 2.0*Mp;
    for (int a=0; a<3; a++) {
	// M(a,b) a,b=1,2,3
	for (int b=0; b<3; b++) {
	    if (a == b) {
		K(vxdof[a], vxdof[b]) = M2;
		K(vydof[a], vydof[b]) = M2;
		K(pdof[a],pdof[b]) = Mp2;
	    } else {
		K(vxdof[a], vxdof[b]) = M;
		K(vydof[a], vydof[b]) = M;
		K(pdof[a],pdof[b]) = Mp;
	    }
	}
	// M(4,a)
	K(vxdof[3],vxdof[a]) = M3;
	K(vydof[3],vydof[a]) = M3;

	// M(a,4)
	K(vxdof[a],vxdof[3]) = M3;
	K(vydof[a],vydof[3]) = M3;
    }

    // M(4,4)
    K(vxdof[3],vxdof[3]) = M4;
    K(vydof[3],vydof[3]) = M4;

    return K;
}

const Matrix&
PFEMElement2DQuasi::getDamp()
{

    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // Km
    double km = mu*thickness/(6.*J);
    double kb = 27.0*mu*thickness/(40.*J);

    // Kp
    double kp = ops_Dt*kappa*thickness/(2.*J);
    double kpb = 81.0*ops_Dt*kappa*thickness/(40.*J);

    // G
    double g = thickness/6.0;
    double gb = -9.0*thickness/40.0;

    // K(a,b) a,b = 1,2,3
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {

	    // K
	    K(vxdof[a], vxdof[b]) += km*(4*cc[a]*cc[b]+3*dd[a]*dd[b]); // Kxx
	    K(vxdof[a], vydof[b]) += km*(3*dd[a]*cc[b]-2*cc[a]*dd[b]); // Kxy
	    K(vydof[a], vxdof[b]) += km*(3*cc[a]*dd[b]-2*dd[a]*cc[b]); // Kyx
	    K(vydof[a], vydof[b]) += km*(3*cc[a]*cc[b]+4*dd[a]*dd[b]); // Kyy

	    // Kp
	    if (tangent) {
		K(vxdof[a], vxdof[b]) += kp*cc[a]*cc[b]; // Kxx
		K(vxdof[a], vydof[b]) += kp*cc[a]*dd[b]; // Kxy
		K(vydof[a], vxdof[b]) += kp*dd[a]*cc[b]; // Kyx
		K(vydof[a], vydof[b]) += kp*dd[a]*dd[b]; // Kyy
	    }

	    // -G 
	    K(vxdof[a], pdof[b]) = -g*cc[a];
	    K(vydof[a], pdof[b]) = -g*dd[a];

	    // Gt
	    K(pdof[b], vxdof[a]) = g*cc[a];
	    K(pdof[b], vydof[a]) = g*dd[a];
	}

	// -G(4,a)
	K(vxdof[3], pdof[a]) = -gb*cc[a];
	K(vydof[3], pdof[a]) = -gb*dd[a];

	// Gt(a,4)
	K(pdof[a], vxdof[3]) = gb*cc[a];
	K(pdof[a], vydof[3]) = gb*dd[a];
    }

    // K44
    double cc2=0.0, dd2=0.0, cd2=0.0;
    for (int a=0; a<3; a++) {
	cc2 += cc[a]*cc[a];
	dd2 += dd[a]*dd[a];
	cd2 += cc[a]*dd[a];
    }

    K(vxdof[3], vxdof[3]) += kb*(4*cc2+3*dd2); // Kxx
    K(vxdof[3], vydof[3]) += kb*cd2; // Kxy
    K(vydof[3], vxdof[3]) += kb*cd2; // Kyx
    K(vydof[3], vydof[3]) += kb*(3*cc2+4*dd2); // Kyy

    // Kp44
    if (tangent) {
	K(vxdof[3], vxdof[3]) += kpb*cc2; // Kxx
	K(vxdof[3], vydof[3]) += kpb*cd2; // Kxy
	K(vydof[3], vxdof[3]) += kpb*cd2; // Kyx
	K(vydof[3], vydof[3]) += kpb*dd2; // Kyy
    }

    return K;
}

const Matrix&
PFEMElement2DQuasi::getTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

const Matrix&
PFEMElement2DQuasi::getInitialStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement2DQuasi::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2DQuasi::getResistingForce()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2DQuasi::getResistingForceIncInertia()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    // F
    double f = rho*J*thickness/6.0;
    for (int a=0; a<3; a++) {
	P(vxdof[a]) = f*b1;
	P(vydof[a]) = f*b2;
    }

    // F4
    double fb = 9.0*rho*J*thickness/40.0;
    P(vxdof[3]) = fb*b1;
    P(vydof[3]) = fb*b2;
    
    // vdot, v
    Vector vdot(ndf), v(ndf);
    for(int a=0; a<4; a++) {
	const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& accel = nodes[2*a]->getTrialAccel();
	
	vdot(vxdof[a]) = accel(0);
	vdot(vydof[a]) = accel(1);

	v(vxdof[a]) = vel(0);
	v(vydof[a]) = vel(1);

	if (a < 3) {
	    const Vector& p = nodes[2*a+1]->getTrialVel();
	    const Vector& pdot = nodes[2*a+1]->getTrialAccel();
	    
	    vdot(pdof[a]) = pdot(0);
	    v(pdof[a]) = p(0);
	}
    }

    // -r = M*vdot+K*v-F
    tangent = false;
    P.addMatrixVector(-1.0, this->getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, this->getDamp(), v, 1.0);
    tangent = true;

    return P;
}


const char*
PFEMElement2DQuasi::getClassType()const
{
    return "PFEMElement2DQuasi";
}

int
PFEMElement2DQuasi::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PFEMElement2DQuasi::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement2DQuasi::setDomain(Domain *theDomain)
{
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    ndf = 0;
    int eletag = this->getTag();
    Vector bcrds(2);
    for(int i=0; i<3; i++) {

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement2DQuasi - setDomain() "<<eletag<<"\n ";
            return;
        }
	if(nodes[2*i]->getNumberDOF() < 2) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" ndf < 2 ";
            opserr<<"in PFEMElement2DQuasi - setDomain() "<<eletag<<"\n ";
            return;
        }
	const Vector& crds = nodes[2*i]->getCrds();
	for (int j=0; j<2; j++) {
	    bcrds(j) += crds(j);
	}
	vxdof[i] = ndf;
	vydof[i] = ndf+1;
        ndf += nodes[2*i]->getNumberDOF();

        // get pc
        int pndf = 1;
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] != 0) {
            thePCs[i]->setDomain(theDomain);
        } else {
            thePCs[i] = new Pressure_Constraint(ntags(2*i), pndf);
            if(thePCs[i] == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMElement2DQuasi::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement2DQuasi::setDomain "<<eletag<<"\n";
                delete thePCs[i];
                thePCs[i] = 0;
                return;
            }
        }

        // connect
        thePCs[i]->connect(eletag);

        // get pressure node
        nodes[2*i+1] = thePCs[i]->getPressureNode();
        if(nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in PFEMElement2DQuasi - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
	pdof[i] = ndf;
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    // get a node tag
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = theNodes();
    ntags[6] = 0;
    if (theNode != 0) {
	ntags[6] = theNode->getTag();
    }
    ntags[6]--;

    // create bubble node
    bcrds /= 3.0;
    nodes[6] = new Node(ntags[6],2,bcrds(0),bcrds(1));
    if (nodes[6] == 0) {
	opserr<<"WARNING: run out of memory in creating node\n";
	return;
    }
    if (theDomain->addNode(nodes[6]) == false) {
	opserr<<"WARNING: failed to add node to domain\n";
	delete nodes[6];
	nodes[6] = 0;
    }
    vxdof[3] = ndf++;
    vydof[3] = ndf++;

}

void
PFEMElement2DQuasi::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2DQuasi: "<<this->getTag()<<endln;
}

int
PFEMElement2DQuasi::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
    return 0;
}

int
PFEMElement2DQuasi::setParameter(const char **argv, int argc, 
					Parameter &parameter)
{
    return 0;
}

int
PFEMElement2DQuasi::updateParameter(int passparameterID, Information &info)
{
    return 0;
}

int
PFEMElement2DQuasi::activateParameter(int passedParameterID)
{
    return 0;
}

const Matrix&
PFEMElement2DQuasi::getDampSensitivity(int gradNumber)
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();
    return K;
}

const Matrix&
PFEMElement2DQuasi::getMassSensitivity(int gradNumber)
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();


    return K;
}

const Vector&
PFEMElement2DQuasi::getResistingForceSensitivity(int gradNumber)
{
    // resize P
    P.resize(ndf);
    P.Zero();
    return P;    
}

int
PFEMElement2DQuasi::commitSensitivity(int gradNumber, int numGrads)
{
    return 0;
}

