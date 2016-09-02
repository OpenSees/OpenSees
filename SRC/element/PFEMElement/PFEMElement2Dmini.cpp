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

// $Revision: 1.0 $
// $Date: 2014/10/1 9:45:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2Dmini.h,v $

// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2Dmini.

#include "PFEMElement2Dmini.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Renderer.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Matrix PFEMElement2Dmini::K;
Vector PFEMElement2Dmini::P;


void* OPS_PFEMElement2DMini()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 8) {
	opserr<<"WARNING: insufficient number of arguments\n";
	return 0;
    }

    // tag, nd1, nd2, nd3
    numdata = 4;
    int idata[4];
    if(OPS_GetIntInput(&numdata,idata)<0) return 0;

    // rho, mu, b1, b2, (thinknes,kappa)
    numdata = OPS_GetNumRemainingInputArgs();
    if(numdata > 6) numdata = 6;
    double data[6] = {0,0,0,0,1.0,1e9};
    if(OPS_GetDoubleInput(&numdata,data) < 0) return 0;

    return new PFEMElement2Dmini(idata[0],idata[1],idata[2],idata[3],
				 data[0],data[1],data[2],data[3],data[4],data[5]);
}


// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2Dmini::PFEMElement2Dmini()
    :Element(0, ELE_TAG_PFEMElement2Dmini), ntags(6), nodes(6,static_cast<Node*>(0)),
     thePCs(3,static_cast<Pressure_Constraint*>(0)), rho(0), mu(0), body(2),
     cc(3), dd(3), J(0.0), Jn(0.0), vdof(6), pdof(3), thickness(1.0),
     kappa(0.0), ndf(0), lumped(false), checkJ(false)
{
}

// for object
PFEMElement2Dmini::PFEMElement2Dmini(
    int tag, int nd1, int nd2, int nd3, 
    double r, double m, double b1,
    double b2, double thk, double ka, bool lmpd, bool chk)
    :Element(tag, ELE_TAG_PFEMElement2Dmini), ntags(6), nodes(6,static_cast<Node*>(0)),
     thePCs(3,static_cast<Pressure_Constraint*>(0)), rho(r), mu(m), body(2),
     cc(3), dd(3), J(0.0), Jn(0.0), vdof(6), pdof(3), thickness(thk),
     kappa(ka), ndf(0), lumped(lmpd), checkJ(chk)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    for(int i=0;i<3;i++) {
        ntags(2*i+1) = ntags(2*i);
    }
    body(0) = b1;
    body(1) = b2;
    if(ops_Dt == 0) ops_Dt = 1.0;
}


PFEMElement2Dmini::~PFEMElement2Dmini()
{
    for(int i=0; i<3; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}



int
PFEMElement2Dmini::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement2Dmini::getExternalNodes()
{
    return ntags;
}

Node **
PFEMElement2Dmini::getNodePtrs(void)
{
    return &nodes[0];
}

int
PFEMElement2Dmini::getNumDOF()
{
    return ndf;
}

int
PFEMElement2Dmini::revertToLastCommit()
{
    return 0;
}

int PFEMElement2Dmini::commitState()
{
    Jn = J;
    return Element::commitState();
}

int
PFEMElement2Dmini::update()
{
    // get nodal coordinates
    Vector x(3), y(3), xn(3), yn(3);
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
        const Vector& disp = nodes[2*i]->getTrialDisp();
	const Vector& dispn = nodes[2*i]->getDisp();
        x(i) = coord(0) + disp(0);
        y(i) = coord(1) + disp(1);
	xn(i) = coord(0) + dispn(0);
        yn(i) = coord(1) + dispn(1);
    }

    // get c and d
    cc(0) = y(1)-y(2);
    dd(0) = x(2)-x(1);
    cc(1) = y(2)-y(0);
    dd(1) = x(0)-x(2);
    cc(2) = y(0)-y(1);
    dd(2) = x(1)-x(0);

    // get Jacobi
    J = cc(0)*dd(1)-dd(0)*cc(1);
    Jn = J;

    if(checkJ && J<=0) {
    	opserr<<"WARNING: element "<<this->getTag()<<" Jacobian determinant ";
    	opserr<<J<<" <= 0\n";
	opserr<<"J = "<<J<<"\n";
	opserr<<"x = "<<x;
	opserr<<"y = "<<y;
	opserr<<"Jn = "<<Jn<<"\n";
	opserr<<"xn = "<<xn;
	opserr<<"yn = "<<yn;
    	return -1;
    }

    return 0;
}

const Matrix&
PFEMElement2Dmini::getMass()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M
    double m = rho*Jn*thickness/6.0;
    for(int i=0; i<vdof.Size(); i++) {
	K(vdof(i),vdof(i)) = m;
    }

    // Mp
    double mp = J*thickness/kappa/24.0;
    if (kappa <= 0) {
	mp = 0.0;
    }
    
    for(int i=0; i<pdof.Size(); i++) {
	for (int j=0; j<pdof.Size(); j++) {
	    if (i==j) {
		K(pdof(i),pdof(j)) = mp*2;
	    } else {
		K(pdof(i),pdof(j)) = mp;
	    }
	    
	}
    }

    return K;
}

const Matrix&
PFEMElement2Dmini::getDamp()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // invMb
    Matrix invMb(2,2);
    double mb = 9.*rho*Jn*thickness/(40.*ops_Dt);
    double kb = 27.*mu*thickness/(40.*J);
    for(int i=0; i<2; i++) {
	invMb(i,i) = mb;
    }
    if(kb > 0) {
	double cc2 = 0., dd2 = 0., cd2 = 0.;
	for(int a=0; a<3; a++) {
	    cc2 += cc(a)*cc(a);
	    dd2 += dd(a)*dd(a);
	    cd2 += cc(a)*dd(a);
	}
	invMb(0,0) += kb*(4*cc2+3*dd2); // Kxx
	invMb(0,1) += kb*cd2; // Kxy
	invMb(1,0) += kb*cd2; // Kyx
	invMb(1,1) += kb*(3*cc2+4*dd2); // Kyy
    }
    invMb(0,1) = invMb(0,0);
    invMb(0,0) = invMb(1,1);
    invMb(1,1) = invMb(0,1);
    invMb(0,1) = -invMb(1,0);
    invMb(1,0) = -invMb(1,0);
    invMb /= invMb(0,0)*invMb(1,1)-invMb(0,1)*invMb(1,0);

    // G, S
    Matrix S(3,3);
    Matrix Gb(2,3);
    double gb = -9.*thickness/40.0;
    for(int a=0; a<1; a++) {
	for(int b=0; b<3; b++) {
	    Gb(2*a,b) = cc(b)*gb;
	    Gb(2*a+1,b) = dd(b)*gb;
	}
    }
    S.addMatrixTripleProduct(0.0, Gb, invMb, 1.0);
    
    // K, G, Gt, S
    double k = mu*thickness/(6.*J);
    double g = thickness/6.0;
    for(int a=0; a<3; a++) {
	for(int b=0; b<3; b++) {

	    // K
	    if(!lumped) {
		K(vdof(2*a),vdof(2*b)) = k*(4*cc(a)*cc(b)+3*dd(a)*dd(b)); // Kxx
		K(vdof(2*a),vdof(2*b+1)) = k*(3*dd(a)*cc(b)-2*cc(a)*dd(b)); // Kxy
		K(vdof(2*a+1),vdof(2*b)) = k*(3*cc(a)*dd(b)-2*dd(a)*cc(b)); // Kyx
		K(vdof(2*a+1),vdof(2*b+1)) = k*(3*cc(a)*cc(b)+4*dd(a)*dd(b)); // Kyy
	    }

	    // -G
	    K(vdof(2*a),pdof(b)) = -cc(a)*g;
	    K(vdof(2*a+1),pdof(b)) = -dd(a)*g;

	    // Gt
	    K(pdof(b),vdof(2*a)) = -K(vdof(2*a),pdof(b));
	    K(pdof(b),vdof(2*a+1)) = -K(vdof(2*a+1),pdof(b));

	    // S
	    K(pdof(a), pdof(b)) = S(a,b);
	}
    }

    return K;
}

const Matrix&
PFEMElement2Dmini::getTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();
    
    return K;
}


const Matrix&
PFEMElement2Dmini::getInitialStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement2Dmini::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2Dmini::getResistingForce()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2Dmini::getResistingForceIncInertia()
{
    // resize p
    P.resize(ndf);
    P.Zero();
    
    // get velocity, accleration
    Vector v(ndf), vdot(ndf);
    for(int a=0; a<3; a++) {

        const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& accel = nodes[2*a]->getTrialAccel();
	for(int i=0; i<2; i++) {
	    v(vdof(2*a+i)) = vel(i);
	    vdot(vdof(2*a+i)) = accel(i);
	}
	const Vector& pressure = nodes[2*a+1]->getTrialVel();
	const Vector& pressuredot = nodes[2*a+1]->getTrialAccel();
	v(pdof(a)) = pressure(0);
	vdot(pdof(a)) = pressuredot(0);
    }

    // M*vdot+K*v
    P.addMatrixVector(0.0, getMass(), vdot, 1.0);
    bool l = lumped;
    lumped = false;
    P.addMatrixVector(1.0, getDamp(), v, 1.0);
    lumped = l;

    // invMb
    Matrix invMb(2,2);
    double mb = 9.*rho*Jn*thickness/(40.*ops_Dt);
    double kb = 27.*mu*thickness/(40.*J);
    for(int i=0; i<2; i++) {
	invMb(i,i) = mb;
    }
    if(kb > 0) {
	double cc2 = 0., dd2 = 0., cd2 = 0.;
	for(int a=0; a<3; a++) {
	    cc2 += cc(a)*cc(a);
	    dd2 += dd(a)*dd(a);
	    cd2 += cc(a)*dd(a);
	}
	invMb(0,0) += kb*(4*cc2+3*dd2); // Kxx
	invMb(0,1) += kb*cd2; // Kxy
	invMb(1,0) += kb*cd2; // Kyx
	invMb(1,1) += kb*(3*cc2+4*dd2); // Kyy
    }
    invMb(0,1) = invMb(0,0);
    invMb(0,0) = invMb(1,1);
    invMb(1,1) = invMb(0,1);
    invMb(0,1) = -invMb(1,0);
    invMb(1,0) = -invMb(1,0);
    invMb /= invMb(0,0)*invMb(1,1)-invMb(0,1)*invMb(1,0);

    // F-M*vdot-K*v
    Vector Fb(2), Fp(3);
    Matrix Gb(2,3);
    double f = rho*Jn*thickness/6.;
    double fb = 9.*rho*Jn*thickness/40.;
    double gb = -9.*thickness/40.0;
    for(int a=0; a<1; a++) {
	for(int b=0; b<3; b++) {
	    Gb(2*a,b) = cc(b)*gb;
	    Gb(2*a+1,b) = dd(b)*gb;
	}
    }
    for(int i=0; i<2; i++) {
	Fb(i) = fb*body(i);
    }
    Fp.addMatrixTransposeVector(0.0, Gb, invMb*Fb, -1.0);
    
    for(int a=0; a<3; a++) {
	for(int i=0; i<2; i++) {
	    P(vdof(2*a+i)) -= f*body(i);
	}
	P(pdof(a)) -= Fp(a);
    }


    
    return P;
}


const char*
PFEMElement2Dmini::getClassType()const
{
    return "PFEMElement2Dmini";
}

int
PFEMElement2Dmini::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PFEMElement2Dmini::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement2Dmini::setDomain(Domain *theDomain)
{
    this->DomainComponent::setDomain(theDomain);
    if(theDomain == 0) {return;}

    Vector x(3), y(3);

    ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<3; i++) {

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement2Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
        if(nodes[2*i]->getNumberDOF() < 2) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" ndf < 2 ";
            opserr<<"in PFEMElement2Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
	vdof(2*i) = ndf;
	vdof(2*i+1) = ndf+1;
        ndf += nodes[2*i]->getNumberDOF();

	// material coordinates
	const Vector& coord = nodes[2*i]->getCrds();
	const Vector& disp = nodes[2*i]->getDisp();
	x(i) = coord(0)+disp(0);
	y(i) = coord(1)+disp(1);

        // get pc
        int pndf = 1;
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] != 0) {
            thePCs[i]->setDomain(theDomain);
        } else {
            thePCs[i] = new Pressure_Constraint(ntags(2*i), pndf);
            if(thePCs[i] == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMElement2Dmini::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement2Dmini::setDomain "<<eletag<<"\n";
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
            opserr<<"in PFEMElement2Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
	pdof(i) = ndf;
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    // Jacobian
    Jn = (y(1)-y(2))*(x(0)-x(2))-(x(2)-x(1))*(y(2)-y(0));
}

void
PFEMElement2Dmini::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2Dmini: "<<this->getTag()<<endln;
}

int
PFEMElement2Dmini::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    return 0;
}
