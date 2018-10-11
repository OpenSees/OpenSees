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
// $Date: 2015/3/10 10:23:06 $
                                                                        
// Written: Minjie Zhu
//
// Description: This file contains the class definition for PFEMElement3Dmini.

#include "PFEMElement3Dmini.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Renderer.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// #define debug3dmini

Matrix PFEMElement3Dmini::K;
Vector PFEMElement3Dmini::P;


// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement3Dmini::PFEMElement3Dmini()
    :Element(0, ELE_TAG_PFEMElement3Dmini), ntags(8), nodes(8,static_cast<Node*>(0)),
     thePCs(4,static_cast<Pressure_Constraint*>(0)), rho(0), mu(0), body(3),
     J(0.0), Jn(0.0), vdof(12), pdof(4),
     kappa(0.0), ndf(0), lumped(false), checkJ(false), dNdx()
{
}

// for object
PFEMElement3Dmini::PFEMElement3Dmini(
    int tag, int nd1, int nd2, int nd3, int nd4,
    double r, double m, double b1,
    double b2, double b3, double ka, bool lmpd, bool chk)
    :Element(tag, ELE_TAG_PFEMElement3Dmini), ntags(8), nodes(8,static_cast<Node*>(0)),
     thePCs(4,static_cast<Pressure_Constraint*>(0)), rho(r), mu(m), body(3),
     J(0.0), Jn(0.0), vdof(12), pdof(4),
     kappa(ka), ndf(0), lumped(lmpd), checkJ(chk), dNdx()
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3; ntags(6)=nd4;
    for(int i=0;i<ntags.Size()/2;i++) {
        ntags(2*i+1) = ntags(2*i);
    }
    body(0) = b1;
    body(1) = b2;
    body(2) = b3;
}


PFEMElement3Dmini::~PFEMElement3Dmini()
{
    for(int i=0; i<(int)thePCs.size(); i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}



int
PFEMElement3Dmini::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement3Dmini::getExternalNodes()
{
    return ntags;
}

Node **
PFEMElement3Dmini::getNodePtrs(void)
{
    return &nodes[0];
}

int
PFEMElement3Dmini::getNumDOF()
{
    return ndf;
}

int
PFEMElement3Dmini::revertToLastCommit()
{
    return 0;
}

int PFEMElement3Dmini::commitState()
{
    Jn = J;
    return Element::commitState();
}

int
PFEMElement3Dmini::update()
{
    int numnodes = ntags.Size()/2;
    int ndm = body.Size();
    
    // get nodal coordinates
    Vector x(numnodes), y(numnodes), z(numnodes);
    for(int i=0; i<numnodes; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
        const Vector& disp = nodes[2*i]->getTrialDisp();
	x(i) = coord(0)+disp(0);
	y(i) = coord(1)+disp(1);
	z(i) = coord(2)+disp(2);
    }

    // get Jacobi
    J = det(x(1)-x(0),y(1)-y(0),z(1)-z(0),
	     x(2)-x(0),y(2)-y(0),z(2)-z(0),
	     x(3)-x(0),y(3)-y(0),z(3)-z(0));

    if(checkJ && J<=0) {
    	opserr<<"WARNING: element "<<this->getTag()<<" Jacobian determinant ";
    	opserr<<J<<" <= 0\n";
    	opserr<<"J = "<<J<<"\n";
    	opserr<<"Jn = "<<Jn<<"\n";
    	opserr<<"x = "<<x;
    	opserr<<", y = "<<y;
	opserr<<", z = "<<z<<"\n";
    	return -1;
    }

    // deformation tensor
    Matrix F1(ndm,ndm), F(ndm,ndm);
    dNdx.resize(ndm,numnodes); dNdx.Zero();
    for(int i=0; i<ndm; i++) {
	F1(i,0) = x(i+1)-x(0);
	F1(i,1) = y(i+1)-y(0);
	F1(i,2) = z(i+1)-z(0);
    }
    F1.Invert(F);

    // get dNdx
    for(int i=0; i<ndm; i++) {
	for(int a=1; a<numnodes; a++) {
	    dNdx(i,a) = F(i,a-1)*J;
	    dNdx(i,0) -= F(i,a-1)*J;
	}
    }
    
#ifdef debug3dmini
    Vector xn(numnodes), yn(numnodes), zn(numnodes);
    for(int i=0; i<numnodes; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
	const Vector& dispn = nodes[2*i]->getDisp();
	xn(i) = coord(0)+dispn(0);
	yn(i) = coord(1)+dispn(1);
	zn(i) = coord(2)+dispn(2);
    }
    Matrix Fn(3,3), invFn(3,3);
    for(int i=0; i<3; i++) {
	Fn(i,0) = xn(i+1)-xn(0);
	Fn(i,1) = yn(i+1)-yn(0);
	Fn(i,2) = zn(i+1)-zn(0);
    }
    Fn.Invert(invFn);
    opserr<<"Fkj = "<<F;
    opserr<<"J = "<<J<<"\n";
    opserr<<"FiA = "<<F1*invFn;
    opserr<<"Jn = "<<Jn<<"\n";
    exit(0);
#endif

    return 0;
}

const Matrix&
PFEMElement3Dmini::getMass()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M
    double m = rho*J/24.0;
    for(int i=0; i<vdof.Size(); i++) {
    	K(vdof(i),vdof(i)) = m;
    }

    // Mp
    if(kappa > 0) {
    	double mp = J/kappa/24.0;
    	for(int i=0; i<pdof.Size(); i++) {
    	    K(pdof(i),pdof(i)) = mp;
    	}
    }

    return K;
}

const Matrix&
PFEMElement3Dmini::getDamp()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    int numnodes = ntags.Size()/2;
    int ndm = body.Size();

    // invMb
    Matrix invMb(ndm,ndm), Mb(ndm,ndm);
    double mb = 16.*rho*Jn/(315.*ops_Dt);
    double kb = 1024.*mu/(2835.*J);
    for(int i=0; i<ndm; i++) {
    	Mb(i,i) = mb;
    }
    if(mu > 0) {
	for(int i=0; i<ndm; i++) {
	    for(int l=0; l<ndm; l++) {
		Mb(i,i) += kb*kbb(l,l);
	    }
	    for(int j=0; j<ndm; j++) {
		Mb(i,j) += kb*kbb(j,i);
		Mb(i,j) -= 2.0*kb*kbb(i,j)/3.0;
	    }
	}
    }
    Mb.Invert(invMb);

    // Gb, S
    Matrix S(numnodes,numnodes);
    Matrix Gb(ndm,numnodes);
    double gb = -16./315.0;
    for(int i=0; i<ndm; i++) {
    	for(int b=0; b<numnodes; b++) {
	    Gb(i,b) = dNdx(i,b)*gb;
    	}
    }
    S.addMatrixTripleProduct(0.0, Gb, invMb, 1.0);

    // K, G, Gt, S
    double k = mu/(6.*J);
    double g = 1.0/24.0;
    for(int a=0; a<numnodes; a++) {
    	for(int b=0; b<numnodes; b++) {

    	    // K
	    if(mu > 0) {
		if(lumped) {
		    
		} else {
		    for(int i=0; i<ndm; i++) {
		    	for(int l=0; l<ndm; l++) {
		    	    K(vdof(ndm*a+i),vdof(ndm*b+i)) += k*dNdx(l,a)*dNdx(l,b);
		    	}
		    	for(int j=0; j<ndm; j++) {
		    	    K(vdof(ndm*a+i),vdof(ndm*b+j)) += k*dNdx(j,a)*dNdx(i,b);
		    	    K(vdof(ndm*a+i),vdof(ndm*b+j)) -= 2.0*k*dNdx(i,a)*dNdx(j,b)/3.0;
		    	}
		    }
		}
	    }

    	    // -G and Gt
	    for(int i=0; i<ndm; i++) {
		K(vdof(ndm*a+i),pdof(b)) = -dNdx(i,a)*g;
		K(pdof(b),vdof(ndm*a+i)) = -K(vdof(ndm*a+i),pdof(b));
	    }

    	    // S
    	    K(pdof(a), pdof(b)) = S(a,b);
    	}
    }

    return K;
}

const Matrix&
PFEMElement3Dmini::getTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();
    
    return K;
}


const Matrix&
PFEMElement3Dmini::getInitialStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement3Dmini::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement3Dmini::getResistingForce()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement3Dmini::getResistingForceIncInertia()
{
    // resize p
    P.resize(ndf);
    P.Zero();

    int numnodes = ntags.Size()/2;
    int ndm = body.Size();
    
    // get velocity, acceleration
    Vector v(ndf), vdot(ndf);
    for(int a=0; a<numnodes; a++) {

        const Vector& vel = nodes[2*a]->getTrialVel();
    	const Vector& accel = nodes[2*a]->getTrialAccel();
    	for(int i=0; i<ndm; i++) {
    	    v(vdof(ndm*a+i)) = vel(i);
    	    vdot(vdof(ndm*a+i)) = accel(i);
    	}
    	const Vector& pressure = nodes[2*a+1]->getTrialVel();
    	v(pdof(a)) = pressure(0);
    }

    // M*vdot+K*v
    P.addMatrixVector(0.0, getMass(), vdot, 1.0);
    bool l = lumped;
    lumped = false;
    P.addMatrixVector(1.0, getDamp(), v, 1.0);
    lumped = l;

    // invMb
    Matrix invMb(ndm,ndm), Mb(ndm,ndm);
    double mb = 16.*rho*Jn/(315.*ops_Dt);
    double kb = 1024.*mu/(2835.*J);
    for(int i=0; i<ndm; i++) {
    	Mb(i,i) = mb;
    }
    if(mu > 0) {
	for(int i=0; i<ndm; i++) {
	    for(int l=0; l<ndm; l++) {
		Mb(i,i) += kb*kbb(l,l);
	    }
	    for(int j=0; j<ndm; j++) {
		Mb(i,j) += kb*kbb(j,i);
		Mb(i,j) -= 2.0*kb*kbb(i,j)/3.0;
	    }
	}
    }
    Mb.Invert(invMb);

    // Gb
    Matrix Gb(ndm,numnodes);
    double gb = -16./315.0;
    for(int i=0; i<ndm; i++) {
    	for(int b=0; b<numnodes; b++) {
	    Gb(i,b) = dNdx(i,b)*gb;
    	}
    }

    // F-M*vdot-K*v
    Vector Fb(ndm), Fp(numnodes);
    double f = rho*Jn/24.;
    double fb = 16.*rho*Jn/315.;
    for(int i=0; i<ndm; i++) {
    	Fb(i) = fb*body(i);
    }
    Fp.addMatrixTransposeVector(0.0, Gb, invMb*Fb, -1.0);
    
    for(int a=0; a<numnodes; a++) {
    	for(int i=0; i<ndm; i++) {
    	    P(vdof(ndm*a+i)) -= f*body(i);
    	}
    	P(pdof(a)) -= Fp(a);
    }
    
    return P;
}


const char*
PFEMElement3Dmini::getClassType()const
{
    return "PFEMElement3Dmini";
}

int
PFEMElement3Dmini::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PFEMElement3Dmini::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement3Dmini::setDomain(Domain *theDomain)
{
    this->DomainComponent::setDomain(theDomain);
    if(theDomain == 0) {return;}

    int numnodes = ntags.Size()/2;
    
    Vector x(numnodes), y(numnodes), z(numnodes);

    ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<numnodes; i++) {

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement3Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
	const Vector& coord = nodes[2*i]->getCrds();
	int ndm = coord.Size();
	if(ndm != 3) {
	    opserr<<"WARNING: node "<<ntags(2*i)<<" ndm != 3 ";
            opserr<<"in PFEMElement3Dmini - setDomain() "<<eletag<<"\n ";
            return;
	}
        if(nodes[2*i]->getNumberDOF() < ndm) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" ndf < ndm ";
            opserr<<"in PFEMElement3Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
	
	for(int j=0; j<ndm; j++) {
	    vdof(ndm*i+j) = ndf+j;
	}
        ndf += nodes[2*i]->getNumberDOF();

	// material coordinates
	const Vector& disp = nodes[2*i]->getDisp();
	x(i) = coord(0)+disp(0);
	y(i) = coord(1)+disp(1);
	z(i) = coord(2)+disp(2);

        // get pc
        int pndf = 1;
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] != 0) {
            thePCs[i]->setDomain(theDomain);
        } else {
            thePCs[i] = new Pressure_Constraint(ntags(2*i), pndf);
            if(thePCs[i] == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMElement3Dmini::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement3Dmini::setDomain "<<eletag<<"\n";
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
            opserr<<"in PFEMElement3Dmini - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
	pdof(i) = ndf;
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    // Jacobian
    Jn = det(x(1)-x(0),y(1)-y(0),z(1)-z(0),
	     x(2)-x(0),y(2)-y(0),z(2)-z(0),
	     x(3)-x(0),y(3)-y(0),z(3)-z(0));

}

void
PFEMElement3Dmini::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement3Dmini: "<<this->getTag()<<endln;
}

int
PFEMElement3Dmini::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    return 0;
}
