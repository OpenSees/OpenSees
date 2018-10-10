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
#include <map>
#include <cmath>

Matrix PFEMElement2Dmini::K;
Vector PFEMElement2Dmini::P;
bool PFEMElement2Dmini::dispon = true;

void* OPS_PFEMElement2Dmini(const ID &info)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
    	opserr << "WARNING: domain is not created\n";
    	return 0;
    }

    int idata[5];
    double data[6] = {0,0,0,0,1.0,-1};
    int numdata;

    // regular element, not in a mesh, get tags
    if (info.Size() == 0) {
	numdata = OPS_GetNumRemainingInputArgs();
	if(numdata < 5) {
	    opserr<<"WARNING: insufficient number of arguments: tag, nd1, nd2, nd3, nd4\n";
	    return 0;
	}

	// tag, nd1, nd2, nd3, nd4
	numdata = 5;
	if(OPS_GetIntInput(&numdata,idata)<0) {
	    opserr << "WARNING: failed to get tags\n";
	    return 0;
	}
    }

    // regular element, or save data
    if (info.Size()==0 || info(0)==1) {
	if(OPS_GetNumRemainingInputArgs() < 4) {
	    opserr<<"insufficient arguments: rho, mu, b1, b2, (thinknes,kappa)\n";
	    return 0;
	}

	// rho, mu, b1, b2, (thinknes,kappa)
	numdata = OPS_GetNumRemainingInputArgs();
	if(numdata > 6) numdata = 6;
	if(OPS_GetDoubleInput(&numdata,data) < 0) {
	    opserr << "WARNING: failed to get fluid properties\n";
	    return 0;
	}
    }

    // save/load data for different mesh
    static std::map<int, Vector> meshdata;
    if (info.Size()>0 && info(0)==1) {
	if (info.Size() < 2) {
	    opserr << "WARNING: need info -- inmesh, meshtag\n";
	    return 0;
	}

	// save the data for a mesh
	Vector& mdata = meshdata[info(1)];
	mdata.resize(6);
	for (int i=0; i<6; ++i) {
	    mdata(i) = data[i];
	}
	return &meshdata;

    } else if (info.Size()>0 && info(0)==2) {
	if (info.Size() < 7) {
	    opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2, nd3, nd4\n";
	    return 0;
	}

	// get the data for a mesh
	Vector& mdata = meshdata[info(1)];
	if (mdata.Size() < 6) return 0;

	idata[0] = info(2);
	for (int i=0; i<4; ++i) {
	    idata[i+1] = info(3+i);
	}

	for (int i=0; i<6; ++i) {
	    data[i] = mdata(i);
	}

    }

    return new PFEMElement2Dmini(idata[0],idata[1],idata[2],idata[3],idata[4],
				 data[0],data[1],data[2],data[3],data[4],
				 data[5]);
}



// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2Dmini::PFEMElement2Dmini()
    :Element(0, ELE_TAG_PFEMElement2Dmini), ntags(7), nodes(7,static_cast<Node*>(0)),
     thePCs(4,static_cast<Pressure_Constraint*>(0)), rho(0), mu(0), bx(0.0), by(0.0),
     thk(1.0),ka(-1),J(0.0),bb(3),cc(3),vxdof(4),vydof(4),pdof(3),ndf(0),bnode(0)
{
}

// for object
PFEMElement2Dmini::PFEMElement2Dmini(int tag, int nd1, int nd2, int nd3, int nd4,
				     double r, double m, double b1,
				     double b2, double t, double k)
    :Element(tag, ELE_TAG_PFEMElement2Dmini), ntags(7), nodes(7,static_cast<Node*>(0)),
     thePCs(4,static_cast<Pressure_Constraint*>(0)), rho(r), mu(m), bx(b1), by(b2),
     thk(t),ka(k),J(0.0),bb(3),cc(3),vxdof(4),vydof(4),pdof(3),ndf(0),bnode(nd4)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    ntags(1)=nd1; ntags(3)=nd2; ntags(5)=nd3;
    ntags(6)=nd1;
}


PFEMElement2Dmini::~PFEMElement2Dmini()
{
    for(int i=0; i<(int)thePCs.size(); i++) {
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
	if (thePCs[3] != 0) {
	    domain->removePressure_Constraint(thePCs[3]->getTag());
	    delete thePCs[3];
	    thePCs[3] = 0;
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
    if (!dispon) {
	if (updateJacobian() < 0) return -1;
    }
    return Element::commitState();
}

int
PFEMElement2Dmini::update()
{
    if (dispon) {
	return updateJacobian();
    }

    return 0;
}

int
PFEMElement2Dmini::updateJacobian()
{
    // get nodal coordinates
    Vector x(3), y(3);
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
        const Vector& disp = nodes[2*i]->getTrialDisp();
        x(i) = coord(0) + disp(0);
        y(i) = coord(1) + disp(1);
    }

    // get c and d
    bb(0) = y(1)-y(2);
    cc(0) = x(2)-x(1);
    bb(1) = y(2)-y(0);
    cc(1) = x(0)-x(2);
    bb(2) = y(0)-y(1);
    cc(2) = x(1)-x(0);

    // get Jacobi
    J = bb(0)*cc(1)-cc(0)*bb(1);

    if(fabs(J)<=1e-15) {
    	opserr<<"WARNING: element "<<this->getTag()<<" Jacobian determinant ";
    	opserr<<J<<" <= 0\n";
	opserr<<"J = "<<J<<"\n";
	opserr<<"x = "<<x;
	opserr<<"y = "<<y;
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

    return K;
}

const Matrix&
PFEMElement2Dmini::getDamp()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // get matrices
    Matrix M1, M2, G, L;
    getM(M1);
    getK(M2);
    getG(G);
    getL(L);

    // set matrix
    for(int a=0; a<vxdof.Size(); a++) {
	for(int b=0; b<vxdof.Size(); b++) {

	    // K
	    K(vxdof(a),vxdof(b)) = M1(2*a,2*b)/ops_Dt + M2(2*a,2*b);
	    K(vxdof(a),vydof(b)) = M1(2*a,2*b+1)/ops_Dt + M2(2*a,2*b+1);
	    K(vydof(a),vxdof(b)) = M1(2*a+1,2*b)/ops_Dt + M2(2*a+1,2*b);
	    K(vydof(a),vydof(b)) = M1(2*a+1,2*b+1)/ops_Dt + M2(2*a+1,2*b+1);
	}

	for(int b=0; b<pdof.Size(); b++) {
	    // G
	    K(vxdof(a),pdof(b)) = G(2*a,b);
	    K(vydof(a),pdof(b)) = G(2*a+1,b);

	    // Gt
	    K(pdof(b),vxdof(a)) = G(2*a,b);
	    K(pdof(b),vydof(a)) = G(2*a+1,b);
	}
    }

    // L
    for(int a=0; a<pdof.Size(); a++) {
	for(int b=0; b<pdof.Size(); b++) {
	    // L
	    K(pdof(a),pdof(b)) = ops_Dt*L(a,b)/rho;
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
PFEMElement2Dmini::addInertiaLoadToUnbalance(const Vector &abbel)
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

    // get velocity and pressure
    Vector vk(vxdof.Size()*2), vn(vxdof.Size()*2), pk(pdof.Size());
    for(int a=0; a<vxdof.Size(); a++) {

        const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& veln = nodes[2*a]->getVel();
	vk(2*a) = vel(0);
	vk(2*a+1) = vel(1);
	vn(2*a) = veln(0);
	vn(2*a+1) = veln(1);
    }
    for (int a=0; a<pdof.Size(); ++a) {
	const Vector& pressure = nodes[2*a+1]->getTrialVel();
	pk(a) = pressure(0);
    }

    // get rhs
    Vector F, Fp;
    getF(F);
    getFp(Fp);

    // get matrices
    Matrix M1, M2, G;
    getM(M1);
    getK(M2);
    getG(G);

    // add forces
    F.addMatrixVector(1.0, M1, vk, -1.0/ops_Dt);
    F.addMatrixVector(1.0, M1, vn, 1.0/ops_Dt);
    F.addMatrixVector(1.0, G, pk, -1.0);
    F.addMatrixVector(1.0, M2, vk, -1.0);

    Fp.addMatrixTransposeVector(-1.0, G, vk, 1.0);

    // add to vector
    for(int a=0; a<vxdof.Size(); a++) {
	P(vxdof(a)) = -F(2*a);
	P(vydof(a)) = -F(2*a+1);
    }
    for(int a=0; a<pdof.Size(); a++) {
	P(pdof(a)) = -Fp(a);
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

    ndf = 0;
    int eletag = this->getTag();
    Vector bcrds(2);
    for(int i=0; i<pdof.Size(); i++) {

	// set velocity ndf
	vxdof(i) = ndf;
	vydof(i) = ndf+1;

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
        ndf += nodes[2*i]->getNumberDOF();

	// get crds
	const Vector& crds = nodes[2*i]->getCrds();
	for (int j=0; j<2; j++) {
	    bcrds(j) += crds(j);
	}

	// set pressure dof
	pdof(i) = ndf;

        // get pc
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] == 0) {
	    opserr << "WARNING: failed to get PC -- PFEMElement2Dmini\n";
	    return;
        }
	thePCs[i]->setDomain(theDomain);

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
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    // create bubble node
    ntags[6] = bnode;
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

    // create pc for bubble node
    thePCs[3] = theDomain->getPressure_Constraint(ntags[6]);
    if (thePCs[3] != 0) {
	opserr<<"WARNING: pc for bubble node already exists\n";
	return;
    }
    thePCs[3] = new Pressure_Constraint(ntags[6], 0.0);
    if (thePCs[3] == 0) {
	opserr<<"WARNING: failed to create PC for buble node\n";
	return;
    }
    if (theDomain->addPressure_Constraint(thePCs[3]) == false) {
	opserr<<"WARNING: failed to add PC to domain\n";
	delete thePCs[3];
	return;
    }
    thePCs[3]->setDomain(theDomain);
    thePCs[3]->connect(eletag);

    // Jacobian
    if (!dispon) {
	updateJacobian();
    }
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

void
PFEMElement2Dmini::getM(Matrix& M)
{
    M.resize(vxdof.Size()*2,vxdof.Size()*2);
    M.Zero();
    double m = rho*thk*J/6.0;
    double mb = 27*rho*thk*J/120.0;
    for (int a=0; a<vxdof.Size()-1; ++a) {
	M(2*a,2*a) = m;
	M(2*a+1,2*a+1) = m;
    }
    M(6,6) = mb;
    M(7,7) = mb;
}

void
PFEMElement2Dmini::getK(Matrix& K)
{
    K.resize(vxdof.Size()*2,vxdof.Size()*2);
    K.Zero();
}

void
PFEMElement2Dmini::getG(Matrix& G)
{
    G.resize(vxdof.Size()*2,pdof.Size());
    G.Zero();
    for (int b=0; b<pdof.Size(); ++b) {
	for (int a=0; a<vxdof.Size()-1; ++a) {
	    G(2*a,b) = bb(b)*thk/6.0;
	    G(2*a+1,b) = cc(b)*thk/6.0;
	}
	G(6,b) = bb(b)*thk*27/120.0;
	G(7,b) = cc(b)*thk*27/120.0;
    }
}

void
PFEMElement2Dmini::getL(Matrix& L)
{
    // Laplace matrix
    L.resize(pdof.Size(),pdof.Size());
    L.Zero();

    for (int a=0; a<pdof.Size(); ++a) {
	for (int b=0; b<pdof.Size(); ++b) {
	    L(a,b) = thk*(bb(a)*bb(b)+cc(a)*cc(b))/(2*J);
	}
    }
}

void
PFEMElement2Dmini::getF(Vector& F)
{
    F.resize(vxdof.Size()*2);
    for (int a=0; a<vxdof.Size()-1; ++a) {
	F(2*a) = rho*J*thk*bx/6.0;
	F(2*a+1) = rho*J*thk*by/6.0;
    }
    F(6) = rho*J*thk*bx*27/120.0;
    F(7) = rho*J*thk*by*27/120.0;
}

void
PFEMElement2Dmini::getFp(Vector& Fp)
{
    Fp.resize(pdof.Size());
    Fp.Zero();
}
