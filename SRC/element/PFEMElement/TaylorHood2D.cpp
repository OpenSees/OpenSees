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
// Description: P2+P-1 element, quadratic velocity, discontinuous linear pressure
//

// for FEM_ObjectBroker, recvSelf must invoke

#include "TaylorHood2D.h"
#include <elementAPI.h>
#include "TriGaussPoints.h"
#include <Node.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <math.h>
#include "HigherOrder.h"

Matrix TaylorHood2D::mat;
Vector TaylorHood2D::vec;

void* OPS_TaylorHood2D()
{
    int num = OPS_GetNumRemainingInputArgs();
    if(num < 4) {
	opserr<<"WARNING: insufficient number of arguments -- TaylorHood2D:";
	opserr<<"tag, nd1, nd2, nd3, (rho, mu, b1, b2, thk, kappa)\n";
	return 0;
    }

    // tag, nd1, nd2, nd3
    num = 4;
    int idata[4];
    if(OPS_GetIntInput(&num,idata)<0) {
	opserr << "WARNING: failed to read integers -- TaylorHood2D\n";
	return 0;
    }

    // (rho, mu, b1, b2, thinkness, kappa) 
    num = OPS_GetNumRemainingInputArgs();
    if(num > 6) num = 6;
    double data[6] = {1000.0,1e-3,0,-9.81,1.0,2.15e9};
    if(OPS_GetDoubleInput(&num,data) < 0) {
	opserr << "WARNING: failed to read doubles -- TaylorHood2D\n";
	return 0;
    }

    return new TaylorHood2D(idata[0],idata[1],idata[2],idata[3],
			    data[0],data[1],data[2],data[3],
			    data[4],data[5]);
}


TaylorHood2D::TaylorHood2D()
    :Element(0, ELE_TAG_TaylorHood2D),
     ntags(9), nodes(9), pcs(3), rho(0), mu(0), b1(0), b2(0),
     thk(0), kappa(0), ndf(0), vxdof(6), vydof(6),
     pdof(3)
{
}

TaylorHood2D::TaylorHood2D(int tag, int nd1, int nd2, int nd3,
			   double r, double m,
			   double bx, double by, 
			   double t, double ka)
    :Element(tag, ELE_TAG_TaylorHood2D),
     ntags(9), nodes(9), pcs(3), rho(r), mu(m), b1(bx), b2(by),
     thk(t), kappa(0.0), ndf(0), vxdof(6), vydof(6),
     pdof(3)
{
    ntags[0] = nd1;
    ntags[1] = nd2;
    ntags[2] = nd3;
    ntags[3] = nd3;
    ntags[4] = nd3;
    ntags[5] = nd3;
    ntags[6] = nd3;
    ntags[7] = nd3;
    ntags[8] = nd3;

    if (ka <= 0) {
	ka = 2.15e9;
    }
    kappa = sqrt(ka);
}

// destructor
TaylorHood2D::~TaylorHood2D()
{
    HigherOrder& ho = OPS_getHigherOrder();
    int eletag = this->getTag();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;
    
    for (int i=0; i<vxdof.Size()/2; ++i) {
	
	// opposite edge
	int j = i + 1;
	int k = j + 1;
	if (j > 2) j-=3;
	if (k > 2) k-= 3;

	// check if edge node exists
	VInt edge(2);
	edge[0] = ntags(j);
	edge[1] = ntags(k);

	// remove element from ho
	// and check if no element is
	// associated with the edge
	bool empty = ho.removeEle(edge,eletag);
	if (empty) {
	    // delete mid node
	    domain->removeNode(ntags(i+3));
	    delete nodes[i+3];
	}
    }
}

// methods dealing with nodes and number of external dof
int
TaylorHood2D::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID &
TaylorHood2D::getExternalNodes()
{
    return ntags;
}

Node **
TaylorHood2D::getNodePtrs()
{
    return &nodes[0];
}

int
TaylorHood2D::getNumDOF()
{
    return ndf;
}

// public methods to set the state of the element
int
TaylorHood2D::revertToLastCommit()
{
    return 0;
}

int
TaylorHood2D::update()
{
    return 0;
}

int
TaylorHood2D::commitState()
{
    return Element::commitState();
}

// public methods to obtain stiffness, mass, damping and residual information
const Matrix &
TaylorHood2D::getTangentStiff()
{
    mat.resize(ndf,ndf);
    mat.Zero();
    return mat;
}

const Matrix &
TaylorHood2D::getInitialStiff()
{
    mat.resize(ndf,ndf);
    mat.Zero();
    return mat;
}

const Matrix &
TaylorHood2D::getMass()
{
    mat.resize(ndf,ndf);
    mat.Zero();

    // current coordinates
    Matrix X;
    if (getX(X) < 0) {
	opserr << "WARING: failed to get current coordinate ";
	opserr << "-- TaylorHood2D::getMass\n";
	return mat;
    }

    // gauss points
    TriGaussPoints gauss;
    VDouble xpts, ypts, wts;
    int np = 6;
    gauss(np,xpts,ypts,wts);

    // integration
    for(unsigned int i=0; i<wts.size(); ++i) {

	// area coordinates
	double l2 = xpts[i];
	double l3 = ypts[i];
	double l1 = 1.0 - l2 - l3;

	// shape functions
	Vector Nv, Np;
	Matrix dNv;
	vshape(Nv,l1,l2,l3);
	pshape(Np,l1,l2,l3);
	derishape(dNv,l1,l2,l3);

	// Jacobian matrix
	Matrix J;
	Jacobian(J,X,dNv);
	double detJ = J(0,0)*J(1,1)-J(1,0)*J(0,1);
	
	// M
	for (int a=0; a<Nv.Size(); ++a) {
	    for (int b=0; b<Nv.Size(); ++b) {

		double m = wts[i]*Mab(a,b,detJ,Nv);
		mat(vxdof(a), vxdof(b)) += m;
		mat(vydof(a), vydof(b)) += m;
	    }
	}

	// Mp must be lumped
	for (int a=0; a<Np.Size(); ++a) {
	    for (int b=0; b<Np.Size(); ++b) {

		double m = wts[i]*Mpab(a,b,detJ,Np);
		mat(pdof(a), pdof(a)) += m;
	    }
	}
    }

    mat *= 0.5;
    
    return mat;
}


const Matrix &
TaylorHood2D::getDamp()
{
    mat.resize(ndf,ndf);
    mat.Zero();

    // current coordinates
    Matrix X;
    if (getX(X) < 0) {
	opserr << "WARING: failed to get current coordinate ";
	opserr << "-- TaylorHood2D::getDamp\n";
	return mat;
    }

    // gauss points
    TriGaussPoints gauss;
    VDouble xpts, ypts, wts;
    int np = 8;
    gauss(np,xpts,ypts,wts);

    // integration
    for(unsigned int i=0; i<wts.size(); ++i) {

	// area coordinates
	double l2 = xpts[i];
	double l3 = ypts[i];
	double l1 = 1.0 - l2 - l3;

	// shape functions
	Vector Nv, Np;
	Matrix dNv;
	vshape(Nv,l1,l2,l3);
	pshape(Np,l1,l2,l3);
	derishape(dNv,l1,l2,l3);

	// Jacobian matrix
	Matrix J;
	Jacobian(J,X,dNv);
	double detJ = J(0,0)*J(1,1)-J(1,0)*J(0,1);
	Matrix invJ = J;
	invJ(0,0) = J(1,1);
	invJ(1,1) = J(0,0);
	invJ(0,1) = -J(0,1);
	invJ(1,0) = -J(1,0);
	if (fabs(detJ) < 1e-15) {
	    opserr << "WARNING: detJ < 1e-15";
	    opserr << "-- TaylorHood2D::getDamp\n";
	}
	invJ /= detJ;

	// dNdx and dNdy
	Matrix dNdx = invJ*dNv;

	// K and -G and Gt
	for (int a=0; a<Nv.Size(); ++a) {

	    // K
	    for (int b=0; b<Nv.Size(); ++b) {

		double k11 = wts[i]*K11(a,b,detJ,dNdx);
		double k12 = wts[i]*K12(a,b,detJ,dNdx);
		double k21 = wts[i]*K21(a,b,detJ,dNdx);
		double k22 = wts[i]*K22(a,b,detJ,dNdx);
		
		mat(vxdof(a), vxdof(b)) += k11;
		mat(vxdof(a), vydof(b)) += k12;
		mat(vydof(a), vxdof(b)) += k21;
		mat(vydof(a), vydof(b)) += k22;
	    }

	    // -G and Gt
	    for (int b=0; b<Np.Size(); ++b) {

		double g1 = wts[i]*Gab(a,b,detJ,dNdx,0, Np);
		double g2 = wts[i]*Gab(a,b,detJ,dNdx,1, Np);
		
		mat(vxdof(a), pdof(b)) -= g1;
		mat(vydof(a), pdof(b)) -= g2;

		mat(pdof(b), vxdof(a)) += kappa*g1;
		mat(pdof(b), vydof(a)) += kappa*g2;
	    }
	}

    }

    mat *= 0.5;

    
    return mat;
}

// methods for applying loads
int
TaylorHood2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

// methods for obtaining resisting force (force includes elemental loads)
const Vector &
TaylorHood2D::getResistingForce()
{
    vec.resize(ndf);
    vec.Zero();
    return vec;
}

const Vector &
TaylorHood2D::getResistingForceIncInertia()
{
    vec.resize(ndf);
    vec.Zero();

    // current coordinates
    Matrix X;
    if (getX(X) < 0) {
	opserr << "WARING: failed to get current coordinate ";
	opserr << "-- TaylorHood2D::getResistingForceIncInertia\n";
	return vec;
    }

    // gauss points
    TriGaussPoints gauss;
    VDouble xpts, ypts, wts;
    int np = 3;
    gauss(np,xpts,ypts,wts);

    // integration
    for(unsigned int i=0; i<wts.size(); ++i) {

	// area coordinates
	double l2 = xpts[i];
	double l3 = ypts[i];
	double l1 = 1.0 - l2 - l3;

	// shape functions
	Vector Nv;
	Matrix dNv;
	vshape(Nv,l1,l2,l3);
	derishape(dNv,l1,l2,l3);

	// Jacobian matrix
	Matrix J;
	Jacobian(J,X,dNv);
	double detJ = J(0,0)*J(1,1)-J(1,0)*J(0,1);

	// F
	for (int a=0; a<Nv.Size(); ++a) {

	    double f1 = wts[i]*Fa(a,detJ,Nv,b1);
	    double f2 = wts[i]*Fa(a,detJ,Nv,b2);

	    vec(vxdof(a)) += f1;
	    vec(vydof(a)) += f2;
	}

    }
    
    vec *= -0.5;

    // vdot, v
    Vector vdot(ndf), v(ndf);
    for(int a=0; a<vxdof.Size(); a++) {
	const Vector& vel = nodes[a]->getTrialVel();
	const Vector& accel = nodes[a]->getTrialAccel();
	
	v(vxdof(a)) = vel(0);
	v(vydof(a)) = vel(1);

	vdot(vxdof(a)) = accel(0);
	vdot(vydof(a)) = accel(1);
    }

    for(int a=0; a<pdof.Size(); a++) {
	const Vector& p = nodes[a+6]->getTrialVel();
	const Vector& pdot = nodes[a+6]->getTrialAccel();
	    
	v(pdof(a)) = p(0);
	vdot(pdof(a)) = pdot(0);
    }

    // -r = M*vdot+K*v-F
    vec.addMatrixVector(1.0, this->getMass(), vdot, 1.0);
    vec.addMatrixVector(1.0, this->getDamp(), v, 1.0);
    
    return vec;
}

// MovableObject
const char *
TaylorHood2D::getClassType() const
{
    return "TaylorHood2D";
}

int
TaylorHood2D::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
TaylorHood2D::recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
    return 0;
}

// DomainComponent
void
TaylorHood2D::setDomain(Domain *theDomain)
{

    // set domain
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) return;

    // ndf
    ndf = 0;
    int eletag = this->getTag();

    // get corner nodes
    for (int i=0; i<vxdof.Size()/2; ++i) {
	
	// get corner node
    	nodes[i] = theDomain->getNode(ntags(i));
    	if (nodes[i] == 0) {
    	    opserr << "WARNING: node "<<ntags(i);
    	    opserr << "does not exist -- TaylorHood2D\n";
    	    return;
    	}
    	if (nodes[i]->getNumberDOF() < 2) {
    	    opserr << "WARNING: node "<<ntags(i);
    	    opserr << "ndf < 2 -- TaylorHood2D\n";
    	    return;
    	}
    	if (nodes[i]->getCrds().Size() < 2) {
    	    opserr << "WARNING: node "<<ntags(i);
    	    opserr << "ndm < 2 -- TaylorHood2D\n";
    	    return;
    	}

	// ndf
	vxdof(i) = ndf;
	vydof(i) = ndf+1;
	ndf += nodes[i]->getNumberDOF();
    }

    // get mid nodes
    HigherOrder& ho = OPS_getHigherOrder();
    //BackgroundUtil util;
    for (int i=0; i<vxdof.Size()/2; ++i) {
	
	// opposite edge
	int j = i + 1;
	int k = j + 1;
	if (j > 2) j-=3;
	if (k > 2) k-= 3;

	// check if edge node exists
	VInt edge(2);
	edge[0] = ntags(j);
	edge[1] = ntags(k);
	VInt& mid = ho(edge);
	if (mid.size() == 1) {
	    nodes[i+3] = theDomain->getNode(mid[0]);
	    if (nodes[i+3] == 0) {
		opserr << "WARNING: node "<<mid[0];
		opserr << "does not exist -- TaylorHood2D\n";
		return;
	    }
	    if (nodes[i+3]->getNumberDOF() < 2) {
		opserr << "WARNING: node "<<mid[0];
		opserr << "ndf < 2 -- TaylorHood2D\n";
		return;
	    }
	    if (nodes[i+3]->getCrds().Size() < 2) {
		opserr << "WARNING: node "<<mid[0];
		opserr << "ndm < 2 -- TaylorHood2D\n";
		return;
	    }
	    ntags[i+3] = nodes[i+3]->getTag();

	} else if (mid.size() > 1) {
	    opserr << "WARNING: edge ("<<edge[0]<<" , "<<edge[1];
	    opserr << ") is larger than 2nd order\n";
	    return;
	} else {
	    nodes[i+3] = 0;
	}
	
	// new edge node tag
	if (nodes[i+3] == 0) {
	    //ntags[i+3] = util.findNodeTag();

	    // edge node states
	    Vector ncrds(2), ndisp(2), nvel(2), naccel(2);

	    // corner node states
	    const Vector& crdsj = nodes[j]->getCrds();
	    const Vector& dispj = nodes[j]->getTrialDisp();
	    const Vector& velj = nodes[j]->getTrialVel();
	    const Vector& accelj = nodes[j]->getTrialAccel();

	    const Vector& crdsk = nodes[k]->getCrds();
	    const Vector& dispk = nodes[k]->getTrialDisp();
	    const Vector& velk = nodes[k]->getTrialVel();
	    const Vector& accelk = nodes[k]->getTrialAccel();

	    for (int m=0; m<2; ++m) {
		ncrds(m) = (crdsj(m)+crdsk(m))/2.0;
		ndisp(m) = (dispj(m)+dispk(m))/2.0;
		nvel(m) = (velj(m)+velk(m))/2.0;
		naccel(m) = (accelj(m)+accelk(m))/2.0;
	    }

	    // create the edge node
	    nodes[i+3] = new Node(ntags[i+3], 2, ncrds[0], ncrds[1]);
	    if (nodes[i+3] == 0) {
		opserr << "WARNING: run out of memory -- TaylorHood2D\n";
		return;
	    }
	    if (theDomain->addNode(nodes[i+3]) == false) {
		opserr << "WARNING: failed to add node to domain";
		opserr << " -- TaylorHood2D\n";
		delete nodes[i+3];
		return;
	    }
	    nodes[i+3]->setTrialDisp(ndisp);
	    nodes[i+3]->setTrialVel(nvel);
	    nodes[i+3]->setTrialAccel(naccel);
	    nodes[i+3]->commitState();

	    // add the edge node to ho
	    mid.push_back(ntags[i+3]);
	}

	// ndf
	vxdof(i+3) = ndf;
	vydof(i+3) = ndf+1;
	ndf += nodes[i+3]->getNumberDOF();

	// add element to ho
	ho.addEle(edge, eletag);
    }

    // get pc
    for (int i=0; i<pdof.Size(); ++i) {
	int pndf = 1;
	pcs[i] = theDomain->getPressure_Constraint(ntags(i));
	if(pcs[i] != 0) {
	    pcs[i]->setDomain(theDomain);
	} else {
	    pcs[i] = new Pressure_Constraint(ntags(i), pndf);
	    if(pcs[i] == 0) {
		opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
		opserr<<"TaylorHood2D::setDomain "<<eletag<<"\n";
		return;
	    }
	    if(theDomain->addPressure_Constraint(pcs[i]) == false) {
		opserr<<"WARNING: failed to add Pressure_Constraint to domain";
		opserr<<" -- TaylorHood2D::setDomain "<<eletag<<"\n";
		delete pcs[i];
		pcs[i] = 0;
		return;
	    }
	}

	// connect
	pcs[i]->connect(eletag);

	// get pressure node
	nodes[i+6] = pcs[i]->getPressureNode();
	if(nodes[i+6] == 0) {
	    opserr<<"WARNING: pressure node does not exist ";
	    opserr<<"in TaylorHood2D - setDomain() "<<eletag<<"\n ";
	    return;
	}
	ntags(i+6) = nodes[i+6]->getTag();
	pdof(i) = ndf;
	ndf += nodes[i+6]->getNumberDOF();
    }

}

// TaggedObject
void
TaylorHood2D::Print(OPS_Stream &s, int flag)
{
    s << this->getClassType() <<"\n";
    s << "tag : " << this->getTag() << "\n";
    s << "nodes : ";
    for (int i=0; i<ntags.Size(); ++i) {
	s << ntags(i)<<" ";
    }
    s << "\n";
    s << "rho : " << rho << "\n";
    s << "mu : " << mu << "\n";
    s << "b1 : " << b1 << "\n";
    s << "b2 : " << b2 << "\n";
    s << "thickness : " << thk << "\n";
    s << "kappa : " << kappa*kappa << "\n";
}

int
TaylorHood2D::displaySelf(Renderer &, int mode, float fact,
			  const char **displayModes,
			  int numModes)
{
    return 0;
}

// get current coordinates
int
TaylorHood2D::getX(Matrix& x)
{
    x.resize(vxdof.Size(), 2);

    for (int i=0; i<x.noRows(); ++i) {
	if (nodes[i] == 0) return -1;
	const Vector& crds = nodes[i]->getCrds();
	const Vector& disp = nodes[i]->getTrialDisp();
	if (crds.Size() < 2 || disp.Size() < 2) return -1;
	x(i,0) = crds(0)+disp(0);
	x(i,1) = crds(1)+disp(1);
    }

    return 0;
}

// Jacobian matrix
inline void
TaylorHood2D::Jacobian(Matrix& J, const Matrix& X,
		       const Matrix& dNdL)
{
    J = dNdL*X;
}

// shape functions
void
TaylorHood2D::vshape(Vector& N, double l1, double l2, double l3)
{
    N.resize(vxdof.Size());

    N(0) = l1*(2*l1-1);
    N(1) = l2*(2*l2-1);
    N(2) = l3*(2*l3-1);
    N(3) = 4*l2*l3;
    N(4) = 4*l1*l3;
    N(5) = 4*l1*l2;
}

void
TaylorHood2D::pshape(Vector& N, double l1, double l2, double l3)
{
    N.resize(pdof.Size());

    N(0) = l1;
    N(1) = l2;
    N(2) = l3;
}

void
TaylorHood2D::derishape(Matrix& dN, double l1, double l2, double l3)
{
    dN.resize(2,vxdof.Size());

    dN(0,0) = 1-4*l1;
    dN(0,1) = 4*l2-1;
    dN(0,2) = 0.0;
    dN(0,3) = 4*l3;
    dN(0,4) = -4*l3;
    dN(0,5) = 4*l1 - 4*l2;

    dN(1,0) = 1-4*l1;
    dN(1,1) = 0.0;
    dN(1,2) = 4*l3-1;
    dN(1,3) = 4*l2;
    dN(1,4) = 4*l1-4*l3;
    dN(1,5) = -4*l2;
}

// matrices
inline double
TaylorHood2D::Mab(int a, int b, double detJ, const Vector& N)
{
    return rho*detJ*N(a)*N(b);
}

inline double
TaylorHood2D::K11(int a, int b, double detJ,
		  const Matrix& dNdx)
{
    return mu*detJ*(4*dNdx(0,a)*dNdx(0,b)+3*dNdx(1,a)*dNdx(1,b))/3.0;
}

inline double
TaylorHood2D::K12(int a, int b, double detJ,
		  const Matrix& dNdx)
{
    return mu*detJ*(3*dNdx(1,a)*dNdx(0,b)-2*dNdx(0,a)*dNdx(1,b))/3.0;
}

inline double
TaylorHood2D::K21(int a, int b, double detJ,
		  const Matrix& dNdx)
{
    return mu*detJ*(3*dNdx(0,a)*dNdx(1,b)-2*dNdx(1,a)*dNdx(0,b))/3.0;
}

inline double
TaylorHood2D::K22(int a, int b, double detJ,
		  const Matrix& dNdx)
{
    return mu*detJ*(4*dNdx(1,a)*dNdx(1,b)+3*dNdx(0,a)*dNdx(0,b))/3.0;
}

inline double
TaylorHood2D::Gab(int a, int b, double detJ,
		  const Matrix& dNdx, int i,
		  const Vector& Np)
{
    return detJ*dNdx(i,a)*Np(b);
}

inline double
TaylorHood2D::Fa(int a, double detJ, const Vector& N,
		 double bx)
{
    return rho*detJ*bx*N(a);
}

inline double
TaylorHood2D::Mpab(int a, int b, double detJ, const Vector& Np)
{
    return detJ*Np(a)*Np(b)/kappa;
}

