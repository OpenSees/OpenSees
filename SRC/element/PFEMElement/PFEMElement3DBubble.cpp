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

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: This file contains the class definition for PFEMElement3DBubble.

#include "PFEMElement3DBubble.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Renderer.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <cmath>
#include <ElementIter.h>
#include <iostream>
#include <map>

Matrix PFEMElement3DBubble::K;
Vector PFEMElement3DBubble::P;

bool PFEMElement3DBubble::dispon = true;

void* OPS_PFEMElement3DBubble(const ID &info)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
    	opserr << "WARNING: domain is not created\n";
    	return 0;
    }

    int idata[5];
    double data[6] = {0,0,0,0,0,-1};
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
	if(OPS_GetNumRemainingInputArgs() < 5) {
	    opserr<<"insufficient arguments: rho, mu, b1, b2, b3, (kappa)\n";
	    return 0;
	}

	// rho, mu, b1, b2, b3, (kappa)
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

    } else if (info.Size()>0 && info(0)==3) {
        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }

        // get the data for a mesh
        Vector& mdata = meshdata[info(1)];
        return &mdata;
    }

    return new PFEMElement3DBubble(idata[0],idata[1],idata[2],idata[3],idata[4],
				   data[0],data[1],data[2],data[3],data[4],data[5]);
}

// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement3DBubble::PFEMElement3DBubble()
    :Element(0, ELE_TAG_PFEMElement3DBubble), ntags(),
     nodes(), thePCs(),
     rho(0), mu(0), bx(0), by(0), bz(0), J(0.0), numDOFs(),
     kappa(-1), parameterID(0),dNdx(),dNdy(),dNdz(),
     M(), D(), F(), Fp()
{
}

// for object
PFEMElement3DBubble::PFEMElement3DBubble(int tag, int nd1, int nd2, int nd3, int nd4,
                                         double r, double m, double b1, double b2,
                                         double b3, double ka)
    :Element(tag, ELE_TAG_PFEMElement3DBubble), ntags(8),
     nodes(8), thePCs(4),
     rho(r), mu(m), bx(b1), by(b2), bz(b3), J(0.0), numDOFs(),
     kappa(ka), parameterID(0), dNdx(4),dNdy(4),dNdz(4),
     M(), D(), F(), Fp()
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3; ntags(6)=nd4;
    ntags(1)=nd1; ntags(3)=nd2; ntags(5)=nd3; ntags(7)=nd4;
}

PFEMElement3DBubble::~PFEMElement3DBubble()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;
    for (unsigned int i=0; i<thePCs.size(); ++i) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}

int
PFEMElement3DBubble::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement3DBubble::getExternalNodes()
{
    if(numDOFs.Size() == 0) return numDOFs;
    return ntags;
}

Node **
PFEMElement3DBubble::getNodePtrs(void)
{
    return &nodes[0];
}

int
PFEMElement3DBubble::getNumDOF()
{
    if(numDOFs.Size() == 0) return 0;
    return numDOFs(numDOFs.Size()-1);
}

int
PFEMElement3DBubble::revertToLastCommit()
{
    return 0;
}

int PFEMElement3DBubble::commitState()
{
//    if (!dispon) {
//	if (updateJacobi() < 0) {
//	    opserr << "WARNING: failed to update Jacobi -- Bubble3D::commitState\n";
//	    return -1;
//	}
//    }
    return Element::commitState();
}

int
PFEMElement3DBubble::updateJacobi()
{
    Matrix Jmat(4,4), Jfact(4,4);
    Jmat(0,0) = 1.0; Jmat(0,1) = 1.0; Jmat(0,2) = 1.0; Jmat(0,3) = 1.0;
    for (int i=0; i<Jmat.noCols(); ++i) {
	const Vector& coord = nodes[2*i]->getCrds();
	const Vector& disp = nodes[2*i]->getTrialDisp();
	for (int j=0; j<coord.Size(); ++j) {
	    Jmat(j+1,i) = coord(j)+disp(j);
	}
    }
    
    cofactor(Jmat,Jfact);

    J = 0.0;
    for (int i=0; i<Jfact.noRows(); ++i) {
	J += Jfact(i,0);
    }

    // devoid element and disconnect pcs
    if (J < 1e-15) {
	opserr<<"J < 1e-15\n";
	return -1;
    }

    dNdx.resize(Jfact.noRows());
    dNdy.resize(Jfact.noRows());
    dNdz.resize(Jfact.noRows());

    for (int i=0; i<Jfact.noRows(); ++i) {
	if (J > 0) {
	    dNdx[i] = Jfact(i,1)/J;
	    dNdy[i] = Jfact(i,2)/J;
	    dNdz[i] = Jfact(i,3)/J;
	} else {
	    dNdx[i] = 0.0;
	    dNdy[i] = 0.0;
	    dNdz[i] = 0.0;
	}
    }

    return 0;
}

int PFEMElement3DBubble::updateMatrix()
{
    int ndf = this->getNumDOF();
    M.resize(ndf, ndf);
    M.Zero();
    D.resize(ndf, ndf);
    D.Zero();
    F.resize(12);
    F.Zero();
    Fp.resize(4);
    Fp.Zero();

    double m = getM();
    double mp = getMp();

    // mass
    for(int a=0; a<(int)thePCs.size(); a++) {
        M(numDOFs(2*a), numDOFs(2*a)) = m;          // Mxd
        M(numDOFs(2*a)+1, numDOFs(2*a)+1) = m;      // Myd
        M(numDOFs(2*a)+2, numDOFs(2*a)+2) = m;      // Mzd
        M(numDOFs(2*a+1), numDOFs(2*a+1)) = mp;      // Mpd
    }

    // damp
    Matrix G,L;
    getG(G);
    getL(L);

    // other matrices
    for(int a=0; a<(int)thePCs.size(); ++a) {
        for(int b=0; b<(int)thePCs.size(); ++b) {

            // Gt
            D(numDOFs(2*a+1), numDOFs(2*b)) = G(3*b,a);   // GxT
            D(numDOFs(2*a+1), numDOFs(2*b)+1) = G(3*b+1,a); // GyT
            D(numDOFs(2*a+1), numDOFs(2*b)+2) = G(3*b+2,a); // GzT

            // G
            D(numDOFs(2*a), numDOFs(2*b+1)) = -G(3*a,b);   // -Gx
            D(numDOFs(2*a)+1, numDOFs(2*b+1)) = -G(3*a+1,b); // -Gy
            D(numDOFs(2*a)+2, numDOFs(2*b+1)) = -G(3*a+2,b); // -Gz

            // L
            D(numDOFs(2*a+1), numDOFs(2*b+1)) = L(a,b);   // bubble
        }
    }

    // force
    getFp(Fp);
    getF(F);

    return 0;
}

int
PFEMElement3DBubble::update()
{
    if (dispon) {
	return updateJacobi();
    }

    return 0;
}

const Matrix&
PFEMElement3DBubble::getMass()
{
    return M;
}

const Matrix&
PFEMElement3DBubble::getDamp()
{
    return D;
}

const Matrix&
PFEMElement3DBubble::getTangentStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}


const Matrix&
PFEMElement3DBubble::getInitialStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement3DBubble::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement3DBubble::getResistingForce()
{

    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement3DBubble::getResistingForceIncInertia()
{
    if (!dispon) {
        if (M.noCols() == 0) {
            updateMatrix();
        }
    }

    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    if (J == 0) {
	return P;
    }

    // get velocity, accleration
    Vector v(ndf), vdot(ndf);
    for(int i=0; i<(int)thePCs.size(); i++) {
        const Vector& accel = nodes[2*i]->getTrialAccel();
        vdot(numDOFs(2*i)) = accel(0);
        vdot(numDOFs(2*i)+1) = accel(1);
        vdot(numDOFs(2*i)+2) = accel(2);

        const Vector& accel2 = nodes[2*i+1]->getTrialAccel();  // pressure
        vdot(numDOFs(2*i+1)) = accel2(0);

        const Vector& vel = nodes[2*i]->getTrialVel();
        v(numDOFs(2*i)) = vel(0);
        v(numDOFs(2*i)+1) = vel(1);
        v(numDOFs(2*i)+2) = vel(2);

        const Vector& vel2 = nodes[2*i+1]->getTrialVel();   // pressure
        v(numDOFs(2*i+1)) = vel2(0);

    }

    // internal force
    P.addMatrixVector(1.0, getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, getDamp(), v, 1.0);

    // external force
    for(int i=0; i<(int)thePCs.size(); i++) {
        P(numDOFs(2*i)) -= F(3*i);
        P(numDOFs(2*i)+1) -= F(3*i+1);
        P(numDOFs(2*i)+2) -= F(3*i+2);
        P(numDOFs(2*i+1)) -= Fp(i);
    }

    //opserr<<"F = "<<F;
    return P;
}


const char*
PFEMElement3DBubble::getClassType()const
{
    return "PFEMElement3DBubble";
}

int
PFEMElement3DBubble::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PFEMElement3DBubble::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement3DBubble::setDomain(Domain *theDomain)
{
    numDOFs.resize(ntags.Size()+1);
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    numDOFs.Zero();
    int ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<(int)thePCs.size(); i++) {

        // set ndf
        numDOFs(2*i) = ndf;

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement3DBubble - setDomain() "<<eletag<<"\n ";
            return;
        }
	if (nodes[2*i]->getNumberDOF() < 3) {
	    opserr<<"WARNING: node "<<ntags(2*i)<<" ndm < 3 ";
            opserr<<"in PFEMElement3DBubble - setDomain() "<<eletag<<"\n ";
            return;
	}
        ndf += nodes[2*i]->getNumberDOF();

        // set ndf
        numDOFs(2*i+1) = ndf;

        // get pc
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] == 0) {
	    opserr << "WARNING: failed to get PC -- PFEMElement3DBubble\n";
	    return;
        }
	thePCs[i]->setDomain(theDomain);

        // connect
        thePCs[i]->connect(eletag);

        // get pressure node
	nodes[2*i+1] = thePCs[i]->getPressureNode();
        if (nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in PFEMElement3DBubble - setDomain() "<<eletag<<"\n ";
            return;
        }
	ntags(2*i+1) = nodes[2*i+1]->getTag();
        ndf += nodes[2*i+1]->getNumberDOF();
    }
    numDOFs(numDOFs.Size()-1) = ndf;

    if (!dispon) {
	if (updateJacobi() < 0) {
	    opserr << "WARNING: failed to update Jacobi -- Bubble3D::setDomain\n";
	    return;
	}
    }
}

void
PFEMElement3DBubble::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement3DBubble: "<<this->getTag()<<endln;
}

int
PFEMElement3DBubble::displaySelf(Renderer &theViewer, int displayMode, float fact,
				 const char **displayModes, int numModes)
{
    return 0;
}

double
PFEMElement3DBubble::getM() const
{
    return rho*J/24.0;
}

double
PFEMElement3DBubble::getMp() const
{
    if(kappa <= 0) return 0.0;
    return J/(kappa*24.0);
}

double
PFEMElement3DBubble::getinvMbub() const
{
    return 5040.*ops_Dt/(256.*rho*J);
}

void
PFEMElement3DBubble::getK(Matrix& k) const
{
    
    double J2 = mu*J/18.0;

    k.resize(12,12);
    k.Zero();

    if (mu <= 0) {
	return;
    }

    // other matrices
    for(int a=0; a<4; ++a) {
        for(int b=0; b<4; ++b) {
            k(3*a, 3*b) += J2*(4*dNdx[a]*dNdx[b] + 3*dNdy[a]*dNdy[b] + 3*dNdz[a]*dNdz[b]);
            k(3*a, 3*b+1) += J2*(3*dNdy[a]*dNdx[b]-2*dNdx[a]*dNdy[b]);
            k(3*a, 3*b+2) += J2*(3*dNdz[a]*dNdx[b]-2*dNdx[a]*dNdz[b]);
	    
            k(3*a+1, 3*b) += J2*(3*dNdx[a]*dNdy[b]-2*dNdy[a]*dNdx[b]);
	    k(3*a+1, 3*b+1) += J2*(4*dNdy[a]*dNdy[b] + 3*dNdx[a]*dNdx[b] + 3*dNdz[a]*dNdz[b]);
            k(3*a+1, 3*b+2) += J2*(3*dNdz[a]*dNdy[b] - 2*dNdy[a]*dNdz[b]);

	    k(3*a+2, 3*b) += J2*(3*dNdx[a]*dNdz[b] - 2*dNdz[a]*dNdx[b]);
	    k(3*a+2, 3*b+1) += J2*(3*dNdy[a]*dNdz[b] - 2*dNdz[a]*dNdy[b]);
	    k(3*a+2, 3*b+2) += J2*(4*dNdz[a]*dNdz[b] + 3*dNdy[a]*dNdy[b] + 3*dNdx[a]*dNdx[b]);
        }
    }
}

void
PFEMElement3DBubble::getGbub(Matrix& gbub) const
{
    double g = -256.0/5040.0;
    gbub.resize(3,4);
    for(int b=0; b<4; b++) {
	gbub(0,b) = g*dNdx[b]*J;
	gbub(1,b) = g*dNdy[b]*J;
	gbub(2,b) = g*dNdz[b]*J;
    }
}

void
PFEMElement3DBubble::getF(Vector& f) const
{
    f.resize(12);
    f.Zero();

    // external force
    for(int a=0; a<4; a++) {
        f(3*a) = bx;
        f(3*a+1) = by;
	f(3*a+2) = bz;
    }
    f *= rho*J/24.0;

    // velocity
    if(mu > 0) {
        Vector v(12);
        for(int a=0; a<4; a++) {
            const Vector& vel = nodes[2*a]->getVel();
            v(3*a) = vel(0);
            v(3*a+1) = vel(1);
	    v(3*a+2) = vel(2);
        }
        Matrix k;
        getK(k);
        f.addMatrixVector(1.0, k, v, -1.0);
    }
}

void
PFEMElement3DBubble::getFbub(Vector& fbub) const
{
    fbub.resize(3);

    // external force
    fbub(0) = rho*J*bx*256/5040.0;
    fbub(1) = rho*J*by*256/5040.0;
    fbub(2) = rho*J*bz*256/5040.0;
}


void
PFEMElement3DBubble::getG(Matrix& g) const
{
    g.resize(12,4);
    double gval = 1.0/24.0*J;
    for (int a=0; a<4; ++a) {
	for (int b=0; b<4; ++b) {
	    g(3*a,b) = gval*dNdx[a];
	    g(3*a+1,b) = gval*dNdy[a];
	    g(3*a+2,b) = gval*dNdz[a];
	}
    }
}

void
PFEMElement3DBubble::getL(Matrix& l) const
{
    Matrix Gbub;
    getGbub(Gbub);
    double invMbub = getinvMbub();

    l.resize(4,4);
    l.addMatrixTransposeProduct(0.0, Gbub, Gbub, invMbub);
}

void
PFEMElement3DBubble::getFp(Vector& fp) const
{
    Matrix Gbub;
    getGbub(Gbub);
    double invMbub = getinvMbub();
    Vector Fbub;
    getFbub(Fbub);

    // bubble force
    fp.resize(4);
    fp.Zero();
    fp.addMatrixTransposeVector(0.0, Gbub, Fbub, -invMbub);
}

int
PFEMElement3DBubble::setParameter(const char **argv, int argc,
                                  Parameter &parameter)
{
    if (argc < 1)
        return -1;

    // viscocity of the fluid
    if (strcmp(argv[0],"mu") == 0) {
        parameter.setValue(mu);
        return parameter.addObject(1, this);
    }
    // Mass densitity of the
    if (strcmp(argv[0],"rho") == 0) {
        parameter.setValue(rho);
        return parameter.addObject(2, this);
    }
    // body acceleration of the fluid
    if (strcmp(argv[0],"bx") == 0) {
        parameter.setValue(bx);
        return parameter.addObject(3, this);
    }
    // body acceleration of the
    if (strcmp(argv[0],"by") == 0) {
        parameter.setValue(by);
        return parameter.addObject(4, this);
    }
    // body acceleration of the
    if (strcmp(argv[0],"bz") == 0) {
        parameter.setValue(bz);
        return parameter.addObject(5, this);
    }

    return -1;
}

int
PFEMElement3DBubble::updateParameter (int passparameterID, Information &info)
{
    switch (passparameterID) {
    case 1:
        mu = info.theDouble;
        return 0;
    case 2:
        rho = info.theDouble;
        return 0;
    case 3:
        bx = info.theDouble;
        return 0;
    case 4:
        by = info.theDouble;
        return 0;
    case 5:
	bz = info.theDouble;
    default:
        return -1;
    }
}

int
PFEMElement3DBubble::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}

double
PFEMElement3DBubble::det(const Matrix& m)
{
    if (m.noRows() != m.noCols()) {
	return 0.0;
    }
    int N = m.noRows();
    if (N == 1) {
	return m(0,0);
    }
    if (N == 2) {
	return m(0,0)*m(1,1)-m(0,1)*m(1,0);
    }
    if (N == 3) {
	return m(0,0)*m(1,1)*m(2,2)+m(0,1)*m(1,2)*m(2,0)+m(0,2)*m(1,0)*m(2,1)
	    -m(0,2)*m(1,1)*m(2,0)-m(0,1)*m(1,0)*m(2,2)-m(0,0)*m(1,2)*m(2,1);
    }

    double res = 0.0;
    bool minus = false;
    for (int i=0; i<N; ++i) {
	Matrix subm(N-1,N-1);
	for (int j=0; j<N-1; ++j) {
	    for (int k=0; k<N-1; ++k) {
		if (k >= i) {
		    subm(j,k) = m(j+1,k+1);
		} else {
		    subm(j,k) = m(j+1,k);
		}
		
	    }
	}
	if (minus) {
	    res -= m(0,i)*det(subm);
	} else {
	    res += m(0,i)*det(subm);
	}
	minus = !minus;
    }
    
    return res;
}

void
PFEMElement3DBubble::cofactor(const Matrix& mat, Matrix& mfact)
{
    if (mat.noRows() != mat.noCols()) {
	return;
    }
    int N = mat.noRows();
    mfact.resize(N,N);
    mfact.Zero();

    for (int i=0; i<N; ++i) {
	for (int j=0; j<N; ++j) {

	    Matrix subm(N-1,N-1);
	    for (int m=0; m<N-1; ++m) {
		for (int n=0; n<N-1; ++n) {
		    int ii = m;
		    int jj = n;
		    if (m >= i) {
			ii = m+1;
		    }
		    if (n >= j) {
			jj = n+1;
		    }
		    subm(m,n) = mat(ii,jj);
		}
	    }

	    bool minus = true;
	    if ((i+j) % 2 == 0) {
		minus = false;
	    }

	    if (minus) {
		mfact(j,i) = -det(subm);		
	    } else {
		mfact(j,i) = det(subm);
	    }
	}
    }
}
