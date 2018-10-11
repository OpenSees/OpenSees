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
                                                                        
// $Revision: 1.00 $
// $Date: 2012/01/12 11:27:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2D.h,v $
                                                                        
// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2D.

#include "PFEMElement2D.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Renderer.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <cmath>
#include <ElementIter.h>

Matrix PFEMElement2D::K;
Vector PFEMElement2D::P;

void* OPS_PFEMElement2D()
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

    // rho, mu, b1, b2, (thickness, kappa, lumped)
    numdata = OPS_GetNumRemainingInputArgs();
    if(numdata > 7) numdata = 7;
    double data[7] = {0,0,0,0,1.0,-1,1};
    if(OPS_GetDoubleInput(&numdata,data) < 0) return 0;

    return new PFEMElement2D(idata[0],idata[1],idata[2],idata[3],
			     data[0],data[1],data[2],data[3],data[4],data[5],data[6]);
}

int OPS_PFEMElement2D(Domain& theDomain, const ID& elenodes, ID& eletags)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 4) {
	opserr<<"WARNING: insufficient number of arguments\n";
	return 0;
    }

    // rho, mu, b1, b2, (thickness, kappa, lumped)
    numdata = OPS_GetNumRemainingInputArgs();
    if(numdata > 7) numdata = 7;
    double data[7] = {0,0,0,0,1.0,-1,1};
    if(OPS_GetDoubleInput(&numdata,data) < 0) return 0;

    // create elements
    ElementIter& theEles = theDomain.getElements();
    Element* theEle = theEles();
    int currTag = 0;
    if (theEle != 0) {
	currTag = theEle->getTag();
    }

    eletags.resize(elenodes.Size()/3);

    for (int i=0; i<eletags.Size(); i++) {
	theEle = new PFEMElement2D(--currTag,elenodes(3*i),elenodes(3*i+1),elenodes(3*i+2),
				   data[0],data[1],data[2],data[3],data[4],data[5],data[6]);
	if (theEle == 0) {
	    opserr<<"WARNING: run out of memory for creating element\n";
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
PFEMElement2D::PFEMElement2D()
    :Element(0, ELE_TAG_PFEMElement2D), ntags(6), 
     rho(0), mu(0), b1(0), b2(0), thickness(1.0), kappa(-1),
     lumped(true), ndf(0), M(0.0), Mp(0.0), Km(6,6), S(3,3),
     Gx(3), Gy(3), F(2), Fp(3)
{
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        thePCs[i] = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	pdof[i] = 0;
    }
}

// for object
PFEMElement2D::PFEMElement2D(int tag, int nd1, int nd2, int nd3,
                             double r, double m, double bx, double by, double thk,
			     double ka, bool l)
    :Element(tag, ELE_TAG_PFEMElement2D), ntags(6),
     rho(r), mu(m), b1(bx), b2(by), thickness(thk), kappa(ka),
     lumped(l), ndf(0), M(0.0), Mp(0.0), Km(6,6), S(3,3),
     Gx(3), Gy(3), F(2), Fp(3)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        ntags(2*i+1) = ntags(2*i);
        thePCs[i] = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	pdof[i] = 0;
    }
}


PFEMElement2D::~PFEMElement2D()
{
    for(int i=0; i<3; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}

int
PFEMElement2D::update()
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
    double cc[3], dd[3];
    cc[0] = y[1]-y[2];
    dd[0] = x[2]-x[1];
    cc[1] = y[2]-y[0];
    dd[1] = x[0]-x[2];
    cc[2] = y[0]-y[1];
    dd[2] = x[1]-x[0];

    // get Jacobi
    double J = cc[0]*dd[1]-dd[0]*cc[1];

    if(fabs(J)<1e-15) {
        opserr<<"WARNING: element area is nearly zero";
        opserr<<" -- PFEMElement2D::update\n";
	for (int i=0; i<3; i++) {
	    opserr<<"node "<<ntags[2*i]<<"\n";
	    opserr<<"x = "<<x[i]<<" , y = "<<y[i]<<"\n";
	}

        return -1;
    }

    // get M
    M = rho*J*thickness/6.0;
    Mp = (kappa<=0? 0.0 : J*thickness/kappa/24.0);
    double Mb = 9.*rho*J*thickness/40.0;

    // get Km
    Km.Zero();
    double fact = mu*thickness/(6.*J);
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {
	    Km(2*a,2*b) = fact*(4*cc[a]*cc[b]+3*dd[a]*dd[b]); // Kxx
	    Km(2*a,2*b+1) = fact*(3*dd[a]*cc[b]-2*cc[a]*dd[b]); // Kxy
	    Km(2*a+1,2*b) = fact*(3*cc[a]*dd[b]-2*dd[a]*cc[b]); // Kyx
	    Km(2*a+1,2*b+1) = fact*(3*cc[a]*cc[b]+4*dd[a]*dd[b]); // Kyy
	}
    }

    // get Kb
    Matrix Kb(2,2);
    fact = 27.*mu*thickness/(40.*J);
    double cc2 = 0., dd2 = 0., cd2 = 0.;
    for(int a=0; a<3; a++) {
	cc2 += cc[a]*cc[a];
	dd2 += dd[a]*dd[a];
	cd2 += cc[a]*dd[a];
    }
    Kb(0,0) = fact*(4*cc2+3*dd2); // Kxx
    Kb(0,1) = fact*cd2; // Kxy
    Kb(1,0) = fact*cd2; // Kyx
    Kb(1,1) = fact*(3*cc2+4*dd2); // Kyy

    // get Gx and Gy
    Gx.Zero(); Gy.Zero();
    fact = thickness/6.0;
    for (int a=0; a<3; a++) {
	Gx(a) = cc[a]*fact;
	Gy(a) = dd[a]*fact;
    }

    // get Gb
    Matrix Gb(2,3);
    fact = -9.*thickness/40.0;
    for (int a=0; a<3; a++) {
	Gb(0,a) = cc[a]*fact;
	Gb(1,a) = dd[a]*fact;
    }

    // get S
    S.Zero();
    if (ops_Dt > 0) {
	Kb(0,0) += Mb/ops_Dt;
	Kb(1,1) += Mb/ops_Dt;
    }
    if (Kb(0,0)!=0 && Kb(1,1)!=0) {
	this->inverse(Kb);
    }
    S.addMatrixTripleProduct(0.0, Gb, Kb, 1);

    // get F
    F.Zero();
    fact = rho*J*thickness/6.0;
    F(0) = fact*b1;
    F(1) = fact*b2;

    // get Fb
    Vector Fb(2);
    fact = 9.*rho*J*thickness/40.;
    Fb(0) = fact*b1;
    Fb(1) = fact*b2;

    // get Fp
    Fp.Zero();
    Fp.addMatrixTransposeVector(0.0, Gb, Kb*Fb, -1);

    return 0;
}

const Matrix&
PFEMElement2D::getMass()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M
    for (int a=0; a<3; a++) {
	K(vxdof[a], vxdof[a]) = M;
	K(vydof[a], vydof[a]) = M;
    }

    // Mp
    for(int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {
	    if (a == b) {
		K(pdof[a],pdof[b]) = Mp*2;
	    } else {
		K(pdof[a],pdof[b]) = Mp;
	    }
	}
    }

    return K;
}

const Matrix&
PFEMElement2D::getDamp()
{

    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // K, G, Gt
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {

	    // K
	    if (lumped == false) {
		K(vxdof[a], vxdof[b]) += Km(2*a,2*b);
		K(vxdof[a], vydof[b]) += Km(2*a,2*b+1);
		K(vydof[a], vxdof[b]) += Km(2*a+1,2*b);
		K(vydof[a], vydof[b]) += Km(2*a+1,2*b+1);
	    }

	    // -G
	    K(vxdof[a], pdof[b]) = -Gx[a];
	    K(vydof[a], pdof[b]) = -Gy[a];

	    // Gt
	    K(pdof[b], vxdof[a]) = Gx[a];
	    K(pdof[b], vydof[a]) = Gy[a];

	    // S
	    K(pdof[a], pdof[b]) = S(a,b);
	}

    }


    return K;
}

const Matrix&
PFEMElement2D::getTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}


const Matrix& 
PFEMElement2D::getInitialStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int 
PFEMElement2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2D::getResistingForce()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2D::getResistingForceIncInertia()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    // get v and p
    Vector v(ndf), vdot(ndf);
    for(int a=0; a<3; a++) {

        const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& accel = nodes[2*a]->getTrialAccel();
	v(vxdof[a]) = vel(0);
	v(vydof[a]) = vel(1);
	vdot(vxdof[a]) = accel(0);
	vdot(vydof[a]) = accel(1);

	const Vector& pressure = nodes[2*a+1]->getTrialVel();
	const Vector& pressuredot = nodes[2*a+1]->getTrialAccel();
	v(pdof[a]) = pressure(0);
	vdot(pdof[a]) = pressuredot(0);

	// F, Fp
	P(vxdof[a]) = F(0);
	P(vydof[a]) = F(1);
	P(pdof[a]) = Fp[a];
    }

    // M*vdot+K*v
    P.addMatrixVector(-1.0, getMass(), vdot, 1.0);
    bool l = lumped;
    lumped = false;
    P.addMatrixVector(1.0, getDamp(), v, 1.0);
    lumped = l;

    return P;
}


const char*
PFEMElement2D::getClassType()const
{
    return "PFEMElement2D";
}

int
PFEMElement2D::sendSelf(int commitTag, Channel &theChannel)
{
    // int res = 0;
    // int dataTag = this->getDbTag();

    // // send vector
    // static Vector data(25);
    // data(0) = this->getTag();
    // data(1) = rho;
    // data(2) = mu;
    // data(3) = b1;
    // data(4) = b2;
    // for(int i=0; i<3; i++) {
    //     data(5+i) = dNdx[i];
    //     data(8+i) = dNdy[i];
    // }
    // data(11) = J;
    // for(int i=0; i<6; i++) {
    //     data(12+i) = ntags(i);
    //     data(18+i) = numDOFs(i);
    // }
    // data(24) = numDOFs(6);

    // res = theChannel.sendVector(dataTag, commitTag, data);
    // if(res < 0) {
    //     opserr<<"WARNING: PFEMElement2D::sendSelf - "<<this->getTag()<<" failed to send vector\n";
    //     return -1;
    // }


    return 0;
}

int
PFEMElement2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // int res;
    // int dataTag = this->getDbTag();

    // // receive vector
    // static Vector data(12);
    // res = theChannel.recvVector(dataTag, commitTag, data);
    // if(res < 0) {
    //     opserr<<"WARNING: PFEMElement2D::recvSelf - failed to receive vector\n";
    //     return -1;
    // }
    // this->setTag((int)data(0));
    // rho = data(1);
    // mu = data(2);
    // bx = data(3);
    // by = data(4);
    // for(int i=0; i<3; i++) {
    //     dNdx[i] = data(5+i);
    //     dNdy[i] = data(8+i);
    // }
    // J = data(11);
    // for(int i=0; i<6; i++) {
    //     ntags(i) = (int)data(12+i);
    //     numDOFs(i) = (int)data(18+i);
    // }
    // numDOFs(6) = (int)data(24);

    return 0;
}

void
PFEMElement2D::setDomain(Domain *theDomain)
{
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<3; i++) {

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement2D - setDomain() "<<eletag<<"\n ";
            return;
        }
	if(nodes[2*i]->getNumberDOF() < 2) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" ndf < 2 ";
            opserr<<"in PFEMElement2D - setDomain() "<<eletag<<"\n ";
            return;
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
                opserr<<"PFEMElement2D::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement2D::setDomain "<<eletag<<"\n";
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
            opserr<<"in PFEMElement2D - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
	pdof[i] = ndf;
        ndf += nodes[2*i+1]->getNumberDOF();
    }

}

void
PFEMElement2D::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2D: "<<this->getTag()<<endln;
}

int
PFEMElement2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes) 
{

  static Vector values(3);

    // now  determine the end points of the Tri31 based on
    // the display factor (a measure of the distorted image)
    // store this information in 3 3d vectors v1 through v3
    const Vector &end1Crd = nodes[0]->getCrds();
    const Vector &end2Crd = nodes[2]->getCrds();	
    const Vector &end3Crd = nodes[4]->getCrds();	

    const int numnodes = 3;
    static Matrix coords(numnodes,3);

    if (displayMode >= 0) {  

        const Vector &end1Disp = nodes[0]->getDisp();
        const Vector &end2Disp = nodes[2]->getDisp();
        const Vector &end3Disp = nodes[4]->getDisp();

        for (int i = 0; i < 2; i++) {
            coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
            coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
            coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
        }
    } else {
        int mode = displayMode * -1;
        const Matrix &eigen1 = nodes[0]->getEigenvectors();
        const Matrix &eigen2 = nodes[2]->getEigenvectors();
        const Matrix &eigen3 = nodes[4]->getEigenvectors();
        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 2; i++) {
                coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
                coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
            }    
        } else {
            for (int i = 0; i < 2; i++) {
                coords(0,i) = end1Crd(i);
                coords(1,i) = end2Crd(i);
                coords(2,i) = end3Crd(i);
            }    
        }
    }
    
    int error = 0;

    // finally we draw the element using drawPolygon
    error += theViewer.drawPolygon (coords, values);

    return error;
}

void
PFEMElement2D::inverse(Matrix& mat) const
{
    double tmp = mat(0,0);
    mat(0,0) = mat(1,1);
    mat(1,1) = tmp;
    mat(0,1) = -mat(0,1);
    mat(1,0) = -mat(1,0);

    tmp = mat(0,0)*mat(1,1)-mat(1,0)*mat(0,1);
    
    mat /= tmp;
}
