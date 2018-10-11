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
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2DCompressible.h,v $

// Written: Minjie Zhu (zhum@oregontate.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2DCompressible.

#include "PFEMElement2DCompressible.h"
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
#include <map>

Matrix PFEMElement2DCompressible::K;
Vector PFEMElement2DCompressible::P;
bool PFEMElement2DCompressible::dispon = true;

void* OPS_PFEMElement2DCompressible(const ID &info)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
    	opserr << "WARNING: domain is not created\n";
    	return 0;
    }

    int idata[5];
    double data[6] = {0,0,0,0,1.0,2.15e9};
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

    // regular element, or in a mesh
    if (info.Size()==0 || info(0)==1) {
	if(OPS_GetNumRemainingInputArgs() < 4) {
	    opserr<<"insufficient arguments: rho, mu, b1, b2, (thinknes,kappa)\n";
	    return 0;
	}

	// rho, mu, b1, b2, (thickness, kappa)
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

    return new PFEMElement2DCompressible(idata[0],idata[1],idata[2],idata[3],idata[4],
					 data[0],data[1],data[2],data[3],data[4],
					 data[5]);
}

// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2DCompressible::PFEMElement2DCompressible()
    :Element(0, ELE_TAG_PFEMElement2DCompressible), ntags(7),
     rho(0), mu(0), b1(0), b2(0), thickness(1.0), kappa(2.15e9),
     ndf(0), J(0.0), parameterID(0), bubblenode(0)
{
    for(int i=0;i<4;i++)
    {
	nodes[2*i] = 0;
	ntags(2*i) = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	pdof[i] = 0;
	thePCs[i] = 0;
	if (i < 3) {
	    nodes[2*i+1] = 0;
	    ntags(2*i+1) = 0;
	    cc[i] = 0;
	    dd[i] = 0;
	}

    }
}

// for object
PFEMElement2DCompressible::PFEMElement2DCompressible(
    int tag, int nd1, int nd2, int nd3, int nd4,
    double r, double m, double bx, double by, double thk,
    double ka)
    :Element(tag, ELE_TAG_PFEMElement2DCompressible), ntags(7),
     rho(r), mu(m), b1(bx), b2(by), thickness(thk), kappa(ka),
     ndf(0), J(0.0), parameterID(0), bubblenode(nd4)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    ntags(1)=nd1; ntags(3)=nd2; ntags(5)=nd3;
    ntags(6)=nd3;
    for(int i=0;i<4;i++)
    {
	nodes[2*i] = 0;
	vxdof[i] = 0;
	vydof[i] = 0;
	pdof[i] = 0;
	thePCs[i] = 0;
	if (i < 3) {
	    nodes[2*i+1] = 0;
	    cc[i] = 0;
	    dd[i] = 0;
	}
    }
    if (kappa <= 0.0) {
	kappa = 2.15e9;
    }
}


PFEMElement2DCompressible::~PFEMElement2DCompressible()
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
	if (thePCs[3] != 0) {
	    domain->removePressure_Constraint(thePCs[3]->getTag());
	    delete thePCs[3];
	    thePCs[3] = 0;
	}
    }
}

int
PFEMElement2DCompressible::commitState()
{
    // Vector bdisp(2);
    // for(int a=0; a<3; a++) {
    //     const Vector& disp = nodes[2*a]->getTrialDisp();
    // 	bdisp(0) += disp(0);
    // 	bdisp(1) += disp(1);
    // }
    // bdisp /= 3.0;
    // nodes[6]->setTrialDisp(bdisp);
    // nodes[6]->commitState();
    if (!dispon) {
	if (updateJacobian() < 0) return -1;
    }
    return Element::commitState();
}

int
PFEMElement2DCompressible::update()
{
    if (dispon) {
	return updateJacobian();
    }

    return 0;
}

int
PFEMElement2DCompressible::updateJacobian()
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
        opserr<<"WARNING: element area is negative";
        opserr<<" -- PFEMElement2DCompressible::update\n";
        return -1;
    }

    return 0;
}

const Matrix&
PFEMElement2DCompressible::getMass()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M
    // double M = rho*J*thickness/24.0;
    // double M2 = 2.0*M;
    // double M3 = 3.0*rho*J*thickness/40.0;
    // double M4 = 81.0*rho*J*thickness/560.0;
    double m = rho*J*thickness/6.0;
    double m4 = 27.0*rho*J*thickness/120.0;
    for (int a=0; a<3; a++) {
	K(vxdof[a], vxdof[a]) = m;
	K(vydof[a], vydof[a]) = m;
	// // M(a,b) a,b=1,2,3
	// for (int b=0; b<3; b++) {
	//     if (a == b) {
	// 	K(vxdof[a], vxdof[b]) = M2;
	// 	K(vydof[a], vydof[b]) = M2;
	//     } else {
	// 	K(vxdof[a], vxdof[b]) = M;
	// 	K(vydof[a], vydof[b]) = M;
	//     }
	// }
	// // M(4,a)
	// K(vxdof[3],vxdof[a]) = M3;
	// K(vydof[3],vydof[a]) = M3;

	// // M(a,4)
	// K(vxdof[a],vxdof[3]) = M3;
	// K(vydof[a],vydof[3]) = M3;
    }

    // M(4,4)
    K(vxdof[3],vxdof[3]) = m4;
    K(vydof[3],vydof[3]) = m4;

    // Mp
    double Mp = J*thickness/6.0/kappa;
    for(int a=0; a<3; a++) {
	K(pdof[a],pdof[a]) = Mp;
    }

    return K;
}

const Matrix&
PFEMElement2DCompressible::getDamp()
{

    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // Km
    double km = mu*thickness/(6.*J);
    double kb = 729.0*mu*thickness/(1080.*J);

    // G
    double g = thickness/6.0;
    double gb = -27.0*thickness/120.0;
    double gt = thickness/6.0;
    double gtb = -27.0*thickness/120.0;

    // K(a,b) a,b = 1,2,3
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {

	    // K
	    K(vxdof[a], vxdof[b]) += km*(4*cc[a]*cc[b]+3*dd[a]*dd[b]); // Kxx
	    K(vxdof[a], vydof[b]) += km*(3*dd[a]*cc[b]-2*cc[a]*dd[b]); // Kxy
	    K(vydof[a], vxdof[b]) += km*(3*cc[a]*dd[b]-2*dd[a]*cc[b]); // Kyx
	    K(vydof[a], vydof[b]) += km*(3*cc[a]*cc[b]+4*dd[a]*dd[b]); // Kyy

	    // -G
	    K(vxdof[a], pdof[b]) = -g*cc[a];
	    K(vydof[a], pdof[b]) = -g*dd[a];

	    // Gt
	    K(pdof[b], vxdof[a]) = gt*cc[a];
	    K(pdof[b], vydof[a]) = gt*dd[a];
	}

	// -G(4,a)
	K(vxdof[3], pdof[a]) = -gb*cc[a];
	K(vydof[3], pdof[a]) = -gb*dd[a];

	// Gt(a,4)
	K(pdof[a], vxdof[3]) = gtb*cc[a];
	K(pdof[a], vydof[3]) = gtb*dd[a];
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

    return K;
}

const Matrix&
PFEMElement2DCompressible::getTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

const Matrix&
PFEMElement2DCompressible::getGeometricTangentStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // get accel, vel and pressure
    Vector vdot(8), v(8), pdot(3), p(3);
    for(int a=0; a<4; a++) {
	const Vector& vel = nodes[2*a]->getTrialVel();
	const Vector& accel = nodes[2*a]->getTrialAccel();

	vdot(2*a) = accel(0);
	vdot(2*a+1) = accel(1);

	v(2*a) = vel(0);
	v(2*a+1) = vel(1);

    }
    for(int a=0; a<3; a++) {
	const Vector& pre = nodes[2*a+1]->getTrialVel();
	const Vector& predot = nodes[2*a+1]->getTrialAccel();

	pdot(a) = predot(0);
	p(a) = pre(0);
    }

    Matrix dM, dG, dF, dGt, dMp, dK;
    getdM(vdot, dM);
    getdK(v,dK);
    getdG(p,v,dG,dGt);
    getdF(dF);
    getdMp(pdot,dMp);

    for (int b=0; b<3; b++) {
	for (int a=0; a<4; a++) {
	    K(vxdof[a],vxdof[b]) += dM(2*a,2*b)+dK(2*a,2*b)-dG(2*a,2*b)-dF(2*a,2*b);
	    K(vxdof[a],vydof[b]) += dM(2*a,2*b+1)+dK(2*a,2*b+1)-dG(2*a,2*b+1)-dF(2*a,2*b+1);
	    K(vydof[a],vxdof[b]) += dM(2*a+1,2*b)+dK(2*a+1,2*b)-dG(2*a+1,2*b)-dF(2*a+1,2*b);
	    K(vydof[a],vydof[b]) += dM(2*a+1,2*b+1)+dK(2*a+1,2*b+1)-dG(2*a+1,2*b+1)-dF(2*a+1,2*b+1);
	}
	for (int a=0; a<3; a++) {
	    K(pdof[a],vxdof[b]) += dMp(a,2*b)+dGt(a,2*b);
	    K(pdof[a],vydof[b]) += dMp(a,2*b+1)+dGt(a,2*b+1);
	}
    }


/*
    // K44
    double cc2=0.0, dd2=0.0, cd2=0.0;
    for (int a=0; a<3; a++) {
	cc2 += cc[a]*cc[a];
	dd2 += dd[a]*dd[a];
	cd2 += cc[a]*dd[a];
    }

    // M*vdot-Fext+dKdA
    Vector mvdot(8);
    mvdot(0) = vdot(0)/12.0+vdot(2)/24.0+vdot(4)/24.0+vdot(6)*3.0/40.0-b1/6.0;
    mvdot(1) = vdot(1)/12.0+vdot(3)/24.0+vdot(5)/24.0+vdot(7)*3.0/40.0-b2/6.0;
    mvdot(2) = vdot(0)/24.0+vdot(2)/12.0+vdot(4)/24.0+vdot(6)*3.0/40.0-b1/6.0;
    mvdot(3) = vdot(1)/24.0+vdot(3)/12.0+vdot(5)/24.0+vdot(7)*3.0/40.0-b2/6.0;
    mvdot(4) = vdot(0)/24.0+vdot(2)/24.0+vdot(4)/12.0+vdot(6)*3.0/40.0-b1/6.0;
    mvdot(5) = vdot(1)/24.0+vdot(3)/24.0+vdot(5)/12.0+vdot(7)*3.0/40.0-b2/6.0;
    mvdot(6) = vdot(0)*3.0/40.0+vdot(2)*3.0/40.0+vdot(4)*3.0/40.0+vdot(6)*207.0/506.0-9.0*b1/40.0;
    mvdot(7) = vdot(1)*3.0/40.0+vdot(3)*3.0/40.0+vdot(5)*3.0/40.0+vdot(7)*207.0/506.0-9.0*b2/40.0;
    mvdot *= rho;

    for (int b=0; b<3; b++) {
	for (int a=0; a<3; a++) {
	    mvdot(2*a) += -mu/(6.0*J*J)*((4*cc[a]*cc[b]+3*dd[a]*dd[b])*v(2*b)+(3*dd[a]*cc[b]-2*cc[a]*dd[b])*v(2*b+1));
	    mvdot(2*a+1) += -mu/(6.0*J*J)*((3*cc[a]*dd[b]-2*dd[a]*cc[b])*v(2*b)+(4*dd[a]*dd[b]+3*cc[a]*cc[b])*v(2*b+1));
	}
    }
    mvdot(6) += -27.0*mu/(40.0*J*J)*((4*cc2+3*dd2)*v(6)+(cd2)*v(7));
    mvdot(7) += -27.0*mu/(40.0*J*J)*((cd2)*v(6)+(4*dd2+3*cc2)*v(7));

    // (M*vdot-Fext)*dAdu
    for (int b=0; b<3; b++) {
	for (int a=0; a<4; a++) {
	    K(vxdof[a], vxdof[b]) += mvdot(2*a)*cc[b];
	    K(vxdof[a], vydof[b]) += mvdot(2*a)*dd[b];
	    K(vydof[a], vxdof[b]) += mvdot(2*a+1)*cc[b];
	    K(vydof[a], vydof[b]) += mvdot(2*a+1)*dd[b];
	}
    }

    // dKdu1
    double k = mu/(6.0*J);
    double kb = 27.0*mu/(40.0*J);
    Matrix dK(8,8);
    dK(0,2) = dK(2,0) = 3*dd[0]*k;
    dK(0,3) = dK(3,0) = -2*cc[0]*k;
    dK(0,4) = dK(4,0) = -3*dd[0]*k;
    dK(0,5) = dK(5,0) = 2*cc[0]*k;
    dK(1,2) = dK(2,1) = 3*cc[0]*k;
    dK(1,3) = dK(3,1) = 4*dd[0]*k;
    dK(1,4) = dK(4,1) = -3*cc[0]*k;
    dK(1,5) = dK(5,1) = -4*dd[0]*k;
    dK(2,2) = 6*dd[1]*k;
    dK(2,3) = dK(3,2) = cc[1]*k;
    dK(2,4) = dK(4,2) = (3*dd[2]-3*dd[1])*k;
    dK(2,5) = dK(5,2) = (3*cc[2]+2*cc[1])*k;
    dK(3,3) = 8*dd[1]*k;
    dK(3,4) = dK(4,3) = (-3*cc[1]-2*cc[2])*k;
    dK(3,5) = dK(5,3) = (4*dd[2]-4*dd[1])*k;
    dK(4,4) = -6*dd[2]*k;
    dK(4,5) = dK(5,4) = -cc[2]*k;
    dK(5,5) = -8*dd[2]*k;
    dK(6,6) = 6*kb*(dd[1]-dd[2]);
    dK(6,7) = dK(7,6) = kb*(cc[1]-cc[2]);
    dK(7,7) = kb*8*(dd[1]-dd[2]);

    Vector kv(8);
    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],0) += kv(2*a);
	K(vydof[a],0) += kv(2*a+1);
    }

    // dKdu2
    dK.Zero();
    dK(0,2) = dK(2,0) = -4*cc[0]*k;
    dK(0,3) = dK(3,0) = -3*dd[0]*k;
    dK(0,4) = dK(4,0) = 4*cc[0]*k;
    dK(0,5) = dK(5,0) = 3*dd[0]*k;
    dK(1,2) = dK(2,1) = 2*dd[0]*k;
    dK(1,3) = dK(3,1) = -3*cc[0]*k;
    dK(1,4) = dK(4,1) = -2*dd[0]*k;
    dK(1,5) = dK(5,1) = 3*cc[0]*k;
    dK(2,2) = -8*cc[1]*k;
    dK(2,3) = dK(3,2) = -dd[1]*k;
    dK(2,4) = dK(4,2) = (4*cc[1]-4*cc[2])*k;
    dK(2,5) = dK(5,2) = (3*dd[1]+2*dd[2])*k;
    dK(3,3) = -6*cc[1]*k;
    dK(3,4) = dK(4,3) = (-3*dd[2]-2*dd[1])*k;
    dK(3,5) = dK(5,3) = (3*cc[1]-3*cc[2])*k;
    dK(4,4) = 8*cc[2]*k;
    dK(4,5) = dK(5,4) = dd[2]*k;
    dK(5,5) = 6*cc[2]*k;
    dK(6,6) = 8*kb*(cc[2]-cc[1]);
    dK(6,7) = dK(7,6) = kb*(dd[2]-dd[1]);
    dK(7,7) = kb*6*(cc[2]-cc[1]);

    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],1) += kv(2*a);
	K(vydof[a],1) += kv(2*a+1);
    }

    // dKdu3
    dK.Zero();
    dK(0,0) = -6*k*dd[0];
    dK(0,1) = dK(1,0) = -k*cc[0];
    dK(0,2) = dK(2,0) = -3*dd[1]*k;
    dK(0,3) = dK(3,0) = -3*cc[1]*k;
    dK(0,4) = dK(4,0) = k*(3*dd[0]-3*dd[2]);
    dK(0,5) = dK(5,0) = k*(-3*cc[2]-2*cc[0]);
    dK(1,1) = -8*k*dd[0];
    dK(1,2) = dK(2,1) = 2*cc[1]*k;
    dK(1,3) = dK(3,1) = -4*dd[1]*k;
    dK(1,4) = dK(4,1) = k*(3*cc[0]+2*cc[2]);
    dK(1,5) = dK(5,1) = k*(4*dd[0]-4*dd[2]);
    dK(2,4) = dK(4,2) = 3*k*dd[1];
    dK(2,5) = dK(5,2) = -2*k*cc[1];
    dK(3,4) = dK(4,3) = 3*k*cc[1];
    dK(3,5) = dK(5,3) = 4*k*dd[1];
    dK(4,4) = 6*k*dd[2];
    dK(4,5) = dK(5,4) = k*cc[2];
    dK(5,5) = 8*k*dd[2];
    dK(6,6) = kb*6*(dd[2]-dd[0]);
    dK(6,7) = dK(7,6) = kb*(cc[2]-cc[0]);
    dK(7,7) = 8*kb*(dd[2]-dd[0]);

    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],2) += kv(2*a);
	K(vydof[a],2) += kv(2*a+1);
    }

    // dKdu4
    dK.Zero();
    dK(0,0) = 8*k*cc[0];
    dK(0,1) = dK(1,0) = k*dd[0];
    dK(0,2) = dK(2,0) = 4*k*cc[1];
    dK(0,3) = dK(3,0) = -2*k*dd[1];
    dK(0,4) = dK(4,0) = 4*k*cc[2];
    dK(0,5) = dK(5,0) = k*(-3*dd[0]-2*dd[2]);
    dK(1,1) = 6*k*cc[0];
    dK(1,2) = dK(2,1) = 3*k*dd[1];
    dK(1,3) = dK(3,1) = 3*k*cc[1];
    dK(1,4) = dK(4,1) = k*(3*dd[2]+2*dd[0]);
    dK(1,5) = dK(5,1) = k*(3*cc[2]-3*cc[0]);
    dK(2,4) = dK(4,2) = -4*k*cc[1];
    dK(2,5) = dK(5,2) = -3*k*dd[1];
    dK(3,4) = dK(4,3) = 3*k*dd[1];
    dK(3,5) = dK(5,3) = -3*k*cc[1];
    dK(4,4) = -8*k*cc[2];
    dK(4,5) = dK(5,4) = -k*dd[2];
    dK(5,5) = -6*k*cc[2];
    dK(6,6) = kb*8*(cc[0]-cc[2]);
    dK(6,7) = kb*(dd[0]-dd[2]);
    dK(7,7) = kb*6*(cc[0]-cc[2]);

    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],3) += kv(2*a);
	K(vydof[a],3) += kv(2*a+1);
    }

    // dKdu5
    dK.Zero();
    dK(0,0) = 6*k*dd[0];
    dK(0,1) = dK(1,0) = k*cc[0];
    dK(0,2) = dK(2,0) = k*(3*dd[1]-3*dd[0]);
    dK(0,3) = dK(3,0) = k*(3*cc[1]+2*cc[0]);
    dK(0,4) = dK(4,0) = k*3*dd[2];
    dK(0,5) = dK(5,0) = k*3*cc[2];
    dK(1,1) = 8*k*dd[0];
    dK(1,2) = dK(2,1) = k*(-3*cc[0]-2*cc[1]);
    dK(1,3) = dK(3,1) = k*(4*dd[1]-4*dd[0]);
    dK(1,4) = dK(4,1) = -2*k*cc[2];
    dK(1,5) = dK(5,1) = 4*k*dd[2];
    dK(2,2) = -6*k*dd[1];
    dK(2,3) = dK(3,2) = -k*cc[1];
    dK(2,4) = dK(4,2) = -3*k*dd[2];
    dK(2,5) = dK(5,2) = -3*k*cc[2];
    dK(3,3) = -8*k*dd[1];
    dK(3,4) = dK(4,3) = 2*k*cc[2];
    dK(3,5) = dK(5,3) = -4*k*dd[2];
    dK(6,6) = kb*6*(dd[0]-dd[1]);
    dK(6,7) = kb*(cc[0]-cc[1]);
    dK(7,7) = kb*8*(dd[0]-dd[1]);

    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],4) += kv(2*a);
	K(vydof[a],4) += kv(2*a+1);
    }

    // dKdu6
    dK.Zero();
    dK(0,0) = -8*k*cc[0];
    dK(0,1) = dK(1,0) = -k*dd[0];
    dK(0,2) = dK(2,0) = k*(4*cc[0]-4*cc[1]);
    dK(0,3) = dK(3,0) = k*(3*dd[0]+2*dd[1]);
    dK(0,4) = dK(4,0) = -4*k*cc[2];
    dK(0,5) = dK(5,0) = 2*k*dd[2];
    dK(1,1) = -6*k*cc[0];
    dK(1,2) = dK(2,1) = k*(-3*dd[1]-2*dd[0]);
    dK(1,3) = dK(3,1) = k*(3*cc[0]-3*cc[1]);
    dK(1,4) = dK(4,1) = -3*k*dd[2];
    dK(1,5) = dK(5,1) = -3*k*cc[2];
    dK(2,2) = 8*k*cc[1];
    dK(2,3) = dK(3,2) = k*dd[1];
    dK(2,4) = dK(4,2) = 4*k*cc[2];
    dK(2,5) = dK(5,2) = -2*k*dd[2];
    dK(3,3) = 6*k*cc[1];
    dK(3,4) = dK(4,3) = 3*k*dd[2];
    dK(3,5) = dK(5,3) = 3*k*cc[2];
    dK(6,6) = kb*8*(cc[1]-cc[0]);
    dK(6,7) = kb*(dd[1]-dd[0]);
    dK(7,7) = kb*6*(cc[1]-cc[0]);

    kv.addMatrixVector(0.0,dK,v,1.0);
    for (int a=0; a<4; a++) {
	K(vxdof[a],5) += kv(2*a);
	K(vydof[a],5) += kv(2*a+1);
    }

    // sum of p
    double sp = (p(0)+p(1)+p(2))/6.0;

    // -dGdu*p
    for (int b=0; b<3; b++) {
	int b1 = b+1; if (b1 > 2) b1-=3;
	int b2 = b+2; if (b2 > 2) b2-=3;
	K(vydof[b1],vxdof[b]) -= sp;
	K(vydof[b2],vxdof[b]) -= -sp;
	K(vydof[3],vxdof[b]) -= (p(b1)-p(b2))*9.0/(-40.0);

	K(vxdof[b2],vydof[b]) -= sp;
	K(vxdof[b1],vydof[b]) -= -sp;
	K(vxdof[3],vydof[b]) -= (p(b2)-p(b1))*9.0*(-40.0);
    }

    // Mp*pdot*dAdu + dGtdu*v
    for (int b=0; b<3; b++) {
	int b1 = b+1; if (b1 > 2) b1-=3;
	int b2 = b+2; if (b2 > 2) b2-=3;
	for (int a=0; a<3; a++) {
	    K(pdof[a], vxdof[b]) += pdot(a)*cc[b]/6.0 + (v(2*b1+1)-v(2*b2+1))/6.0*kappa;
	    K(pdof[a], vydof[b]) += pdot(a)*dd[b]/6.0 + (v(2*b2)-v(2*b1))/6.0*kappa;
	}
	K(pdof[b1], vxdof[b]) += v(7)*9.0/(-40.0)*kappa;
	K(pdof[b2], vxdof[b]) -= v(7)*9.0/(-40.0)*kappa;
	K(pdof[b1], vydof[b]) -= v(6)*9.0/(-40.0)*kappa;
	K(pdof[b2], vydof[b]) += v(6)*9.0/(-40.0)*kappa;
    }


    K *= thickness;
*/
    return K;
}

const Matrix&
PFEMElement2DCompressible::getInitialStiff()
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement2DCompressible::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2DCompressible::getResistingForce()
{
    // resize P
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2DCompressible::getResistingForceIncInertia()
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
    double fb = 27.0*rho*J*thickness/120.0;
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
    P.addMatrixVector(-1.0, this->getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, this->getDamp(), v, 1.0);

    return P;
}


const char*
PFEMElement2DCompressible::getClassType()const
{
    return "PFEMElement2DCompressible";
}

int
PFEMElement2DCompressible::sendSelf(int commitTag, Channel &theChannel)
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
    //     opserr<<"WARNING: PFEMElement2DCompressible::sendSelf - "<<this->getTag()<<" failed to send vector\n";
    //     return -1;
    // }


    return 0;
}

int
PFEMElement2DCompressible::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // int res;
    // int dataTag = this->getDbTag();

    // // receive vector
    // static Vector data(12);
    // res = theChannel.recvVector(dataTag, commitTag, data);
    // if(res < 0) {
    //     opserr<<"WARNING: PFEMElement2DCompressible::recvSelf - failed to receive vector\n";
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
PFEMElement2DCompressible::setDomain(Domain *theDomain)
{
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    ndf = 0;
    int eletag = this->getTag();
    Vector bcrds(2);
    for(int i=0; i<3; i++) {

	// set ndf
	vxdof[i] = ndf;
	vydof[i] = ndf+1;

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement2DCompressible - setDomain() "<<eletag<<"\n ";
            return;
        }
	if(nodes[2*i]->getNumberDOF() < 2) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" ndf < 2 ";
            opserr<<"in PFEMElement2DCompressible - setDomain() "<<eletag<<"\n ";
            return;
        }
	const Vector& crds = nodes[2*i]->getCrds();
	for (int j=0; j<2; j++) {
	    bcrds(j) += crds(j);
	}
        ndf += nodes[2*i]->getNumberDOF();

	// set ndf
	pdof[i] = ndf;

        // get pc
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
	if(thePCs[i] == 0) {
	    opserr << "WARNING: failed to get PC -- PFEMElement2DBubble\n";
	    return;
        }
	thePCs[i]->setDomain(theDomain);

        // connect
        thePCs[i]->connect(eletag);

        // get pressure node
        nodes[2*i+1] = thePCs[i]->getPressureNode();
        if(nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in PFEMElement2DCompressible - setDomain() "<<eletag<<"\n ";
            return;
        }
	ntags(2*i+1) = nodes[2*i+1]->getTag();
        ndf += nodes[2*i+1]->getNumberDOF();
    }

    // create bubble node
    ntags[6] = bubblenode;
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

    if (!dispon) {
	updateJacobian();
    }
}

void
PFEMElement2DCompressible::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2DCompressible: "<<this->getTag()<<endln;
}

int
PFEMElement2DCompressible::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
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

int
PFEMElement2DCompressible::setParameter(const char **argv, int argc,
					Parameter &parameter)
{
    if (argc < 1)
        return -1;

    if (strcmp(argv[0],"mu") == 0) {
        parameter.setValue(mu);
        return parameter.addObject(1, this);
    } else if (strcmp(argv[0],"rho") == 0) {
        parameter.setValue(rho);
        return parameter.addObject(2, this);
    } else if (strcmp(argv[0],"bx") == 0) {
        parameter.setValue(b1);
        return parameter.addObject(3, this);
    } else if (strcmp(argv[0],"by") == 0) {
        parameter.setValue(b2);
        return parameter.addObject(4, this);
    } else if (strcmp(argv[0],"thickness") == 0) {
	parameter.setValue(thickness);
	return parameter.addObject(5, this);
    }

    return -1;
}

int
PFEMElement2DCompressible::updateParameter(int passparameterID, Information &info)
{
    switch (passparameterID) {
    case 1:
        mu = info.theDouble;
        return 0;
    case 2:
        rho = info.theDouble;
        return 0;
    case 3:
        b1 = info.theDouble;
        return 0;
    case 4:
        b2 = info.theDouble;
        return 0;
    case 5:
	thickness = info.theDouble;
	return 0;

    default:
        return -1;
    }
}

int
PFEMElement2DCompressible::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}

const Matrix&
PFEMElement2DCompressible::getDampSensitivity(int gradNumber)
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // Km and G
    double km=0.0, kb=0.0, g=0.0, gb=0.0, gt=0.0, gtb=0.0;
    if (parameterID == 1) {

	// mu
	km = thickness/(6.*J);
	kb = 27.0*thickness/(40.*J);

    } else if (parameterID == 6) {
	// thickness
	km = mu/(6.*J);
	kb = 27.0*mu/(40.*J);
	g = 1.0/6.0;
	gb = -9.0/40.0;
	gt = 1/6.0;
	gtb = -9.0/40.0;
    }

    // K(a,b) a,b = 1,2,3
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {

	    // K
	    K(vxdof[a], vxdof[b]) += km*(4*cc[a]*cc[b]+3*dd[a]*dd[b]); // Kxx
	    K(vxdof[a], vydof[b]) += km*(3*dd[a]*cc[b]-2*cc[a]*dd[b]); // Kxy
	    K(vydof[a], vxdof[b]) += km*(3*cc[a]*dd[b]-2*dd[a]*cc[b]); // Kyx
	    K(vydof[a], vydof[b]) += km*(3*cc[a]*cc[b]+4*dd[a]*dd[b]); // Kyy

	    // -G
	    K(vxdof[a], pdof[b]) = -g*cc[a];
	    K(vydof[a], pdof[b]) = -g*dd[a];

	    // Gt
	    K(pdof[b], vxdof[a]) = gt*cc[a];
	    K(pdof[b], vydof[a]) = gt*dd[a];
	}

	// -G(4,a)
	K(vxdof[3], pdof[a]) = -gb*cc[a];
	K(vydof[3], pdof[a]) = -gb*dd[a];

	// Gt(a,4)
	K(pdof[a], vxdof[3]) = gtb*cc[a];
	K(pdof[a], vydof[3]) = gtb*dd[a];
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

    return K;
}

const Matrix&
PFEMElement2DCompressible::getMassSensitivity(int gradNumber)
{
    // resize K
    K.resize(ndf, ndf);
    K.Zero();

    // M, Mp
    double M = 0.0;
    double M2 = 0.0;
    double M3 = 0.0;
    double M4 = 0.0;
    double Mp = 0.0;

    if (parameterID == 2) {
	// rho
	M = J*thickness/24.0;
	M2 = 2.0*M;
	M3 = 3.0*J*thickness/40.0;
	M4 = 207.0*J*thickness/506.0;

    } else if (parameterID == 6) {
	// thickness
	M =  rho*J/24.0;
	M2 = 2.0*M;
	M3 = 3.0*rho*J/40.0;
	M4 = 207.0*rho*J/506.0;
	Mp = J/6.0/kappa;
    }

    for (int a=0; a<3; a++) {
	// M(a,b) a,b=1,2,3
	for (int b=0; b<3; b++) {
	    if (a == b) {
		K(vxdof[a], vxdof[b]) = M2;
		K(vydof[a], vydof[b]) = M2;
	    } else {
		K(vxdof[a], vxdof[b]) = M;
		K(vydof[a], vydof[b]) = M;
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

    // Mp
    for(int a=0; a<3; a++) {
	K(pdof[a],pdof[a]) = Mp;
    }

    return K;
}

const Vector&
PFEMElement2DCompressible::getResistingForceSensitivity(int gradNumber)
{
    // resize P
    P.resize(ndf);
    P.Zero();

    // F
    double fx = 0.0;
    double fy = 0.0;
    double fbx = 0.0;
    double fby = 0.0;

    if (parameterID == 2) {
	// rho
	fx = J*thickness*b1/6.0;
	fy = J*thickness*b2/6.0;
	fbx = 9.0*J*thickness*b1/40.0;
	fby = 9.0*J*thickness*b2/40.0;

    } else if (parameterID == 6) {
	// thickness
	fx = rho*J*b1/6.0;
	fy = rho*J*b2/6.0;
	fbx = 9.0*rho*J*b1/40.0;
	fby = 9.0*rho*J*b2/40.0;

    } else if (parameterID == 3) {
	// bx
	fx = rho*J*thickness/6.0;
	fbx = 9.0*rho*J*thickness/40.0;

    } else if (parameterID == 4) {
	// by
	fy = rho*J*thickness/6.0;
	fby = 9.0*rho*J*thickness/40.0;
    }


    for (int a=0; a<3; a++) {
	P(vxdof[a]) = -fx;
	P(vydof[a]) = -fy;
    }

    // F4
    P(vxdof[3]) = -fbx;
    P(vydof[3]) = -fby;

    return P;
}

int
PFEMElement2DCompressible::commitSensitivity(int gradNumber, int numGrads)
{
    return 0;
}

// geometric sensitivity
void
PFEMElement2DCompressible::getdM(const Vector& vdot, Matrix& dm) const
{
    dm.resize(8,8);
    dm.Zero();

    // M
    double M = rho*thickness/24.0;
    double M2 = 2.0*M;
    double M3 = 3.0*rho*thickness/40.0;
    double M4 = 207.0*rho*thickness/506.0;
    for (int a=0; a<3; a++) {
	// M(a,b) a,b=1,2,3
	for (int b=0; b<3; b++) {
	    if (a == b) {
		dm(2*a, 2*b) = M2;
		dm(2*a+1, 2*b+1) = M2;
	    } else {
		dm(2*a,2*b) = M;
		dm(2*a+1,2*b+1) = M;
	    }
	}
	// M(4,a)
	dm(6,2*a) = M3;
	dm(7,2*a+1) = M3;

	// M(a,4)
	dm(2*a,6) = M3;
	dm(2*a+1,7) = M3;
    }

    // M(4,4)
    dm(6,6) = M4;
    dm(7,7) = M4;

    // dMdA*vdot
    Vector v = dm*vdot;
    dm.Zero();

    for(int a=0; a<8; a++) {
        for(int b=0; b<3; b++) {
            dm(a,2*b) = v(a)*cc[b];
	    dm(a,2*b+1) = v(a)*dd[b];
        }
    }
}

void
PFEMElement2DCompressible::getdK(const Vector& v, Matrix& dk) const
{
    dk.resize(8,8);
    dk.Zero();

    Vector kv[7];

    // dKdA
    double k = -mu*thickness/(6.*J*J);
    double kb = -27.0*mu*thickness/(40.*J*J);
    for (int a=0; a<3; a++) {
	for (int b=0; b<3; b++) {
	    dk(2*a,2*b) += k*(4*cc[a]*cc[b]+3*dd[a]*dd[b]); // Kxx
	    dk(2*a,2*b+1) += k*(3*dd[a]*cc[b]-2*cc[a]*dd[b]); // Kxy
	    dk(2*a+1,2*b) += k*(3*cc[a]*dd[b]-2*dd[a]*cc[b]); // Kyx
	    dk(2*a+1,2*b+1) += k*(3*cc[a]*cc[b]+4*dd[a]*dd[b]); // Kyy
	}
    }
    double cc2=0.0, dd2=0.0, cd2=0.0;
    for (int a=0; a<3; a++) {
	cc2 += cc[a]*cc[a];
	dd2 += dd[a]*dd[a];
	cd2 += cc[a]*dd[a];
    }
    dk(6,6) += kb*(4*cc2+3*dd2); // Kxx
    dk(6,7) += kb*cd2; // Kxy
    dk(7,6) += kb*cd2; // Kyx
    dk(7,7) += kb*(3*cc2+4*dd2); // Kyy
    kv[6] = dk*v;

    // dKdu
    k = mu*thickness/(6.0*J);
    kb = 27.0*mu*thickness/(40.0*J);

    // dKdu1
    dk.Zero();
    dk(0,2) = dk(2,0) = 3*dd[0]*k;
    dk(0,3) = dk(3,0) = -2*cc[0]*k;
    dk(0,4) = dk(4,0) = -3*dd[0]*k;
    dk(0,5) = dk(5,0) = 2*cc[0]*k;
    dk(1,2) = dk(2,1) = 3*cc[0]*k;
    dk(1,3) = dk(3,1) = 4*dd[0]*k;
    dk(1,4) = dk(4,1) = -3*cc[0]*k;
    dk(1,5) = dk(5,1) = -4*dd[0]*k;
    dk(2,2) = 6*dd[1]*k;
    dk(2,3) = dk(3,2) = cc[1]*k;
    dk(2,4) = dk(4,2) = (3*dd[2]-3*dd[1])*k;
    dk(2,5) = dk(5,2) = (3*cc[2]+2*cc[1])*k;
    dk(3,3) = 8*dd[1]*k;
    dk(3,4) = dk(4,3) = (-3*cc[1]-2*cc[2])*k;
    dk(3,5) = dk(5,3) = (4*dd[2]-4*dd[1])*k;
    dk(4,4) = -6*dd[2]*k;
    dk(4,5) = dk(5,4) = -cc[2]*k;
    dk(5,5) = -8*dd[2]*k;
    dk(6,6) = 6*kb*(dd[1]-dd[2]);
    dk(6,7) = dk(7,6) = kb*(cc[1]-cc[2]);
    dk(7,7) = kb*8*(dd[1]-dd[2]);
    kv[0] = dk*v;

    // dkdu2
    dk.Zero();
    dk(0,2) = dk(2,0) = -4*cc[0]*k;
    dk(0,3) = dk(3,0) = -3*dd[0]*k;
    dk(0,4) = dk(4,0) = 4*cc[0]*k;
    dk(0,5) = dk(5,0) = 3*dd[0]*k;
    dk(1,2) = dk(2,1) = 2*dd[0]*k;
    dk(1,3) = dk(3,1) = -3*cc[0]*k;
    dk(1,4) = dk(4,1) = -2*dd[0]*k;
    dk(1,5) = dk(5,1) = 3*cc[0]*k;
    dk(2,2) = -8*cc[1]*k;
    dk(2,3) = dk(3,2) = -dd[1]*k;
    dk(2,4) = dk(4,2) = (4*cc[1]-4*cc[2])*k;
    dk(2,5) = dk(5,2) = (3*dd[1]+2*dd[2])*k;
    dk(3,3) = -6*cc[1]*k;
    dk(3,4) = dk(4,3) = (-3*dd[2]-2*dd[1])*k;
    dk(3,5) = dk(5,3) = (3*cc[1]-3*cc[2])*k;
    dk(4,4) = 8*cc[2]*k;
    dk(4,5) = dk(5,4) = dd[2]*k;
    dk(5,5) = 6*cc[2]*k;
    dk(6,6) = 8*kb*(cc[2]-cc[1]);
    dk(6,7) = dk(7,6) = kb*(dd[2]-dd[1]);
    dk(7,7) = kb*6*(cc[2]-cc[1]);
    kv[1] = dk*v;

    // dkdu3
    dk.Zero();
    dk(0,0) = -6*k*dd[0];
    dk(0,1) = dk(1,0) = -k*cc[0];
    dk(0,2) = dk(2,0) = -3*dd[1]*k;
    dk(0,3) = dk(3,0) = -3*cc[1]*k;
    dk(0,4) = dk(4,0) = k*(3*dd[0]-3*dd[2]);
    dk(0,5) = dk(5,0) = k*(-3*cc[2]-2*cc[0]);
    dk(1,1) = -8*k*dd[0];
    dk(1,2) = dk(2,1) = 2*cc[1]*k;
    dk(1,3) = dk(3,1) = -4*dd[1]*k;
    dk(1,4) = dk(4,1) = k*(3*cc[0]+2*cc[2]);
    dk(1,5) = dk(5,1) = k*(4*dd[0]-4*dd[2]);
    dk(2,4) = dk(4,2) = 3*k*dd[1];
    dk(2,5) = dk(5,2) = -2*k*cc[1];
    dk(3,4) = dk(4,3) = 3*k*cc[1];
    dk(3,5) = dk(5,3) = 4*k*dd[1];
    dk(4,4) = 6*k*dd[2];
    dk(4,5) = dk(5,4) = k*cc[2];
    dk(5,5) = 8*k*dd[2];
    dk(6,6) = kb*6*(dd[2]-dd[0]);
    dk(6,7) = dk(7,6) = kb*(cc[2]-cc[0]);
    dk(7,7) = 8*kb*(dd[2]-dd[0]);
    kv[2] = dk*v;

    // dkdu4
    dk.Zero();
    dk(0,0) = 8*k*cc[0];
    dk(0,1) = dk(1,0) = k*dd[0];
    dk(0,2) = dk(2,0) = 4*k*cc[1];
    dk(0,3) = dk(3,0) = -2*k*dd[1];
    dk(0,4) = dk(4,0) = (4*cc[2]-4*cc[0])*k;
    dk(0,5) = dk(5,0) = k*(-3*dd[0]-2*dd[2]);
    dk(1,1) = 6*k*cc[0];
    dk(1,2) = dk(2,1) = 3*k*dd[1];
    dk(1,3) = dk(3,1) = 3*k*cc[1];
    dk(1,4) = dk(4,1) = k*(3*dd[2]+2*dd[0]);
    dk(1,5) = dk(5,1) = k*(3*cc[2]-3*cc[0]);
    dk(2,4) = dk(4,2) = -4*k*cc[1];
    dk(2,5) = dk(5,2) = -3*k*dd[1];
    dk(3,4) = dk(4,3) = 2*k*dd[1];
    dk(3,5) = dk(5,3) = -3*k*cc[1];
    dk(4,4) = -8*k*cc[2];
    dk(4,5) = dk(5,4) = -k*dd[2];
    dk(5,5) = -6*k*cc[2];
    dk(6,6) = kb*8*(cc[0]-cc[2]);
    dk(6,7) = kb*(dd[0]-dd[2]);
    dk(7,7) = kb*6*(cc[0]-cc[2]);
    kv[3] = dk*v;

    // dkdu5
    dk.Zero();
    dk(0,0) = 6*k*dd[0];
    dk(0,1) = dk(1,0) = k*cc[0];
    dk(0,2) = dk(2,0) = k*(3*dd[1]-3*dd[0]);
    dk(0,3) = dk(3,0) = k*(3*cc[1]+2*cc[0]);
    dk(0,4) = dk(4,0) = k*3*dd[2];
    dk(0,5) = dk(5,0) = k*3*cc[2];
    dk(1,1) = 8*k*dd[0];
    dk(1,2) = dk(2,1) = k*(-3*cc[0]-2*cc[1]);
    dk(1,3) = dk(3,1) = k*(4*dd[1]-4*dd[0]);
    dk(1,4) = dk(4,1) = -2*k*cc[2];
    dk(1,5) = dk(5,1) = 4*k*dd[2];
    dk(2,2) = -6*k*dd[1];
    dk(2,3) = dk(3,2) = -k*cc[1];
    dk(2,4) = dk(4,2) = -3*k*dd[2];
    dk(2,5) = dk(5,2) = -3*k*cc[2];
    dk(3,3) = -8*k*dd[1];
    dk(3,4) = dk(4,3) = 2*k*cc[2];
    dk(3,5) = dk(5,3) = -4*k*dd[2];
    dk(6,6) = kb*6*(dd[0]-dd[1]);
    dk(6,7) = dk(7,6) = kb*(cc[0]-cc[1]);
    dk(7,7) = kb*8*(dd[0]-dd[1]);
    kv[4] = dk*v;

    // dkdu6
    dk.Zero();
    dk(0,0) = -8*k*cc[0];
    dk(0,1) = dk(1,0) = -k*dd[0];
    dk(0,2) = dk(2,0) = k*(4*cc[0]-4*cc[1]);
    dk(0,3) = dk(3,0) = k*(3*dd[0]+2*dd[1]);
    dk(0,4) = dk(4,0) = -4*k*cc[2];
    dk(0,5) = dk(5,0) = 2*k*dd[2];
    dk(1,1) = -6*k*cc[0];
    dk(1,2) = dk(2,1) = k*(-3*dd[1]-2*dd[0]);
    dk(1,3) = dk(3,1) = k*(3*cc[0]-3*cc[1]);
    dk(1,4) = dk(4,1) = -3*k*dd[2];
    dk(1,5) = dk(5,1) = -3*k*cc[2];
    dk(2,2) = 8*k*cc[1];
    dk(2,3) = dk(3,2) = k*dd[1];
    dk(2,4) = dk(4,2) = 4*k*cc[2];
    dk(2,5) = dk(5,2) = -2*k*dd[2];
    dk(3,3) = 6*k*cc[1];
    dk(3,4) = dk(4,3) = 3*k*dd[2];
    dk(3,5) = dk(5,3) = 3*k*cc[2];
    dk(6,6) = kb*8*(cc[1]-cc[0]);
    dk(6,7) = dk(7,6) = kb*(dd[1]-dd[0]);
    dk(7,7) = kb*6*(cc[1]-cc[0]);
    kv[5] = dk*v;

    // dAdu
    dk.Zero();
    for(int a=0; a<8; a++) {
        for(int b=0; b<3; b++) {
            dk(a,2*b) = kv[6](a)*cc[b]+kv[2*b](a);
	    dk(a,2*b+1) = kv[6](a)*dd[b]+kv[2*b+1](a);
        }
    }
}

void
PFEMElement2DCompressible::getdG(const Vector& p, const Vector& v, Matrix& dg, Matrix& dgt) const
{
    // set C
    Matrix C(6,6);
    double l[12] = {3,5,4,2,5,1,0,4,1,3,2,0};
    for(int b=0; b<6; b++) {
        C(l[2*b],b) = 1.0;
        C(l[2*b+1],b) = -1.0;
    }

    // dG
    dg.resize(8,8);
    dg.Zero();
    double sump = 0.0;
    for(int i=0; i<p.Size(); i++) {
        sump += p(i);
    }
    sump *= thickness/6.0;
    for(int b=0; b<6; b++) {
        dg(l[2*b],b) = sump;
        dg(l[2*b+1],b) = -sump;
    }

    // dGb
    double factor = -9.0*thickness/40.0;
    for(int a=0; a<2; a++) {
        for(int b=0; b<6; b++) {
            for(int i=0; i<3; i++) {
                dg(a+6,b) += C(2*i+a,b)*p(i)*factor;
            }
        }
    }

    // dGt
    dgt.resize(3,8);
    dgt.Zero();
    Matrix vt3(3,6);
    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            vt3(a,b) = v(b);
        }
    }
    Matrix temp = vt3*C;
    temp *= thickness/6.0;
    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            dgt(a,b) = temp(a,b);
        }
    }

    // dGbt
    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            for(int i=0; i<2; i++) {
                dgt(a,b) += C(2*a+i,b)*v(i+6)*factor;
            }
        }
    }
}

void
PFEMElement2DCompressible::getdF(Matrix& df) const
{
    df.resize(8,8);
    df.Zero();

    // dF
    double factor = rho*thickness/6.0;
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            df(2*a,2*b) = b1*cc[b]*factor;
	    df(2*a,2*b+1) = b1*dd[b]*factor;
            df(2*a+1,2*b) = b2*cc[b]*factor;
	    df(2*a+1,2*b+1) = b2*dd[b]*factor;
        }
    }

    // dFb
    factor = rho*thickness*9.0/40.0;
    for(int b=0; b<3; b++) {
        df(6,2*b) = b1*cc[b]*factor;
        df(7,2*b) = b2*cc[b]*factor;
	df(6,2*b+1) = b1*dd[b]*factor;
        df(7,2*b+1) = b2*dd[b]*factor;
    }
}

void
PFEMElement2DCompressible::getdMp(const Vector& pdot, Matrix& dmp) const
{
    dmp.resize(3,8);
    dmp.Zero();

    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            dmp(a,2*b) = pdot(a)*cc[b];
	    dmp(a,2*b+1) = pdot(a)*dd[b];
        }
    }
    dmp *= thickness/6.0/kappa;
}
