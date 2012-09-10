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
#include <Node.h>

Matrix PFEMElement2D::K(15,15);
Vector PFEMElement2D::P(15);

PFEMElement2D::PFEMElement2D(int tag, int nd1, int nd2, int nd3,
                             double r, double m, double b1, double b2)
    :Element(tag, ELE_TAG_PFEMElement2D), ntags(3),
     rho(r), mu(m), bx(b1), by(b2)
{
    ntags[0]=nd1; ntags[1]=nd2; ntags[2]=nd3;
    for(int i=0;i<3;i++)
    {
        nodes[i] = 0;
    }
}


PFEMElement2D::~PFEMElement2D()
{

}



int
PFEMElement2D::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement2D::getExternalNodes()
{
    return ntags;
}

Node **
PFEMElement2D::getNodePtrs(void) 
{
    return nodes;
}

int
PFEMElement2D::getNumDOF()
{
    return P.Size();
}

int
PFEMElement2D::revertToLastCommit()
{
    return 0;
}

int PFEMElement2D::commitState()
{
    return this->Element::commitState();
}

int
PFEMElement2D::update()
{
    return 0;
}

const Matrix&
PFEMElement2D::getMass()
{
    K.Zero();

    // get nodal coordinates 
    double x[3], y[3];
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[i]->getCrds();
        const Vector& disp = nodes[i]->getTrialDisp();
        x[i] = coord[0] + disp[0];
        y[i] = coord[1] + disp[1];
    }

    // get Jacobi
    double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    double J2 = J/2.;

    // get derivatives
    double pts[3][3] = {{0.5,0.5,0}, {0.5,0,0.5}, {0,0.5,0.5}};
    
    // lumped mass 
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            double m = rho*(pts[0][a]*pts[0][b]+pts[1][a]*pts[1][b]+pts[2][a]*pts[2][b])/3.*J2;
            K(5*a,5*a) += m;          // Mxd
            K(5*a+1,5*a+1) += m;      // Myd
        }
    }

    return K;
}


const Matrix&
PFEMElement2D::getDamp()
{
    this->getDampWithK();
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            K(5*a+1, 5*b+1) = K(5*a, 5*b) = 0.0; // no K
        }
    }


    return K;
}

const Matrix&
PFEMElement2D::getDampWithK()
{
    K.Zero();

    // get nodal coordinates 
    double x[3], y[3];
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[i]->getCrds();
        const Vector& disp = nodes[i]->getTrialDisp();
        x[i] = coord[0] + disp[0];
        y[i] = coord[1] + disp[1];
    }

    // get Jacobi
    double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    double J2 = J/2.;

    // get derivatives
    double dNdx[3] = {(y[1]-y[2])/J, (y[2]-y[0])/J, (y[0]-y[1])/J};
    double dNdy[3] = {(x[2]-x[1])/J, (x[0]-x[2])/J, (x[1]-x[0])/J};

    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            K(5*a+1, 5*b+1) = K(5*a, 5*b) = mu*J2*(dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]); // K
            K(5*a+2, 5*b) = dNdx[b]/3.*J2; // GxT
            K(5*a+2, 5*b+1) = dNdy[b]/3.*J2; // GyT
        }
    }


    return K;
}

const Matrix&
PFEMElement2D::getTangentStiff()
{
    K.Zero();

    // get nodal coordinates 
    double x[3], y[3];
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[i]->getCrds();
        const Vector& disp = nodes[i]->getTrialDisp();
        x[i] = coord[0] + disp[0];
        y[i] = coord[1] + disp[1];
    }

    // get Jacobi
    double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    double J2 = J/2.;
    double tau = 1./(rho/ops_Dt+8*mu/(3*4*J2));

    // get derivatives
    double dNdx[3] = {(y[1]-y[2])/J, (y[2]-y[0])/J, (y[0]-y[1])/J};
    double dNdy[3] = {(x[2]-x[1])/J, (x[0]-x[2])/J, (x[1]-x[0])/J};
    double pts[3][3] = {{0.5,0.5,0}, {0.5,0,0.5}, {0,0.5,0.5}};

    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            K(5*a, 5*b+2) = -dNdx[a]/3.*J2; // -Gx
            K(5*a+1, 5*b+2) = -dNdy[a]/3.*J2; // -Gy
            K(5*a+2, 5*b+2) = tau*J2*(dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]); // L
            K(5*b+3, 5*a+2) = K(5*a+2, 5*b+3) = tau*dNdx[a]/3.*J2;   // QxT, Qx
            K(5*b+4, 5*a+2) = K(5*a+2, 5*b+4) = tau*dNdy[a]/3.*J2;   // QyT ,Qy
            double m = tau*(pts[0][a]*pts[0][b]+pts[1][a]*pts[1][b]+pts[2][a]*pts[2][b])/3.*J2;
            K(5*a+3, 5*a+3) += m;  // Mhatxd
            K(5*a+4, 5*a+4) += m;  // Mhatyd
        }
    }

    return K;
}


const Matrix& 
PFEMElement2D::getInitialStiff()
{
    return getTangentStiff();
}

int 
PFEMElement2D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2D::getResistingForce()
{
    P.Zero();
    return P;
}

const Vector&
PFEMElement2D::getResistingForceIncInertia()
{
    P.Zero();

    double x[3], y[3];
    Vector a(15), v(15), u(15);
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[i]->getCrds();
        const Vector& accel = nodes[i]->getTrialAccel();
        const Vector& vel = nodes[i]->getTrialVel();
        const Vector& disp = nodes[i]->getTrialDisp();
        x[i] = coord[0] + disp[0];
        y[i] = coord[1] + disp[1];
        for(int j=0; j<5; j++) {
            a[5*i+j] = accel[j];
            v[5*i+j] = vel[j];
            u[5*i+j] = disp[j];
        }
    }

    P.addMatrixVector(1.0, this->getMass(), a, 1.0);
    P.addMatrixVector(1.0, this->getDampWithK(), v, 1.0);
    P.addMatrixVector(1.0, this->getTangentStiff(), u, 1.0);
 
    // get Jacobi
    double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    double J2 = J/2.;

    // -F
    for(int i=0; i<3; i++) {
        P(5*i) -= rho*bx/3.*J2;
        P(5*i+1) -= rho*by/3.*J2;
    }


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
    return 0;
}

int
PFEMElement2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement2D::setDomain(Domain *domain)
{
    this->DomainComponent::setDomain(domain);

    if(domain == 0)
    {
        opserr<<"WARINING: domain is not available -- PFEMElement2D::setDomain\n";
        return;
    }

    int eletag = this->getTag();
    for(int i=0; i<3; i++) {
        nodes[i] = domain->getNode(ntags[i]);
        if(nodes[i] == 0) {
            opserr<<"WARNING: node "<<ntags[i]<<" does not exist ";
            opserr<<"in PFEMElement2D - setDomain() "<<eletag<<"\n";
            continue;
        }
        if(nodes[i]->getNumberDOF() < 5) {
            opserr<<"WARNING: wrong number of dof for node "<<ntags[i]<<" ";
            opserr<<"in PFEMElement2D - setDomain() "<<eletag<<"\n";
        }

    }

}

void
PFEMElement2D::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2D: "<<this->getTag()<<endln;
    s << "     Node tags  :  ";
    for(int i=0;i<3;i++)
    {
        s<<ntags(i)<<" ";
    }
    s << endln;
    s << "     Nodal ndf  :  ";
    for(int i=0;i<3;i++)
    {
        s<<nodes[i]->getNumberDOF()<<" ";
    }
    s << endln;
}
