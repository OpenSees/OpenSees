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
// $Date: 2013/01/29 13:41:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement3D.h,v $
                                                                        
// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2013
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement3D.

#include "PFEMElement3D.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <cmath>
#include <stdlib.h>

Matrix PFEMElement3D::K;
Vector PFEMElement3D::P;

PFEMElement3D::PFEMElement3D(int tag, int nd1, int nd2, int nd3, int nd4,
                             double r, double m, double b1, double b2, double b3)
    :Element(tag, ELE_TAG_PFEMElement3D), ntags(8), 
     rho(r), mu(m), bx(b1), by(b2), bz(b3), J(0.0), numDOFs()
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3; ntags(6)=nd4;
    for(int i=0;i<4;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        ntags(2*i+1) = ntags(2*i);
        thePCs[i] = 0;
        dNdx[i] = 0.0;
        dNdy[i] = 0.0;
        dNdz[i] = 0.0;
    }
}


PFEMElement3D::~PFEMElement3D()
{
    for(int i=0; i<4; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}



int
PFEMElement3D::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement3D::getExternalNodes()
{
    if(numDOFs.Size()==0) return numDOFs;
    return ntags;
}

Node **
PFEMElement3D::getNodePtrs(void) 
{
    return nodes;
}

int
PFEMElement3D::getNumDOF()
{
    return numDOFs(numDOFs.Size()-1);
}

int
PFEMElement3D::revertToLastCommit()
{
    return 0;
}

int PFEMElement3D::commitState()
{
    return Element::commitState();
}

int
PFEMElement3D::update()
{
    //  get nodal coordinates 
    Matrix A(4,4);
    for(int i=0; i<4; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
        const Vector& disp = nodes[2*i]->getTrialDisp();
        A(i,0) = 1.0;
        for(int j=1; j<4; j++) {
            A(i,j) = coord(j-1) + disp(j-1);
        }
    }

    // get jacobi
    Matrix coef(4,4), sub(3,3);
    Vector jacobi(4);
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {

            // sub matrix
            for(int k=0; k<4; k++) {
                if(k==i) continue;
                int k1 = k;
                if(k>i) k1--;
                for(int l=0; l<4; l++) {
                    if(l==j) continue;
                    int l1 = l;
                    if(l>j) l1--;
                    sub(k1,l1) = A(k,l);
                }
            }

            // cofactor
            int sign = 1;
            if((i+j)%2 == 1) sign = -1;
            coef(i,j) = sign*det33(sub);
            jacobi(i) += coef(i,j) * A(i,j);
        }
    }

    // opserr<<"J = "<<jacobi;
    J = jacobi(0);

    // if(fabs(J) <= 1e-6) {
    //     opserr<<"WARNING: zero Jacobian for element ";
    //     opserr<<J<<" --PFEMElement3D::update\n";
    //     //return -1;
    // }

    // get derivatives
    for(int i=0; i<4; i++) {
        dNdx[i] = coef(i,1)/J;
        dNdy[i] = coef(i,2)/J;
        dNdz[i] = coef(i,3)/J;
        //opserr<<"dNdx = "<<dNdx[i]<<"\n";
    }

    return 0;
}

const Matrix&
PFEMElement3D::getMass()
{

    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    // mass 
    for(int a=0; a<4; a++) {
        double m = rho*J/24.0;
        K(numDOFs(2*a), numDOFs(2*a)) = m;          // Mxd
        K(numDOFs(2*a)+1, numDOFs(2*a)+1) = m;      // Myd
        K(numDOFs(2*a)+2, numDOFs(2*a)+2) = m;      // Mzd
    }

    return K;
}

const Matrix&
PFEMElement3D::getDamp()
{

    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    double V = fabs(J/6.0);
    double h = exp(log(V)/3.0);     // cubic root of V
    double tau = 1./(rho/ops_Dt+8*mu/(3*h*h));

    for(int a=0; a<4; a++) {
        for(int b=0; b<4; b++) {

            K(numDOFs(2*a), numDOFs(2*b+1)) = -dNdx[a]/24.*J;   // -Gx
            K(numDOFs(2*a)+1, numDOFs(2*b+1)) = -dNdy[a]/24.*J; // -Gy
            K(numDOFs(2*a)+2, numDOFs(2*b+1)) = -dNdz[a]/24.*J; // -Gz

            K(numDOFs(2*a+1), numDOFs(2*b)) = dNdx[b]/24.*J;    // GxT
            K(numDOFs(2*a+1), numDOFs(2*b)+1) = dNdy[b]/24.*J;  // GyT
            K(numDOFs(2*a+1), numDOFs(2*b)+2) = dNdz[b]/24.*J;  // GzT

            K(numDOFs(2*a+1), numDOFs(2*b+1)) 
                = tau*J/6.0*(dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b]+dNdz[a]*dNdz[b]); // L

            K(numDOFs(2*a+1), numDOFs(2*b+1)+1) = tau*dNdx[a]/24.*J;   // Qx
            K(numDOFs(2*a+1), numDOFs(2*b+1)+2) = tau*dNdy[a]/24.*J;   // Qy
            K(numDOFs(2*a+1), numDOFs(2*b+1)+3) = tau*dNdz[a]/24.*J;   // Qz

            K(numDOFs(2*a+1)+1, numDOFs(2*b+1)) = tau*dNdx[b]/24.*J;   // QxT
            K(numDOFs(2*a+1)+2, numDOFs(2*b+1)) = tau*dNdy[b]/24.*J;   // QyT
            K(numDOFs(2*a+1)+3, numDOFs(2*b+1)) = tau*dNdz[b]/24.*J;   // QzT

        }
        double m = tau*J/24.0;
        K(numDOFs(2*a+1)+1, numDOFs(2*a+1)+1) = m; // Mhatxd
        K(numDOFs(2*a+1)+2, numDOFs(2*a+1)+2) = m; // Mhatyd
        K(numDOFs(2*a+1)+3, numDOFs(2*a+1)+3) = m; // Mhatzd
    }

    return K;
}

const Matrix&
PFEMElement3D::getDampWithK()
{

    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    double V = fabs(J/6.0);
    double h = exp(log(V)/3.0);     // cubic root of V
    double tau = 1./(rho/ops_Dt+8*mu/(3*h*h));

    for(int a=0; a<4; a++) {
        for(int b=0; b<4; b++) {

            // K(1,1:3)
            K(numDOFs(2*a), numDOFs(2*b)) 
                = mu*J/6.0*(2*dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b]+dNdz[a]*dNdz[b]);
            K(numDOFs(2*a), numDOFs(2*b)+1) = mu*J/6.0*dNdy[a]*dNdx[b];
            K(numDOFs(2*a), numDOFs(2*b)+2) = mu*J/6.0*dNdz[a]*dNdx[b];

            // K(2,1:3)
            K(numDOFs(2*a)+1, numDOFs(2*b)) = mu*J/6.0*dNdx[a]*dNdy[b];
            K(numDOFs(2*a)+1, numDOFs(2*b)+1) 
                = mu*J/6.0*(2*dNdy[a]*dNdy[b]+dNdx[a]*dNdx[b]+dNdz[a]*dNdz[b]);
            K(numDOFs(2*a)+1, numDOFs(2*b)+2) = mu*J/6.0*dNdz[a]*dNdy[b];

            // K(3,1:3)
            K(numDOFs(2*a)+2, numDOFs(2*b)) = mu*J/6.0*dNdx[a]*dNdz[b];
            K(numDOFs(2*a)+2, numDOFs(2*b)+1) = mu*J/6.0*dNdy[a]*dNdz[b];
            K(numDOFs(2*a)+2, numDOFs(2*b)+2) 
                = mu*J/6.0*(2*dNdz[a]*dNdz[b]+dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b]);

            K(numDOFs(2*a), numDOFs(2*b+1)) = -dNdx[a]/24.*J;   // -Gx
            K(numDOFs(2*a)+1, numDOFs(2*b+1)) = -dNdy[a]/24.*J; // -Gy
            K(numDOFs(2*a)+2, numDOFs(2*b+1)) = -dNdz[a]/24.*J; // -Gz

            K(numDOFs(2*a+1), numDOFs(2*b)) = dNdx[b]/24.*J;    // GxT
            K(numDOFs(2*a+1), numDOFs(2*b)+1) = dNdy[b]/24.*J;  // GyT
            K(numDOFs(2*a+1), numDOFs(2*b)+2) = dNdz[b]/24.*J;  // GzT

            K(numDOFs(2*a+1), numDOFs(2*b+1)) 
                = tau*J/6.0*(dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b]+dNdz[a]*dNdz[b]); // L

            K(numDOFs(2*a+1), numDOFs(2*b+1)+1) = tau*dNdx[a]/24.*J;   // Qx
            K(numDOFs(2*a+1), numDOFs(2*b+1)+2) = tau*dNdy[a]/24.*J;   // Qy
            K(numDOFs(2*a+1), numDOFs(2*b+1)+3) = tau*dNdz[a]/24.*J;   // Qz

            K(numDOFs(2*a+1)+1, numDOFs(2*b+1)) = tau*dNdx[b]/24.*J;   // QxT
            K(numDOFs(2*a+1)+2, numDOFs(2*b+1)) = tau*dNdy[b]/24.*J;   // QyT
            K(numDOFs(2*a+1)+3, numDOFs(2*b+1)) = tau*dNdz[b]/24.*J;   // QzT

        }
        double m = tau*J/24.0;
        K(numDOFs(2*a+1)+1, numDOFs(2*a+1)+1) = m; // Mhatxd
        K(numDOFs(2*a+1)+2, numDOFs(2*a+1)+2) = m; // Mhatyd
        K(numDOFs(2*a+1)+3, numDOFs(2*a+1)+3) = m; // Mhatzd
    }

    return K;
}

const Matrix&
PFEMElement3D::getTangentStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}


const Matrix& 
PFEMElement3D::getInitialStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int 
PFEMElement3D::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement3D::getResistingForce()
{

    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement3D::getResistingForceIncInertia()
{
    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    // get velocity, acceleration
    Vector v(ndf), vdot(ndf);
    for(int i=0; i<4; i++) {
        const Vector& accel = nodes[2*i]->getTrialAccel();
        vdot(numDOFs(2*i)) = accel(0);
        vdot(numDOFs(2*i)+1) = accel(1);
        vdot(numDOFs(2*i)+2) = accel(2);

        const Vector& vel = nodes[2*i]->getTrialVel();
        v(numDOFs(2*i)) = vel(0);
        v(numDOFs(2*i)+1) = vel(1);
        v(numDOFs(2*i)+2) = vel(2);

        const Vector& vel2 = nodes[2*i+1]->getTrialVel();
        v(numDOFs(2*i+1)) = vel2(0);
        v(numDOFs(2*i+1)+1) = vel2(1);
        v(numDOFs(2*i+1)+2) = vel2(2);
        v(numDOFs(2*i+1)+3) = vel2(3);
    }

    double d = v.Norm();
    if(d!=d) opserr<<"v "<<this->getTag()<<"\n";
    d = vdot.Norm();
    if(d!=d) opserr<<"vdot "<<this->getTag()<<"\n";

    P.addMatrixVector(1.0, getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, getDampWithK(), v, 1.0);

    Vector ones(ndf); 
    for(int i=0;i<ndf;i++) ones(i) = 1.0;
    Vector Q1(ndf), Q2(ndf);
    Q1.addMatrixVector(1.0, getMass(), ones, 1.0);
    Q2.addMatrixVector(1.0, getDampWithK(), ones, 1.0);
    d = Q1.Norm();
    if(d!=d) opserr<<"mass "<<this->getTag()<<"\n";
    d = Q2.Norm();
    if(d!=d) {
        opserr<<"damp "<<this->getTag()<<"\n";
        exit(1);
    }

    // get Jacobi
    for(int i=0; i<4; i++) {
        P(numDOFs(2*i)) -= rho*bx/24.*J;
        P(numDOFs(2*i)+1) -= rho*by/24.*J;
        P(numDOFs(2*i)+2) -= rho*bz/24.*J;
    }


    return P;
}


const char*
PFEMElement3D::getClassType()const
{
    return "PFEMElement3D";
}

int
PFEMElement3D::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PFEMElement3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PFEMElement3D::setDomain(Domain *theDomain)
{
    numDOFs.resize(9);
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    numDOFs.Zero();
    int ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<4; i++) {

        // set ndf
        numDOFs(2*i) = ndf;

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement3D - setDomain() "<<eletag<<"\n ";
            return;
        }
        ndf += nodes[2*i]->getNumberDOF();
 
        // set ndf
        numDOFs(2*i+1) = ndf;

        // get pc 
        int pndf = 3;
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] == 0) {
            thePCs[i] = new Pressure_Constraint(ntags(2*i), pndf);
            if(thePCs[i] == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMElement3D::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement3D::setDomain "<<eletag<<"\n";
                delete thePCs[i];
                thePCs[i] = 0;
                return;
            }
        }

        // connect
        thePCs[i]->connect(eletag);

        // get pressure node
        // ntags(2*i+1) = thePCs[i]->getPressureNode();
        // nodes[2*i+1] = theDomain->getNode(ntags(2*i+1));
        nodes[2*i+1] = thePCs[i]->getPressureNode();
        if(nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in PFEMElement3D - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
        ndf += nodes[2*i+1]->getNumberDOF();
    }
    numDOFs(numDOFs.Size()-1) = ndf;
}

void
PFEMElement3D::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement3D: "<<this->getTag()<<endln;
}


double
PFEMElement3D::det33(const Matrix& A)
{
    return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
        -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
        +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
}

int 
PFEMElement3D::displaySelf(Renderer &, int mode, float fact, const char **displayModes, int numModes)
{
  return 0;
}
