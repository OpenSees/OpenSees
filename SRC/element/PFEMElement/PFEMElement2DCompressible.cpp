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
                                                                        
// Written: Minjie Zhu (zhum@engr.orst.edu)
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
#include <PFEMAnalysis.h>

Matrix PFEMElement2DCompressible::K;
Vector PFEMElement2DCompressible::P;


// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2DCompressible::PFEMElement2DCompressible()
    :Element(0, ELE_TAG_PFEMElement2DCompressible), ntags(6), 
     rho(0), mu(0), bx(0), by(0), J(0.0), 
     numDOFs(),thickness(1.0), kappa(0.0)
{
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        thePCs[i] = 0;
        dNdx[i] = 0.0;
        dNdy[i] = 0.0;
    }
}

// for object
PFEMElement2DCompressible::PFEMElement2DCompressible(int tag, int nd1, int nd2, int nd3,
                             double r, double m, double b1, 
                             double b2, double thk, double ka)
    :Element(tag, ELE_TAG_PFEMElement2DCompressible), ntags(6), 
     rho(r), mu(m), bx(b1), by(b2), J(0.0), 
     numDOFs(), thickness(thk), kappa(ka)
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        ntags(2*i+1) = ntags(2*i);
        thePCs[i] = 0;
        dNdx[i] = 0.0;
        dNdy[i] = 0.0;
    }
}


PFEMElement2DCompressible::~PFEMElement2DCompressible()
{
    for(int i=0; i<3; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}



int
PFEMElement2DCompressible::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement2DCompressible::getExternalNodes()
{
    if(numDOFs.Size() == 0) return numDOFs;
    return ntags;
}

Node **
PFEMElement2DCompressible::getNodePtrs(void) 
{
    return nodes;
}

int
PFEMElement2DCompressible::getNumDOF()
{
    if(numDOFs.Size() == 0) return 0;
    return numDOFs(numDOFs.Size()-1);
}

int
PFEMElement2DCompressible::revertToLastCommit()
{
    return 0;
}

int PFEMElement2DCompressible::commitState()
{
    
    return Element::commitState();
}

int
PFEMElement2DCompressible::update()
{
    // get nodal coordinates 
    double x[3], y[3];
    for(int i=0; i<3; i++) {
        const Vector& coord = nodes[2*i]->getCrds();
        const Vector& disp = nodes[2*i]->getTrialDisp();
        x[i] = coord[0] + disp[0];
        y[i] = coord[1] + disp[1];
    }

    // get Jacobi
    J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    J *= thickness;

    // get derivatives
    double dndx[3] = {(y[1]-y[2])/J, (y[2]-y[0])/J, (y[0]-y[1])/J};
    double dndy[3] = {(x[2]-x[1])/J, (x[0]-x[2])/J, (x[1]-x[0])/J};
    for(int i=0; i<3; i++) {
        dNdx[i] = dndx[i];
        dNdy[i] = dndy[i];
    }

    return 0;
}

const Matrix&
PFEMElement2DCompressible::getMass()
{

    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    double J2 = J/2.;
    
    // mass 
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            double m = rho*J2/12.0;
            if(a == b) m *= 2.0;
            K(numDOFs(2*a), numDOFs(2*b)) = m;          // Mx
            K(numDOFs(2*a)+1, numDOFs(2*b)+1) = m;      // My

            m = J2/12.0/kappa;
            if(a == b) m *= 2.0;
            K(numDOFs(2*a+1), numDOFs(2*a+1)) += m;   // Mp
        }
    }

    return K;
}

const Matrix&
PFEMElement2DCompressible::getDamp()
{

    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    double J2 = J/2.;
    double lambda = -2*mu/3.0;

    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {

            K(numDOFs(2*a), numDOFs(2*b)) += mu*J2*(2*dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]); // K1
            K(numDOFs(2*a), numDOFs(2*b)+1) += mu*J2*dNdy[a]*dNdx[b]; // K1
            K(numDOFs(2*a)+1, numDOFs(2*b)) += mu*J2*dNdx[a]*dNdy[b]; // K1
            K(numDOFs(2*a)+1, numDOFs(2*b)+1) += mu*J2*(2*dNdy[a]*dNdy[b] + dNdx[a]*dNdx[b]); // K1

            K(numDOFs(2*a), numDOFs(2*b)) += lambda*J2*dNdx[a]*dNdx[b]; // K2
            K(numDOFs(2*a), numDOFs(2*b)+1) += lambda*J2*dNdx[a]*dNdy[b]; // K2
            K(numDOFs(2*a)+1, numDOFs(2*b)) += lambda*J2*dNdy[a]*dNdx[b]; // K2
            K(numDOFs(2*a)+1, numDOFs(2*b)+1) += lambda*J2*dNdy[a]*dNdy[b]; // K2

            K(numDOFs(2*a), numDOFs(2*b+1)) += -dNdx[a]/3.*J2;   // -Gx
            K(numDOFs(2*a)+1, numDOFs(2*b+1)) += -dNdy[a]/3.*J2; // -Gy
            K(numDOFs(2*a+1), numDOFs(2*b)) += dNdx[b]/3.*J2;   // GxT
            K(numDOFs(2*a+1), numDOFs(2*b)+1) += dNdy[b]/3.*J2; // GyT

            // K(numDOFs(2*a+1), numDOFs(2*b+1)) += Lfact*tau*J2*(dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]); // L

            // K(numDOFs(2*a+1), numDOFs(2*b+1)+1) = tau*dNdx[a]/3.*J2;   // Qx
            // K(numDOFs(2*a+1), numDOFs(2*b+1)+2) = tau*dNdy[a]/3.*J2;   // Qy

            // K(numDOFs(2*a+1)+1, numDOFs(2*b+1)) = tau*dNdx[b]/3.*J2;   // QxT
            // K(numDOFs(2*a+1)+2, numDOFs(2*b+1)) = tau*dNdy[b]/3.*J2;   // QyT
        }
        // double m = tau*J2/3.0;
        // K(numDOFs(2*a+1)+1, numDOFs(2*a+1)+1) = m; // Mhatxd
        // K(numDOFs(2*a+1)+2, numDOFs(2*a+1)+2) = m; // Mhatyd
    }

    return K;
}

const Matrix&
PFEMElement2DCompressible::getTangentStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();


    return K;
}


const Matrix& 
PFEMElement2DCompressible::getInitialStiff()
{
    // resize K
    int ndf = this->getNumDOF();
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
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2DCompressible::getResistingForceIncInertia()
{

    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    // get velocity, accleration
    Vector v(ndf), vdot(ndf);
    for(int i=0; i<3; i++) {
        const Vector& accel = nodes[2*i]->getTrialAccel();
        vdot(numDOFs(2*i)) = accel(0);
        vdot(numDOFs(2*i)+1) = accel(1);

        const Vector& accel2 = nodes[2*i+1]->getTrialAccel();
        vdot(numDOFs(2*i+1)) = accel2(0);

        const Vector& vel = nodes[2*i]->getTrialVel();
        v(numDOFs(2*i)) = vel(0);
        v(numDOFs(2*i)+1) = vel(1);

        const Vector& vel2 = nodes[2*i+1]->getTrialVel();
        v(numDOFs(2*i+1)) = vel2(0);
    }

    // Ma+k1v-Gp
    // Gt*v+Mp*p
    P.addMatrixVector(1.0, getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, getDamp(), v, 1.0);

    // f
    double J2 = J/2.0;
    for(int i=0; i<3; i++) {
        P(numDOFs(2*i)) -= rho*bx/3.*J2;
        P(numDOFs(2*i)+1) -= rho*by/3.*J2;
    }


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
    int res = 0;
    int dataTag = this->getDbTag();

    // send vector
    static Vector data(25);
    data(0) = this->getTag();
    data(1) = rho;
    data(2) = mu;
    data(3) = bx;
    data(4) = by;
    for(int i=0; i<3; i++) {
        data(5+i) = dNdx[i];
        data(8+i) = dNdy[i];
    }
    data(11) = J;
    for(int i=0; i<6; i++) {
        data(12+i) = ntags(i);
        data(18+i) = numDOFs(i);
    }
    data(24) = numDOFs(6);

    res = theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
        opserr<<"WARNING: PFEMElement2DCompressible::sendSelf - "<<this->getTag()<<" failed to send vector\n";
        return -1;
    }


    return 0;
}

int
PFEMElement2DCompressible::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res;
    int dataTag = this->getDbTag();

    // receive vector
    static Vector data(12);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
        opserr<<"WARNING: PFEMElement2DCompressible::recvSelf - failed to receive vector\n";
        return -1;
    }
    this->setTag((int)data(0));
    rho = data(1);
    mu = data(2);
    bx = data(3);
    by = data(4);
    for(int i=0; i<3; i++) {
        dNdx[i] = data(5+i);
        dNdy[i] = data(8+i);
    }
    J = data(11);
    for(int i=0; i<6; i++) {
        ntags(i) = (int)data(12+i);
        numDOFs(i) = (int)data(18+i);
    }
    numDOFs(6) = (int)data(24);

    return 0;
}

void
PFEMElement2DCompressible::setDomain(Domain *theDomain)
{
    numDOFs.resize(7);
    this->DomainComponent::setDomain(theDomain);

    if(theDomain == 0) {
        return;
    }

    numDOFs.Zero();
    int ndf = 0;
    int eletag = this->getTag();
    for(int i=0; i<3; i++) {

        // set ndf
        numDOFs(2*i) = ndf;

        // get node
        nodes[2*i] = theDomain->getNode(ntags(2*i));
        if(nodes[2*i] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i)<<" does not exist ";
            opserr<<"in PFEMElement2DCompressible - setDomain() "<<eletag<<"\n ";
            return;
        }
        ndf += nodes[2*i]->getNumberDOF();
 
        // set ndf
        numDOFs(2*i+1) = ndf;

        // get pc 
        int pndf = 1;
        thePCs[i] = theDomain->getPressure_Constraint(ntags(2*i));
        if(thePCs[i] != 0) {
            thePCs[i]->setDomain(theDomain);
        } else {
            thePCs[i] = new Pressure_Constraint(ntags(2*i), by, pndf);
            if(thePCs[i] == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMElement2DCompressible::setDomain "<<eletag<<"\n";
                return;
            }
            if(theDomain->addPressure_Constraint(thePCs[i]) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMElement2DCompressible::setDomain "<<eletag<<"\n";
                delete thePCs[i];
                thePCs[i] = 0;
                return;
            }
        }

        // set gravity
        thePCs[i]->setGravity(by);

        // connect
        thePCs[i]->connect(eletag);

        // get pressure node
        ntags(2*i+1) = thePCs[i]->getPressureNode();
        nodes[2*i+1] = theDomain->getNode(ntags(2*i+1));
        if(nodes[2*i+1] == 0) {
            opserr<<"WARNING: node "<<ntags(2*i+1)<<" does not exist ";
            opserr<<"in PFEMElement2DCompressible - setDomain() "<<eletag<<"\n ";
            return;
        }
        ndf += nodes[2*i+1]->getNumberDOF();
    }
    numDOFs(numDOFs.Size()-1) = ndf;
}

void
PFEMElement2DCompressible::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2DCompressible: "<<this->getTag()<<endln;
}

int
PFEMElement2DCompressible::displaySelf(Renderer &theViewer, int displayMode, float fact)                                          
{

    // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

    /*
  static Vector values(numgp);

    for (int j=0; j<numgp; j++) values(j) = 0.0;

    if (displayMode < numgp && displayMode > 0) {
		for (int i=0; i<numgp; i++) {
			const Vector &stress = theMaterial[i]->getStress();
	        values(i) = stress(displayMode-1);
	    }
    }
    */

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
