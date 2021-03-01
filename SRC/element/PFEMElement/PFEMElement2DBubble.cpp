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
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2DBubble.h,v $

// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2DBubble.

#include "PFEMElement2DBubble.h"
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

Matrix PFEMElement2DBubble::K;
Vector PFEMElement2DBubble::P;
Matrix PFEMElement2DBubble::C;


bool PFEMElement2DBubble::dispon = true;

void* OPS_PFEMElement2DBubble(const ID &info)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return 0;
    }

    int idata[4];
    double data[6] = {0,0,0,0,1.0,-1};
    int numdata;

    // regular element, not in a mesh, get tags
    if (info.Size() == 0) {
        numdata = OPS_GetNumRemainingInputArgs();
        if(numdata < 4) {
            opserr<<"WARNING: insufficient number of arguments: tag, nd1, nd2, nd3\n";
            return 0;
        }

        // tag, nd1, nd2, nd3
        numdata = 4;
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
        if (info.Size() < 6) {
            opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2, nd3\n";
            return 0;
        }

        // get the data for a mesh
        Vector& mdata = meshdata[info(1)];
        if (mdata.Size() < 6) return 0;

        idata[0] = info(2);
        for (int i=0; i<3; ++i) {
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

    return new PFEMElement2DBubble(idata[0],idata[1],idata[2],idata[3],
                                   data[0],data[1],data[2],data[3],data[4],data[5]);
}

// for FEM_ObjectBroker, recvSelf must invoke
PFEMElement2DBubble::PFEMElement2DBubble()
        :Element(0, ELE_TAG_PFEMElement2DBubble), ntags(6),
         rho(0), mu(0), bx(0), by(0), J(0.0), dJ(6),
         numDOFs(),thickness(1.0), kappa(-1), parameterID(0),
         M(), D(), F(), Fp()
{
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        thePCs[i] = 0;
    }

    setC();
}

// for object
PFEMElement2DBubble::PFEMElement2DBubble(int tag, int nd1, int nd2, int nd3,
                                         double r, double m, double b1, double b2,
                                         double thk, double ka)
        :Element(tag, ELE_TAG_PFEMElement2DBubble), ntags(6),
         rho(r), mu(m), bx(b1), by(b2), J(0.0), dJ(6), numDOFs(),
         thickness(thk), kappa(ka), parameterID(0),
         M(), D(), F(), Fp()
{
    ntags(0)=nd1; ntags(2)=nd2; ntags(4)=nd3;
    ntags(1)=nd1; ntags(3)=nd2; ntags(5)=nd3;
    for(int i=0;i<3;i++)
    {
        nodes[2*i] = 0;
        nodes[2*i+1] = 0;
        thePCs[i] = 0;
    }

    setC();
}


PFEMElement2DBubble::~PFEMElement2DBubble()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;
    for(int i=0; i<3; i++) {
        if(thePCs[i] != 0) {
            thePCs[i]->disconnect(this->getTag());
        }
    }
}

int
PFEMElement2DBubble::getNumExternalNodes() const
{
    return ntags.Size();
}

const ID&
PFEMElement2DBubble::getExternalNodes()
{
    if(numDOFs.Size() == 0) return numDOFs;
    return ntags;
}

Node **
PFEMElement2DBubble::getNodePtrs(void)
{
    return nodes;
}

int
PFEMElement2DBubble::getNumDOF()
{
    if(numDOFs.Size() == 0) return 0;
    return numDOFs(numDOFs.Size()-1);
}

int
PFEMElement2DBubble::revertToLastCommit()
{
    return 0;
}

int PFEMElement2DBubble::commitState()
{
    return Element::commitState();
}

int
PFEMElement2DBubble::updateMatrix()
{
    int ndf = getNumDOF();
    M.resize(ndf, ndf);
    M.Zero();
    D.resize(ndf, ndf);
    D.Zero();
    F.resize(6);
    F.Zero();
    Fp.resize(3);
    Fp.Zero();

    // mass
    double m = getM();
    double mp = getMp();
    for(int a=0; a<3; a++) {
        M(numDOFs(2*a), numDOFs(2*a)) = m;          // Mxd
        M(numDOFs(2*a)+1, numDOFs(2*a)+1) = m;      // Myd

        for(int b=0; b<3; b++) {
            if(a == b) {
                M(numDOFs(2*a+1), numDOFs(2*b+1)) = 2*mp;   // Mp
            } else {
                M(numDOFs(2*a+1), numDOFs(2*b+1)) = mp;   // Mp
            }
        }
    }

    // damp
    Vector G(6);
    getG(G);
    Matrix L(3,3);
    getL(L);
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            D(numDOFs(2*a+1), numDOFs(2*b)) = G(2*b);   // GxT
            D(numDOFs(2*a+1), numDOFs(2*b)+1) = G(2*b+1); // GyT

            D(numDOFs(2*a), numDOFs(2*b+1)) = -G(2*a);   // -Gx
            D(numDOFs(2*a)+1, numDOFs(2*b+1)) = -G(2*a+1); // -Gy

            D(numDOFs(2*a+1), numDOFs(2*b+1)) = L(a,b);   // bubble
        }
    }

    // force
    getFp(Fp);
    getF(F);

    return 0;
}

int
PFEMElement2DBubble::update()
{
    if (dispon) {
        setJ();
    }

    // this is the trick to check negative jacobian
    if((kappa==-2 && J<0) || (kappa!=-2 && fabs(J)<1e-15)) {
        opserr<<"WARING: element "<<this->getTag()<<" area is "<<J<<"\n";
        for (int i=0; i<3; i++) {
            opserr << "node "<<nodes[2*i]->getTag()<<": \n";
            opserr << "coordinates - "<<nodes[2*i]->getCrds();
            opserr << "displacement - "<<nodes[2*i]->getTrialDisp();
        }
        opserr<<" -- PFEMElement2DBubble::update\n";
        return -1;
    }

    if (dispon) {
        setdJ();
        updateMatrix();
    }

    return 0;
}

const Matrix&
PFEMElement2DBubble::getMass()
{
    return M;
}

const Matrix&
PFEMElement2DBubble::getDamp()
{
    return D;
}

const Matrix&
PFEMElement2DBubble::getTangentStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}


const Matrix&
PFEMElement2DBubble::getInitialStiff()
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

int
PFEMElement2DBubble::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector&
PFEMElement2DBubble::getResistingForce()
{

    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    return P;
}

const Vector&
PFEMElement2DBubble::getResistingForceIncInertia()
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

    // get velocity, acceleration
    Vector v(ndf), vdot(ndf);
    for(int i=0; i<3; i++) {
        const Vector& accel = nodes[2*i]->getTrialAccel();
        vdot(numDOFs(2*i)) = accel(0);
        vdot(numDOFs(2*i)+1) = accel(1);

        const Vector& accel2 = nodes[2*i+1]->getTrialAccel();  // pressure
        vdot(numDOFs(2*i+1)) = accel2(0);

        const Vector& vel = nodes[2*i]->getTrialVel();
        v(numDOFs(2*i)) = vel(0);
        v(numDOFs(2*i)+1) = vel(1);

        const Vector& vel2 = nodes[2*i+1]->getTrialVel();   // pressure
        v(numDOFs(2*i+1)) = vel2(0);

    }

    // internal force
    P.addMatrixVector(1.0, getMass(), vdot, 1.0);
    P.addMatrixVector(1.0, getDamp(), v, 1.0);

    // external force
    for(int i=0; i<3; i++) {
        P(numDOFs(2*i)) -= F(2*i);
        P(numDOFs(2*i)+1) -= F(2*i+1);
        P(numDOFs(2*i+1)) -= Fp(i);
    }

    return P;
}


const char*
PFEMElement2DBubble::getClassType()const
{
    return "PFEMElement2DBubble";
}

int
PFEMElement2DBubble::sendSelf(int commitTag, Channel &theChannel)
{
    // int res = 0;
    // int dataTag = this->getDbTag();

    // // send vector
    // static Vector data(25);
    // data(0) = this->getTag();
    // data(1) = rho;
    // data(2) = mu;
    // data(3) = bx;
    // data(4) = by;
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
    //     opserr<<"WARNING: PFEMElement2DBubble::sendSelf - "<<this->getTag()<<" failed to send vector\n";
    //     return -1;
    // }


    return 0;
}

int
PFEMElement2DBubble::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // int res;
    // int dataTag = this->getDbTag();

    // // receive vector
    // static Vector data(12);
    // res = theChannel.recvVector(dataTag, commitTag, data);
    // if(res < 0) {
    //     opserr<<"WARNING: PFEMElement2DBubble::recvSelf - failed to receive vector\n";
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
PFEMElement2DBubble::setDomain(Domain *theDomain)
{
    numDOFs.resize(7);
    DomainComponent::setDomain(theDomain);

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
            opserr<<"in PFEMElement2DBubble - setDomain() "<<eletag<<"\n ";
            return;
        }
        ndf += nodes[2*i]->getNumberDOF();

        // set ndf
        numDOFs(2*i+1) = ndf;

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
        if (nodes[2*i+1] == 0) {
            opserr<<"WARNING: pressure node does not exist ";
            opserr<<"in PFEMElement2DBubble - setDomain() "<<eletag<<"\n ";
            return;
        }
        ntags(2*i+1) = nodes[2*i+1]->getTag();
        ndf += nodes[2*i+1]->getNumberDOF();
    }
    numDOFs(numDOFs.Size()-1) = ndf;

    if (!dispon) {
        setJ();
        setdJ();
    }
}

void
PFEMElement2DBubble::Print(OPS_Stream &s, int flag)
{
    s << "PFEMElement2DBubble: "<<this->getTag()<<endln;
}

int
PFEMElement2DBubble::displaySelf(Renderer &theViewer, int displayMode, float fact,
                                 const char **displayModes, int numModes)
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

    // determine the end points of the Tri31 based on
    // the display factor (a measure of the distorted image)
    // store this information in 3 3d vectors v1 through v3
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);

    nodes[0]->getDisplayCrds(v1, fact, displayMode);
    nodes[2]->getDisplayCrds(v2, fact, displayMode);
    nodes[4]->getDisplayCrds(v3, fact, displayMode);

    // point coordinates for polygon
    static Matrix coords(3, 3);
    for (int i = 0; i < 2; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
    }

    // color map
    static Vector values(3);
    values(0) = 0.0;
    values(1) = 0.0;
    values(2) = 0.0;

    // finally we draw the element using drawPolygon
    return theViewer.drawPolygon(coords, values);
}

// C = [0 0 0 1 0 -1]
//     [0 0 -1 0 1 0]
//     [0 -1 0 0 0 1]
//     [1 0 0 0 -1 0]
//     [0 1 0 -1 0 0]
//     [-1 0 1 0 0 0]
void
PFEMElement2DBubble::setC()
{
    if(C.noRows() == 6) return;
    C.resize(6,6); C.Zero();
    double l[12] = {3,5,4,2,5,1,0,4,1,3,2,0};
    for(int b=0; b<6; b++) {
        C(l[2*b],b) = 1.0;
        C(l[2*b+1],b) = -1.0;
    }
}

// J = x2y3−y2x3 + x3y1−x1y3 + x1y2-x2y1
void
PFEMElement2DBubble::setJ()
{
    Vector x(6);
    for(int a=0; a<3; a++) {
        const Vector& coord = nodes[2*a]->getCrds();
        const Vector& disp = nodes[2*a]->getTrialDisp();
        for(int i=0; i<2; i++) {
            x(2*a+i) = coord(i)+disp(i);
        }
    }

    J = x(2)*x(5)-x(3)*x(4)+x(1)*x(4)-x(0)*x(5)+x(0)*x(3)-x(1)*x(2);
}

// dJ =[y2-y3,x3-x2,y3-y1,x1-x3,y1-y2] = [b1,c1,b2,c2,b3,c3]
void
PFEMElement2DBubble::setdJ()
{
    Vector x(6);
    for(int a=0; a<3; a++) {
        const Vector& coord = nodes[2*a]->getCrds();
        const Vector& disp = nodes[2*a]->getTrialDisp();
        for(int i=0; i<2; i++) {
            x(2*a+i) = coord(i)+disp(i);
        }
    }

    dJ.addMatrixVector(0.0, C, x, 1.0);
}


double
PFEMElement2DBubble::getM() const
{
    return rho*thickness*J/6.0;
}

double
PFEMElement2DBubble::getMp() const
{
    if(kappa <= 0) return 0.0;
    return J*thickness/kappa/24.0;
}

double
PFEMElement2DBubble::getMbub() const
{
    return 27.0*rho*J*thickness/120.0;
}

void
PFEMElement2DBubble::getK(Matrix& k) const
{
    double J2 = mu*J/6.*thickness;

    k.resize(6,6);
    k.Zero();

    if (mu <= 0) return;

    double dNdx[3], dNdy[3];
    for(int a=0; a<3; a++) {
        dNdx[a] = dJ(2*a)/J;
        dNdy[a] = dJ(2*a+1)/J;
    }

    // other matrices
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            k(2*a, 2*b) += J2*(4*dNdx[a]*dNdx[b] + 3*dNdy[a]*dNdy[b]); // K1
            k(2*a, 2*b+1) += J2*(3*dNdy[a]*dNdx[b]-2*dNdx[a]*dNdy[b]); // K1
            k(2*a+1, 2*b) += J2*(3*dNdx[a]*dNdy[b]-2*dNdy[a]*dNdx[b]); // K1
            k(2*a+1, 2*b+1) += J2*(4*dNdy[a]*dNdy[b] + 3*dNdx[a]*dNdx[b]); // K1
        }
    }
}

void
PFEMElement2DBubble::getKbub(Matrix& kbub) const
{
    kbub.resize(2,2);
    kbub.Zero();

    if (mu <= 0) return;

    double b2 = 0, c2 = 0, bc = 0;
    for(int a=0; a<3; a++) {
        b2 += dJ(2*a)*dJ(2*a);
        c2 += dJ(2*a+1)*dJ(2*a+1);
        bc += dJ(2*a)*dJ(2*a+1);
    }

    double factor = 729.0*mu*thickness/(1080.0*J);

    kbub(0,0) = factor*(4*b2+3*c2);
    kbub(0,1) = factor*bc;
    kbub(1,0) = factor*bc;
    kbub(1,1) = factor*(4*c2+3*b2);
}

void
PFEMElement2DBubble::getGbub(Matrix& gbub) const
{
    gbub.resize(2,3);
    for(int a=0; a<2; a++) {
        for(int b=0; b<3; b++) {
            gbub(a,b) = dJ(2*b+a);
        }
    }

    gbub *= -27.0*thickness/120.0;
}



void
PFEMElement2DBubble::getF(Vector& f) const
{
    f.resize(6);
    f.Zero();

    // external force
    for(int a=0; a<3; a++) {
        f(2*a) = bx;
        f(2*a+1) = by;
    }
    f *= rho*thickness*J/6.0;

    // velocity
    if(mu > 0) {
        Vector v(6);
        for(int a=0; a<3; a++) {
            const Vector& vel = nodes[2*a]->getVel();
            v(2*a) = vel(0);
            v(2*a+1) = vel(1);
        }
        Matrix k(6,6);
        getK(k);
        f.addMatrixVector(1.0, k, v, -1.0);
    }
}

void
PFEMElement2DBubble::getFbub(Vector& fbub) const
{
    fbub.resize(2);

    // external force
    fbub(0) = rho*J*thickness*bx*27/120.0;
    fbub(1) = rho*J*thickness*by*27/120.0;
}


void
PFEMElement2DBubble::getG(Vector& g) const
{
    g = dJ;
    g *= thickness/6.0;
}

void
PFEMElement2DBubble::getL(Matrix& l) const
{
    Matrix Gbub(2,3);
    getGbub(Gbub);
    double Mbub = getMbub();

    Matrix Kbub(2,2);
    getKbub(Kbub);

    if (ops_Dt > 0) {
        Kbub(0,0) += Mbub/ops_Dt;
        Kbub(1,1) += Mbub/ops_Dt;
    }

    Matrix invKbub(2,2);
    Kbub.Invert(invKbub);

    l.resize(3,3);
    l.addMatrixTripleProduct(0.0, Gbub, invKbub, 1.0);
}

void
PFEMElement2DBubble::getFp(Vector& fp) const
{
    Matrix Gbub(2,3);
    getGbub(Gbub);
    double Mbub = getMbub();

    Matrix Kbub(2,2);
    getKbub(Kbub);

    if (ops_Dt > 0) {
        Kbub(0,0) += Mbub/ops_Dt;
        Kbub(1,1) += Mbub/ops_Dt;
    }

    Matrix invKbub(2,2);
    Kbub.Invert(invKbub);

    Vector Fbub(2);
    getFbub(Fbub);

    // bubble force
    fp.resize(3);
    fp.Zero();
    fp.addMatrixTransposeVector(0.0, Gbub, invKbub*Fbub, -1.0);
}

int
PFEMElement2DBubble::setParameter(const char **argv, int argc,
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

    return -1;
}

int
PFEMElement2DBubble::updateParameter (int passparameterID, Information &info)
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
        default:
            return -1;
    }
}

int
PFEMElement2DBubble::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}

const Matrix &
PFEMElement2DBubble::getDampSensitivity(int gradNumber)
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    return K;
}

const Matrix &
PFEMElement2DBubble::getMassSensitivity(int gradNumber)
{
    // resize K
    int ndf = this->getNumDOF();
    K.resize(ndf, ndf);
    K.Zero();

    // double dm = getdM();

    // for(int a=0; a<3; a++) {
    //     K(numDOFs(2*a), numDOFs(2*a)) = dm;          // Mxd
    //     K(numDOFs(2*a)+1, numDOFs(2*a)+1) = dm;
    // }

    return K;
}

const Vector &
PFEMElement2DBubble::getResistingForceSensitivity(int gradnumber)
{
    // resize P
    int ndf = this->getNumDOF();
    P.resize(ndf);
    P.Zero();

    Vector dF(6), dFp(3), vdot(6), v(6), p(3), du(6);
    for(int i=0; i<3; i++) {
        const Vector& accel = nodes[2*i]->getTrialAccel();
        vdot(2*i) = accel(0);
        vdot(2*i+1) = accel(1);

        const Vector& vel = nodes[2*i]->getTrialVel();
        v(2*i) = vel(0);
        v(2*i+1) = vel(1);

        const Vector& vel2 = nodes[2*i+1]->getTrialVel();   // pressure
        p(i) = vel2(0);

        du(2*i) = nodes[2*i]->getDispSensitivity(1,gradnumber);
        du(2*i+1) = nodes[2*i]->getDispSensitivity(2,gradnumber);
    }

    // consditional sensitivity
    getdF(dF);
    double dm = getdM();
    dF.addVector(-1.0, vdot, dm);

    getdFp(dFp);
    Matrix dl;
    getdL(dl);
    dFp.addMatrixVector(-1.0, dl, p, 1.0);

    // geometric sensitivity
    Matrix dM, dg, df;
    getdM(vdot, dM);
    getdG(p, dg);
    getdF(df);
    dF.addMatrixVector(1.0, dM, du, 1.0);
    dF.addMatrixVector(1.0, dg, du, -1.0);
    dF.addMatrixVector(1.0, df, du, -1.0);

    Matrix dgt, dL, dfp;
    getdGt(v, dgt);
    getdL(p, dL);
    getdFp(dfp);
    dFp.addMatrixVector(1.0, dgt, du, 1.0);
    dFp.addMatrixVector(1.0, dL, du, 1.0);
    dFp.addMatrixVector(1.0, dfp, du, -1.0);

    // copy
    for(int i=0; i<3; i++) {
        P(numDOFs(2*i)) += dF(2*i);
        P(numDOFs(2*i)+1) += dF(2*i+1);
        P(numDOFs(2*i+1)) += dFp(i);
    }

    return P;
}

int
PFEMElement2DBubble::commitSensitivity(int gradNumber, int numGrads)
{
    return 0;
}

// conditional sensitivity
double
PFEMElement2DBubble::getdM() const
{
    if(parameterID != 2) return 0.0;
    //return (1.0/6.0+3.0/40.0)*J*thickness;
    return thickness*J/6.0;
}

double
PFEMElement2DBubble::getdinvMbub() const
{
    if(parameterID != 2) return 0.0;
    //return -5040.0*ops_Dt/(1863.0*thickness*J*rho*rho);
    return -40.*ops_Dt/(9.*rho*rho*J*thickness);
}

void
PFEMElement2DBubble::getdF(Vector& df) const
{
    df.resize(6);

    // external force
    if(parameterID == 2) {
        for(int a=0; a<3; a++) {
            df(2*a) = bx/6.*J*thickness;
            df(2*a+1) = by/6.*J*thickness;
        }
    } else if(parameterID == 3) {
        for(int a=0; a<3; a++) {
            df(2*a) = rho/6.*J*thickness;
            df(2*a+1) = 0.0;
        }
    } else if(parameterID == 4) {
        for(int a=0; a<3; a++) {
            df(2*a) = 0.0;
            df(2*a+1) = rho/6.*J*thickness;
        }
    }

    // velocity
    if(mu > 0 && parameterID == 1) {
        Vector v(6);
        for(int a=0; a<3; a++) {
            const Vector& vel = nodes[2*a]->getTrialVel();
            v(2*a) = vel(0);
            v(2*a+1) = vel(1);
        }
        Matrix dk(6,6);
        getdK(dk);
        df.addMatrixVector(1.0, dk, v, -1.0);
    }
}

void
PFEMElement2DBubble::getdK(Matrix& dk) const
{
    double J2 = J/2.*thickness;
    // double lambda = -2*mu/3.0;

    dk.resize(6,6);
    dk.Zero();

    double dNdx[3], dNdy[3];
    for(int a=0; a<3; a++) {
        dNdx[a] = dJ(2*a)/J;
        dNdy[a] = dJ(2*a+1)/J;
    }

    // other matrices
    for(int a=0; a<3; a++) {
        for(int b=0; b<3; b++) {
            dk(2*a, 2*b) += J2*(2*dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]); // K1
            dk(2*a, 2*b+1) += J2*dNdy[a]*dNdx[b]; // K1
            dk(2*a+1, 2*b) += J2*dNdx[a]*dNdy[b]; // K1
            dk(2*a+1, 2*b+1) += J2*(2*dNdy[a]*dNdy[b] + dNdx[a]*dNdx[b]); // K1

            // K(2*a, 2*b) += lambda*J2*dNdx[a]*dNdx[b]; // K2
            // K(2*a, 2*b+1) += lambda*J2*dNdx[a]*dNdy[b]; // K2
            // K(2*a+1, 2*b) += lambda*J2*dNdy[a]*dNdx[b]; // K2
            // K(2*a+1, 2*b+1) += lambda*J2*dNdy[a]*dNdy[b]; // K2
        }
    }
}

void
PFEMElement2DBubble::getdFbub(Vector& dfb) const
{
    dfb.resize(2);

    // external force
    if(parameterID == 2) {
        dfb(0) = J*thickness*bx*27/120.0;
        dfb(1) = J*thickness*by*27/120.0;
    } else if(parameterID == 3) {
        dfb(0) = rho*thickness*J*27/120.0;
        dfb(1) = 0.0;
    } else if(parameterID == 4) {
        dfb(0) = 0.0;
        dfb(1) = rho*thickness*J*27/120.0;
    }

    // velocity
    if(mu > 0) {
        // Vector vb(2);
        // getVb(vb);
        // Matrix kb(2,2);
        // getKbub(kb);
        // fbub.addMatrixVector(1.0, kb, vb, -1.0);
    }
}

void
PFEMElement2DBubble::getdFp(Vector& dfp) const
{
    Matrix gb(2,3);
    getGbub(gb);

    double invmb = ops_Dt/getMbub();

    Vector fb(2);
    getFbub(fb);

    double dinvmb = getdinvMbub();

    Vector dfb(2);
    getdFbub(dfb);

    // bubble force
    dfp.resize(3);
    dfp.Zero();
    dfp.addMatrixTransposeVector(0.0, gb, fb, -dinvmb);
    dfp.addMatrixTransposeVector(1.0, gb, dfb, -invmb);
}

void
PFEMElement2DBubble::getdL(Matrix& dl) const
{
    Matrix gb;
    getGbub(gb);
    dl.resize(3,3);
    dl.addMatrixTransposeProduct(0.0,gb,gb,getdinvMbub());
}

// geometric sensitivity
void
PFEMElement2DBubble::getdM(const Vector& vdot, Matrix& dm) const
{
    dm.resize(6,6);
    dm.Zero();

    for(int a=0; a<6; a++) {
        for(int b=0; b<6; b++) {
            dm(a,b) = vdot(a)*dJ(b);
        }
    }
    //dm *= (1.0/6.0+3.0/40.0)*rho*thickness;
    dm *= rho*thickness/6.0;
}

void
PFEMElement2DBubble::getdinvMbub(const Vector& vb, Matrix& dmb) const {

    dmb.resize(2,6);
    dmb.Zero();

    for(int a=0; a<2; a++) {
        for(int b=0; b<6; b++) {
            dmb(a,b) = vb(a)*dJ(b);
        }
    }
    // dmb *= -5040.0*ops_Dt/(1863.0*thickness*rho*J*J);
    dmb *= -40.*ops_Dt/(9.*rho*J*J*thickness);
}

void
PFEMElement2DBubble::getdK(const Vector& v, Matrix& dk) const
{
    dk.resize(6,6);
    getK(dk);
    dk *= -1.0/J;

    Vector kv = dk*v;

    dk.Zero();
    for(int a=0; a<6; a++) {
        for(int b=0; b<6; b++) {
            dk(a,b) = kv(a)*dJ(b);
        }
    }
}

void
PFEMElement2DBubble::getdF(Matrix& df) const {

    df.resize(6,6);
    df.Zero();

    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            df(2*a,b) = bx*dJ(b);
            df(2*a+1,b) = by*dJ(b);
        }
    }

    df *= rho*thickness/6.0;

    // velocity
    if(mu > 0) {
        Vector v(6);
        for(int a=0; a<3; a++) {
            const Vector& vel = nodes[2*a]->getTrialVel();
            v(2*a) = vel(0);
            v(2*a+1) = vel(1);
        }
        Matrix dk(6,6);
        getdK(v,dk);
        df -= dk;
    }
}

void
PFEMElement2DBubble::getdFbub(Matrix& dfb) const {

    dfb.resize(2,6);
    dfb.Zero();

    for(int b=0; b<6; b++) {
        dfb(0,b) = bx*dJ(b);
        dfb(1,b) = by*dJ(b);
    }

    dfb *= rho*thickness*27.0/120.0;
}

void
PFEMElement2DBubble::getdG(const Vector& p, Matrix& dg) const {

    dg = C;
    double sump = 0.0;
    for(int i=0; i<p.Size(); i++) {
        sump += p(i);
    }
    dg *= thickness*sump/6.0;
}

void
PFEMElement2DBubble::getdGt(const Vector& v, Matrix& dgt) const {

    Matrix vt3(3,6);
    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            vt3(a,b) = v(b);
        }
    }
    dgt = vt3*C;
    dgt *= thickness/6.0;
}

void
PFEMElement2DBubble::getdGb(const Vector& p, Matrix& dgb) const {

    dgb.resize(2,6);
    dgb.Zero();

    for(int a=0; a<2; a++) {
        for(int b=0; b<6; b++) {
            for(int i=0; i<3; i++) {
                dgb(a,b) += C(2*i+a,b)*p(i);
            }
        }
    }

    dgb *= -27.0*thickness/120.0;
}

void
PFEMElement2DBubble::getdGbt(const Vector& vb, Matrix& dgbt) const {

    dgbt.resize(3,6);
    dgbt.Zero();

    for(int a=0; a<3; a++) {
        for(int b=0; b<6; b++) {
            for(int i=0; i<vb.Size(); i++) {
                dgbt(a,b) += C(2*a+i,b)*vb(i);
            }
        }
    }
    dgbt *= -27.0*thickness/120.0;
}

void
PFEMElement2DBubble::getdFp(Matrix& dfp) const {

    Matrix gb(2,3);
    getGbub(gb);

    Vector fb(2);
    getFbub(fb);

    double invmb = ops_Dt/getMbub();

    getdGbt(fb*invmb, dfp);

    Matrix dmb(2,6);
    getdinvMbub(fb,dmb);

    dfp.addMatrixTransposeProduct(-1.0, gb, dmb, -1.0);

    Matrix dfb(2,6);
    getdFbub(dfb);
    dfp.addMatrixTransposeProduct(1.0, gb, dfb, -invmb);
}

void
PFEMElement2DBubble::getdL(const Vector& p, Matrix& dl) const {
    Matrix gb(2,3);
    getGbub(gb);

    double invmb = ops_Dt/getMbub();

    getdGbt(gb*p*invmb, dl);

    Matrix dmb(2,6);
    getdinvMbub(gb*p, dmb);

    dl.addMatrixTransposeProduct(1.0,gb,dmb,1.0);

    Matrix dgb(2,6);
    getdGb(p,dgb);
    dl.addMatrixTransposeProduct(1.0,gb,dgb,invmb);
}
