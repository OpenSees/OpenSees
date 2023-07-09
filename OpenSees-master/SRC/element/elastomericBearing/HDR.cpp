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

// Written: Manish Kumar (mkumar2@buffalo.edu)
// Created: 02/12
// Revision: A
//
// Description: This file contains the implementation of the
// HDR class.

#include "HDR.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <elementAPI.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <G3Globals.h>
#include <Message.h>
using namespace std;
#include <iostream>

#define PI 3.14159l


// initialize the class wide variables
Matrix HDR::theMatrix(12,12);
Vector HDR::theVector(12);


static int numMyBearing = 0;
static int tag = 0;  // Tag to identify if bearing has failed in buckling
void *OPS_HDR()
{
    // print out a message about who wrote this element & any copyright info wanted
    if (numMyBearing == 0) {
        opserr << "HDR element - Written by Manish Kumar, University at Buffalo, 2012\n";
        numMyBearing++;
    }
    
    Element *theEle = 0;
    
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs == 0) { // parallel processing
        theEle = new HDR();
        return theEle;
    }
    
    if (numArgs !=20 && numArgs !=26 && numArgs !=27 && numArgs !=28
        && numArgs !=29 && numArgs !=30 && numArgs !=31 && numArgs !=32) {
            opserr << "ERROR - HDR incorrect # args provided";
            return theEle;
    }
    
    // get the id and end nodes
    int iData[3];
    double dData[17];
    int numData;
    
    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data\n";
        return 0;
    }
    
    int eleTag = iData[0];
    
    numData = 17;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element properties for element" << eleTag << endln;
        return 0;
    }
    
    // get the orientation vector
    Vector x(0);
    Vector y(3); y(0) = -1.0; y(1) = 0.0; y(2) = 0.0;
    
    // The default values of the parameters
    double kl = 10.0;       // Cavitation parameter
    double phi = 0.5;       // Damage index
    double al = 1.0;        // Strength degradation parameter
    double sDratio = 0.5;   // Shear distance ratio
    double m = 0.0;         // Mass of the bearing
    double tc1 = 0.0;       // Cover thickness
    
    if (numArgs >= 26) {
        double value;
        x.resize(3);
        numData = 1;
        for (int i=0; i<3; i++) {
            if (OPS_GetDoubleInput(&numData, &value) != 0) {
                opserr << "WARNING invalid orientation value for element" << eleTag << endln;
                return 0;
            } else {
                x(i) = value;
            }
        }
        for (int i=0; i<3; i++) {
            if (OPS_GetDoubleInput(&numData, &value) != 0) {
                opserr << "WARNING invalid orientation value for element" << eleTag << endln;
                return 0;
            } else {
                y(i) = value;
            }
        }
        if (numArgs >= 27) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &kl) != 0) {
                opserr << "WARNING error reading element property cavitation parameter for element" << eleTag << endln;
                return 0;
            }
            if (numArgs >= 28) {
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &phi) != 0) {
                    opserr << "WARNING error reading element property damage index for element" << eleTag << endln;
                    return 0;
                }
                if (numArgs >= 29) {
                    numData = 1;
                    if (OPS_GetDoubleInput(&numData, &al) != 0) {
                        opserr << "WARNING error reading element property strength degradation parameter for element" << eleTag << endln;
                        return 0;
                    }
                    if (numArgs >= 30) {
                        numData = 1;
                        if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
                            opserr << "WARNING error reading element property shear distance ratio for element" << eleTag << endln;
                            return 0;
                        }
                        if (numArgs >= 31) {
                            numData = 1;
                            if (OPS_GetDoubleInput(&numData, &m) != 0) {
                                opserr << "WARNING error reading element property mass for element" << eleTag << endln;
                                return 0;
                            }
                            if (numArgs == 32) {
                                numData = 1;
                                if (OPS_GetDoubleInput(&numData, &tc1) != 0) {
                                    opserr << "WARNING error reading element property cover thickness for element" << eleTag << endln;
                                    return 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm == 3) {
        // check space frame problem has 6 dof per node
        if (ndf != 6)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for space problem need 6 - HDR \n"; 
        }
        theEle = new HDR(iData[0], iData[1], iData[2], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], int(dData[6]), dData[7], dData[8], 
            dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], y, x, kl, phi, al, sDratio, m, tc1);
    }
    
    if (theEle == 0) {
        opserr << "WARNING this element only works for 3D problem with 6 dof per node: " << eleTag << endln;
        return 0;
    }
    
    return theEle;
}


HDR::HDR(int tag, int Nd1, int Nd2, double Gr, double kbulk, double Di, double Do, double ts1, double tr, int n1, 
    double _a1, double _a2, double _a3, double _b1, double _b2, double _b3, double _c1, double _c2, double _c3, double _c4, 
    const Vector _y, const Vector _x, double kl, double PhiMax, double al, double sDratio, double m, double tc1)
    : Element(tag, ELE_TAG_HDR), connectedExternalNodes(2),
    a1(_a1), a2(_a2), a3(_a3), b1(_b1), b2(_b2), b3(_b3), c1(_c1), c2(_c2), c3(_c3), c4(_c4), PhiM(PhiMax), ac(al),
    G(Gr), Kbulk(kbulk), x(_x), y(_y), shearDistI(sDratio), mass(m), tc(tc1), D1(Di), D2(Do), L(0.0), n(n1), ts(ts1),
    Fcrn(0.0), ucrn(0.0), Fcrmin(0.0), Fcn(0.0), ucn(0.0), Fmax(0.0), umax(0.0), F2(2),
    ub(6), ubdot(6), qb(6), kb(6,6), ul(12), Tgl(12,12), Tlb(6,12), ubC(6), F2C(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "HDR::HDR() - element: "
            << this->getTag() << " failed to create an ID of size 2\n";
        exit(-1);
    }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
    
    // vertical motion
    A = (PI/4.0)*((D2+tc)*(D2+tc) - D1*D1);                 // Bonded rubber area of bearing
    S = (D2 - D1)/(4*tr);                                   // Shape factor for case with central hole
    Tr = n*tr;                                              // height of rubber in the bearing
    h = Tr + (n-1)*ts;                                      // height of rubber + shims
    double F;
    if (D1 < DBL_EPSILON) {
        F = 1.0;                                            // Dimension modification factor
    } else {
        double r = D2/D1;                                   // Outer to inner diameter ratio
        F = (r*r+1)/((r-1)*(r-1)) + (1+r)/((1-r)*log(r));   // Dimension modification factor
    }
    Ec = 1.0/((1.0/(6*G*S*S*F))+(4.0/3.0)*(1.0/Kbulk));     // Compressive modulus of elastomeric bearing
    double E = 3.0*G;                                       // Elastic modulus
    double I = (PI/64.0)*(pow((D2+tc),4) - pow(D1,4));      // Moment of inertia of bearing
    rg = sqrt(I/A);                                         // Radius of gyration
    Kv0 = A*Ec/Tr;                                          // Pre-cavitation stiffness at zero lateral displacement
    Kv = Kv0;                                               // Pre-cavitation stiffness
    if (kl < DBL_EPSILON) {
        kc = 0.0001;                                        // cavitation parameter
    } else {
        kc = kl;                                            // cavitation parameter
    }
    double Er = Ec/3.0;                                     // Rotation modulus of bearing
    double As = A*h/Tr;                                     // Adjusted shear area of bearing
    double Is = I*h/Tr;                                     // Adjusted moment of inertia of bearing
    double Pe = PI*PI*Er*Is/(h*h);                          // Euler buckling load of bearing
    Fcr = -sqrt(Pe*G*As);                                   // Critical buckling load in compression
    Fcrn = Fcr;                                             // Current value of critical buckling load
    Fcrmin = Fcr;                                           // Initial value of critical buckling load during loading
    ucr = Fcr/Kv0;                                          // Initial value of critical buckling deformation
    ucrn = ucr;                                             // Current value of critical buckling deformation
    Fc = 3.0*G*A;                                           // Cavitation force
    Fcn = Fc;                                               // Initial value of cavitation load
    uc = Fc/Kv0;                                            // Initial cavitation deformation
    ucn = uc;                                               // Initial value of cavitation deformation
    Fmax = Fc;                                              // Initial value of maximum tensile force
    umax = uc;                                              // Initial value of maximum tensile deformation
    
    // torsion
    Kt = G*(2*Is)/h;                                        // torsional stiffness of bearing
    
    // rotation
    Kr = Er*Is/h;                                           // rotational stiffness of bearing
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = a1;
    kbInit(2,2) = a1;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
}


HDR::HDR()
    : Element(0, ELE_TAG_HDR), connectedExternalNodes(2),
    a1(0.0), a2(0.0), a3(0.0), b1(0.0), b2(0.0), b3(0.0), c1(0.0), c2(0.0), c3(0.0), c4(0.0),
    S(0.0), Ec(0.0), Kv0(0.0), Kv(0.0), kc(10.0), PhiM(0.5), ac(1.0), Fcr(0.0), ucr(0.0), Fc(0.0),
    uc(0.0), Kt(0.0), Kr(0.0), G(0.0), Kbulk(0.0), x(0), y(0), shearDistI(0.5), mass(0.0), tc(0.0),
    Tr(0.0), D1(0.0), D2(0.0), L(0.0), h(0.0), rg(0.0), A(0.0), Ar(0.0), n(0.0), ts(0.0), Fcrn(0.0),
    ucrn(0.0), Fcrmin(0.0), Fcn(0.0), ucn(0.0), Fmax(0.0), umax(0.0), DSplus(0.0), DSminus(0.0),
    DS(0.0), DM(0.0), Delta(0.0), F2(2), ub(6), ubdot(6), qb(6), kb(6,6), ul(6), Tgl(12,12), Tlb(6,12),
    ubC(6), DSplusC(0.0), DSminusC(0.0), DSC(0.0), DMC(0.0), F2C(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "HDR::HDR() - "
            <<  "failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}


HDR::~HDR()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
}


int HDR::getNumExternalNodes() const
{
    return 2;
}


const ID& HDR::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** HDR::getNodePtrs()
{
    return theNodes;
}


int HDR::getNumDOF()
{
    return 12;
}


void HDR::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        
        return;
    }
    
    // first set the node pointers
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "WARNING HDR::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING HDR::setDomain() - Nd2: "
                << connectedExternalNodes(1)
                << " does not exist in the model for";
        }
        opserr << " element: " << this->getTag() << endln;
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != 6)  {
        opserr << "HDR::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "HDR::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int HDR::commitState()
{
    int errCode = 0;
    
    double uh = sqrt(ub(1)*ub(1) + ub(2)*ub(2));
    //double uhC = sqrt(ubC(1)*ubC(1) + ubC(2)*ubC(2));
    
    // vertical motion
    Kv = Kv0*(1.0/(1.0+(3.0/(PI*PI))*(uh/rg)*(uh/rg)));
    if (uh > DBL_EPSILON)
        uc = Fc/Kv;
    
    // tension
    if (ub(0) > umax) {
        umax = ub(0);
        Fcn = Fc*(1.0-PhiM*(1.0-exp(-ac*(ub(0)-uc)/uc)));
    }
    
    // compression
    double Delta = 2.0*acos(uh/D2);   // this becomes undefined for uh/D2 > 1.0
    //Ar = (D2*D2/4.0)*(Delta-sin(Delta));
    Ar = ((D2+tc)*(D2+tc) - D1*D1)/4.0*(Delta-sin(Delta));  // A does not include central hole
    if (Ar/A < 0.2 || uh/D2 >= 1.0) {
        Fcrn = 0.2*Fcr;
    } else {
        Fcrn = Fcr*Ar/A;
    }
    
    if (Fcrn > Fcrmin)
        Fcrmin = Fcrn;
    ucrn = Fcrn/Kv;
    
    // horizontal motion
    // commit trial history variables for horizontal direction
    DSplusC = DSplus;
    DSminusC = DSminus;
    DSC = DS;
    DMC = DM;
    F2C = F2;
    ubC = ub;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int HDR::revertToLastCommit()
{
    return 0;
}


int HDR::revertToStart()
{
    int errCode = 0;
    
    // reset trial history variables
    ub.Zero();
    qb.Zero();
    F2.Zero();
    DSplus = 0;
    DSminus = 0;
    DS = 0;
    DM = 0;
    Delta = 0;
    
    // reset committed history variables
    ubC.Zero();
    F2C.Zero();
    DSplusC = 0;
    DSminusC = 0;
    DSC = 0;
    DMC = 0;
    
    // reset stiffness matrix in basic system
    kb = kbInit;
    
    return errCode;
}


int HDR::update()
{
    // get global trial displacements and velocities
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    
    static Vector ug(12), ugdot(12), uldot(12);
    for (int i=0; i<6; i++)  {
        ug(i)   = dsp1(i);  ugdot(i)   = vel1(i);
        ug(i+6) = dsp2(i);  ugdot(i+6) = vel2(i);
    }
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
    
    // 1) get axial force and stiffness in basic x-direction
    ucn = Fcn/Kv;
    Fmax = Fc*(1.0+(1.0/(Tr*kc))*(1.0-exp(-kc*(umax-uc))));
    
    if (ub(0) <= ucrn) {
        kb(0,0) = Kv0/10000.0;
        qb(0) = Fcrmin + kb(0,0)*(ub(0)-ucrn);
    } else if (ub(0) < ucn) {
        kb(0,0) = Kv;
        qb(0) = kb(0,0)*ub(0);
    } else if(ub(0) < umax) {
        kb(0,0) = (Fmax-Fcn)/(umax-ucn);
        qb(0) = Fcn + (Fmax-Fcn)/(umax-ucn)*(ub(0)-ucn);
    } else {
        kb(0,0) = ((Fc/Tr)*exp(-kc*(ub(0)-uc)));
        qb(0) = Fc*(1.0+(1.0/(Tr*kc))*(1.0-exp(-kc*(ub(0)-uc))));
    }
    
    //2) calculate shear forces and stiffnesses in basic y- and z-direction
    Matrix I(2,2), DF_DU(2,2), DF1_DU(2,2), DF2_DU(2,2), dF1_dU(2,2), dF2_dn(2,2), dF2_dMu(2,2), dn_dU(2,2), dMu_dU(2,2), dFhat_dU(2,2);
    Vector U(2), UC(2), Fhat(2), F1(2), Mu(2), n(2), dF1_dKS1(2), dF1_dKM(2), dDS_dU(2), dDM_dU(2), dF2_dR(2), dF2_dDelta(2), dR_dU(2), dDelta_dU(2);
    double dKS1_dDS, dKS2_dDS, dKM_dDM;
    
    U(0) = ub(1), U(1) = ub(2);
    UC(0) = ubC(1), UC(1) = ubC(2);
    double uh = U.Norm();
    double uhC = UC.Norm();
    double dU = (U-UC).Norm();
    if (dU > DBL_EPSILON) {
        n = (U-UC)/dU;
    } else {
        n.Zero();
    }
    I.Zero();
    I(0,0) = I(1,1) = 1.0;
    
    if (uh > DSplusC) {
        DSplus = uh;
        DSminus = DSminusC + (DSplus-DSplusC);
        DS = DSC;
    }
    if (uh < DSminusC) {
        DSplus = DSplusC;
        DSminus = uh;
        DS = DSC - (DSminus-DSminusC);
    }
    if (uh < uhC) {
        DM = DMC + uhC - uh;
    }
    
    double KS1 = exp(-c1*pow(DS,3));
    double KS2 = exp(-c2*pow(DS,3));
    double KS2C = exp(-c2*pow(DSC,3));
    double KM = c3 + (1.0-c3)*exp(-c4*pow(DM,3));
    double R = b1 + KS2*b2*pow(uh,2);
    double RC = b1 + KS2C*b2*pow(uhC,2);
    F1 = KS1*KM*(a1 + a2*pow(uh,2) + a3*pow(uh,4))*U;
    
    Fhat = R*n;
    double DeltaC = (RC*n-F2C).Norm();
    Delta = DeltaC/(1.0 + b3*dU);
    
    double Fnorm = (Fhat-F2C).Norm();
    if (Fnorm > DBL_EPSILON) {
        Mu = (Fhat-F2C)/Fnorm;
    } else {
        Mu.Zero();
    }
    
    F2 = R*n - Delta*Mu;
    
    // forces
    qb(1) = F1(0) + F2(0);
    qb(2) = F1(1) + F2(1);
    
    // consistent tangent matrix
    dF1_dU = KS1*KM*((a1+a2*pow(uh,2)+a3*pow(uh,4))*I+2*(a2+2*a3*pow(uh,2))*U%U);
    dF1_dKS1 = KM*(a1+a2*pow(uh,2)+a3*pow(uh,4))*U;
    dF1_dKM = KS1*(a1+a2*pow(uh,2)+a3*pow(uh,4))*U;
    dKS1_dDS = -3*c1*pow(DS,2)*exp(-c1*pow(DS,3));
    dKS2_dDS = -3*c2*pow(DS,2)*exp(-c2*pow(DS,3));
    dKM_dDM = -3*c4*(1-c4)*pow(DM,2)*exp(-c4*pow(DM,3));
    if (uh < DSminusC && uh > DBL_EPSILON) {
        dDS_dU = -(1.0/uh)*U;
    } else {
        dDS_dU.Zero();
    }
    if (uh < uhC && uh > DBL_EPSILON) {
        dDM_dU = -(1.0/uh)*U;
    } else {
        dDM_dU.Zero();
    }
    dF2_dR = n;
    dF2_dn = R*I;
    dF2_dDelta = -1.0*Mu;
    dF2_dMu = -Delta*I;
    dR_dU = b2*(2*KS2*U+pow(uh,2)*dKS2_dDS*dDS_dU);
    if (dU > DBL_EPSILON) {
        dn_dU = (1.0/dU)*(I-n%n);
    } else {
        dn_dU.Zero();
    }
    dDelta_dU = -1.0*DeltaC*b3*n/pow(1.0+b3*dU,2);
    
    //cout<<"R: "<<R<<" n(0): "<<n(0)<<" n(1): "<<n(1)<<" dR_dU(0): "<<dR_dU(0)<<" dR_dU(1): "<<dR_dU(1)<<" dn_dU(0,0): "<<dn_dU(0,0)<<" dn_dU(0,1): "<<dn_dU(0,1)<<" dn_dU(1,0): "<<dn_dU(1,0)<<" dn_dU(1,1): "<<dn_dU(1,1)<<endln;
    dFhat_dU = n%dR_dU + R*dn_dU;
    if (Fnorm > DBL_EPSILON) {
        dMu_dU = (1.0/Fnorm)*(I-Mu%Mu)*dFhat_dU;
    } else {
        dMu_dU.Zero();
    }
    
    DF1_DU = dF1_dU+dF1_dKS1%(dKS1_dDS*dDS_dU)+dF1_dKM%(dKM_dDM*dDM_dU);
    DF2_DU = dF2_dR%dR_dU+dF2_dn*dn_dU+dF2_dDelta%dDelta_dU+dF2_dMu*dMu_dU;
    DF_DU = DF1_DU + DF2_DU;
    kb(1,1) = DF_DU(0,0);
    kb(1,2) = DF_DU(0,1);
    kb(2,1) = DF_DU(1,0);
    kb(2,2) = DF_DU(1,1);
    
    //cout<<"Delta: "<<Delta<<" DeltaC: "<<DeltaC<<" R: "<<R<<" dn_dU(0,0): "<<dn_dU(0,0)<<" dn_dU(1,1): "<<dn_dU(1,1)<<endln;
    //cout<<"a1: "<<a1<<" a2: "<<a2<<" a3: "<<a3<<" b1: "<<b1<<" b2: "<<b2<<" b3: "<<b3<<" c1: "<<c1<<" c2: "<<c2<<" c3: "<<c3<<" c4: "<<c4<<endln;
    //cout<<"Delta: "<<Delta<<" dKS1_dDS: "<<dKS1_dDS<<" dKS2_dDS: "<<dKS2_dDS<<" dKM_dDM: "<<dKM_dDM<<" DF_DU(0,0): "<<DF_DU(0,0)<<" DF_DU(0,1): "<<DF_DU(0,1)<<" DF_DU(1,0): "<<DF_DU(1,0)<<" DF_DU(1,1): "<<DF_DU(1,1)<<endln;
    
    // 3) get moment and stiffness in basic x-direction
    qb(3) = Kt*ub(3);
    kb(3,3) = Kt;
    
    // 4) get moment and stiffness in basic y-direction
    qb(4) = Kr*ub(4);
    kb(4,4) = Kr;
    
    // 5) get moment and stiffness in basic z-direction
    qb(5) = Kr*ub(5);
    kb(5,5) = Kr;
    
    return 0;
}


const Matrix& HDR::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    double kGeo1 = 0.5*qb(0);
    kl(5,1)   -= kGeo1;
    kl(5,7)   += kGeo1;
    kl(11,1)  -= kGeo1;
    kl(11,7)  += kGeo1;
    kl(4,2)   += kGeo1;
    kl(4,8)   -= kGeo1;
    kl(10,2)  += kGeo1;
    kl(10,8)  -= kGeo1;
    double kGeo2 = kGeo1*shearDistI*L;
    kl(5,5)   += kGeo2;
    kl(11,5)  -= kGeo2;
    kl(4,4)   += kGeo2;
    kl(10,4)  -= kGeo2;
    double kGeo3 = kGeo1*(1.0 - shearDistI)*L;
    kl(5,11)  -= kGeo3;
    kl(11,11) += kGeo3;
    kl(4,10)  -= kGeo3;
    kl(10,10) += kGeo3;
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& HDR::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& HDR::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    return theMatrix;
}


const Matrix& HDR::getMass()
{
    // zero the matrix
    theMatrix.Zero();
    
    // check for quick return
    if (mass == 0.0)  {
        return theMatrix;
    }
    
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theMatrix(i,i)     = m;
        theMatrix(i+6,i+6) = m;
    }
    
    return theMatrix;
}


void HDR::zeroLoad()
{
    theLoad.Zero();
}


int HDR::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"HDR::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int HDR::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
        opserr << "HDR::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible.\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+6) -= m * Raccel2(i);
    }
    
    return 0;
}


const Vector& HDR::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(12);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
    
    // add P-Delta moments to local forces
    double kGeo1 = 0.5*qb(0);
    double MpDelta1 = kGeo1*(ul(7)-ul(1));
    ql(5)  += MpDelta1;
    ql(11) += MpDelta1;
    double MpDelta2 = kGeo1*shearDistI*L*ul(5);
    ql(5)  += MpDelta2;
    ql(11) -= MpDelta2;
    double MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(11);
    ql(5)  -= MpDelta3;
    ql(11) += MpDelta3;
    double MpDelta4 = kGeo1*(ul(8)-ul(2));
    ql(4)  -= MpDelta4;
    ql(10) -= MpDelta4;
    double MpDelta5 = kGeo1*shearDistI*L*ul(4);
    ql(4)  += MpDelta5;
    ql(10) -= MpDelta5;
    double MpDelta6 = kGeo1*(1.0 - shearDistI)*L*ul(10);
    ql(4)  -= MpDelta6;
    ql(10) += MpDelta6;
    
    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    return theVector;
}


const Vector& HDR::getResistingForceIncInertia()
{
    this->getResistingForce();
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    // add the damping forces if rayleigh damping
    //if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    //    theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
    // now include the mass portion
    if (mass != 0.0)  {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();
        
        double m = 0.5*mass;
        for (int i=0; i<3; i++)  {
            theVector(i)   += m * accel1(i);
            theVector(i+6) += m * accel2(i);
        }
    }
    
    return theVector;
}


int HDR::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(26);
    data(0)  = this->getTag();
    data(1)  = a1;
    data(2)  = a2;
    data(3)  = a3;
    data(4)  = b1;
    data(5)  = b2;
    data(6)  = b3;
    data(7)  = c1;
    data(8)  = c2;
    data(9)  = c3;
    data(10) = c4;
    data(11) = G;
    data(12) = Kbulk;
    data(13) = D1;
    data(14) = D2;
    data(15) = ts;
    data(16) = Tr;
    data(17) = n;
    data(18) = x.Size();
    data(19) = y.Size();
    data(20) = kc;
    data(21) = PhiM;
    data(22) = ac;
    data(23) = shearDistI;
    data(24) = mass;
    data(25) = tc;
    
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    
    return 0;
}


int HDR::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // receive element parameters
    static Vector data(26);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    
    a1 = data(1);
    a2 = data(2);
    a3 = data(3);
    b1 = data(4);
    b2 = data(5);
    b3 = data(6);
    c1 = data(7);
    c2 = data(8);
    c3 = data(9);
    c4 = data(10);
    
    G = data(11);
    Kbulk = data(12);
    D1 = data(13);
    D2 = data(14);
    ts = data(15);
    Tr = data(16);
    n = data(17);
    kc = data(20);
    PhiM = data(21);
    ac = data(22);
    shearDistI = data(23);
    mass = data(24);
    tc = data(25);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive remaining data
    if ((int)data(18) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(19) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    
    // vertical motion
    A = (PI/4.0)*((D2+tc)*(D2+tc) - D1*D1);
    S = (D2 - D1)/(4*Tr/n);
    h = Tr + (n-1)*ts;
    double F;
    if (D1 < DBL_EPSILON) {
        F = 1.0;
    } else {
        double r = D2/D1;
        F = (r*r+1)/((r-1)*(r-1)) + (1+r)/((1-r)*log(r));
    }
    Ec = 1.0/((1.0/(6*G*S*S*F))+(4.0/3.0)*(1.0/Kbulk));
    double E = 3.0*G;
    double I = (PI/64.0)*(pow((D2+tc),4) - pow(D1,4));
    rg = sqrt(I/A);
    Kv0 = A*Ec/Tr;
    Kv = Kv0;
    double Er = Ec/3.0;
    double As = A*h/Tr;
    double Is = I*h/Tr;
    double Pe = PI*PI*Er*Is/(h*h);
    Fcr = -sqrt(Pe*G*As);
    Fcrn = Fcr;
    Fcrmin = Fcr;
    ucr = Fcr/Kv0;
    ucrn = ucr;
    Fc = 3.0*G*A;
    Fcn = Fc;
    uc = Fc/Kv0;
    ucn = uc;
    Fmax = Fc;
    umax = uc;
    
    // torsion
    Kt = G*(2*Is)/h;
    
    // rotation
    Kr = Er*Is/h;
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = a1;
    kbInit(2,2) = a1;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
    
    return 0;
}


int HDR::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void HDR::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "************************************************************" << endln;
        s << "Element: " << this->getTag();
        s << "  type: HDR  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "************************************************************" << endln;
        s << "GEOMETRIC PROPERTIES" << endln;
        s << "D1: "<< D1 << " D2: " << D2 << " L: " << L << " Tr: " << Tr << " n: " << n << " A: " << A << endln;
        s << "MATERIAL PROPERTIES" << endln;
        s << "G: " << G << " kc: " << kc << " ac: " << ac << " PhiM: " << PhiM << " shearDistI: " << shearDistI << " mass: " << mass << endln;
        s << "MECHANICAL PROPERTIES: HORIZONTAL MOTION" << endln;
        s << "a1: " << a1 << " a2: " << a2 << " a3: " << a3 <<endln;
        s << "b1: " << b1 << " b2: " << b2 << " b3: " << b3 <<endln;
        s << "c1: " << c1 << " c2: " << c2 << " c3: " << c3 << " c4: " << c4 << endln;
        s << "MECHANICAL PROPERTIES: VERTICAL MOTION" << endln;
        s << "Ec: " << Ec << " Kv0: " << Kv0 << " Kv: " << Kv << " uc: " << uc << " Fcr: " << Fcr << " ucr: " << ucr << " Fcn: " << Fcn << " umax: " << umax << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
        s << "************************************************************" << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"HDR\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"D1\": " << D1 << ", ";
        s << "\"D2\": " << D2 << ", ";
        s << "\"L\": " << L << ", ";
        s << "\"Tr\": " << Tr << ", ";
        s << "\"n\": " << n << ", ";
        s << "\"A\": " << A << ", ";
        s << "\"G\": " << G << ", ";
        s << "\"kc\": " << kc << ", ";
        s << "\"ac\": " << ac << ", ";
        s << "\"PhiM\": " << PhiM << ", ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"mass\": " << mass << ", ";
        s << "\"a1\": " << a1 << ", ";
        s << "\"a2\": " << a2 << ", ";
        s << "\"a3\": " << a3 << ", ";
        s << "\"b1\": " << b1 << ", ";
        s << "\"b2\": " << b2 << ", ";
        s << "\"b3\": " << b3 << ", ";
        s << "\"c1\": " << c1 << ", ";
        s << "\"c2\": " << c2 << ", ";
        s << "\"c3\": " << c3 << ", ";
        s << "\"c4\": " << c4 << "}";
    }
}


Response* HDR::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","HDR");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Pz_1");
        output.tag("ResponseType","Mx_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","Mx_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
        output.tag("ResponseType","N_ 1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Tz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 2, theVector);
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
        output.tag("ResponseType","qb1");
        output.tag("ResponseType","qb2");
        output.tag("ResponseType","qb3");
        output.tag("ResponseType","qb4");
        output.tag("ResponseType","qb5");
        output.tag("ResponseType","qb6");
        
        theResponse = new ElementResponse(this, 3, Vector(6));
    }
    // local displacements
    else if (strcmp(argv[0],"localDisplacement") == 0 ||
        strcmp(argv[0],"localDisplacements") == 0)
    {
        output.tag("ResponseType","ux_1");
        output.tag("ResponseType","uy_1");
        output.tag("ResponseType","uz_1");
        output.tag("ResponseType","rx_1");
        output.tag("ResponseType","ry_1");
        output.tag("ResponseType","rz_1");
        output.tag("ResponseType","ux_2");
        output.tag("ResponseType","uy_2");
        output.tag("ResponseType","uz_2");
        output.tag("ResponseType","rx_2");
        output.tag("ResponseType","ry_2");
        output.tag("ResponseType","rz_2");
        
        theResponse = new ElementResponse(this, 4, theVector);
    }
    // basic displacements
    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
        output.tag("ResponseType","ub1");
        output.tag("ResponseType","ub2");
        output.tag("ResponseType","ub3");
        output.tag("ResponseType","ub4");
        output.tag("ResponseType","ub5");
        output.tag("ResponseType","ub6");
        
        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // hysteretic and damage evolution parameters
    else if (strcmp(argv[0],"Evolution") == 0 || strcmp(argv[0],"evolution") == 0 || 
        strcmp(argv[0],"EvolutionParameter") == 0 || strcmp(argv[0],"evolutionparameter") == 0 ||
        strcmp(argv[0],"z") == 0)
    {
        output.tag("ResponseType","DSplus");
        output.tag("ResponseType","DSminus");
        output.tag("ResponseType","DS");
        output.tag("ResponseType","DM");
        output.tag("ResponseType","Delta");
        
        theResponse = new ElementResponse(this, 6, Vector(5));
    }
    // basic stiffness
    else if (strcmp(argv[0],"kb") == 0 || strcmp(argv[0],"basicStiff") == 0 ||
        strcmp(argv[0],"basicStiffness") == 0)
    {
        output.tag("ResponseType","kb22");
        output.tag("ResponseType","kb23");
        output.tag("ResponseType","kb32");
        output.tag("ResponseType","kb33");
        
        theResponse = new ElementResponse(this, 7, Vector(4));
    }
    // parameters
    else if (strcmp(argv[0],"param") == 0 || strcmp(argv[0],"Param") == 0 ||
        strcmp(argv[0],"parameters") == 0 || strcmp(argv[0],"Parameters") == 0)
    {
        output.tag("ResponseType","Fcn");
        output.tag("ResponseType","Fcrn");
        output.tag("ResponseType","Kv");
        
        theResponse = new ElementResponse(this, 8, Vector(3));
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int HDR::getResponse(int responseID, Information &eleInfo)
{
    double kGeo1, MpDelta1, MpDelta2, MpDelta3, MpDelta4, MpDelta5, MpDelta6;
    Vector evol(5), kbVec(4), Param(3);
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta moments
        kGeo1 = 0.5*qb(0);
        MpDelta1 = kGeo1*(ul(7)-ul(1));
        theVector(5)  += MpDelta1;
        theVector(11) += MpDelta1;
        MpDelta2 = kGeo1*shearDistI*L*ul(5);
        theVector(5)  += MpDelta2;
        theVector(11) -= MpDelta2;
        MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(11);
        theVector(5)  -= MpDelta3;
        theVector(11) += MpDelta3;
        MpDelta4 = kGeo1*(ul(8)-ul(2));
        theVector(4)  -= MpDelta4;
        theVector(10) -= MpDelta4;
        MpDelta5 = kGeo1*shearDistI*L*ul(4);
        theVector(4)  += MpDelta5;
        theVector(10) -= MpDelta5;
        MpDelta6 = kGeo1*(1.0 - shearDistI)*L*ul(10);
        theVector(4)  -= MpDelta6;
        theVector(10) += MpDelta6;
        return eleInfo.setVector(theVector);
        
    case 3:  // basic forces
        return eleInfo.setVector(qb);
        
    case 4:  // local displacements
        return eleInfo.setVector(ul);
        
    case 5:  // basic displacements
        return eleInfo.setVector(ub);
        
    case 6:  // hysteretic, scragging, and Mullins damage parameters
        evol(0) = DSplus;
        evol(1) = DSminus;
        evol(2) = DS;
        evol(3) = DM;
        evol(4) = Delta; 
        return eleInfo.setVector(evol);
        
    case 7:  // basic stiffness
        kbVec(0) = kb(1,1); kbVec(1) = kb(1,2);
        kbVec(2) = kb(2,1); kbVec(3) = kb(2,2);
        return eleInfo.setVector(kbVec);
        
    case 8:  // parameters that vary with time
        Param(0) = Fcn;
        Param(1) = Fcrn;
        Param(2) = Kv;
        return eleInfo.setVector(Param);
        
    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void HDR::setUp()
{
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    if (L > DBL_EPSILON)  {
        if (x.Size() == 0)  {
            x.resize(3);
            x = xp;
        } /*else  {
          opserr << "WARNING ElastomericX::setUp() - "
          << "element: " << this->getTag() << endln
          << "ignoring nodes and using specified "
          << "local x vector to determine orientation\n";
          }*/
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "HDR::setUp() - "
            << "element: " << this->getTag() << endln
            << "incorrect dimension of orientation vectors\n";
        exit(-1);
    }
    
    // establish orientation of element for the transformation matrix
    // z = x cross y
    Vector z(3);
    z(0) = x(1)*y(2) - x(2)*y(1);
    z(1) = x(2)*y(0) - x(0)*y(2);
    z(2) = x(0)*y(1) - x(1)*y(0);
    
    // y = z cross x
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);
    
    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();
    
    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0)  {
        opserr << "HDR::setUp() - "
            << "element: " << this->getTag() << endln
            << "invalid orientation vectors\n";
        exit(-1);
    }
    
    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = x(0)/xn;
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = x(1)/xn;
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = x(2)/xn;
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = y(0)/yn;
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = y(1)/yn;
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = y(2)/yn;
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = z(0)/zn;
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = z(1)/zn;
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = z(2)/zn;
    
    // create transformation matrix from local to basic system (linear)
    Tlb.Zero();
    Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
    Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;
    Tlb(1,5) = -shearDistI*L;
    Tlb(1,11) = -(1.0 - shearDistI)*L;
    Tlb(2,4) = -Tlb(1,5);
    Tlb(2,10) = -Tlb(1,11);
}


double HDR::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
