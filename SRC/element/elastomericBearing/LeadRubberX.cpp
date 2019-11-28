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
// Credits: This element extends the formulation of elastomericBearing element written by Andreas Schellenberg 
// Created: 02/12
//
// Description: This file contains the implementation of the
// LeadRubberX class.

#include "LeadRubberX.h"

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
Matrix LeadRubberX::theMatrix(12,12);
Vector LeadRubberX::theVector(12);


static int numMyBearing = 0;
static int tag = 0;  // Tag to identify if bearing has failed in buckling
void *OPS_LeadRubberX()
{
    // print out a message about who wrote this element & any copyright info wanted
    if (numMyBearing == 0) {
        opserr << "LeadRubberX element - Written by Manish Kumar, University at Buffalo, 2012\n";
        numMyBearing++;
    }
    
    Element *theEle = 0;
    
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs == 0) { // parallel processing
        theEle = new LeadRubberX();
        return theEle;
    }
    
    if (numArgs !=12 && numArgs !=18 && numArgs !=19 && numArgs !=20
        && numArgs !=24 && numArgs !=25 && numArgs !=29 && numArgs !=30
        && numArgs !=31 && numArgs !=32 && numArgs !=33 && numArgs !=34) {
            opserr << "ERROR - LeadRubberX incorrect # args provided";
            return theEle;
    }
    
    // get the id and end nodes
    int iData[3];
    double dData[9];
    int numData;
    
    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data\n";
        return 0;
    }
    
    int eleTag = iData[0];
    
    numData = 9;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element properties for element" << eleTag << endln;
        return 0;
    }
    
    // get the orientation vector
    Vector x(0);
    Vector y(3); y(0) = -1.0; y(1) = 0.0; y(2) = 0.0;
    
    // get the tags of the properties
    int tag1 = 0;           // Cavitation and post-cavitation
    int tag2 = 0;           // Buckling load variation
    int tag3 = 0;           // Horizontal stiffness variation
    int tag4 = 0;           // Vertical stiffness variation
    int tag5 = 0;           // Strength degradation due to heating 
    
    // The default values of the parameters
    double kl = 10.0;       // Cavitation parameter
    double phi = 0.5;       // Damage index
    double al = 1.0;        // Strength degradation parameter
    double sDratio = 0.5;   // Shear distance ratio
    double m = 0.0;         // Mass of the bearing
    double cd1 = 0.0;       // Viscous damping parameter
    double tc1 = 0.0;       // Cover thickness
    double qL1 = 11200.0;   // Density of lead in kg/m3
    double cL1 = 130.0;     // Specific heat of lead in N-m/kg oC
    double kS1 = 50.0;      // Thermal conductivity of steel m2/s
    double aS1 = 1.41e-05;  // Diffusitivity 
    
    if (numArgs >= 18) {
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
        if (numArgs >= 19) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &kl) != 0) {
                opserr << "WARNING error reading element property cavitation parameter for element" << eleTag << endln;
                return 0;
            }
            if (numArgs >= 20) {
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &phi) != 0) {
                    opserr << "WARNING error reading element property damage index for element" << eleTag << endln;
                    return 0;
                }
                if (numArgs >= 21) {
                    numData = 1;
                    if (OPS_GetDoubleInput(&numData, &al) != 0) {
                        opserr << "WARNING error reading element property strength degradation parameter for element" << eleTag << endln;
                        return 0;
                    }
                    if (numArgs >= 22) {
                        numData = 1;
                        if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
                            opserr << "WARNING error reading element property shear distance ratio for element" << eleTag << endln;
                            return 0;
                        }
                        if (numArgs >= 23) {
                            numData = 1;
                            if (OPS_GetDoubleInput(&numData, &m) != 0) {
                                opserr << "WARNING error reading element property mass for element" << eleTag << endln;
                                return 0;
                            }
                            if (numArgs >= 24) {
                                numData = 1;
                                if (OPS_GetDoubleInput(&numData, &cd1) != 0) {
                                    opserr << "WARNING error reading element property viscous damping parameter for element" << eleTag << endln;
                                    return 0;
                                }
                                if (numArgs >= 25) {
                                    numData = 1;
                                    if (OPS_GetDoubleInput(&numData, &tc1) != 0) {
                                        opserr << "WARNING error reading element property cover thickness for element" << eleTag << endln;
                                        return 0;
                                    }
                                    if (numArgs >= 29) {
                                        numData = 1;
                                        if (OPS_GetDoubleInput(&numData, &qL1) != 0) {
                                            opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                            return 0;
                                        }
                                        if (OPS_GetDoubleInput(&numData, &cL1) != 0) {
                                            opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                            return 0;
                                        }
                                        if (OPS_GetDoubleInput(&numData, &kS1) != 0) {
                                            opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                            return 0;
                                        }
                                        if (OPS_GetDoubleInput(&numData, &aS1) != 0) {
                                            opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                            return 0;
                                        }
                                        if (numArgs >= 30) {
                                            numData = 1;
                                            if (OPS_GetIntInput(&numData, &tag1) != 0) {
                                                opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                                return 0;
                                            }
                                            if (numArgs >= 31) {
                                                numData = 1;
                                                if (OPS_GetIntInput(&numData, &tag2) != 0) {
                                                    opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                                    return 0;
                                                }
                                                if (numArgs >= 32) {
                                                    numData = 1;
                                                    if (OPS_GetIntInput(&numData, &tag3) != 0) {
                                                        opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                                        return 0;
                                                    }
                                                    if (numArgs >= 33) {
                                                        numData = 1;
                                                        if (OPS_GetIntInput(&numData, &tag4) != 0) {
                                                            opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                                            return 0;
                                                        }
                                                        if (numArgs == 34) {
                                                            numData = 1;
                                                            if (OPS_GetIntInput(&numData, &tag5) != 0) {
                                                                opserr << "WARNING error reading element properties for element" << eleTag << endln;
                                                                return 0;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
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
            opserr << ", for space problem need 6 - LeadRubberX \n"; 
        }
        theEle = new LeadRubberX(iData[0], iData[1], iData[2], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], y, x, kl, phi, al, sDratio, m, cd1, tc1, qL1, cL1, kS1, aS1, tag1, tag2, tag3, tag4, tag5);
    }
    
    if (theEle == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }
    
    return theEle;
}


LeadRubberX::LeadRubberX(int eleTag, int Nd1, int Nd2, double qd, double alpha1, double Gr, double kbulk, double Di, double Do,
    double ts1, double tr, double n1, const Vector _y, const Vector _x, double kl, double PhiMax, double al, double sDratio,
    double m, double cd1, double tc1, double qL2, double cL2, double kS2, double aS2, int tg1, int tg2, int tg3, int tg4, int tg5)
    : Element(eleTag, ELE_TAG_LeadRubberX), connectedExternalNodes(2),
    qYield0(qd), alpha(alpha1), cd(cd1), TL_trial(0.0), TL_commit(0.0), qL(qL2), cL(cL2),
    kS(kS2), aS(aS2), PhiM(PhiMax), ac(al), tCurrent(0.0), tCommit(0.0), G(Gr), Kbulk(kbulk),
    x(_x), y(_y), tag1(tg1), tag2(tg2), tag3(tg3), tag4(tg4), tag5(tg5),
    shearDistI(sDratio), mass(m), tc(tc1), D1(Di), D2(Do), L(0.0), n(n1), ts(ts1),
    Fcrn(0.0), ucrn(0.0), Fcrmin(0.0), Fcn(0.0), ucn(0.0), Fmax(0.0), umax(0.0),
    ub(6), ubdot(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "LeadRubberX::LeadRubberX() - element: "
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
    S = (D2*D2 - D1*D1)/(4*D2*tr);                          // Shape factor for case with lead core
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
    double Is = I*h/Tr;                                     // Adjusted moment of intertia of bearing
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
    
    // horizontal motion
    qYield = qYield0;                                       // This yield stress changes with time
    ke = G*A/Tr;                                            // Stiffness of elastic component (due to rubber)
    k0 = (1.0/alpha - 1.0)*ke;                              // Initial stiffness of hysteretic component (due to lead)
    
    // torsion
    Kt = G*(2*Is)/h;                                        // torsional stiffness of bearing
    
    // rotation
    Kr = Er*Is/h;                                           // rotational stiffness of bearing
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = k0 + ke;
    kbInit(2,2) = k0 + ke;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
}


LeadRubberX::LeadRubberX()
    : Element(0, ELE_TAG_LeadRubberX), connectedExternalNodes(2),
    k0(0.0), qYield0(0.0), qYield(0.0), alpha(0.0), ke(0.0), cd(0.0),
    TL_trial(0.0), TL_commit(0.0), qL(0.0), cL(0.0), kS(0.0), aS(0.0),
    S(0.0), Ec(0.0), Kv0(0.0), Kv(0.0), kc(10.0), PhiM(0.5), ac(1.0), Fcr(0.0),
    ucr(0.0), Fc(0.0), uc(0.0), tCurrent(0.0), tCommit(0.0), Kt(0.0), Kr(0.0),
    G(0.0), Kbulk(0.0), x(0), y(0), tag1(0), tag2(0), tag3(0), tag4(0), tag5(0),
    shearDistI(0.5), mass(0.0), tc(0.0), Tr(0.0), D1(0.0), D2(0.0), L(0.0), h(0.0),
    rg(0.0), A(0.0), Ar(0.0), n(0.0), ts(0.0), Fcrn(0.0), ucrn(0.0), Fcrmin(0.0),
    Fcn(0.0), ucn(0.0), Fmax(0.0), umax(0.0),
    ub(6), ubdot(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "LeadRubberX::LeadRubberX() - "
            <<  "failed to create an ID of size 2\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
}


LeadRubberX::~LeadRubberX()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
}


int LeadRubberX::getNumExternalNodes() const
{
    return 2;
}


const ID& LeadRubberX::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** LeadRubberX::getNodePtrs()
{
    return theNodes;
}


int LeadRubberX::getNumDOF()
{
    return 12;
}


void LeadRubberX::setDomain(Domain *theDomain)
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
            opserr << "WARNING LeadRubberX::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING LeadRubberX::setDomain() - Nd2: "
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
        opserr << "LeadRubberX::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "LeadRubberX::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int LeadRubberX::commitState()
{
    int errCode = 0;
    
    double uh = sqrt(ub(1)*ub(1) + ub(2)*ub(2));
    
    // vertical motion
    if (tag4 == 1) {
        Kv = Kv0*(1.0/(1.0+(3.0/(PI*PI))*(uh/rg)*(uh/rg)));
        if (uh > DBL_EPSILON)
            uc = Fc/Kv;
    }
    
    // tension
    if (tag1 == 1) {
        if (ub(0) > umax) {
            umax = ub(0);
            Fcn = Fc*(1.0-PhiM*(1.0-exp(-ac*(ub(0)-uc)/uc)));
        }
    }
    
    // compression
    if (tag2 == 1) {
        double Delta = 2.0*acos(uh/D2);   // this becomes undefined for uh/D2 > 1.0
        //Ar = (D2*D2/4.0)*(Delta-sin(Delta));
        Ar = ((D2+tc)*(D2+tc) - D1*D1)/4.0*(Delta-sin(Delta));  // A does not include lead core
        if (Ar/A < 0.2 || uh/D2 >= 1.0) {
            Fcrn = 0.2*Fcr;
        } else {
            Fcrn = Fcr*Ar/A;
        }
        
        if (Fcrn > Fcrmin)
            Fcrmin = Fcrn;
        ucrn = Fcrn/Kv;
    }
    
    // horizontal motion
    if (tag3 == 1) {
        ke = (G*A/Tr)*(1.0-pow(qb(0)/Fcrn,2));
        //if (ke < 0) {
        //    ke = 0.01*(G*A/Tr);  // a fraction of ke to avoid convergence issues
        //    opserr << "WARNING LeadRubberX::commitState() - Negative horizontal stiffness\n";
        //}
    }
    
    // lead core heating
    TL_commit = TL_trial;
    tCommit = (this->getDomain())->getCurrentTime();
    if (tag5 == 1) {
        qYield = qYield0*exp(-0.0069*TL_commit);
    }
    
    // commit trial history variables for horizontal direction
    ubC = ub;
    zC = z;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int LeadRubberX::revertToLastCommit()
{
    return 0;
}


int LeadRubberX::revertToStart()
{
    int errCode = 0;
    
    // reset trial history variables
    ub.Zero();
    z.Zero();
    qb.Zero();
    
    // reset committed history variables
    ubC.Zero();
    zC.Zero();
    
    // reset tangent of hysteretic evolution parameters
    dzdu(0,0) = dzdu(1,1) = k0/qYield;
    dzdu(1,0) = dzdu(0,1) = 0.0;
    
    // reset stiffness matrix in basic system
    kb = kbInit;
    
    return errCode;
}


int LeadRubberX::update()
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
    
    // get the current temperature of the lead core
    double v = sqrt(ubdot(1)*ubdot(1) + ubdot(2)*ubdot(2));
    TL_trial = getCurrentTemp(qYield, TL_commit, v);
    
    // 1) get axial force and stiffness in basic x-direction
    double qTrial = Kv*ub(0);
    ucn = Fcn/Kv;
    Fmax = Fc*(1.0+(1.0/(Tr*kc))*(1.0-exp(-kc*(umax-uc))));
    
    if (ub(0) <= ucrn) {
        if (tag2 == 1) {
            kb(0,0) = Kv0/10000.0;
            qb(0) = Fcrmin + kb(0,0)*(ub(0)-ucrn);
            // tag=1;
        } else {
            kb(0,0) = Kv;
            qb(0) = kb(0,0)*ub(0);
        }
    }
    
    if (ub(0) > ucrn) {
        if (tag1 == 1) {
            if (ub(0) <= ucn) {
                kb(0,0) = Kv;
                qb(0) = kb(0,0)*ub(0);
            } else if (ub(0) < umax) {
                kb(0,0) = (Fmax-Fcn)/(umax-ucn);
                qb(0) = Fcn + (Fmax-Fcn)/(umax-ucn)*(ub(0)-ucn);
            } else {
                kb(0,0) = (Fc/Tr)*exp(-kc*(ub(0)-uc));
                qb(0) = Fc*(1.0+(1.0/(Tr*kc))*(1.0-exp(-kc*(ub(0)-uc))));
            }
        } else {
            kb(0,0) = Kv;
            qb(0) = kb(0,0)*ub(0);
        }
    }
    
    //2) calculate shear forces and stiffnesses in basic y- and z-direction
    // get displacement increments (trial-committed)
    Vector delta_ub = ub - ubC;
    if (sqrt(pow(delta_ub(1),2)+pow(delta_ub(2),2)) > 0.0)  {
        
        // get yield displacement
        double uy = qYield/k0;
        
        // calculate hysteretic evolution parameter z using Newton-Raphson
        int iter = 0;
        int maxIter = 100;
        double tol = 1E-8;
        double beta = 0.1; // note here beta and gamma are as per Nagarajaih(1991), which is opposite to Park et al.(1986)
        double gamma = 0.9;
        double zNrm, tmp1, tmp2, tmp3;
        Vector f(2), delta_z(2);
        Matrix Df(2,2);
        do  {
            zNrm = z.Norm();
            if (zNrm == 0.0)  // check because of negative exponents
                zNrm = DBL_EPSILON;
            tmp1 = beta + gamma*sgn(z(0)*delta_ub(1));
            tmp2 = beta + gamma*sgn(z(1)*delta_ub(2));
            tmp3 = z(0)*delta_ub(1)*tmp1 + z(1)*delta_ub(2)*tmp2;
            
            // function and derivative
            f(0) = z(0) - zC(0) - 1.0/uy*(delta_ub(1) - z(0)*tmp3);
            f(1) = z(1) - zC(1) - 1.0/uy*(delta_ub(2) - z(1)*tmp3);
            
            Df(0,0) = 1.0 + (1.0/uy)*(2*z(0)*delta_ub(1)*tmp1+z(1)*delta_ub(2)*tmp2);
            Df(1,0) = (tmp1/uy)*z(1)*delta_ub(1);
            Df(0,1) = (tmp2/uy)*z(0)*delta_ub(2);
            Df(1,1) = 1.0 + (1.0/uy)*(z(0)*delta_ub(1)*tmp1+2*z(1)*delta_ub(2)*tmp2);
            
            // issue warning if diagonal of derivative Df is zero
            if ((fabs(Df(0,0)) <= DBL_EPSILON) || (fabs(Df(1,1)) <= DBL_EPSILON))  {
                opserr << "WARNING: LeadRubberX::update() - "
                    << "zero Jacobian in Newton-Raphson scheme for hysteretic "
                    << "evolution parameter z.\n";
                return -1;
            }
            
            // advance one step
            // delta_z = f/Df; either write a function to do matrix division or use the solution below
            delta_z(0) = (f(0)*Df(1,1)-f(1)*Df(0,1))/(Df(0,0)*Df(1,1)-Df(0,1)*Df(1,0));
            delta_z(1) = (f(0)*Df(1,0)-f(1)*Df(0,0))/(Df(0,1)*Df(1,0)-Df(0,0)*Df(1,1));
            z -= delta_z;
            iter++;
        } while ((delta_z.Norm() >= tol) && (iter < maxIter));
        
        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter)   {
            opserr << "WARNING: LeadRubberX::update() - "
                << "did not find the hysteretic evolution parameters z after "
                << iter << " iterations and norm: " << delta_z.Norm() << endln;
            return -2;
        }
        
        // get derivative of hysteretic evolution parameter
        double du1du2 = delta_ub(1)/delta_ub(2);
        double du2du1 = delta_ub(2)/delta_ub(1);
        if (delta_ub(1)*delta_ub(2) == 0)  {
            du1du2 = 0.0;
            du2du1 = 0.0;
        }
        
        dzdu(0,0) = (1.0/uy)*(1.0 - z(0)*(z(0)*tmp1+z(1)*tmp2*du2du1));
        dzdu(0,1) = (1.0/uy)*(du1du2-z(0)*(z(0)*tmp1*du1du2+z(1)*tmp2));
        dzdu(1,0) = (1.0/uy)*(du2du1-z(1)*(z(0)*tmp1+z(1)*tmp2*du2du1));
        dzdu(1,1) = (1.0/uy)*(1.0 - z(1)*(z(0)*tmp1*du1du2+z(1)*tmp2));
        
        tCurrent = (this->getDomain())->getCurrentTime();
        double dT = tCurrent - tCommit;
        
        // set shear force
        qb(1) = cd*ubdot(1) + qYield*z(0) + ke*ub(1);
        qb(2) = cd*ubdot(2) + qYield*z(1) + ke*ub(2);
        // set tangent stiffness
        kb(1,1) = cd/dT + qYield*dzdu(0,0) + ke;
        kb(1,2) = qYield*dzdu(0,1);
        kb(2,1) = qYield*dzdu(1,0);
        kb(2,2) = cd/dT + qYield*dzdu(1,1) + ke;
    }
    
    // 3) get moment and stiffness in basic x-direction
    qb(3) = Kt*ub(3);
    kb(3,3) = Kt;
    
    // 4) get moment and stiffness in basic y-direction
    qb(4) = Kr*ub(4);
    kb(4,4) = Kr;
    
    // 5) get moment and stiffness in basic z-direction
    qb(5) = Kr*ub(5);
    kb(5,5) = Kr;
    
    /* if buckling 
    if (tag==1) {
        tag = 0;
        return -1;  // return any negative integer
    }*/
    
    return 0;
}


const Matrix& LeadRubberX::getTangentStiff()
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


const Matrix& LeadRubberX::getInitialStiff()
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


const Matrix& LeadRubberX::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    return theMatrix;
}


const Matrix& LeadRubberX::getMass()
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


void LeadRubberX::zeroLoad()
{
    theLoad.Zero();
}


int LeadRubberX::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"LeadRubberX::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << endln;
    
    return -1;
}


int LeadRubberX::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
        opserr << "LeadRubberX::addInertiaLoadToUnbalance() - "
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


const Vector& LeadRubberX::getResistingForce()
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


const Vector& LeadRubberX::getResistingForceIncInertia()
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


int LeadRubberX::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(28);
    data(0)  = this->getTag();
    data(1)  = qYield0;
    data(2)  = alpha;
    data(3)  = G;
    data(4)  = Kbulk;
    data(5)  = D1;
    data(6)  = D2;
    data(7)  = ts;
    data(8)  = Tr;
    data(9)  = n;
    data(10) = x.Size();
    data(11) = y.Size();
    data(12) = kc;
    data(13) = PhiM;
    data(14) = ac;
    data(15) = shearDistI;
    data(16) = mass;
    data(17) = cd;
    data(18) = tc;
    data(19) = qL;
    data(20) = cL;
    data(21) = kS;
    data(22) = aS;
    data(23) = tag1;
    data(24) = tag2;
    data(25) = tag3;
    data(26) = tag4;
    data(27) = tag5;
    
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


int LeadRubberX::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // receive element parameters
    static Vector data(28);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    
    qYield0 = data(1);
    alpha = data(2);
    G = data(3);
    Kbulk = data(4);
    D1 = data(5);
    D2 = data(6);
    ts = data(7);
    Tr = data(8);
    n = data(9);
    kc = data(12);
    PhiM = data(13);
    ac = data(14);
    shearDistI = data(15);
    mass = data(16);
    cd = data(17);
    tc = data(18);
    qL = data(19);
    cL = data(20);
    kS = data(21);
    aS = data(22);
    tag1 = (int)data(23);
    tag2 = (int)data(24);
    tag3 = (int)data(25);
    tag4 = (int)data(26);
    tag5 = (int)data(27);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
    // receive remaining data
    if ((int)data(10) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(11) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    
    // vertical motion
    A = (PI/4.0)*((D2+tc)*(D2+tc) - D1*D1);
    S = (D2*D2 - D1*D1)/(4*D2*Tr/n);
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
    
    // horizontal motion
    qYield = qYield0;
    ke = G*A/Tr;
    k0 = (1.0/alpha - 1.0)*ke;
    
    // torsion
    Kt = G*(2*Is)/h;
    
    // rotation
    Kr = Er*Is/h;
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = k0 + ke;
    kbInit(2,2) = k0 + ke;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
    
    return 0;
}


int LeadRubberX::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    
    static Vector v1(3);
    static Vector v2(3);
    
    if (displayMode >= 0)  {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i=0; i<3; i++)  {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v2(i) = end2Crd(i) + end2Disp(i)*fact;
        }
    } else  {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        
        if (eigen1.noCols() >= mode)  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;
            }
        } else  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i);
                v2(i) = end2Crd(i);
            }
        }
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag(), 0);
}


void LeadRubberX::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "************************************************************" << endln;
        s << "Element: " << this->getTag();
        s << "  type: LeadRubberX  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "************************************************************" << endln;
        s << "GEOMETRIC PROPERTIES" << endln;
        s << "D1: " << D1 << " D2: " << D2 << " L: " << L << " Tr: " << Tr << " S: " << S << " A: " << A << endln;
        s << "MATERIAL PROPERTIES" << endln;
        s << "G: " << G << " kc: " << kc << " ac: " << ac << " PhiM: " << PhiM << " shearDistI: " << shearDistI << " mass: " << mass << endln;
        s << " qL: " << qL << " cL: " << cL << " kS: " << kS << " aS: " << aS << endln;
        s << "MECHANICAL PROPERTIES: HORIZONTAL MOTION" << endln;
        s << "k0: " << k0 << " ke: " << ke << " qYield: " << qYield << " DeltaT: " << TL_commit << " Fcrmin: " << Fcrmin << endln;
        s << "MECHANICAL PROPERTIES: VERTICAL MOTION" << endln;
        s << "Kv: "<< Kv << " Fc: " << Fc << " Fcr: " << Fcr << " Fcn: " << Fcn << " umax: " << umax << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
        s << "************************************************************" << endln;
        //s <<"time: " << tCommit <<" ke0: " << G*A/Tr  <<" ke: " << ke <<" Fcr: "<< Fcr << " Fcrmin: "<< Fcrmin <<" Kv0: "<< Kv0<<" Kv: "<< Kv << " qYield0: " << qYield0 << " qYield: " << qYield << " DeltaT: " << TL_commit << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"LeadRubberX\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"D1\": " << D1 << ", ";
        s << "\"D2\": " << D2 << ", ";
        s << "\"L\": " << L << ", ";
        s << "\"Tr\": " << Tr << ", ";
        s << "\"S\": " << S << ", ";
        s << "\"A\": " << A << ", ";
        s << "\"G\": " << G << ", ";
        s << "\"kc\": " << kc << ", ";
        s << "\"ac\": " << ac << ", ";
        s << "\"PhiM\": " << PhiM << ", ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"mass\": " << mass << ", ";
        s << "\"qL\": " << qL << ", ";
        s << "\"cL\": " << cL << ", ";
        s << "\"kS\": " << kS << ", ";
        s << "\"aS\": " << aS << "}";
    }
}


Response* LeadRubberX::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","LeadRubberX");
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
    // hysteretic evolution parameter
    else if (strcmp(argv[0],"hystereticParameter") == 0 || strcmp(argv[0],"hystParameter") == 0 ||
        strcmp(argv[0],"hystereticParam") == 0 || strcmp(argv[0],"hystParam") == 0 ||
        strcmp(argv[0],"z") == 0)
    {
        output.tag("ResponseType","z1");
        output.tag("ResponseType","z2");
        
        theResponse = new ElementResponse(this, 6, Vector(2));
    }
    // dzdu
    else if (strcmp(argv[0],"dzdu") == 0)
    {
        output.tag("ResponseType","dz1du1");
        output.tag("ResponseType","dz1du2");
        output.tag("ResponseType","dz2du1");
        output.tag("ResponseType","dz2du2");
        
        theResponse = new ElementResponse(this, 7, Vector(4));
    }
    // basic stiffness
    else if (strcmp(argv[0],"kb") == 0 || strcmp(argv[0],"basicStiff") == 0 ||
        strcmp(argv[0],"basicStiffness") == 0)
    {
        output.tag("ResponseType","kb22");
        output.tag("ResponseType","kb23");
        output.tag("ResponseType","kb32");
        output.tag("ResponseType","kb33");
        
        theResponse = new ElementResponse(this, 8, Vector(4));
    }
    // parameters
    else if (strcmp(argv[0],"param") == 0 || strcmp(argv[0],"Param") == 0 ||
        strcmp(argv[0],"parameters") == 0 || strcmp(argv[0],"Parameters") == 0)
    {
        output.tag("ResponseType","Fcn");
        output.tag("ResponseType","Fcrn");
        output.tag("ResponseType","Kv");
        output.tag("ResponseType","ke");
        output.tag("ResponseType","DeltaT");
        output.tag("ResponseType","qYield");
        
        theResponse = new ElementResponse(this, 9, Vector(6));
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int LeadRubberX::getResponse(int responseID, Information &eleInfo)
{
    double kGeo1, MpDelta1, MpDelta2, MpDelta3, MpDelta4, MpDelta5, MpDelta6;
    Vector dzduVec(4), kbVec(4), Param(6);
    
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
        
    case 6:  // hysteretic evolution parameter
        return eleInfo.setVector(z);
        
    case 7:  // dzdu
        dzduVec(0) = dzdu(0,0); dzduVec(1) = dzdu(0,1);
        dzduVec(2) = dzdu(1,0); dzduVec(3) = dzdu(1,1);
        return eleInfo.setVector(dzduVec);
        
    case 8:  // basic stiffness
        kbVec(0) = kb(1,1); kbVec(1) = kb(1,2);
        kbVec(2) = kb(2,1); kbVec(3) = kb(2,2);
        return eleInfo.setVector(kbVec);
        
    case 9:  // parameters that vary with time
        Param(0) = Fcn;
        Param(1) = Fcrn;
        Param(2) = Kv;
        Param(3) = ke;
        Param(4) = TL_commit;
        Param(5) = qYield;
        return eleInfo.setVector(Param);
        
    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void LeadRubberX::setUp()
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
          opserr << "WARNING LeadRubberX::setUp() - "
          << "element: " << this->getTag() << endln
          << "ignoring nodes and using specified "
          << "local x vector to determine orientation\n";
          }*/
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "LeadRubberX::setUp() - "
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
        opserr << "LeadRubberX::setUp() - "
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


double LeadRubberX::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}


double LeadRubberX::getCurrentTemp(double qYield, double TL_commit, double v)
{
    // lead core heating
    tCurrent = (this->getDomain())->getCurrentTime();
    if (tCurrent < tCommit) {
        tCommit = 0.0;
    }
    
    double a = D1/2;
    double dT = tCurrent - tCommit;
    // if (dT>1) dT = 0;
    double ALead = PI*pow(a,2);
    double tau = (aS*tCurrent)/(pow(a,2));
    double F;
    if (tau < 0.6) {
        F = 2.0*sqrt(tau/PI)-(tau/PI)*(2.0-(tau/4.0)-pow(tau/4.0,2)-(15.0/4.0)*(pow(tau/4.0,3)));
    } else {
        F = 8.0/(3.0*PI)-(1.0/(2.0*sqrt(PI*tau)))*(1.0-(1.0/(12.0*tau))+(1.0/(6.0*pow(4.0*tau,2)))-(1.0/(12.0*pow(4.0*tau,3))));
    }
    double deltaT1 = (dT/(qL*cL*h))*((qYield*v*zC.Norm())/ALead-(kS*TL_commit/a)*(1.0/F+1.274*((n-1)*ts/a)*pow(tau,-1.0/3.0))); 
    if (deltaT1 <= 0.0) {
        deltaT1 = 0.0;
    }
    
    // use improved euler method to obtain final temperature
    double TL_trial1 = TL_commit + deltaT1;
    double tCurrent2 = tCurrent + dT;
    tau = (aS*tCurrent2)/(pow(a,2));
    if (tau < 0.6) {
        F = 2.0*sqrt(tau/PI)-(tau/PI)*(2.0-(tau/4.0)-pow(tau/4.0,2)-(15.0/4.0)*(pow(tau/4.0,3)));
    } else {
        F = 8.0/(3.0*PI)-(1.0/(2.0*sqrt(PI*tau)))*(1.0-(1.0/(12.0*tau))+(1.0/(6.0*pow(4.0*tau,2)))-(1.0/(12.0*pow(4.0*tau,3))));
    }
    double deltaT2 = (dT/(qL*cL*h))*((qYield*v*zC.Norm())/ALead-(kS*TL_trial1/a)*(1.0/F+1.274*((n-1)*ts/a)*pow(tau,-1.0/3.0))); 
    if (deltaT2 <= 0.0) {
        deltaT2 = 0.0;
    }
    
    double TL_trial = TL_commit + 0.5*(deltaT1+deltaT2);
    
    return TL_trial;
}
