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
// $URL$

// TripleFrictionPendulumX element
// Written by Hyunmyung Kim (hkim59@buffalo.edu) and Michael C. Constantinou (constan1@buffalo.edu)
// Created: 2021/11
// Description: This command is used to construct the TripleFrictionPendulumX element (Kim and Constantinou, 2022,2023) object, 
// which is an extension of the TripleFrictionPendulum element (Dao et al., 2013) capable of accounting for heating effects on 
// the frictional behavior of triple Friction Pendulum isolators. The horizontal behavior of the element is achieved by the series model, 
// which consists of properly combined hysteretic/frictional and multidirectional gap elements. 

// Three main modifications in the TripleFrictionPendulumX element include : 
// 1) computation of the displacement and velocity histories at each of the four sliding interfaces of the isolator, 
// 2) computation of the temperature history at each sliding interface, and 
// 3) accounting for the dependency of the coefficient of friction on the instantaneous temperature at each sliding interface.



#include "TripleFrictionPendulumX.h"
#include <MovableObject.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>
#include <elementAPI.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <G3Globals.h>
#include <Message.h>
using namespace std;
#include <iostream>
#include <Vector.h>
#include <SingleFPSimple3d.h>


// initialize the class wide variables
Matrix TripleFrictionPendulumX::eleK(12, 12);
Matrix TripleFrictionPendulumX::eleKinit(12, 12);
Matrix TripleFrictionPendulumX::eleM(12, 12);
Matrix TripleFrictionPendulumX::eleD(12, 12);
Vector TripleFrictionPendulumX::eleR(12);

static int numTripleFrictionPendulumX = 0;


void* OPS_TripleFrictionPendulumX()
{
    if (numTripleFrictionPendulumX == 0) {
        numTripleFrictionPendulumX++;
        opserr << "TripleFrictionPendulumX \n";
    }

    // get the id and end nodes 
    int iData[11]; 
    double dData[26];
    int numData, eleTag;

    numData = 11;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data";
        return 0;
    }
    eleTag = iData[0];

    
    // get the uniaxial materials
    UniaxialMaterial* theMaterials[4];
    for (int i = 0; i < 4; i++) {
        theMaterials[i] = OPS_getUniaxialMaterial(iData[4 + i]); 
        if (theMaterials[i] == 0) {
            opserr << "WARNING uniaxial material not found\n";
            opserr << "uniaxialMaterial: " << iData[4 + i] << endln;
            opserr << "TripleFrictionPendulumX element: " << eleTag << endln;
            return 0;
        }
    }

    

    numData = 26;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element" << eleTag << endln;
        return 0;
    }

    // create the element and add it to the Domain
    Element* theTripleFrictionPendulumX = new TripleFrictionPendulumX(eleTag, iData[1], iData[2], iData[3], theMaterials, iData[8], iData[9], iData[10], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20], dData[21], dData[22], dData[23], dData[24], dData[25]);

    if (theTripleFrictionPendulumX == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }

    return theTripleFrictionPendulumX;
}


// typical constructor
TripleFrictionPendulumX::TripleFrictionPendulumX(int tag,
    int Nd1, int Nd2, int tg,
    UniaxialMaterial** materials,
    int Pfactor, int Tfactor, int Vfactor,
    double muref1,
    double muref2,
    double muref3,
    double l1,
    double l2,
    double l3,
    double ubar1,
    double ubar2,
    double ubar3,
    double b1,
    double b2,
    double b3,
    double w,
    double uy,
    double kvt,
    double minFv,
    double tol,
    double refp1,
    double refp2,
    double refp3,
    double diffu,
    double conduct,
    double t0,
    double velPara,
    double tempPara,
    double PVunit
    
    )
    : Element(tag, ELE_TAG_TripleFrictionPendulumX), externalNodes(2), tag1(tg), kpFactor(Pfactor), kTFactor(Tfactor), kvFactor(Vfactor),
    Mu_ref1(muref1), Mu_ref2(muref2), Mu_ref3(muref3), L1(l1), L2(l2), L3(l3), Ubar1(ubar1), Ubar2(ubar2), Ubar3(ubar3), B1(b1), B2(b2), B3(b3),
    W(w), Uy(uy), Kvt(kvt), MinFv(minFv), TOL(tol), refPressure1(refp1), refPressure2(refp2), refPressure3(refp3), Diffusivity(diffu), Conductivity(conduct), Temperature0(t0), rateParam(velPara), tempParam(tempPara), unit(PVunit), Niter(20),
    K(2, 2), Kpr(2, 2), f(2), fpr(2),
    k12(2, 2), k12pr(2, 2), k34(2, 2), k34pr(2, 2), k56(2, 2), k56pr(2, 2),
    d1(2), d1pr(2), d3(2), d3pr(2), d5(2), d5pr(2), v1(2), v3(2), v5(2),
    ep1(2), ep1pr(2), ep3(2), ep3pr(2), ep5(2), ep5pr(2),
    q1(2), q1pr(2), q3(2), q3pr(2), q5(2), q5pr(2),
    ep1tmp(2), ep3tmp(2), ep5tmp(2), q1tmp(2), q3tmp(2), q5tmp(2)
    , d1ppr(2), d3ppr(2), d5ppr(2)
    , DomainTime(1), DomainHeatFlux1(1), DomainHeatFlux2(1), DomainHeatFlux3(1), DomainTimeTemp(0), DomainHeatFluxTemp1(0), DomainHeatFluxTemp2(0), DomainHeatFluxTemp3(0), kpFTemp1(1), kpFTemp2(1), kpFTemp3(1), kTFTemp1(1), kTFTemp2(1), kTFTemp3(1), kvFTemp1(1), kvFTemp2(1), kvFTemp3(1), iCountTime(0), TemperatureCenter1(1), TemperatureCenter2(1), TemperatureCenter3(1), HeatFluxCenter1(1), HeatFluxCenter2(1), HeatFluxCenter3(1)
    
{
    // fill in the ID containing external node info with node id's
    if (externalNodes.Size() != 2) {
        opserr << "FATAL TripleFrictionPendulumX::TripleFrictionPendulumX() - out of memory, could not create an ID of size 2\n";
        exit(-1);
    }
    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;
    theNodes[0] = 0;
    theNodes[1] = 0;

    // check material input
    if (materials == 0) {
        opserr << "TripleFrictionPendulumX::TripleFrictionPendulumX() - "
            << "null material array passed.\n";
        exit(-1);
    }

    // get copies of the uniaxial materials
    for (int i = 0; i < 4; i++) {
        if (materials[i] == 0) {
            opserr << "TripleFrictionPendulumX::TripleFrictionPendulumX() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0) {
            opserr << "TripleFrictionPendulumX::TripleFrictionPendulumX() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }

    // initialize constants
    v1Fact = 0.5;
    v3Fact = L2 / (L2 - L1);
    v5Fact = L3 / (L3 - L1);

    Gap2 = 2 * Ubar1 + Ubar2 + Ubar3 + B1 / 2 - Ubar2 * (1 - L1 / L2) - Ubar3 * (1 - L1 / L3); 
    Gap4 = Ubar2 * (1 - L1 / L2);
    Gap6 = Ubar3 * (1 - L1 / L3);

    //Unit conversion for pressure to be used in the pressure factor computation
    p_Unit_Convert; //To convert to MPa
    if (unit == 1) { p_Unit_Convert = 0.000001; }
    if (unit == 2) { p_Unit_Convert = 0.001; }
    if (unit == 3) { p_Unit_Convert = 1.0; }
    if (unit == 4) { p_Unit_Convert = 1000.0; }
    if (unit == 5) { p_Unit_Convert = 0.006894; }
    if (unit == 6) { p_Unit_Convert = 6.894; }
    if (unit == 7) { p_Unit_Convert = 0.00004788; }
    if (unit == 8) { p_Unit_Convert = 0.04788; }

    v_Unit_Convert; //To convert velocity 
    if (unit == 1) { v_Unit_Convert = 1.0; }
    if (unit == 2) { v_Unit_Convert = 1.0; }
    if (unit == 3) { v_Unit_Convert = 0.001; }
    if (unit == 4) { v_Unit_Convert = 0.001; }
    if (unit == 5) { v_Unit_Convert = 0.0254; }
    if (unit == 6) { v_Unit_Convert = 0.0254; }
    if (unit == 7) { v_Unit_Convert = 0.3048; }
    if (unit == 8) { v_Unit_Convert = 0.3048; }
    
    // initialize other variables
    this->revertToStart();
}


// constructor which should be invoked by an FE_ObjectBroker only
TripleFrictionPendulumX::TripleFrictionPendulumX()
    : Element(0, ELE_TAG_TripleFrictionPendulumX), externalNodes(2), tag1(0), kpFactor(0), kTFactor(0), kvFactor(0),
    Mu_ref1(0.0), Mu_ref2(0.0), Mu_ref3(0.0), L1(0.0), L2(0.0), L3(0.0), Ubar1(0.0), Ubar2(0.0), Ubar3(0.0), B1(0.0), B2(0.0), B3(0.0),
    W(0.0), Uy(0.0), Kvt(0.0), MinFv(0.0), TOL(1E-6), refPressure1(0.0), refPressure2(0.0), refPressure3(0.0), Diffusivity(0.0), Conductivity(0.0), Temperature0(0.0), rateParam(0.0), tempParam(0.0), unit(0.0), Niter(20),
    K(2, 2), Kpr(2, 2), f(2), fpr(2),
    k12(2, 2), k12pr(2, 2), k34(2, 2), k34pr(2, 2), k56(2, 2), k56pr(2, 2),
    d1(2), d1pr(2), d3(2), d3pr(2), d5(2), d5pr(2), v1(2), v3(2), v5(2),
    ep1(2), ep1pr(2), ep3(2), ep3pr(2), ep5(2), ep5pr(2),
    q1(2), q1pr(2), q3(2), q3pr(2), q5(2), q5pr(2),
    ep1tmp(2), ep3tmp(2), ep5tmp(2), q1tmp(2), q3tmp(2), q5tmp(2)
    , d1ppr(2), d3ppr(2), d5ppr(2)
    , DomainTime(1), DomainHeatFlux1(1), DomainHeatFlux2(1), DomainHeatFlux3(1), DomainTimeTemp(0), DomainHeatFluxTemp1(0), DomainHeatFluxTemp2(0), DomainHeatFluxTemp3(0), kpFTemp1(1), kpFTemp2(1), kpFTemp3(1), kTFTemp1(1), kTFTemp2(1), kTFTemp3(1), kvFTemp1(1), kvFTemp2(1), kvFTemp3(1), iCountTime(0), TemperatureCenter1(1), TemperatureCenter2(1), TemperatureCenter3(1), HeatFluxCenter1(1), HeatFluxCenter2(1), HeatFluxCenter3(1)
{

    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;

    // set material pointers to NULL
    for (int i = 0; i < 4; i++)
        theMaterials[i] = 0;
}


//  destructor - provided to clean up any memory
TripleFrictionPendulumX::~TripleFrictionPendulumX()
{
    // clean up all the material objects
    for (int i = 0; i < 4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int TripleFrictionPendulumX::getNumExternalNodes() const
{
    return 2;
}


const ID& TripleFrictionPendulumX::getExternalNodes()
{
    return externalNodes;
}


Node** TripleFrictionPendulumX::getNodePtrs()
{
    return theNodes;
}


int TripleFrictionPendulumX::getNumDOF()
{
    return 12;
}


// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void TripleFrictionPendulumX::setDomain(Domain* theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        opserr << "Domain does not exist" << endln;
        exit(0);
    }

    // first ensure nodes exist in Domain and set the node pointers
    Node* end1Ptr, * end2Ptr;
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);
    if (end1Ptr == 0) {
        opserr << "WARNING TripleFrictionPendulumX::setDomain() - at TripleFrictionPendulumX " << this->getTag() << " node " <<
            Nd1 << "  does not exist in domain\n";

        return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {
        opserr << "WARNING TripleFrictionPendulumX::setDomain() - at TripleFrictionPendulumX " << this->getTag() << " node " <<
            Nd2 << "  does not exist in domain\n";

        return;  // don't go any further - otherwise segemntation fault
    }
    theNodes[0] = end1Ptr;
    theNodes[1] = end2Ptr;

    // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();

    if ((dofNd1 != 6) || (dofNd2 != 6)) {
        opserr << "TripleFrictionPendulumX::setDomain(): 6 dof required at nodes\n";
        return;
    }
}


int TripleFrictionPendulumX::commitState()
{
    int errCode = 0;

    // commit material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->commitState();

    // commit the base class
    errCode += this->Element::commitState();

    // commit other history variables
    Wpr = Wcr;
    Fy1pr = Fy1; Fy3pr = Fy3; Fy5pr = Fy5;
    Kpr = K;
    fpr = f;
    k12pr = k12; k34pr = k34; k56pr = k56;

    ep1pr = ep1tmp; ep3pr = ep3tmp; ep5pr = ep5tmp;
    q1pr = q1tmp; q3pr = q3tmp; q5pr = q5tmp;
    
    // Variables for displacement and velocity histories (Approach 1)
    fx = eleR(6), fy = eleR(7);
    forceSlope_x_stored = forceSlope_x; // x-direction slope sign  
    forceSlope_y_stored = forceSlope_y; // y-direction slope sign  

    u1_stored2 = u1_stored; u4_stored2 = u4_stored; u2_stored2 = u2_stored; u3_stored2 = u3_stored; // x-direction displacement (for 3point)
    u1y_stored2 = u1y_stored; u4y_stored2 = u4y_stored; u2y_stored2 = u2y_stored; u3y_stored2 = u3y_stored; // y-direction displacement (for 3point)
    
    u1_stored = u1; u4_stored = u4; u2_stored = u2; u3_stored = u3; // x-direction displacement
    u1y_stored = u1y; u4y_stored = u4y; u2y_stored = u2y; u3y_stored = u3y; // y-direction displacement
    
    // Series model inner surface
    u23yy_storedpr = u23yy_stored; u23xx_storedpr = u23xx_stored;
    u23yy_stored = u23yy; u23xx_stored = u23xx;

    // Displacement and Velocity histories
    d1ppr = d1pr; d3ppr = d3pr; d5ppr = d5pr; // Displacement history in last 2 time step (3 point formula)
    d1pr = d1; d3pr = d3; d5pr = d5;
    D1prAvg = d1pr.Norm(); D3prAvg = d3pr.Norm(); D5prAvg = d5pr.Norm(); // Approach 1
    v23sumpr = v23sum; u23sumpr = u23sum; // Approach 1
    disp1pr = (u23t / 2); disp2pr = u1t; disp3pr = u4t; // Approach 2
    vel1pr = (v23t / 2); vel2pr = v1t; vel3pr = v4t; // Approach 2

    Mu_Adj1ppr = Mu_Adj1pr; Mu_Adj2ppr = Mu_Adj2pr; Mu_Adj3ppr = Mu_Adj3pr; // COF (2 steps before)
    Mu_Adj1pr = Mu_Adj1; Mu_Adj2pr = Mu_Adj2; Mu_Adj3pr = Mu_Adj3;
    
    Dxpppr = Dxppr; Dypppr = Dyppr;
    Dxppr = Dxpr; Dyppr = Dypr;
    Dxpr = Dx; Dypr = Dy;

    resultantVpr = sqrt(((Dxppr - Dxpppr) / ops_Dt) * ((Dxppr - Dxpppr) / ops_Dt) + ((Dyppr - Dypppr) / ops_Dt) * ((Dyppr - Dypppr) / ops_Dt));
    resultantV = sqrt(((Dxpr - Dxppr) / ops_Dt) * ((Dxpr - Dxppr) / ops_Dt) + ((Dypr - Dyppr) / ops_Dt) * ((Dypr - Dyppr) / ops_Dt));
    resultantV_AVG = (resultantVpr + resultantV) / 2;
    return 0;
}


int TripleFrictionPendulumX::revertToLastCommit()
{
    int errCode = 0;

    // revert material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->revertToLastCommit();

    return 0;
}


int TripleFrictionPendulumX::revertToStart()
{
    int errCode = 0;
    Vector tmp1(2), tmp2(2), tmp3(2);

    Vel1Avg = Vel3Avg = Vel5Avg = 0.0;
    Fy1pr = Fy3pr = Fy5pr = 0.0;
    Wpr = Wcr = Wavg = W;
        
    // initialize friction coefficients
    Temperature_Surface1 = Temperature_Surface2 = Temperature_Surface3 = Temperature0;
    trialP1 = refPressure1; trialP2 = refPressure2; trialP3 = refPressure3;
    trialVel1 = trialVel2 = trialVel3 = 0.0;

    //temperature factor
    if (fabs(kTFactor - 1.0) <= 0.00001) { 
        kTF1 = 0.789 * (pow(0.7, (Temperature_Surface1 / 50.0)) + 0.4);
        kTF2 = 0.789 * (pow(0.7, (Temperature_Surface2 / 50.0)) + 0.4);
        kTF3 = 0.789 * (pow(0.7, (Temperature_Surface3 / 50.0)) + 0.4);
    }
    else { kTF1 = 1.0; kTF2 = 1.0; kTF3 = 1.0; }
    //Pressure factor
    if (fabs(kpFactor - 1.0) <= 0.00001) {
        kpF1 = pow(0.7, (0.02 * p_Unit_Convert * (trialP1 - refPressure1))); 
        kpF2 = pow(0.7, (0.02 * p_Unit_Convert * (trialP2 - refPressure2))); 
        kpF3 = pow(0.7, (0.02 * p_Unit_Convert * (trialP3 - refPressure3))); 
    }
    else { kpF1 = 1.0; kpF2 = 1.0; kpF3 = 1.0; }
    //Velocity factor
    if (fabs(kvFactor - 1.0) <= 0.00001) { 
        kvF1 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel1);
        kvF2 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel2);
        kvF3 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel3);
    }
    else { kvF1 = 1.0; kvF2 = 1.0; kvF3 = 1.0; }
         
    Fy1 = Mu_ref1 * kpF1 * kTF1 * kvF1;
    Fy3 = Mu_ref2 * kpF2 * kTF2 * kvF2;
    Fy5 = Mu_ref3 * kpF3 * kTF3 * kvF3;

    E1 = E2 = 3.0 * Fy1 / Uy;
    E3 = E4 = 3.0 * Fy1 / Uy; 
    E5 = E6 = 3.0 * Fy1 / Uy; 

    /*E1 = E2 = 3.0 * Fy1 / Uy;
    E3 = E4 = 3.0 * Fy3 / Uy; 
    E5 = E6 = 3.0 * Fy5 / Uy; */

    double E1p = 1.0 / (2.0 * L1);
    double E3p = 1.0 / (L2 - L1);
    double E5p = 1.0 / (L3 - L1);

    H1 = E1 * E1p / (E1 - E1p);
    H3 = E3 * E3p / (E3 - E3p);
    H5 = E5 * E5p / (E5 - E5p);
    
    
    // revert material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->revertToStart();

    // initialize remaining variables
    Fvert = 0.0;
    Kvert = theMaterials[0]->getInitialTangent();
    TorqX = 0.0;
    KrotX = theMaterials[2]->getInitialTangent();
    TorqY = 0.0;
    KrotY = theMaterials[3]->getInitialTangent();
    TorqZ = 0.0;
    KrotZ = theMaterials[1]->getInitialTangent();
    Hisolator = 0.0;
    Dx = Dy = Dz = 0.0;

    // reset history variables
    d1pr.Zero();
    d3pr.Zero();
    d5pr.Zero();
    ep1pr.Zero();
    ep3pr.Zero();
    ep5pr.Zero();
    q1pr.Zero();
    q3pr.Zero();
    q5pr.Zero();
    fpr.Zero();

    d1ppr.Zero();
    d3ppr.Zero();
    d5ppr.Zero();

    BidirectionalPlastic(k12pr, tmp1, tmp2, tmp3, Fy1, E1, H1, ep1pr, q1pr, d1pr);
    BidirectionalPlastic(k34pr, tmp1, tmp2, tmp3, Fy3, E3, H3, ep3pr, q3pr, d3pr);
    BidirectionalPlastic(k56pr, tmp1, tmp2, tmp3, Fy5, E5, H5, ep5pr, q5pr, d5pr);
    StiffnessForm(Kpr, k12pr, k34pr, k56pr);

    return errCode;
}


int TripleFrictionPendulumX::update()
{
    // get current time
    //Current domain time
    double tCurrent = (this->getDomain())->getCurrentTime();
      
    const Vector& duNd1 = theNodes[0]->getIncrDisp();
    const Vector& duNd2 = theNodes[1]->getIncrDisp();
    const Vector& utrialNd1 = theNodes[0]->getTrialDisp();
    const Vector& utrialNd2 = theNodes[1]->getTrialDisp();
    const Vector& vtrialNd1 = theNodes[0]->getTrialVel();
    const Vector& vtrialNd2 = theNodes[1]->getTrialVel();
    const Vector& uNd1 = theNodes[0]->getDisp();
    const Vector& uNd2 = theNodes[1]->getDisp();
    const Vector& end1Crd = theNodes[0]->getCrds();
    const Vector& end2Crd = theNodes[1]->getCrds();

    const Vector& vNd1 = theNodes[0]->getVel();
    const Vector& vNd2 = theNodes[1]->getVel();

    Vector u(2);
    u(0) = uNd2(0) - uNd1(0);  // converged displacement from previous step
    u(1) = uNd2(1) - uNd1(1);
    
    Vector utrial(2);
    Dx = utrial(0) = utrialNd2(0) - utrialNd1(0);  // trial displacement (target displacement)
    Dy = utrial(1) = utrialNd2(1) - utrialNd1(1);
    Dz = utrialNd2(2) - utrialNd1(2);
    
    Vector vtrial(2);
    Vx = vtrial(0) = vtrialNd1(0) - vtrialNd2(0);
    Vy = vtrial(1) = vtrialNd1(1) - vtrialNd2(1);
    Vz = vtrial(2) = vtrialNd1(2) - vtrialNd2(2);
    
    Vector r(3), rdot(3);
    r(0) = utrialNd2(3) - utrialNd1(3);
    r(1) = utrialNd2(4) - utrialNd1(4);
    r(2) = utrialNd2(5) - utrialNd1(5);
    rdot(0) = vtrialNd2(3) - vtrialNd1(3);
    rdot(1) = vtrialNd2(4) - vtrialNd1(4);
    rdot(2) = vtrialNd2(5) - vtrialNd1(5);

    Vector dusub(2);
    dusub(0) = duNd2(0) - duNd1(0);  // incremental displacement
    dusub(1) = duNd2(1) - duNd1(1);

    // isolator height
    Hisolator = end2Crd(2) - end1Crd(2);


    // 1) get axial force and stiffness in vertical direction
    double DzOld = theMaterials[0]->getStrain();
    double Vz = vtrialNd2(2) - vtrialNd1(2);
    theMaterials[0]->setTrialStrain(Dz, Vz);
    Fvert = theMaterials[0]->getStress();
    Kvert = theMaterials[0]->getTangent();

    // vertical force for computing friction
    if (Fvert >= 0.0) {
        Kvert = theMaterials[0]->getInitialTangent();
        if (Fvert > Kvert * DBL_EPSILON) {
            theMaterials[0]->setTrialStrain(DzOld, 0.0);
            Kvert = Kvt;
        }
        Fvert = -MinFv;
    }
    Wcr = -Fvert;


    // 2) calculate shear forces and stiffnesses in horizontal direction
    double Tol = dusub.Norm() * TOL;
    K = Kpr;
    f = fpr;
    k12 = k12pr; k34 = k34pr; k56 = k56pr;
    d1 = d1pr; d3 = d3pr; d5 = d5pr;
    ep1 = ep1pr; ep3 = ep3pr; ep5 = ep5pr;
    q1 = q1pr; q3 = q3pr; q5 = q5pr;
    ep1tmp = ep1pr; ep3tmp = ep3pr; ep5tmp = ep5pr;
    q1tmp = q1pr; q3tmp = q3pr; q5tmp = q5pr;
    Vector ErrDisp = dusub;

    // averaging vertical force
    Wavg = (Wpr + Wcr) / 2.0; 
    
    // get coefficients of friction
    double Fy1cr, Fy3cr, Fy5cr, dFy1, dFy3, dFy5;
    


    // Partitioning Approaches: Displacement and Velocity Histories of each surfaces
    
    // Regime limit
    // Displacement limit
    u_star = 2 * (Fy3 - Fy1) * L1;  // Regime 1 (Displacement)
    u_star2 = u_star + (Fy5 - Fy3)*(L2 + L1);  // Regime 2 (Displacement)
    udr1 = u_star2 + Ubar2 * (1 + L3 / L2) - (Fy5 - Fy3) * (L2 + L3);  // Regime 3 (Displacement)
    udr4 = udr1 + ((Ubar3 / L3 + Fy5) - (Ubar2 / L2 + Fy3)) * (L1 + L3); // Regime 4(TFP) displacement limit state //  Table 4-1 in MCEER 08-0007
    
    // Force limit
    F_f1 = Wavg * Fy3;
    F_f2 = Wavg * Fy1;
    F_f4 = Wavg * Fy5;
    F_dr1 = Wavg * (Ubar2 / L2 + Fy3);
    F_dr4 = Wavg * (Ubar3 / L3 + Fy5);
    

    // Partitioning Approach 2 (X-direction)

    // 1. Save sign for first cycle    
    if (fabs(fx) <= F_f2) { // Before regime 1
        changeSignX = forceSlope_x;
    }
    
    // 2. First Cycle
    if (unloading_x == 0 && loading_x == 0 && forceSlope_x == changeSignX) { 
        sign_fx = sgn(forceSlope_x); // Keep loading direction for next cycles
        // Regime 1 
        if (sign_fx*fx <= F_f1) { 
            u1 = u1_stored; u4 = u4_stored; u2 = (fx - sign_fx * F_f2) / Wavg * L1; u3 = (fx - sign_fx * F_f2) / Wavg * L1;
            if (sign_fx * fx <= F_f2) { u2 = u3 = 0; }
            // Save force and displacement histories for unloading phase 
            F_tr = fx; u1_tr = u1; u2_tr = u2; u3_tr = u3; u4_tr = u4;
        }
        // Regime 2
        else if (sign_fx*fx > F_f1 && sign_fx * fx <= F_f4) {
            u1 = (fx - sign_fx * F_f1) / Wavg * L2; u4 = u4_stored; u2 = u2_stored; u3 = (fx - sign_fx * F_f2) / Wavg * L1;
            // Save histories for unloading phase 
            F_tr = fx; u1_tr = u1; u2_tr = u2; u3_tr = u3; u4_tr = u4;
        }
        // Regime 3
        else if (sign_fx * fx >= F_f4 && sign_fx * fx <= F_dr1) {
            u1 = (fx - sign_fx * F_f1) / Wavg * L2; u4 = (fx - sign_fx * F_f4) / Wavg * L3; u2 = u2_stored; u3 = u3_stored;
            // Save histories for unloading phase 
            F_tr = fx; u1_tr = u1; u2_tr = u2; u3_tr = u3; u4_tr = u4;
        }
        // Regime 4
        else if (sign_fx * fx >= F_dr1 && sign_fx * fx <= F_dr4) {
            u1 = u1_stored; u4 = (fx - sign_fx * F_f4) / Wavg * L3; u2 = ((fx - sign_fx * F_f2) / Wavg - sign_fx * Ubar2 / L2) * L1; u3 = u3_stored;
            // Save histories for unloading phase 
            F_tr = fx; u1_tr = u1; u2_tr = u2; u3_tr = u3; u4_tr = u4;
        }
        // Regime 5
        else if (sign_fx * fx >= F_dr4) {
            u1 = u1_stored; u4 = u4_stored; u2 = ((fx - sign_fx * F_f2) / Wavg - sign_fx * Ubar2 / L2) * L1; u3 = ((fx - sign_fx * F_f2) / Wavg - sign_fx * Ubar3 / L3) * L1; 
            if (fabs(u2) >= Ubar1) { u2 = sign_fx * Ubar1; }
            if (fabs(u3) >= Ubar1) { u3 = sign_fx * Ubar1; }
            // Save histories for unloading phase 
            F_tr = fx; u1_tr = u1; u2_tr = u2; u3_tr = u3; u4_tr = u4;
        }
    }

    // 2. Change loading phase - from bottom to top of the force-displacement loop
    else if (forceSlope_x == 1){ 
        loading_x = 1; // tag for closing first loop 

        // Setting the starting points depending on the direction of first cycle
        if (sign_fx == -1) { 
            F_ref = F_tr; u1_ref = u1_tr; u4_ref = u4_tr; u2_ref = u2_tr; u3_ref = u3_tr;
        }
        else {
            F_ref = F_tr_u; u1_ref = u1_tr_u; u4_ref = u4_tr_u; u2_ref = u2_tr_u; u3_ref = u3_tr_u;
        }
        
        // Starting Regime (5)
        if (F_ref < -F_dr4) { // If F_ref requires loop for regime 5
            // Regime 0 (just unloading)
            if (fx <= F_ref + 2 * F_f2) {
                u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 1 
            else if (fx > F_ref + 2 * F_f2 && fx <= -F_dr1 + 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 2
            else if (fx > -F_dr1 + 2 * F_f1 && fx <= -F_dr4 + 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                u2 = u2_stored; u4 = u4_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 3
            else if (fx > -F_dr4 + 2 * F_f4 && fx <= F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 4
            else if (fx > F_dr1 && fx <= F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                else { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                u3 = u3_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) {F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 5
            else if (fx >= F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
        } // end of starting regime 5
        
        // Starting Regime (4)
       else if (F_ref >= -F_dr4 && F_ref <= -F_dr1) { // If F_ref requires loop for regime 4
            // Regime 0 (just unloading)
            if (fx <= F_ref + 2 * F_f2) { u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 1 
            else if (fx > F_ref + 2 * F_f2 && fx <= -F_dr1 + 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 2
            else if (fx > -F_dr1 + 2 * F_f1 && fx <= F_ref + 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                u2 = u2_stored; u4 = u4_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 3
            else if (fx > F_ref + 2 * F_f4 && fx <= F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 4
            else if (fx > F_dr1 && fx <= F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                else { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                u3 = u3_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 5
            else if (fx >= F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
        } // end of starting regime 4    

        // Starting Regime (3, 2, 1)
       else if (F_ref > -F_dr1) { // If F_ref requires loop for regime 3 or lower
            // Regime 0 (just unloading)
            if (fx <= F_ref + 2 * F_f2) {
                u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 1 
            else if (fx > F_ref + 2 * F_f2 && fx <= F_ref + 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 2
            else if (fx > F_ref + 2 * F_f1 && fx <= F_ref + 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                u2 = u2_stored; u4 = u4_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 3
            else if (fx > F_ref + 2 * F_f4 && fx <= F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                else { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 4
            else if (fx > F_dr1 && fx <= F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                else { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                u3 = u3_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
            // Regime 5
            else if (fx >= F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                else { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
                else { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
            }
        } // end of starting regime 3,2,1   
    }
    
    // 3. Change loading phase - from top to bottom of the force-displacement loop
    else if (forceSlope_x == -1) {
        unloading_x = 1; // tag for closing first loop 
        // Setting the starting points depending on the direction of first cycle
        if (sign_fx == -1) { F_ref = F_tr_u; u1_ref = u1_tr_u; u4_ref = u4_tr_u; u2_ref = u2_tr_u; u3_ref = u3_tr_u; }
        else { F_ref = F_tr; u1_ref = u1_tr; u4_ref = u4_tr; u2_ref = u2_tr; u3_ref = u3_tr; }                          
        // Starting Regime (5)
        if (F_ref > F_dr4) {
            // Regime 0 (just unloading)
            if (fx >= F_ref - 2 * F_f2) {
                u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fx < F_ref - 2 * F_f2 && fx >= F_dr1 - 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 2
            else if (fx < F_dr1 - 2 * F_f1 && fx >= F_dr4 - 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                u4 = u4_stored; u2 = u2_stored;
                // displacement limit states for inner surface
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 3
            else if (fx <= F_dr4 - 2 * F_f4 && fx >= -F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 4
            else if (fx < -F_dr1 && fx >= -F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                else { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                u3 = u3_stored;
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 5
            else if (fx <= -F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
        } // end of starting regime 5

        // Starting Regime (4)
        else if (F_ref >= F_dr1 && F_ref < F_dr4) {
            // Regime 0 (just unloading)
            if (fx >= F_ref - 2 * F_f2) {
                u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fx < F_ref - 2 * F_f2 && fx >= F_dr1 - 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 2
            else if (fx < F_dr1 - 2 * F_f1 && fx >= F_ref - 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                u4 = u4_stored; u2 = u2_stored;
                // displacement limit states for inner surface
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 3
            else if (fx <= F_ref - 2 * F_f4 && fx >= -F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1)*Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4)*Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 4
            else if (fx < -F_dr1 && fx >= -F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                else { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                u3 = u3_stored;
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 5
            else if (fx <= -F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
        } // end of starting regime 4

        // Starting Regime (3,2,1)
       else if (F_ref < F_dr1) {
            // Regime 0 (just unloading)
            if (fx >= F_ref - 2 * F_f2) {
                u1 = u1_stored; u4 = u4_stored; u2 = u2_stored; u3 = u3_stored;
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fx < F_ref - 2 * F_f2 && fx >= F_ref - 2 * F_f1) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 2
            else if (fx < F_ref - 2 * F_f1 && fx >= F_ref - 2 * F_f4) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                u4 = u4_stored; u2 = u2_stored;
                // displacement limit states for inner surface
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 3
            else if (fx <= F_ref - 2 * F_f4 && fx >= -F_dr1) {
                if (sign_fx == -1) { u1 = ((fx - F_tr) / Wavg * L2 + u1_tr); u4 = ((fx - F_tr) / Wavg * L3 + u4_tr); }
                else { u1 = ((fx - F_tr_u) / Wavg * L2 + u1_tr_u); u4 = ((fx - F_tr_u) / Wavg * L3 + u4_tr_u); }
                u2 = u2_stored; u3 = u3_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1) >= Ubar2) { u1 = sgn(u1) * Ubar2; }
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 4
            else if (fx < -F_dr1 && fx >= -F_dr4) {
                u1 = u1_stored;
                if (sign_fx == -1) { u4 = (fx - F_tr) / Wavg * L3 + u4_tr; u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); }
                else { u4 = (fx - F_tr_u) / Wavg * L3 + u4_tr_u; u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); }
                u3 = u3_stored;
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4) >= Ubar3) { u4 = sgn(u4) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
            // Regime 5
            else if (fx <= -F_dr4) {
                u1 = u1_stored; u4 = u4_stored;
                if (sign_fx == -1) { u2 = ((fx - F_tr) / Wavg * L1 + u2_tr); u3 = ((fx - F_tr) / Wavg * L1 + u3_tr); }
                else { u2 = ((fx - F_tr_u) / Wavg * L1 + u2_tr_u); u3 = ((fx - F_tr_u) / Wavg * L1 + u3_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2) >= Ubar1) { u2 = sgn(u2) * Ubar1; }
                if (fabs(u3) >= Ubar1) { u3 = sgn(u3) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fx == -1) { F_tr = fx; u1_tr = u1; u4_tr = u4; u2_tr = u2; u3_tr = u3; }
                else { F_tr_u = fx; u1_tr_u = u1; u4_tr_u = u4; u2_tr_u = u2; u3_tr_u = u3; }
            }
        } // end of starting regime 3,2,1
    }
    
    // Histories of inner surface for series model element 1
    u23 = u2 + u3;
    v2x = (u2 - u2_stored2) * (1.0 / (2*ops_Dt));
    v3x = (u3 - u3_stored2) * (1.0 / (2*ops_Dt));
    v23 = v2x + v3x;

    // Velocity for outer surfaces
    v1x = (u1 - u1_stored2) * (1.0 / (2*ops_Dt));
    v4x = (u4 - u4_stored2) * (1.0 / (2*ops_Dt));

    // Partitioning Approach 2 (Y-direction)
    
    // 1. Save sign for first cycle    
    if (fabs(fy) <= F_f2) { // Before regime 1
        changeSignY = forceSlope_y;
    }

    // 2. First Cycle
    if (unloading_y == 0 && loading_y == 0 && forceSlope_y == changeSignY) {
        sign_fy = sgn(forceSlope_y); // Keep loading direction for next cycles
        // Regime 1 
        if (sign_fy * fy <= F_f1) {
            u1y = u1y_stored; u4y = u4y_stored; u2y = (fy - sign_fy * F_f2) / Wavg * L1; u3y = (fy - sign_fy * F_f2) / Wavg * L1;
            if (sign_fy * fy <= F_f2) { u2y = u3y = 0; }
            // Save force and displacement histories for unloading phase 
            Fy_tr = fy; u1y_tr = u1y; u2y_tr = u2y; u3y_tr = u3y; u4y_tr = u4y;
        }
        // Regime 2
        else if (sign_fy * fy > F_f1 && sign_fy * fy <= F_f4) {
            u1y = (fy - sign_fy * F_f1) / Wavg * L2; u4y = u4y_stored; u2y = u2y_stored; u3y = (fy - sign_fy * F_f2) / Wavg * L1;
            // Save histories for unloading phase 
            Fy_tr = fy; u1y_tr = u1y; u2y_tr = u2y; u3y_tr = u3y; u4y_tr = u4y;
        }
        // Regime 3
        else if (sign_fy * fy >= F_f4 && sign_fy * fy <= F_dr1) {
            u1y = (fy - sign_fy * F_f1) / Wavg * L2; u4y = (fy - sign_fy * F_f4) / Wavg * L3; u2y = u2y_stored; u3y = u3y_stored;
            // Save histories for unloading phase 
            Fy_tr = fy; u1y_tr = u1y; u2y_tr = u2y; u3y_tr = u3y; u4y_tr = u4y;
        }
        // Regime 4
        else if (sign_fy * fy >= F_dr1 && sign_fy * fy <= F_dr4) {
            u1y = u1y_stored; u4y = (fy - sign_fy * F_f4) / Wavg * L3; u2y = ((fy - sign_fy * F_f2) / Wavg - sign_fy * Ubar2 / L2) * L1; u3y = u3y_stored;
            // Save histories for unloading phase 
            Fy_tr = fy; u1y_tr = u1y; u2y_tr = u2y; u3y_tr = u3y; u4y_tr = u4y;
        }
        // Regime 5
        else if (sign_fy * fy >= F_dr4) {
            u1y = u1y_stored; u4y = u4y_stored; u2y = ((fy - sign_fy * F_f2) / Wavg - sign_fy * Ubar2 / L2) * L1; u3y = ((fy - sign_fy * F_f2) / Wavg - sign_fy * Ubar3 / L3) * L1;
            if (fabs(u2y) >= Ubar1) { u2y = sign_fy * Ubar1; }
            if (fabs(u3y) >= Ubar1) { u3y = sign_fy * Ubar1; }
            // Save histories for unloading phase 
            Fy_tr = fy; u1y_tr = u1y; u2y_tr = u2y; u3y_tr = u3y; u4y_tr = u4y;
        }
    }

    // 2. Change loading phase - from bottom to top of the force-displacement loop
    else if (forceSlope_y == 1) {
        loading_y = 1; // tag for closing first loop 

        // Setting the starting points depending on the direction of first cycle
        if (sign_fy == -1) { Fy_ref = Fy_tr; u1y_ref = u1y_tr; u4y_ref = u4y_tr; u2y_ref = u2y_tr; u3y_ref = u3y_tr; }
        else { Fy_ref = Fy_tr_u; u1y_ref = u1y_tr_u; u4y_ref = u4y_tr_u; u2y_ref = u2y_tr_u; u3y_ref = u3y_tr_u; }

        // Starting Regime (5)
        if (Fy_ref < -F_dr4) { // If F_ref requires loop for regime 5
            // Regime 0 (just unloading)
            if (fy <= Fy_ref + 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 1 
            else if (fy > Fy_ref + 2 * F_f2 && fy <= -F_dr1 + 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 2
            else if (fy > -F_dr1 + 2 * F_f1 && fy <= -F_dr4 + 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                u2y = u2y_stored; u4y = u4y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y;  }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 3
            else if (fy > -F_dr4 + 2 * F_f4 && fy <= F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 4
            else if (fy > F_dr1 && fy <= F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u; u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                else { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 5
            else if (fy >= F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
        } // end of starting regime 5

        // Starting Regime (4)
        else if (Fy_ref >= -F_dr4 && Fy_ref <= -F_dr1) { // If F_ref requires loop for regime 4
             // Regime 0 (just unloading)
            if (fy <= Fy_ref + 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 1 
            else if (fy> Fy_ref + 2 * F_f2 && fy <= -F_dr1 + 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 2
            else if (fy > -F_dr1 + 2 * F_f1 && fy <= Fy_ref + 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                u2y = u2y_stored; u4y = u4y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 3
            else if (fy > Fy_ref + 2 * F_f4 && fy <= F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 4
            else if (fy > F_dr1 && fy <= F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u; u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                else { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 5
            else if (fy >= F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
        } // end of starting regime 4    

        // Starting Regime (3, 2, 1)
        else if (Fy_ref > -F_dr1) { // If F_ref requires loop for regime 3 or lower
             // Regime 0 (just unloading)
            if (fy <= Fy_ref + 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 1 
            else if (fy > Fy_ref + 2 * F_f2 && fy <= Fy_ref + 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 2
            else if (fy > Fy_ref + 2 * F_f1 && fy <= Fy_ref + 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                u2y = u2y_stored; u4y = u4y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 3
            else if (fy > Fy_ref + 2 * F_f4 && fy <= F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                else { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 4
            else if (fy > F_dr1 && fy <= F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u; u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                else { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
            // Regime 5
            else if (fy >= F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                else { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr);  u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
                else { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
            }
        } // end of starting regime 3,2,1   
    }

    // 3. Change loading phase - from top to bottom of the force-displacement loop
    else if (forceSlope_y == -1) {
        unloading_y = 1; // tag for closing first loop 

        // Setting the starting points depending on the direction of first cycle
        if (sign_fy == -1) { Fy_ref = Fy_tr_u; u1y_ref = u1y_tr_u; u4y_ref = u4y_tr_u; u2y_ref = u2y_tr_u; u3y_ref = u3y_tr_u; }
        else { Fy_ref = Fy_tr; u1y_ref = u1y_tr; u4y_ref = u4y_tr; u2y_ref = u2y_tr; u3y_ref = u3y_tr; }

        // Starting Regime (5)
        if (Fy_ref > F_dr4) {
            // Regime 0 (just unloading)
            if (fy >= Fy_ref - 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fy < Fy_ref - 2 * F_f2 && fy >= F_dr1 - 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 2
            else if (fy < F_dr1 - 2 * F_f1 && fy >= F_dr4 - 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                u4y = u4y_stored; u2y = u2y_stored;
                // displacement limit states for inner surface
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 3
            else if (fy <= F_dr4 - 2 * F_f4 && fy >= -F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 4
            else if (fy < -F_dr1 && fy >= -F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                else { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u; u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y;}
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 5
            else if (fy <= -F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
        } // end of starting regime 5

        // Starting Regime (4)
        else if (Fy_ref >= F_dr1 && Fy_ref < F_dr4) {
            // Regime 0 (just unloading) 
            if (fy >= Fy_ref - 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fy < Fy_ref - 2 * F_f2 && fy >= F_dr1 - 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 2
            else if (fy < F_dr1 - 2 * F_f1 && fy >= Fy_ref - 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                u4y = u4y_stored; u2y = u2y_stored;
                // displacement limit states for inner surface
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 3
            else if (fy <= Fy_ref - 2 * F_f4 && fy >= -F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 4
            else if (fy < -F_dr1 && fy >= -F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                else { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u; u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 5
            else if (fy <= -F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
        } // end of starting regime 4

        // Starting Regime (3,2,1)
        else if (Fy_ref < F_dr1) {
            // Regime 0 (just unloading)
            if (fy >= Fy_ref - 2 * F_f2) {
                u1y = u1y_stored; u4y = u4y_stored; u2y = u2y_stored; u3y = u3y_stored;
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 1 // Sliding on inner only (Regime 1 in TFP, force limit)
            else if (fy < Fy_ref - 2 * F_f2 && fy >= Fy_ref - 2 * F_f1) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // displacement limit states for inner surface
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 2
            else if (fy < Fy_ref - 2 * F_f1 && fy >= Fy_ref - 2 * F_f4) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                u4y = u4y_stored;
                u2y = u2y_stored;
                // displacement limit states for inner surface
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 3
            else if (fy <= Fy_ref - 2 * F_f4 && fy >= -F_dr1) {
                if (sign_fy == -1) { u1y = ((fy - Fy_tr) / Wavg * L2 + u1y_tr); u4y = ((fy - Fy_tr) / Wavg * L3 + u4y_tr); }
                else { u1y = ((fy - Fy_tr_u) / Wavg * L2 + u1y_tr_u); u4y = ((fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u); }
                u2y = u2y_stored; u3y = u3y_stored;
                // Limit states for outer surfaces displacement
                if (fabs(u1y) >= Ubar2) { u1y = sgn(u1y) * Ubar2; }
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 4
            else if (fy < -F_dr1 && fy >= -F_dr4) {
                u1y = u1y_stored;
                if (sign_fy == -1) { u4y = (fy - Fy_tr) / Wavg * L3 + u4y_tr; u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); }
                else { u4y = (fy - Fy_tr_u) / Wavg * L3 + u4y_tr_u;  u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); }
                u3y = u3y_stored;
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                // Limit states for outer surfaces displacement
                if (fabs(u4y) >= Ubar3) { u4y = sgn(u4y) * Ubar3; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
            // Regime 5
            else if (fy <= -F_dr4) {
                u1y = u1y_stored; u4y = u4y_stored;
                if (sign_fy == -1) { u2y = ((fy - Fy_tr) / Wavg * L1 + u2y_tr); u3y = ((fy - Fy_tr) / Wavg * L1 + u3y_tr); }
                else { u2y = ((fy - Fy_tr_u) / Wavg * L1 + u2y_tr_u); u3y = ((fy - Fy_tr_u) / Wavg * L1 + u3y_tr_u); }
                // Limit states for inner surfaces displacement
                if (fabs(u2y) >= Ubar1) { u2y = sgn(u2y) * Ubar1; }
                if (fabs(u3y) >= Ubar1) { u3y = sgn(u3y) * Ubar1; }
                // Save histories for unloading phase 
                if (sign_fy == -1) { Fy_tr = fy; u1y_tr = u1y; u4y_tr = u4y; u2y_tr = u2y; u3y_tr = u3y; }
                else { Fy_tr_u = fy; u1y_tr_u = u1y; u4y_tr_u = u4y; u2y_tr_u = u2y; u3y_tr_u = u3y; }
            }
        } // end of starting regime 3,2,1
    }

    // Histories of inner surface for series model element 1
    u23y = u2y + u3y;
    v2y = (u2y - u2y_stored2) * (1.0 / (2*ops_Dt));
    v3y = (u3y - u3y_stored2) * (1.0 / (2*ops_Dt));
    v23y = v2y + v3y;

    // Velocity for outer surfaces
    v1y = (u1y - u1y_stored2) * (1.0 / (2*ops_Dt));
    v4y = (u4y - u4y_stored2) * (1.0 / (2*ops_Dt));

    // Resultant Histories
    u23t = sqrt(u23 * u23 + u23y * u23y);
    u1t = sqrt(u1 * u1 + u1y * u1y);
    u4t = sqrt(u4 * u4 + u4y * u4y);

    v23t = sqrt(v23 * v23 + v23y * v23y);
    v1t = sqrt(v1x * v1x + v1y * v1y);
    v4t = sqrt(v4x * v4x + v4y * v4y);

    

    // Partitioning Approach 1
    // Total displacement in x and y
    uxx = d1(0) + d3(0) + d5(0);
    uyy = d1(1) + d3(1) + d5(1);

    // Displacement history of Inner slider
    u23xx = uxx - (v3Fact * d3(0) + v5Fact * d5(0));
    u23yy = uyy - (v3Fact * d3(1) + v5Fact * d5(1));
    u23sum = sqrt(u23xx * u23xx + u23yy * u23yy);
        
    // 3point velocity
    v23xx = (u23xx - u23xx_storedpr) * (1.0 / (2*ops_Dt));
    v23yy = (u23yy - u23yy_storedpr) * (1.0 / (2*ops_Dt));
    v23sum = sqrt(v23xx * v23xx + v23yy * v23yy);

    
    // Partitioning Approaches (tag1 == 1 for Approach 1, 0 for Approach 2)
    if (tag1 == 1) {
        
        disp1 = v1Fact * u23sumpr; vel1 = v1Fact * v23sumpr; //
        //disp1 = v1Fact * D1prAvg; vel1 = v1Fact * Vel1Avg;
        disp2 = v3Fact * D3prAvg; vel2 = v3Fact * Vel3Avg; 
        disp3 = v5Fact * D5prAvg; vel3 = v5Fact * Vel5Avg;
        
    } else {
               
        disp1 = disp1pr; disp2 = disp2pr; disp3 = disp3pr;
        vel1 = vel1pr; vel2 = vel2pr; vel3 = vel3pr;
    }

    // Calculation for temperature 
    PI = 3.14159;
    trialT = tCurrent;
    rContact1 = B1 / 2; rContact2 = B2 / 2; rContact3 = B3 / 2;
    trialP1 = Wavg / (PI * pow(rContact1, 2)); trialP2 = Wavg / (PI * pow(rContact2, 2)); trialP3 = Wavg / (PI * pow(rContact3, 2));// PRESSURE on the sliding surface
    trialDisp1 = disp1; trialDisp2 = disp2; trialDisp3 = disp3;
    trialVel1 = vel1; trialVel2 = vel2; trialVel3 = vel3;


    if (resultantV_AVG == 0) { trialVel1 = 0; trialVel2 = 0; trialVel3 = 0; } 
   
    //Factors to account for pressure, temperature and velocity dependences of friction
    if (iCountTime == 0) {
        kpF1 = 1.0;	kTF1 = 1.0;	kvF1 = 1.0;
        kpF2 = 1.0;	kTF2 = 1.0;	kvF2 = 1.0;
        kpF3 = 1.0;	kTF3 = 1.0;	kvF3 = 1.0;

        Mu_Adj1 = Mu_ref1 * kpF1 * kTF1 * kvF1;
        Mu_Adj2 = Mu_ref2 * kpF2 * kTF2 * kvF2;
        Mu_Adj3 = Mu_ref3 * kpF3 * kTF3 * kvF3;
    }
    else {
        kpF1 = kpFTemp1(0); kTF1 = kTFTemp1(0); kvF1 = kvFTemp1(0);
        kpF2 = kpFTemp2(0); kTF2 = kTFTemp2(0); kvF2 = kvFTemp2(0);
        kpF3 = kpFTemp3(0); kTF3 = kTFTemp3(0); kvF3 = kvFTemp3(0);

        Mu_Adj1 = Mu_ref1 * kpF1 * kTF1 * kvF1;
        Mu_Adj2 = Mu_ref2 * kpF2 * kTF2 * kvF2;
        Mu_Adj3 = Mu_ref3 * kpF3 * kTF3 * kvF3;
            
    }

    //Compute factors for pressure, temperature and velocity dependencies	
    if (trialT > DomainTime(iCountTime - 1)) {
        //Initiate the DomainTimeVector
        if (iCountTime == 0) {
            DomainTime.resize(2000);
            for (int iTemp = 1; iTemp <= 2000; iTemp++) {
                DomainTime(iTemp - 1) = 0.0;
            }
        }
        DomainTime(iCountTime) = trialT;
        //Resize and repopulate the DomainTime Vector
        for (int iResize = 1; iResize <= 1000; iResize++) { //Vectors can be as long as 2000*1000
            if (iCountTime == (iResize * 2000 - 3)) {
                DomainTimeTemp.resize(iResize * 2000 - 2);
                //Populate DomainTimeTemp 
                for (int iTemp = 1; iTemp <= iResize * 2000 - 2; iTemp++) {
                    DomainTimeTemp(iTemp - 1) = DomainTime(iTemp - 1);
                }
                DomainTime.resize((iResize + 1) * 2000);
                for (int iTemp = 1; iTemp <= (iResize + 1) * 2000; iTemp++) {
                    DomainTime(iTemp - 1) = 0.0;
                }
                //Populate DomainTime 
                for (int iTemp = 1; iTemp <= iResize * 2000 - 2; iTemp++) {
                    DomainTime(iTemp - 1) = DomainTimeTemp(iTemp - 1);
                }
                DomainTimeTemp.resize(0);
                break;
            }
        }
        //Initialize DomainHeatFlux
        if (iCountTime == 0) {
            DomainHeatFlux1.resize(2000);
            DomainHeatFlux2.resize(2000);
            DomainHeatFlux3.resize(2000);
        
            for (int iTemp = 1; iTemp <= 2000; iTemp++) {
                DomainHeatFlux1(iTemp - 1) = 0.0;
                DomainHeatFlux2(iTemp - 1) = 0.0;
                DomainHeatFlux3(iTemp - 1) = 0.0;
            }
        }
        //Resize and repopulate the DomainHeatFlux Vector
        for (int iResize = 1; iResize <= 1000; iResize++) { //Vectors can be as long as 2000*1000
            if (iCountTime == (iResize * 2000 - 3)) {
                DomainHeatFluxTemp1.resize(iResize * 2000 - 2);
                DomainHeatFluxTemp2.resize(iResize * 2000 - 2);
                DomainHeatFluxTemp3.resize(iResize * 2000 - 2);      

                //Populate DomainHeatFluxTemp
                for (int iTemp = 1; iTemp <= iResize * 2000 - 2; iTemp++) {
                    DomainHeatFluxTemp1(iTemp - 1) = DomainHeatFlux1(iTemp - 1);
                    DomainHeatFluxTemp2(iTemp - 1) = DomainHeatFlux2(iTemp - 1);
                    DomainHeatFluxTemp3(iTemp - 1) = DomainHeatFlux3(iTemp - 1);             
                }
                DomainHeatFlux1.resize((iResize + 1) * 2000);
                DomainHeatFlux2.resize((iResize + 1) * 2000);
                DomainHeatFlux3.resize((iResize + 1) * 2000);
       
                for (int iTemp = 1; iTemp <= (iResize + 1) * 2000; iTemp++) {
                    DomainHeatFlux1(iTemp - 1) = 0.0;
                    DomainHeatFlux2(iTemp - 1) = 0.0;
                    DomainHeatFlux3(iTemp - 1) = 0.0;
                }
                //Populate DomainHeatFlux
                for (int iTemp = 1; iTemp <= iResize * 2000 - 2; iTemp++) {
                    DomainHeatFlux1(iTemp - 1) = DomainHeatFluxTemp1(iTemp - 1);
                    DomainHeatFlux2(iTemp - 1) = DomainHeatFluxTemp2(iTemp - 1);
                    DomainHeatFlux3(iTemp - 1) = DomainHeatFluxTemp3(iTemp - 1);
                }
                DomainHeatFluxTemp1.resize(0);
                DomainHeatFluxTemp2.resize(0);
                DomainHeatFluxTemp3.resize(0);
                break;
            }
        }
        // Heat flux at the sliding surface //
        // Series element 1
        if (trialDisp1 < (rContact1 * sqrt(22.0 / (7.0 * 4.0)))) { 
            DomainHeatFlux1(iCountTime) = Mu_Adj1ppr * trialP1 * trialVel1; // Mu_Adj * trialP * trialVel;
        }
        else DomainHeatFlux1(iCountTime) = 0.0;
        HeatFluxCenter1(0) = DomainHeatFlux1(iCountTime);

        // Series element 2
        if (trialDisp2 < (rContact2 * sqrt(22.0 / (7.0 * 4.0)))) { 
            DomainHeatFlux2(iCountTime) = Mu_Adj2ppr * trialP2 * trialVel2; 
        }
        else DomainHeatFlux2(iCountTime) = 0.0;
        HeatFluxCenter2(0) = DomainHeatFlux2(iCountTime);

        // Series element 3
        if (trialDisp3 < (rContact3 * sqrt(22.0 / (7.0 * 4.0)))) { 
            DomainHeatFlux3(iCountTime) = Mu_Adj3ppr * trialP3 * trialVel3; 
        }
        else DomainHeatFlux3(iCountTime) = 0.0;
        HeatFluxCenter3(0) = DomainHeatFlux3(iCountTime);

        // Temperature at the sliding surface //	
        Temperature_Change1 = 0.0;
        Temperature_Change2 = 0.0;
        Temperature_Change3 = 0.0;
 

        // Series element 1
        if (iCountTime > 1 && trialDisp1 > 0.0) { 
            DtAnalysis = 0.000001; 
            tau = 0.000001;
            for (int jTemp = 0; jTemp <= iCountTime; jTemp++) {
                if (iCountTime == 0) {
                    DtAnalysis = DomainTime(iCountTime);
                    tau = DtAnalysis;
                }
                else {
                    DtAnalysis = DomainTime(iCountTime) - DomainTime(iCountTime - 1);
                    tau = DomainTime(jTemp);
                }
                Temperature_Change1 = Temperature_Change1 + sqrt(Diffusivity * 7.0 / 22.0) * DomainHeatFlux1(iCountTime - jTemp) * DtAnalysis / (Conductivity * sqrt(tau));
                
            }
        }

        Temperature_Surface1 = Temperature0 + Temperature_Change1;
        TemperatureCenter1(0) = Temperature_Surface1;

        // Series element 2
        if (iCountTime > 1 && trialDisp2 > 0.0) { 
            DtAnalysis = 0.000001; 
            tau = 0.000001;
            for (int jTemp = 0; jTemp <= iCountTime; jTemp++) {
                if (iCountTime == 0) {
                    DtAnalysis = DomainTime(iCountTime);
                    tau = DtAnalysis;
                }
                else {
                    DtAnalysis = DomainTime(iCountTime) - DomainTime(iCountTime - 1);
                    tau = DomainTime(jTemp);
                }
                Temperature_Change2 = Temperature_Change2 + sqrt(Diffusivity * 7.0 / 22.0) * DomainHeatFlux2(iCountTime - jTemp) * DtAnalysis / (Conductivity * sqrt(tau));
               
            }
        }

        Temperature_Surface2 = Temperature0 + Temperature_Change2;
        TemperatureCenter2(0) = Temperature_Surface2;
   
        // Series element 3
        if (iCountTime > 1 && trialDisp3 > 0.0) { 
            DtAnalysis = 0.000001; 
            tau = 0.000001;
            for (int jTemp = 0; jTemp <= iCountTime; jTemp++) {
                if (iCountTime == 0) {
                    DtAnalysis = DomainTime(iCountTime);
                    tau = DtAnalysis;
                }
                else {
                    DtAnalysis = DomainTime(iCountTime) - DomainTime(iCountTime - 1);
                    tau = DomainTime(jTemp);
                }
                Temperature_Change3 = Temperature_Change3 + sqrt(Diffusivity * 7.0 / 22.0) * DomainHeatFlux3(iCountTime - jTemp) * DtAnalysis / (Conductivity * sqrt(tau));
            }
        }

        Temperature_Surface3 = Temperature0 + Temperature_Change3;
        TemperatureCenter3(0) = Temperature_Surface3;

        //Temperature factor
        if (fabs(kTFactor - 1.0) <= 0.00001) { 

            if (tempParam == 1) { // kTF = 1/2 when T = 200 Celsius
                kTF1 = 0.789 * (pow(0.7, (Temperature_Surface1 / 50.0)) + 0.4);
                kTF2 = 0.789 * (pow(0.7, (Temperature_Surface2 / 50.0)) + 0.4);
                kTF3 = 0.789 * (pow(0.7, (Temperature_Surface3 / 50.0)) + 0.4);
            }
            if (tempParam == 2) { // kTF = 1/3 when T = 200 Celsius
                kTF1 = 0.9696 * pow(0.7, (Temperature_Surface1 * 0.0289)) + 0.2099;
                kTF2 = 0.9696 * pow(0.7, (Temperature_Surface2 * 0.0289)) + 0.2099;
                kTF3 = 0.9696 * pow(0.7, (Temperature_Surface3 * 0.0289)) + 0.2099;
            }
            if (tempParam == 3) { // kTF = 2/3 when T = 200 Celsius
                kTF1 = 0.8383 * pow(0.7, (Temperature_Surface1 * 0.00851)) + 0.2099;
                kTF2 = 0.8383 * pow(0.7, (Temperature_Surface2 * 0.00851)) + 0.2099;
                kTF3 = 0.8383 * pow(0.7, (Temperature_Surface3 * 0.00851)) + 0.2099;
            }


        }
        else { kTF1 = 1.0; kTF2 = 1.0; kTF3 = 1.0; }
        //Pressure factor
        if (fabs(kpFactor - 1.0) <= 0.00001) {
            kpF1 = pow(0.7, (0.02 * p_Unit_Convert * (trialP1 - refPressure1))); 
            kpF2 = pow(0.7, (0.02 * p_Unit_Convert * (trialP2 - refPressure2))); 
            kpF3 = pow(0.7, (0.02 * p_Unit_Convert * (trialP3 - refPressure3))); 
        }
        else { kpF1 = 1.0; kpF2 = 1.0; kpF3 = 1.0; }
        //Velocity factor
        if (fabs(kvFactor - 1.0) <= 0.00001) { 
            kvF1 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel1);
            kvF2 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel2);
            kvF3 = 1 - 0.50 * exp(-rateParam * v_Unit_Convert * trialVel3);
        }
        else { kvF1 = 1.0; kvF2 = 1.0; kvF3 = 1.0; }
        //Coefficients for all the stages of an iteration during a time-step
        kpFTemp1(0) = kpF1; kTFTemp1(0) = kTF1; kvFTemp1(0) = kvF1;
        kpFTemp2(0) = kpF2; kTFTemp2(0) = kTF2; kvFTemp2(0) = kvF2;
        kpFTemp3(0) = kpF3; kTFTemp3(0) = kTF3; kvFTemp3(0) = kvF3;
            
        //Update the counter
        iCountTime = iCountTime + 1;
            
    }


    //When domain time is reset
    if (trialT < DomainTime(iCountTime - 1)) {
        iCountTime = 0;
        DomainTime(iCountTime) = trialT;
        iCountTime = iCountTime + 1;
        kpF1 = 1.0; kTF1 = 1.0; kvF1 = 1.0;
        kpF2 = 1.0; kTF2 = 1.0; kvF2 = 1.0;
        kpF3 = 1.0; kTF3 = 1.0; kvF3 = 1.0;
        Mu_Adj1 = Mu_ref1 * kpF1 * kTF1 * kvF1;
        Mu_Adj2 = Mu_ref2 * kpF2 * kTF2 * kvF2;
        Mu_Adj3 = Mu_ref3 * kpF3 * kTF3 * kvF3;
        kpFTemp1(0) = kpF1; kTFTemp1(0) = kTF1; kvFTemp1(0) = kvF1;
        kpFTemp2(0) = kpF2; kTFTemp2(0) = kTF2; kvFTemp2(0) = kvF2;
        kpFTemp3(0) = kpF3; kTFTemp3(0) = kTF3; kvFTemp3(0) = kvF3;
    }

    //Factors used during iterations to reach the next time-step
    if (fabs(trialT - DomainTime(iCountTime - 1)) < 0.00000001) {
        kpF1 = kpFTemp1(0); kTF1 = kTFTemp1(0); kvF1 = kvFTemp1(0);
        kpF2 = kpFTemp2(0); kTF2 = kTFTemp2(0); kvF2 = kvFTemp2(0);
        kpF3 = kpFTemp3(0); kTF3 = kTFTemp3(0); kvF3 = kvFTemp3(0);
    }
    
    Mu_Adj1 = Mu_ref1 * kpF1 * kTF1 * kvF1;
    Mu_Adj2 = Mu_ref2 * kpF2 * kTF2 * kvF2;
    Mu_Adj3 = Mu_ref3 * kpF3 * kTF3 * kvF3;

    

    Fy1cr = Mu_Adj1; Fy3cr = Mu_Adj2; Fy5cr = Mu_Adj3;
    dFy1 = Fy1cr - Fy1pr; dFy3 = Fy3cr - Fy3pr; dFy5 = Fy5cr - Fy5pr;
    Fy1 = Fy1pr; Fy3 = Fy3pr; Fy5 = Fy5pr;

    int nDiv = 0; int nWhileIter = 0;
    double TolOriginal = Tol;

    while ((nDiv < 10) && (ErrDisp.Norm() > TolOriginal)) {
        Fy1 = Fy1 + dFy1; Fy3 = Fy3 + dFy3; Fy5 = Fy5 + dFy5;
        TFPElement(Conv, ep1tmp, ep3tmp, ep5tmp, q1tmp, q3tmp, q5tmp, K, f, k12, k34, k56, d1, d3, d5, ep1, ep3, ep5, q1, q3, q5, u, dusub, Fy1, Fy3, Fy5, E1, E3, E5, H1, H3, H5, E2, E4, E6, Gap2, Gap4, Gap6, Tol, Niter);
        
        if ((!Conv) && (nDiv < 7)) {
            u(0) = uNd2(0) - uNd1(0);
            u(1) = uNd2(1) - uNd1(1);
            dFy1 /= 2.0; dFy3 /= 2.0; dFy5 /= 2.0;
            Fy1 = Fy1pr; Fy3 = Fy3pr; Fy5 = Fy5pr;
            K = Kpr;
            f = fpr;
            k12 = k12pr; k34 = k34pr; k56 = k56pr;
            d1 = d1pr; d3 = d3pr; d5 = d5pr;
            ep1 = ep1pr; ep3 = ep3pr; ep5 = ep5pr;
            q1 = q1pr; q3 = q3pr; q5 = q5pr;
            dusub /= 2.0;
            nDiv++;
            nWhileIter = 0;
        }

        else {
            if (nWhileIter >= pow(2.0, nDiv)) {
                break;
            }
            ep1 = ep1tmp; ep3 = ep3tmp; ep5 = ep5tmp;
            q1 = q1tmp; q3 = q3tmp; q5 = q5tmp;
            u += dusub;
            ErrDisp(0) = utrial(0) - u(0);
            ErrDisp(1) = utrial(1) - u(1);
            nWhileIter++;
        }
                
        // Three-point formula for velocity calculation 
        v1 = 1.0 / (2 * ops_Dt) * (d1 - d1ppr);
        v3 = 1.0 / (2 * ops_Dt) * (d3 - d3ppr);  
        v5 = 1.0 / (2 * ops_Dt) * (d5 - d5ppr);
                
        Vel1Avg = v1.Norm(); Vel3Avg = v3.Norm(); Vel5Avg = v5.Norm(); 
        Disp1Avg = d1.Norm(); Disp3Avg = d3.Norm(); Disp5Avg = d5.Norm();

        forceSlope_x = sgn((eleR(6) - fx)); forceSlope_y = sgn((eleR(7) - fy)); // Direciton of force-displacement slope for capturing histories
                

    } // end of while


    if (nDiv == 10) {
        if ((!Conv) || (Tol < ErrDisp.Norm())) {
            opserr << "Warning: isolator " << this->getTag() << " has not converged, ErrDisp = " << ErrDisp << endln;
        }
    }

    // 3) get moment and stiffness about vertical direction
    theMaterials[1]->setTrialStrain(r(2), rdot(2));
    TorqZ = theMaterials[1]->getStress();
    KrotZ = theMaterials[1]->getTangent();

    // 4) get moment and stiffness about horizontal 1 direction
    theMaterials[2]->setTrialStrain(r(0), rdot(0));
    TorqX = theMaterials[2]->getStress();
    KrotX = theMaterials[2]->getTangent();

    // 5) get moment and stiffness about horizontal 2 direction
    theMaterials[3]->setTrialStrain(r(1), rdot(1));
    TorqY = theMaterials[3]->getStress();
    KrotY = theMaterials[3]->getTangent();


    return 0;

}    


const Matrix& TripleFrictionPendulumX::getTangentStiff()
{
    Matrix a(2, 12);
    Matrix aT(12, 2);

    a.Zero();
    aT.Zero();
    a(0, 0) = a(1, 1) = -1;
    a(0, 6) = a(1, 7) = 1;
    aT(0, 0) = aT(1, 1) = -1;
    aT(6, 0) = aT(7, 1) = 1;
    eleK = aT * K * a;
    eleK *= -Fvert; 

    eleK(2, 2) = eleK(8, 8) = Kvert;
    eleK(2, 8) = eleK(8, 2) = -Kvert;
    eleK(3, 3) = eleK(9, 9) = KrotX;
    eleK(3, 9) = eleK(9, 3) = -KrotX;
    eleK(4, 4) = eleK(10, 10) = KrotY;
    eleK(4, 10) = eleK(10, 4) = -KrotY;
    eleK(5, 5) = eleK(11, 11) = KrotZ;
    eleK(5, 11) = eleK(11, 5) = -KrotZ;

    return eleK;
}


const Matrix& TripleFrictionPendulumX::getInitialStiff()
{
    Matrix a(2, 12);
    Matrix aT(12, 2);
    Matrix Kinit(2, 2);

    Kinit.Zero();
    Kinit(0, 0) = Kinit(1, 1) = E1 / 3.0;
    a.Zero();
    aT.Zero();
    a(0, 0) = a(1, 1) = -1;
    a(0, 6) = a(1, 7) = 1;
    aT(0, 0) = aT(1, 1) = -1;
    aT(6, 0) = aT(7, 1) = 1;
    eleKinit = aT * Kinit * a;
    eleKinit *= W;

    eleKinit(2, 2) = eleKinit(8, 8) = theMaterials[0]->getInitialTangent();
    eleKinit(2, 8) = eleKinit(8, 2) = -eleKinit(2, 2);
    eleKinit(3, 3) = eleKinit(9, 9) = theMaterials[2]->getInitialTangent();
    eleKinit(3, 9) = eleKinit(9, 3) = -eleKinit(3, 3);
    eleKinit(4, 4) = eleKinit(10, 10) = theMaterials[3]->getInitialTangent();
    eleKinit(4, 10) = eleKinit(10, 4) = -eleKinit(4, 4);
    eleKinit(5, 5) = eleKinit(11, 11) = theMaterials[1]->getInitialTangent();
    eleKinit(5, 11) = eleKinit(11, 5) = -eleKinit(5, 5);

    return eleKinit;
}


const Matrix& TripleFrictionPendulumX::getDamp()
{
    // zero the matrix
    eleD.Zero();

    // add damping tangent from materials
    eleD(2, 2) = eleD(8, 8) = theMaterials[0]->getDampTangent();
    eleD(2, 8) = eleD(8, 2) = -eleD(2, 2);
    eleD(3, 3) = eleD(9, 9) = theMaterials[2]->getDampTangent();
    eleD(3, 9) = eleD(9, 3) = -eleD(3, 3);
    eleD(4, 4) = eleD(10, 10) = theMaterials[3]->getDampTangent();
    eleD(4, 10) = eleD(10, 4) = -eleD(4, 4);
    eleD(5, 5) = eleD(11, 11) = theMaterials[1]->getDampTangent();
    eleD(5, 11) = eleD(11, 5) = -eleD(5, 5);

    return eleD;
}


const Matrix& TripleFrictionPendulumX::getMass()
{
    eleM.Zero();

    return eleM;
}


const Vector& TripleFrictionPendulumX::getResistingForce()
{
    Matrix aT(12, 2);
    aT.Zero();
    aT(0, 0) = aT(1, 1) = -1;
    aT(6, 0) = aT(7, 1) = 1;
    eleR = aT * f;  
    eleR *= -Fvert;
    double Mx, My, Mz;
    Mx = -Fvert * Dy + eleR(7) * Hisolator;
    My = Fvert * Dx - eleR(6) * Hisolator;
    Mz = eleR(6) * Dy - eleR(7) * Dx;
    eleR(3) = eleR(9) = TorqX + Mx / 2;
    eleR(4) = eleR(10) = TorqY + My / 2;
    eleR(5) = eleR(11) = TorqZ + Mz / 2;
    eleR(2) = -Fvert;
    eleR(8) = Fvert;
    
    return eleR;
}


Element* TripleFrictionPendulumX::getCopy()
{
    TripleFrictionPendulumX* theCopy = new TripleFrictionPendulumX(this->getTag(),
        externalNodes(0), externalNodes(1), tag1, theMaterials, kpFactor, kTFactor, kvFactor, Mu_ref1, Mu_ref2, Mu_ref3, L1, L2, L3,
        Ubar1, Ubar2, Ubar3, B1, B2, B3, W, Uy, Kvt, MinFv, TOL, refPressure1, refPressure2, refPressure3, Diffusivity, Conductivity, Temperature0, rateParam, tempParam, unit); //, rContact, kTFactor, diffuse, conduct
    theCopy->Kpr = Kpr;
    theCopy->fpr = fpr;
    theCopy->k12pr = k12pr;
    theCopy->k34pr = k34pr;
    theCopy->k56pr = k56pr;
    theCopy->d1pr = d1pr;
    theCopy->d3pr = d3pr;
    theCopy->d5pr = d5pr;
    theCopy->ep1pr = ep1pr;
    theCopy->ep3pr = ep3pr;
    theCopy->ep5pr = ep5pr;
    theCopy->q1pr = q1pr;
    theCopy->q3pr = q3pr;
    theCopy->q5pr = q5pr;
    theCopy->Wpr = Wpr;
    theCopy->Fy1pr = Fy1pr;
    theCopy->Fy3pr = Fy3pr;
    theCopy->Fy5pr = Fy5pr;

    return theCopy;
}


int TripleFrictionPendulumX::sendSelf(int commitTag, Channel& theChannel)
{
    // send element parameters
    int res;
    int dataTag = this->getDbTag();
    static Vector data(31);
    data(0) = this->getTag();
    data(1) = tag1;
    data(2) = kpFactor;
    data(3) = kTFactor;
    data(4) = kvFactor;
    data(5) = Mu_ref1;
    data(6) = Mu_ref2;
    data(7) = Mu_ref3;
    data(8) = L1;
    data(9) = L2;
    data(10) = L3;
    data(11) = Ubar1;
    data(12) = Ubar2;
    data(13) = Ubar3;
    data(14) = B1;
    data(15) = B2;
    data(16) = B3;
    data(17) = W;
    data(18) = Uy;
    data(19) = Kvt;
    data(20) = MinFv;
    data(21) = TOL;
    data(22) = refPressure1;
    data(23) = refPressure2;
    data(24) = refPressure3;
    data(25) = Diffusivity;
    data(26) = Conductivity;
    data(27) = Temperature0;
    data(28) = rateParam;
    data(29) = tempParam;
    data(30) = unit;
    
    
    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::sendSelf() - failed to send Vector\n";
        return -1;
    }

    // send the two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::sendSelf() - failed to send ID\n";
        return -2;
    }

    
    // send the material class tags
    ID matClassTags(4);
    for (int i = 0; i < 4; i++)
        matClassTags(i) = theMaterials[i]->getClassTag();
    res = theChannel.sendID(dataTag, commitTag, matClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::sendSelf() - failed to send ID\n";
        return -4;
    }

    // send the material models
    for (int i = 0; i < 4; i++)
        theMaterials[i]->sendSelf(commitTag, theChannel);

    return 0;
}


int TripleFrictionPendulumX::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    
    // delete material memory
    for (int i = 0; i < 4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];

    int res;
    int dataTag = this->getDbTag();
    static Vector data(31);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::recvSelf() - failed to receive Vector\n";
        return -1;
    }

    this->setTag((int)data(0));
    tag1 = (int)data(1);
    kpFactor = (int)data(2);
    kTFactor = (int)data(3);
    kvFactor = (int)data(4);
    Mu_ref1 = data(5);
    Mu_ref2 = data(6);
    Mu_ref3 = data(7);
    L1 = data(8);
    L2 = data(9);
    L3 = data(10);
    Ubar1 = data(11);
    Ubar2 = data(12);
    Ubar3 = data(13);
    B1 = data(14);
    B2 = data(15);
    B3 = data(16);
    W = data(17);
    Uy = data(18);
    Kvt = data(19);
    MinFv = data(20);
    TOL = data(21);
    refPressure1 = data(22);
    refPressure2 = data(23);
    refPressure3 = data(24);
    Diffusivity = data(25);
    Conductivity = data(26);
    Temperature0 = data(27);
    rateParam = data(28);
    tempParam = data(29);
    unit = data(30);
    
    

    // receive the two end nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::recvSelf() - failed to receive ID\n";
        return -2;
    }

    // receive the material class tags
    ID matClassTags(4);
    res = theChannel.recvID(dataTag, commitTag, matClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulumX::recvSelf() - failed to receive ID\n";
        return -5;
    }

    // receive the material models
    for (int i = 0; i < 4; i++) {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "TripleFrictionPendulumX::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -6;
        }
        theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    // initialize constants
    v1Fact = 0.5;
    v3Fact = L2 / (L2 - L1);
    v5Fact = L3 / (L3 - L1);

    Gap2 = 2 * Ubar1 + Ubar2 + Ubar3 + B1 / 2 - Ubar2 * (1 - L1 / L2) - Ubar3 * (1 - L1 / L3); 
    Gap4 = Ubar2 * (1 - L1 / L2);
    Gap6 = Ubar3 * (1 - L1 / L3);

    //Unit conversion for pressure to be used in the pressure factor computation
    p_Unit_Convert; //To convert to MPa
    if (unit == 1) { p_Unit_Convert = 0.000001; }
    if (unit == 2) { p_Unit_Convert = 0.001; }
    if (unit == 3) { p_Unit_Convert = 1.0; }
    if (unit == 4) { p_Unit_Convert = 1000.0; }
    if (unit == 5) { p_Unit_Convert = 0.006894; }
    if (unit == 6) { p_Unit_Convert = 6.894; }
    if (unit == 7) { p_Unit_Convert = 0.00004788; }
    if (unit == 8) { p_Unit_Convert = 0.04788; }

    v_Unit_Convert; //To convert velocity 
    if (unit == 1) { v_Unit_Convert = 1.0; }
    if (unit == 2) { v_Unit_Convert = 1.0; }
    if (unit == 3) { v_Unit_Convert = 0.001; }
    if (unit == 4) { v_Unit_Convert = 0.001; }
    if (unit == 5) { v_Unit_Convert = 0.0254; }
    if (unit == 6) { v_Unit_Convert = 0.0254; }
    if (unit == 7) { v_Unit_Convert = 0.3048; }
    if (unit == 8) { v_Unit_Convert = 0.3048; }
        
    // initialize other variables
    this->revertToStart();

    return 0;
}


int TripleFrictionPendulumX::displaySelf(Renderer& theViewer,
    int displayMode, float fact, const char** modes, int numMode)
{
    int errCode = 0;

    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector& end1Crd = theNodes[0]->getCrds();
    const Vector& end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;

    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);

    if (displayMode >= 0) {
        const Vector& end1Disp = theNodes[0]->getDisp();
        const Vector& end2Disp = theNodes[1]->getDisp();

        for (int i = 0; i < 3; i++) {
            v1(i) = end1Crd(i) + end1Disp(i) * fact;
            v3(i) = end2Crd(i) + end2Disp(i) * fact;
        }
        v2(0) = end1Crd(0) + (end2Disp(0) + xp(1) * end2Disp(5) - xp(2) * end2Disp(4)) * fact;
        v2(1) = end1Crd(1) + (end2Disp(1) - xp(0) * end2Disp(5) + xp(2) * end2Disp(3)) * fact;
        v2(2) = end1Crd(2) + (end2Disp(2) + xp(0) * end2Disp(4) - xp(1) * end2Disp(3)) * fact;
    }
    else {
        int mode = displayMode * -1;
        const Matrix& eigen1 = theNodes[0]->getEigenvectors();
        const Matrix& eigen2 = theNodes[1]->getEigenvectors();

        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
                v3(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
            }
            v2(0) = end1Crd(0) + (eigen2(0, mode - 1) + xp(1) * eigen2(5, mode - 1) - xp(2) * eigen2(4, mode - 1)) * fact;
            v2(1) = end1Crd(1) + (eigen2(1, mode - 1) - xp(0) * eigen2(5, mode - 1) + xp(2) * eigen2(3, mode - 1)) * fact;
            v2(2) = end1Crd(2) + (eigen2(2, mode - 1) + xp(0) * eigen2(4, mode - 1) - xp(1) * eigen2(3, mode - 1)) * fact;
        }
        else {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i);
                v2(i) = end1Crd(i);
                v3(i) = end2Crd(i);
            }
        }
    }

    errCode += theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(), 0);
    errCode += theViewer.drawLine(v2, v3, 1.0, 1.0, this->getTag(), 0);

    return errCode;
}


void TripleFrictionPendulumX::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag();
        s << "  type: TripleFrictionPendulumX, iNode: " << externalNodes(0);
        s << ", jNode: " << externalNodes(1) << endln;
        s << "  Materials: " << theMaterials[0]->getTag() << ", ";
        s << theMaterials[1]->getTag() << ", " << theMaterials[2]->getTag();
        s << ", " << theMaterials[3]->getTag() << endln;
        s << "  Mu_ref1: " << Mu_ref1 << ", Mu_ref2: " << Mu_ref2 << ", Mu_ref3: " << Mu_ref3 << endln;
        s << "  L1: " << L1 << ", L2: " << L2 << ", L3: " << L3 << endln;
        s << "  d1: " << Ubar1 << ", d2: " << Ubar2 << ", d3: " << Ubar3 << endln;
        s << "  B1: " << B1 << ", B2: " << B2 << ", B3: " << B3 << endln;
        s << "  uy: " << Uy << ", kvt: " << Kvt << ",  minFv: " << MinFv << endln;
        s << "  refPressure1: " << refPressure1 << ", refPressure2: " << refPressure2 << ",  refPressure3: " << refPressure3 << endln;
        s << "  Diffusivity: " << Diffusivity << ", Conductivity: " << Conductivity << ",  Temperature0: " << Temperature0 << endln;
        s << "  rateParam: " << rateParam << "  tempParam: " << tempParam <<  ", unit: " << unit << endln;
                
        
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"TripleFrictionPendulumX\", ";
        s << "\"nodes\": [" << externalNodes(0) << ", " << externalNodes(1) << "], ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\", \"";
        s << theMaterials[2]->getTag() << "\", \"";
        s << theMaterials[3]->getTag() << "\"], ";
        s << "\"Mu_ref1\": " << Mu_ref1 << ", ";
        s << "\"Mu_ref2\": " << Mu_ref2 << ", ";
        s << "\"Mu_ref3\": " << Mu_ref3 << ", ";
        s << "\"L1\": " << L1 << ", ";
        s << "\"L2\": " << L2 << ", ";
        s << "\"L3\": " << L3 << ", ";
        s << "\"d1\": " << Ubar1 << ", ";
        s << "\"d2\": " << Ubar2 << ", ";
        s << "\"d3\": " << Ubar3 << ", ";
        s << "\"B1\": " << B1 << ", ";
        s << "\"B2\": " << B2 << ", ";
        s << "\"B3\": " << B3 << ", ";
        s << "\"uy\": " << Uy << ", ";
        s << "\"kvt\": " << Kvt << ", ";
        s << "\"minFv\": " << MinFv << ", ";
        s << "\"refPressure1\": " << refPressure1 << ", ";
        s << "\"refPressure2\": " << refPressure2 << ", ";
        s << "\"refPressure3\": " << refPressure3 << ", ";
        s << "\"Diffusivity\": " << Diffusivity << ", ";
        s << "\"Conductivity\": " << Conductivity << ", ";
        s << "\"Temperature0\": " << Temperature0 << ", ";
        s << "\"rateParam\": " << rateParam << ", ";
        s << "\"tempParam\": " << tempParam << ", ";
        s << "\"unit\": " << unit << "}";
    }
}


Response* TripleFrictionPendulumX::setResponse(const char** argv, int argc,
    OPS_Stream& output)
{
    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "TripleFrictionPendulumX");
    output.attr("eleTag", this->getTag());
    output.attr("node1", externalNodes[0]);
    output.attr("node2", externalNodes[1]);

    // global forces
    if (strcmp(argv[0], "force") == 0 ||
        strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 ||
        strcmp(argv[0], "globalForces") == 0)
    {
        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 1, eleR);
    }
    // local forces
    else if (strcmp(argv[0], "localForce") == 0 ||
        strcmp(argv[0], "localForces") == 0)
    {
        output.tag("ResponseType", "N_1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Vy_2");
        output.tag("ResponseType", "Vz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 2, Vector(12));
    }
    // basic forces
    else if (strcmp(argv[0], "basicForce") == 0 ||
        strcmp(argv[0], "basicForces") == 0)
    {
        output.tag("ResponseType", "qb1");
        output.tag("ResponseType", "qb2");
        output.tag("ResponseType", "qb3");
        output.tag("ResponseType", "qb4");
        output.tag("ResponseType", "qb5");
        output.tag("ResponseType", "qb6");

        theResponse = new ElementResponse(this, 3, Vector(6));
    }
    // local displacements
    else if (strcmp(argv[0], "localDisplacement") == 0 ||
        strcmp(argv[0], "localDisplacements") == 0)
    {
        output.tag("ResponseType", "ux_1");
        output.tag("ResponseType", "uy_1");
        output.tag("ResponseType", "uz_1");
        output.tag("ResponseType", "rx_1");
        output.tag("ResponseType", "ry_1");
        output.tag("ResponseType", "rz_1");
        output.tag("ResponseType", "ux_2");
        output.tag("ResponseType", "uy_2");
        output.tag("ResponseType", "uz_2");
        output.tag("ResponseType", "rx_2");
        output.tag("ResponseType", "ry_2");
        output.tag("ResponseType", "rz_2");

        theResponse = new ElementResponse(this, 4, Vector(12));
    }
    // basic displacements
    else if (strcmp(argv[0], "deformation") == 0 ||
        strcmp(argv[0], "deformations") == 0 ||
        strcmp(argv[0], "basicDeformation") == 0 ||
        strcmp(argv[0], "basicDeformations") == 0 ||
        strcmp(argv[0], "basicDisplacement") == 0 ||
        strcmp(argv[0], "basicDisplacements") == 0)
    {
        output.tag("ResponseType", "ub1");
        output.tag("ResponseType", "ub2");
        output.tag("ResponseType", "ub3");
        output.tag("ResponseType", "ub4");
        output.tag("ResponseType", "ub5");
        output.tag("ResponseType", "ub6");

        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // displacement components
    else if (strcmp(argv[0], "compDeformation") == 0 ||
        strcmp(argv[0], "compDeformations") == 0 ||
        strcmp(argv[0], "compDisplacement") == 0 ||
        strcmp(argv[0], "compDisplacements") == 0)
    {
        output.tag("ResponseType", "d1x");
        output.tag("ResponseType", "d3x");
        output.tag("ResponseType", "d5x");
        output.tag("ResponseType", "d1y");
        output.tag("ResponseType", "d3y");
        output.tag("ResponseType", "d5y");
        output.tag("ResponseType", "v1x");
        output.tag("ResponseType", "v3x");
        output.tag("ResponseType", "v5x");
        output.tag("ResponseType", "v1y");
        output.tag("ResponseType", "v3y");
        output.tag("ResponseType", "v5y");
        

        theResponse = new ElementResponse(this, 6, Vector(12));
    }

    // parameters
    else if (strcmp(argv[0],"param") == 0 || strcmp(argv[0],"Param") == 0 ||
        strcmp(argv[0],"parameters") == 0 || strcmp(argv[0],"Parameters") == 0)
    {
        output.tag("ResponseType","Temp2 and 3");
        output.tag("ResponseType","Temp1");
        output.tag("ResponseType","Temp4");
        output.tag("ResponseType","COF 2 and 3");
        output.tag("ResponseType","COF1");
        output.tag("ResponseType","COF4");
        output.tag("ResponseType", "HEATFLUX 2 and 3");
        output.tag("ResponseType", "HEATFLUX1");
        output.tag("ResponseType", "HEATFLUX4");
        output.tag("ResponseType", "kpF 2 and 3");
        output.tag("ResponseType", "kpF 1");
        output.tag("ResponseType", "kpF 4");
        output.tag("ResponseType", "kTF 2 and 3");
        output.tag("ResponseType", "kTF 1");
        output.tag("ResponseType", "kTF 4");
        output.tag("ResponseType", "kVF 2 and 3");
        output.tag("ResponseType", "kVF 1");
        output.tag("ResponseType", "kVF 4");
        
        theResponse = new ElementResponse(this, 7, Vector(18));
    }

    // Surface History
    else if (strcmp(argv[0], "Surface") == 0 || strcmp(argv[0], "SurfaceHistory") == 0 )
    {
    output.tag("ResponseType", "Surface_disp_1x");
    output.tag("ResponseType", "Surface_disp_2x");
    output.tag("ResponseType", "Surface_disp_3x");
    output.tag("ResponseType", "Surface_disp_4x");
    output.tag("ResponseType", "Surface_disp_1y");
    output.tag("ResponseType", "Surface_disp_2y");
    output.tag("ResponseType", "Surface_disp_3y");
    output.tag("ResponseType", "Surface_disp_4y");
    output.tag("ResponseType", "Surface_vel_1x");
    output.tag("ResponseType", "Surface_vel_2x");
    output.tag("ResponseType", "Surface_vel_3x");
    output.tag("ResponseType", "Surface_vel_4x");
    output.tag("ResponseType", "Surface_vel_1y");
    output.tag("ResponseType", "Surface_vel_2y");
    output.tag("ResponseType", "Surface_vel_3y");
    output.tag("ResponseType", "Surface_vel_4y");
    
    theResponse = new ElementResponse(this, 8, Vector(16));
    }

    // material output
    else if (strcmp(argv[0], "material") == 0) {
        if (argc > 2) {
            int matNum = atoi(argv[1]);
            if (matNum >= 1 && matNum <= 4)
                theResponse = theMaterials[matNum - 1]->setResponse(&argv[2], argc - 2, output);
        }
    }
    
    output.endTag(); // ElementOutput

    return theResponse;
}


int TripleFrictionPendulumX::getResponse(int responseID, Information& eleInfo)
{
    Vector locForce(12), locDisp(12), basForce(6), basDisp(6), cmpDisp(12), Param(18), Surf(16); 
    switch (responseID) {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());

    case 2:  // local forces
        this->getResistingForce();
        locForce(0) = eleR(2);
        locForce(1) = eleR(0);
        locForce(2) = eleR(1);
        locForce(3) = eleR(5);
        locForce(4) = eleR(3);
        locForce(5) = eleR(4);
        locForce(6) = eleR(8);
        locForce(7) = eleR(6);
        locForce(8) = eleR(7);
        locForce(9) = eleR(11);
        locForce(10) = eleR(9);
        locForce(11) = eleR(10);
        return eleInfo.setVector(locForce);

    case 3:  // basic forces
        this->getResistingForce();
        basForce(0) = eleR(8);
        basForce(1) = eleR(6);
        basForce(2) = eleR(7);
        basForce(3) = eleR(11);
        basForce(4) = eleR(9);
        basForce(5) = eleR(10);
        return eleInfo.setVector(basForce);

    case 4:  // local displacements
        locDisp.Zero();
        return eleInfo.setVector(locDisp);

    case 5:  // basic displacements
        basDisp(0) = Dz;
        basDisp(1) = Dx;
        basDisp(2) = Dy;
        basDisp(3) = 0.0;
        basDisp(4) = 0.0;
        basDisp(5) = 0.0;
        return eleInfo.setVector(basDisp);

    case 6:  // displacement components 
        cmpDisp(0) = u23xx * v1Fact; 
        cmpDisp(1) = d3(0) * v3Fact;
        cmpDisp(2) = d5(0) * v5Fact; 
        cmpDisp(3) = u23yy * v1Fact; 
        cmpDisp(4) = d3(1) * v3Fact;
        cmpDisp(5) = d5(1) * v5Fact;
        cmpDisp(6) = v23xx * v1Fact;
        cmpDisp(7) = v3(0) * v3Fact;
        cmpDisp(8) = v5(0) * v5Fact;
        cmpDisp(9) = v23yy * v1Fact;
        cmpDisp(10) = v3(1) * v3Fact;
        cmpDisp(11) = v5(1) * v5Fact;
       
        return eleInfo.setVector(cmpDisp);

    case 7:  // Temperature on each sliding surface
        Param(0) = TemperatureCenter1(0);
        Param(1) = TemperatureCenter2(0);
        Param(2) = TemperatureCenter3(0);
        Param(3) = Mu_Adj1;
        Param(4) = Mu_Adj2;
        Param(5) = Mu_Adj3;
        Param(6) = HeatFluxCenter1(0);
        Param(7) = HeatFluxCenter2(0);
        Param(8) = HeatFluxCenter3(0);
        Param(9) = kpF1;
        Param(10) = kpF2;
        Param(11) = kpF3;
        Param(12) = kTF1;
        Param(13) = kTF2;
        Param(14) = kTF3;
        Param(15) = kvF1;
        Param(16) = kvF2;
        Param(17) = kvF3;

        return eleInfo.setVector(Param);
    
    case 8:  // Temperature on each sliding surface
        Surf(0) = u1;
        Surf(1) = u2;
        Surf(2) = u3;
        Surf(3) = u4;
        Surf(4) = u1y;
        Surf(5) = u2y;
        Surf(6) = u3y;
        Surf(7) = u4y;
        Surf(8) = v1x;
        Surf(9) = v2x;
        Surf(10) = v3x;
        Surf(11) = v4x;
        Surf(12) = v1y;
        Surf(13) = v2y;
        Surf(14) = v3y;
        Surf(15) = v4y;

        return eleInfo.setVector(Surf);


    default:
        return -1;
    }
}


void TripleFrictionPendulumX::CircularElasticGap(Matrix &kj, Vector &fj, double Ej, double Gapj, Vector di)
{
    double r = di.Norm();
    if (r == 0) {
        kj.Zero();
        fj.Zero();
    } else {
        double sn = di(1)/r;
        double cs = di(0)/r;
        
        if (r <= Gapj) {
            kj.Zero();
            fj.Zero();
        } else {
            kj(0,0) = Ej*(1 - Gapj/r*sn*sn);
            kj(0,1) = kj(1,0) = Ej*Gapj/r*sn*cs;
            kj(1,1) = Ej*(1 - Gapj/r*cs*cs);
            fj(0) = Ej*(r-Gapj)*cs;
            fj(1) = Ej*(r-Gapj)*sn;
        }
    }
}


void TripleFrictionPendulumX::BidirectionalPlastic(Matrix &ki, Vector &fi, Vector &epitmp, Vector &qitmp, double Fyi, double Ei, double Hi, Vector epi, Vector qi, Vector di)
{                                                   
    // This part was adapted to the source code of Bidirectional section in OpenSees
    Vector xsi;
    Vector ntmp(2);
    double normxsi;
    double fn;
    fi = Ei*(di - epi);  //compute trial stress using elastic tangent(eqn 6.2-7), di = u_n+1, epi = up_n
    xsi = fi - qi;  // predicted stress minus back stress(line 164) // fi = f_n+1 trial, qi = q_n
    normxsi = xsi.Norm(); //xsi_n+1 trial
    fn = normxsi - Fyi;  //yield function (eqn 6.2-9 Dao)
    
    // Elastic step
    if (fn <= 0.0) {
        ki(0,0) = ki(1,1) = Ei; //Elastic tangent modulus, K_T (eqn 6.2-11 Dao)
        ki(1,0) = ki(0,1) = 0.0;
        epitmp = epi; //tmp:"update" n to n+1 (up_n+1)
        qitmp = qi; //q_n+1 (backstress)
    }

    // Plastic step
    else {
        // compute normal vector (6.2-12)
        double n1 = xsi(0)/normxsi;
        double n2 = xsi(1)/normxsi;
        // compute consistency parameter (6.2-13)
        double dlam = fn / (Ei + Hi); // fn = F_n+1 trial

        double A = Ei*Ei/(Ei+Hi);
        double B = Ei*Ei*dlam/normxsi;
        double EB = Ei-B;
        double BA = B-A;
        ki(0,0) = EB + BA*n1*n1; // Elastoplastic tangent modulus (6.2-17, 6.2-24~26)
        ki(1,1) = EB + BA*n2*n2;
        ki(0,1) = ki(1,0) = BA*n1*n2;
        
        n1 = n1*dlam;
        n2 = n2*dlam;
        fi(0) -= Ei*n1; // update f (eqn 6.2-14)
        fi(1) -= Ei*n2;	
        ntmp(0) = n1; 
        ntmp(1) = n2;
        epitmp = epi + ntmp; // update up (eqn 6.2-15)
        qitmp = qi + ntmp*Hi; // update q (eqn 6.2-16)
    }
}


void TripleFrictionPendulumX::Segment(Vector &epitmp, Vector &qitmp, bool &conv, Matrix &kij, Vector &di, Vector epi, Vector qi, Vector f, Vector df, double Fyi, double Ei, double Hi, double Ej, double Gapj, double Tol, int Niter)
{
    Vector dftmp = df;
    Vector dd;
    Matrix ki(2,2);
    Matrix kj(2,2);
    Vector fi(2);
    Vector fj(2);
    Vector fprime(2);
    Matrix invkij(2,2);
    
    kij.Invert(invkij);
    dd = invkij*dftmp;
    int iter = 1;
    epitmp = epi;
    qitmp = qi;
    
    bool enterloop = false;
    
    while (((dd.Norm() > 0.0001*Tol) && (iter <= Niter)) || (!enterloop))
    {
        enterloop = true;
        iter++;
        di = di + dd;
        BidirectionalPlastic(ki, fi, epitmp, qitmp, Fyi, Ei, Hi, epi, qi, di);
        CircularElasticGap(kj, fj, Ej, Gapj, di);
        kij = ki + kj;
        fprime = fi + fj;
        dftmp = f + df - fprime;
        kij.Invert(invkij);
        dd = invkij*dftmp;
    }
    if (iter > Niter)
        conv = false;
    else
        conv = true;
}


void TripleFrictionPendulumX::TFPElement(bool &Conv, Vector &ep1tmp, Vector &ep3tmp, Vector &ep5tmp, Vector &q1tmp, Vector &q3tmp, Vector &q5tmp, Matrix &K, Vector &f, Matrix &k12, Matrix &k34, Matrix &k56, Vector &d1, Vector &d3, Vector &d5, Vector ep1, Vector ep3, Vector ep5, Vector q1, Vector q3, Vector q5, Vector u, Vector dusub, double Fy1, double Fy3, double Fy5, double E1, double E3, double E5, double H1, double H3, double H5, double E2, double E4, double E6, double Gap2, double Gap4, double Gap6, double Tol, int Niter)
{
    Vector df(2);
    Vector du(2);
    du = dusub;
    int iter = 1;
    bool conv = true;
    Vector uprime(2);
    
    Conv = true;
    ep1tmp = ep1;
    ep3tmp = ep3;
    ep5tmp = ep5;
    q1tmp = q1;
    q3tmp = q3;
    q5tmp = q5;
    
    bool enterloop = false;
    while (((du.Norm() > Tol) && (iter <= Niter) && Conv) || (!enterloop))
    {
        enterloop = true;
        iter++;
        df = K*du;
        Segment(ep1tmp, q1tmp, conv, k12, d1, ep1, q1, f, df, Fy1, E1, H1, E2, Gap2, Tol, Niter);
        if (!conv)
        {
            Conv = false;
            break;
        }
        Segment(ep3tmp, q3tmp, conv, k34, d3, ep3, q3, f, df, Fy3, E3, H3, E4, Gap4, Tol, Niter);
        if (!conv)
        {
            Conv = false;
            break;
        }
        Segment(ep5tmp, q5tmp, conv, k56, d5, ep5, q5, f, df, Fy5, E5, H5, E6, Gap6, Tol, Niter);
        if (!conv)
        {
            Conv = false;
            break;
        }
        f = f + df;
        
        uprime(0) = d1(0) + d3(0) + d5(0);
        uprime(1) = d1(1) + d3(1) + d5(1);
        du(0) = u(0) + dusub(0) - uprime(0);
        du(1) = u(1) + dusub(1) - uprime(1);
        StiffnessForm(K, k12, k34, k56);
    }
    if (iter > Niter)
        Conv = false;
}


void TripleFrictionPendulumX::StiffnessForm(Matrix &K, Matrix k12, Matrix k34, Matrix k56)
{
    Matrix K88(8,8);
    Matrix ktt(4,4);
    Matrix Ktmp1(4,4);
    Matrix Ktmp2(4,4);
    Matrix kot(4,4);
    Matrix kto(4,4);
    Matrix invktt(4,4);
    
    K88.Zero();	
    K88(0,0) = k12(0,0);
    K88(0,1) = K88(1,0) = k12(0,1);
    K88(0,4) = K88(4,0) = -k12(0,0);
    K88(0,5) = K88(5,0) = -k12(0,1);
    K88(1,1) = k12(1,1);
    K88(1,4) = K88(4,1) = -k12(0,1);
    K88(1,5) = K88(5,1) = -k12(1,1);
    K88(2,2) = k56(0,0);
    K88(2,3) = K88(3,2) = k56(0,1);
    K88(2,6) = K88(6,2) = -k56(0,0);
    K88(2,7) = K88(7,2) = -k56(0,1);
    K88(3,3) = k56(1,1);
    K88(3,6) = K88(6,3) = -k56(0,1);
    K88(3,7) = K88(7,3) = -k56(1,1);
    K88(4,4) = k12(0,0) + k34(0,0);
    K88(4,5) = K88(5,4) = k12(0,1) + k34(0,1);
    K88(4,6) = K88(6,4) = -k34(0,0);
    K88(4,7) = K88(7,4) = -k34(0,1);
    K88(5,5) = k12(1,1) + k34(1,1);
    K88(5,6) = K88(6,5) = -k34(0,1);
    K88(5,7) = K88(7,5) = -k34(1,1);
    K88(6,6) = k34(0,0) + k56(0,0);
    K88(6,7) = K88(7,6) = k34(0,1) + k56(0,1);
    K88(7,7) = k34(1,1) + k56(1,1);
    
    for (int i=0; i < 4; i++) {
        for (int j=0; j < 4; j++) {
            ktt(i,j) = K88(i+4,j+4);
            kot(i,j) = kto(j,i) = K88(i+4,j);
            Ktmp1(i,j) = K88(i,j);
        }
    }
    invktt.Zero();
    ktt.Invert(invktt);
    Ktmp2 = Ktmp1 - kot*invktt*kto;
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
            K(i,j) = Ktmp2(i+2,j+2);
        }
    }
}


double TripleFrictionPendulumX::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}

