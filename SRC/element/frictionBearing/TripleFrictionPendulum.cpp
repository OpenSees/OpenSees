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

// TripleFrictionPendulum element
// Written by Nhan Dao, nhan.unr@gmail.com

#include "TripleFrictionPendulum.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <FrictionModel.h>
#include <UniaxialMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


// initialize the class wide variables
Matrix TripleFrictionPendulum::eleK(12,12);
Matrix TripleFrictionPendulum::eleKinit(12,12);
Matrix TripleFrictionPendulum::eleM(12,12);
Matrix TripleFrictionPendulum::eleD(12,12);
Vector TripleFrictionPendulum::eleR(12);

static int numTripleFrictionPendulum = 0;


void *OPS_TripleFrictionPendulum()
{
    if (numTripleFrictionPendulum == 0) {
        numTripleFrictionPendulum++;
        opserr << "TripleFrictionPendulum element v2.0.0 - Written by Nhan@unr\n";
    }
    
    // get the id and end nodes 
    int iData[10];
    double dData[11];
    int numData, eleTag;
    
    numData = 10;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid element data";
        return 0;
    }
    eleTag = iData[0];
    
    // get the friction models
    FrictionModel *theFrnMdls[3];
    for (int i=0; i<3; i++)  {
        theFrnMdls[i] = OPS_getFrictionModel(iData[3+i]);
        if (theFrnMdls[i] == 0)  {
            opserr << "WARNING friction model not found\n";
            opserr << "frictionModel: " << iData[3+i] << endln;
            opserr << "TripleFrictionPendulum element: " << eleTag << endln;
            return 0;
        }
    }
    
    // get the uniaxial materials
    UniaxialMaterial *theMaterials[4];
    for (int i=0; i<4; i++)  {
        theMaterials[i] = OPS_getUniaxialMaterial(iData[6+i]);
        if (theMaterials[i] == 0)  {
            opserr << "WARNING uniaxial material not found\n";
            opserr << "uniaxialMaterial: " << iData[6+i] << endln;
            opserr << "TripleFrictionPendulum element: " << eleTag << endln;
            return 0;
        }
    }
    
    numData = 11;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING error reading element" << eleTag << endln;
        return 0;
    }
    
    // create the element and add it to the Domain
    Element *theTripleFrictionPendulum = new TripleFrictionPendulum(eleTag, iData[1], iData[2], theFrnMdls, theMaterials, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10]);
    
    if (theTripleFrictionPendulum == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }
    
    return theTripleFrictionPendulum;
}


// typical constructor
TripleFrictionPendulum::TripleFrictionPendulum(int tag,
    int Nd1, int Nd2,
    FrictionModel **frnmdls, 
    UniaxialMaterial **materials,
    double l1,
    double l2,
    double l3,
    double ubar1,
    double ubar2,
    double ubar3,
    double w,
    double uy,
    double kvt,
    double minFv,
    double tol)
    : Element(tag, ELE_TAG_TripleFrictionPendulum), externalNodes(2),
    L1(l1), L2(l2), L3(l3), Ubar1(ubar1), Ubar2(ubar2), Ubar3(ubar3),
    W(w), Uy(uy), Kvt(kvt), MinFv(minFv), TOL(tol), Niter(20),
    K(2,2), Kpr(2,2), f(2), fpr(2),
    k12(2,2), k12pr(2,2), k34(2,2), k34pr(2,2), k56(2,2), k56pr(2,2),
    d1(2), d1pr(2), d3(2), d3pr(2), d5(2), d5pr(2), v1(2), v3(2), v5(2),
    ep1(2), ep1pr(2), ep3(2), ep3pr(2), ep5(2), ep5pr(2),
    q1(2), q1pr(2), q3(2), q3pr(2), q5(2), q5pr(2),
    ep1tmp(2), ep3tmp(2), ep5tmp(2), q1tmp(2), q3tmp(2), q5tmp(2)
{
    // fill in the ID containing external node info with node id's
    if (externalNodes.Size() != 2) {
        opserr << "FATAL TripleFrictionPendulum::TripleFrictionPendulum() - out of memory, could not create an ID of size 2\n";
        exit(-1);
    }
    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;
    theNodes[0] = 0; 
    theNodes[1] = 0;
    
    // check friction model input
    if (frnmdls == 0)  {
        opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
            << "null friction model array passed.\n";
        exit(-1);
    }
    
    // get copies of the friction models
    for (int i=0; i<3; i++)  {
        if (frnmdls[i] == 0) {
            opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
                "null friction model pointer passed.\n";
            exit(-1);
        }
        theFrnMdls[i] = frnmdls[i]->getCopy();
        if (theFrnMdls[i] == 0) {
            opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
                << "failed to copy friction model.\n";
            exit(-1);
        }
    }
    
    // check material input
    if (materials == 0)  {
        opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
            << "null material array passed.\n";
        exit(-1);
    }
    
    // get copies of the uniaxial materials
    for (int i=0; i<4; i++)  {
        if (materials[i] == 0) {
            opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0) {
            opserr << "TripleFrictionPendulum::TripleFrictionPendulum() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }
    
    // initialize constants
    v1Fact = 0.5;
    v3Fact = L2/(L2 - L1);
    v5Fact = L3/(L3 - L1);
    
    Gap2 = 2.0*(L1/L3*Ubar3 + Ubar1);
    Gap4 = Ubar2*(1 - L1/L2);
    Gap6 = Ubar3*(1 - L1/L3);
    
    // initialize other variables
    this->revertToStart();
}


// constructor which should be invoked by an FE_ObjectBroker only
TripleFrictionPendulum::TripleFrictionPendulum()
    : Element(0, ELE_TAG_TripleFrictionPendulum), externalNodes(2),
    L1(0.0), L2(0.0), L3(0.0), Ubar1(0.0), Ubar2(0.0), Ubar3(0.0),
    W(0.0), Uy(0.0), Kvt(0.0), MinFv(0.0), TOL(1E-6), Niter(20),
    K(2,2), Kpr(2,2), f(2), fpr(2),
    k12(2,2), k12pr(2,2), k34(2,2), k34pr(2,2), k56(2,2), k56pr(2,2),
    d1(2), d1pr(2), d3(2), d3pr(2), d5(2), d5pr(2), v1(2), v3(2), v5(2),
    ep1(2), ep1pr(2), ep3(2), ep3pr(2), ep5(2), ep5pr(2),
    q1(2), q1pr(2), q3(2), q3pr(2), q5(2), q5pr(2),
    ep1tmp(2), ep3tmp(2), ep5tmp(2), q1tmp(2), q3tmp(2), q5tmp(2)
{
    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    // set friction model pointers to NULL
    for (int i=0; i<3; i++)
        theFrnMdls[i] = 0;
    
    // set material pointers to NULL
    for (int i=0; i<4; i++)
        theMaterials[i] = 0;
}


//  destructor - provided to clean up any memory
TripleFrictionPendulum::~TripleFrictionPendulum()
{
    // clean up all the friction model objects
    for (int i=0; i<3; i++)
        if (theFrnMdls[i] != 0)
            delete theFrnMdls[i];
    
    // clean up all the material objects
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int TripleFrictionPendulum::getNumExternalNodes() const
{
    return 2;
}


const ID& TripleFrictionPendulum::getExternalNodes() 
{
    return externalNodes;
}


Node** TripleFrictionPendulum::getNodePtrs() 
{
    return theNodes;
}


int TripleFrictionPendulum::getNumDOF()
{
    return 12;
}


// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void TripleFrictionPendulum::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        opserr << "Domain does not exist" << endln;	
        exit(0);
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    Node *end1Ptr, *end2Ptr;    
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);
    if (end1Ptr == 0) {
        opserr << "WARNING TripleFrictionPendulum::setDomain() - at TripleFrictionPendulum " << this->getTag() << " node " <<
            Nd1 << "  does not exist in domain\n";
    
        return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {
        opserr << "WARNING TripleFrictionPendulum::setDomain() - at TripleFrictionPendulum " << this->getTag() << " node " <<
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
        opserr << "TripleFrictionPendulum::setDomain(): 6 dof required at nodes\n";
        return;
    }
}


int TripleFrictionPendulum::commitState()
{
    int errCode = 0;
    
    // commit friction models
    for (int i=0; i<3; i++)
        errCode += theFrnMdls[i]->commitState();
    
    // commit material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->commitState();
    
    // commit the base class
    errCode += this->Element::commitState();
    
    // commit other history variables
    Wpr = Wcr;
    Fy1pr = Fy1; Fy3pr = Fy3; Fy5pr = Fy5;
    Kpr = K;
    fpr = f;
    k12pr = k12; k34pr = k34; k56pr = k56;
    d1pr = d1; d3pr = d3; d5pr = d5;
    ep1pr = ep1tmp; ep3pr = ep3tmp; ep5pr = ep5tmp;
    q1pr = q1tmp; q3pr = q3tmp; q5pr = q5tmp;
    
    return 0;
}


int TripleFrictionPendulum::revertToLastCommit()
{
    int errCode = 0;
    
    // revert friction models
    for (int i=0; i<3; i++)
        errCode += theFrnMdls[i]->revertToLastCommit();
    
    // revert material models
    for (int i=0; i<4; i++)
        errCode += theMaterials[i]->revertToLastCommit();
    
    return 0;
}


int TripleFrictionPendulum::revertToStart()
{
    int errCode = 0;
    Vector tmp1(2), tmp2(2), tmp3(2);
    
    Vel1Avg = Vel3Avg = Vel5Avg = 0.0;
    Fy1pr = Fy3pr = Fy5pr = 0.0;
    Wpr = Wcr = Wavg = W;
    
    // revert friction models and re-initialize
    for (int i=0; i<3; i++)  {
        errCode += theFrnMdls[i]->revertToStart();
        theFrnMdls[i]->setTrial(Wavg, 0.0);
    }
    Fy1 = theFrnMdls[0]->getFrictionCoeff();
    Fy3 = theFrnMdls[1]->getFrictionCoeff();
    Fy5 = theFrnMdls[2]->getFrictionCoeff();
    
    E1 = E2 = 3.0*Fy1/Uy;
    E3 = E4 = 3.0*Fy1/Uy;
    E5 = E6 = 3.0*Fy1/Uy;
    
    double E1p = 1.0/(2.0*L1);
    double E3p = 1.0/(L2 - L1);
    double E5p = 1.0/(L3 - L1);
    
    H1 = E1*E1p/(E1 - E1p);
    H3 = E3*E3p/(E3 - E3p);
    H5 = E5*E5p/(E5 - E5p);
    
    // revert material models
    for (int i=0; i<4; i++)
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
    
    BidirectionalPlastic(k12pr, tmp1, tmp2, tmp3, Fy1, E1, H1, ep1pr, q1pr, d1pr);
    BidirectionalPlastic(k34pr, tmp1, tmp2, tmp3, Fy3, E3, H3, ep3pr, q3pr, d3pr);
    BidirectionalPlastic(k56pr, tmp1, tmp2, tmp3, Fy5, E5, H5, ep5pr, q5pr, d5pr);
    StiffnessForm(Kpr, k12pr, k34pr, k56pr);
    
    return errCode;
}


int TripleFrictionPendulum::update()
{
    // get current time
    Domain *theDomain = this->getDomain();
    double time = theDomain->getCurrentTime();
    
    const Vector &duNd1 = theNodes[0]->getIncrDisp();
    const Vector &duNd2 = theNodes[1]->getIncrDisp();
    const Vector &utrialNd1 = theNodes[0]->getTrialDisp();
    const Vector &utrialNd2 = theNodes[1]->getTrialDisp();
    const Vector &vtrialNd1 = theNodes[0]->getTrialVel();
    const Vector &vtrialNd2 = theNodes[1]->getTrialVel();
    const Vector &uNd1 = theNodes[0]->getDisp();
    const Vector &uNd2 = theNodes[1]->getDisp();
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    
    Vector u(2);
    u(0) = uNd2(0) - uNd1(0);  // converged displacement from previous step
    u(1) = uNd2(1) - uNd1(1);
    Vector utrial(2);
    Dx = utrial(0) = utrialNd2(0) - utrialNd1(0);  // trial displacement (target displacement)
    Dy = utrial(1) = utrialNd2(1) - utrialNd1(1);
    Dz = utrialNd2(2) - utrialNd1(2);
    
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
        if (Fvert > Kvert*DBL_EPSILON) {
            theMaterials[0]->setTrialStrain(DzOld, 0.0);
            Kvert = Kvt;
        }
        Fvert = -MinFv;
    }
    Wcr = -Fvert;
    
    // 2) calculate shear forces and stiffnesses in horizontal direction
    double Tol = dusub.Norm()*TOL;
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
    Wavg = (Wpr + Wcr)/2.0;
    
    // get coefficients of friction
    double Fy1cr, Fy3cr, Fy5cr, dFy1, dFy3, dFy5;
    theFrnMdls[0]->setTrial(Wavg, v1Fact*Vel1Avg);
    theFrnMdls[1]->setTrial(Wavg, v3Fact*Vel3Avg);
    theFrnMdls[2]->setTrial(Wavg, v5Fact*Vel5Avg);
    Fy1cr = theFrnMdls[0]->getFrictionCoeff();
    Fy3cr = theFrnMdls[1]->getFrictionCoeff();
    Fy5cr = theFrnMdls[2]->getFrictionCoeff();
    
    dFy1 = Fy1cr - Fy1pr; dFy3 = Fy3cr - Fy3pr; dFy5 = Fy5cr - Fy5pr;
    Fy1 = Fy1pr; Fy3 = Fy3pr; Fy5 = Fy5pr;
    
    int nDiv = 0; int nWhileIter = 0;
    double TolOriginal = Tol;
    while ((nDiv < 10) && (ErrDisp.Norm() > TolOriginal)) {
        Fy1 = Fy1 + dFy1; Fy3 = Fy3 + dFy3; Fy5 = Fy5 + dFy5;
        TFPElement(Conv, ep1tmp, ep3tmp, ep5tmp, q1tmp, q3tmp, q5tmp, K, f, k12, k34, k56, d1, d3, d5, ep1, ep3, ep5, q1, q3, q5, u, dusub, Fy1, Fy3, Fy5, E1, E3, E5, H1, H3, H5, E2, E4, E6, Gap2, Gap4, Gap6, Tol, Niter);
        if ((!Conv) && (nDiv < 7)){
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
        } else {
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
        v1 = (1.0/ops_Dt)*(d1 - d1pr);
        v3 = (1.0/ops_Dt)*(d3 - d3pr);
        v5 = (1.0/ops_Dt)*(d5 - d5pr);
        Vel1Avg = v1.Norm();
        Vel3Avg = v3.Norm();
        Vel5Avg = v5.Norm();
    }
    if (nDiv == 10) {
        if ((!Conv) || (Tol < ErrDisp.Norm())) {
            opserr << "Warning: isolator " << this->getTag() << " has not converged, ErrDisp = "<< ErrDisp << endln;
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


const Matrix& TripleFrictionPendulum::getTangentStiff()
{
    Matrix a(2,12);
    Matrix aT(12,2);
    
    a.Zero();
    aT.Zero();
    a(0,0) = a(1,1) = -1;
    a(0,6) = a(1,7) =  1;
    aT(0,0) = aT(1,1) = -1;
    aT(6,0) = aT(7,1) =  1;
    eleK= aT*K*a;
    eleK *= -Fvert;
    
    eleK(2,2) = eleK(8,8) = Kvert;
    eleK(2,8) = eleK(8,2) = -Kvert;
    eleK(3,3) = eleK(9,9) = KrotX;
    eleK(3,9) = eleK(9,3) = -KrotX;
    eleK(4,4) = eleK(10,10) = KrotY;
    eleK(4,10) = eleK(10,4) = -KrotY;
    eleK(5,5) = eleK(11,11) = KrotZ;
    eleK(5,11) = eleK(11,5) = -KrotZ;
    
    return eleK;
}


const Matrix& TripleFrictionPendulum::getInitialStiff()
{
    Matrix a(2,12);
    Matrix aT(12,2);
    Matrix Kinit(2,2);
    
    Kinit.Zero();
    Kinit(0,0) = Kinit(1,1) = E1/3.0;
    a.Zero();
    aT.Zero();
    a(0,0) = a(1,1) = -1;
    a(0,6) = a(1,7) =  1;
    aT(0,0) = aT(1,1) = -1;
    aT(6,0) = aT(7,1) =  1;
    eleKinit = aT*Kinit*a;
    eleKinit *= W;
    
    eleKinit(2,2) = eleKinit(8,8) = theMaterials[0]->getInitialTangent();
    eleKinit(2,8) = eleKinit(8,2) = -eleKinit(2,2);
    eleKinit(3,3) = eleKinit(9,9) = theMaterials[2]->getInitialTangent();
    eleKinit(3,9) = eleKinit(9,3) = -eleKinit(3,3);
    eleKinit(4,4) = eleKinit(10,10) = theMaterials[3]->getInitialTangent();
    eleKinit(4,10) = eleKinit(10,4) = -eleKinit(4,4);
    eleKinit(5,5) = eleKinit(11,11) = theMaterials[1]->getInitialTangent();
    eleKinit(5,11) = eleKinit(11,5) = -eleKinit(5,5);
    
    return eleKinit;
}


const Matrix& TripleFrictionPendulum::getDamp()
{
    // zero the matrix
    eleD.Zero();
    
    // add damping tangent from materials
    eleD(2,2) = eleD(8,8) = theMaterials[0]->getDampTangent();
    eleD(2,8) = eleD(8,2) = -eleD(2,2);
    eleD(3,3) = eleD(9,9) = theMaterials[2]->getDampTangent();
    eleD(3,9) = eleD(9,3) = -eleD(3,3);
    eleD(4,4) = eleD(10,10) = theMaterials[3]->getDampTangent();
    eleD(4,10) = eleD(10,4) = -eleD(4,4);
    eleD(5,5) = eleD(11,11) = theMaterials[1]->getDampTangent();
    eleD(5,11) = eleD(11,5) = -eleD(5,5);
    
    return eleD;
}


const Matrix& TripleFrictionPendulum::getMass()
{
    eleM.Zero();
    
    return eleM;
}


const Vector& TripleFrictionPendulum::getResistingForce()
{
    Matrix aT(12,2);
    aT.Zero();
    aT(0,0) = aT(1,1) = -1;
    aT(6,0) = aT(7,1) = 1;
    eleR = aT*f;
    eleR *= -Fvert;
    double Mx, My, Mz;
    Mx = -Fvert*Dy + eleR(7)*Hisolator;
    My = Fvert*Dx - eleR(6)*Hisolator;
    Mz = eleR(6)*Dy - eleR(7)*Dx;
    eleR(3) = eleR(9) = TorqX + Mx/2;
    eleR(4) = eleR(10)= TorqY + My/2;
    eleR(5) = eleR(11)= TorqZ + Mz/2;
    eleR(2) = -Fvert;
    eleR(8) = Fvert;
    
    return eleR;
}


Element* TripleFrictionPendulum::getCopy()
{
    TripleFrictionPendulum *theCopy = new TripleFrictionPendulum(this->getTag(),
        externalNodes(0), externalNodes(1), theFrnMdls, theMaterials, L1, L2, L3,
        Ubar1, Ubar2, Ubar3, W, Uy, Kvt, MinFv, TOL);
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


int TripleFrictionPendulum::sendSelf(int commitTag, Channel &theChannel)
{
    // send element parameters
    int res;
    int dataTag = this->getDbTag();
    static Vector data(12);
    data(0)  = this->getTag();
    data(1)  = L1;
    data(2)  = L2;
    data(3)  = L3;
    data(4)  = Ubar1;
    data(5)  = Ubar2;
    data(6)  = Ubar3;
    data(7)  = W;
    data(8)  = Uy;
    data(9) = Kvt;
    data(10) = MinFv;
    data(11) = TOL;
    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send Vector\n";
        return -1;
    }
    
    // send the two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send ID\n";
        return -2;
    }
    
    // send the friction model class tags
    ID frnClassTags(3);
    for (int i=0; i<3; i++)
        frnClassTags(i) = theFrnMdls[i]->getClassTag();
    res = theChannel.sendID(dataTag, commitTag, frnClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send ID\n";
        return -3;
    }
    
    // send the friction models
    for (int i=0; i<3; i++)
        theFrnMdls[i]->sendSelf(commitTag, theChannel);
    
    // send the material class tags
    ID matClassTags(4);
    for (int i=0; i<4; i++)
        matClassTags(i) = theMaterials[i]->getClassTag();
    res = theChannel.sendID(dataTag, commitTag, matClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send ID\n";
        return -4;
    }
    
    // send the material models
    for (int i=0; i<4; i++)
        theMaterials[i]->sendSelf(commitTag, theChannel);
    
    return 0;
}


int TripleFrictionPendulum::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // delete friction model memory
    for (int i=0; i<3; i++)
        if (theFrnMdls[i] != 0)
            delete theFrnMdls[i];
    
    // delete material memory
    for (int i=0; i<4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
    
    int res;
    int dataTag = this->getDbTag();
    static Vector data(12);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive Vector\n";
        return -1;
    }
    
    this->setTag((int)data(0));
    L1 = data(1);
    L2 = data(2);
    L3 = data(3);
    Ubar1 = data(4);
    Ubar2 = data(5);
    Ubar3 = data(6);
    W = data(7);
    Uy = data(8);
    Kvt = data(9);
    MinFv = data(10);
    TOL = data(11);
    
    // receive the two end nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive ID\n";
        return -2;
    }
    
    // receive the friction model class tag
    ID frnClassTags(3);
    res = theChannel.recvID(dataTag, commitTag, frnClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive ID\n";
        return -3;
    }
    
    // receive the friction model
    for (int i=0; i<3; i++)  {
        theFrnMdls[i] = theBroker.getNewFrictionModel(frnClassTags(i));
        if (theFrnMdls[i] == 0) {
            opserr << "TripleFrictionPendulum::recvSelf() - "
                << "failed to get blank friction model.\n";
            return -4;
        }
        theFrnMdls[i]->recvSelf(commitTag, theChannel, theBroker);
    }
    
    // receive the material class tags
    ID matClassTags(4);
    res = theChannel.recvID(dataTag, commitTag, matClassTags);
    if (res < 0) {
        opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive ID\n";
        return -5;
    }
    
    // receive the material models
    for (int i=0; i<4; i++)  {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "TripleFrictionPendulum::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -6;
        }
        theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }
    
    // initialize constants
    v1Fact = 0.5;
    v3Fact = L2/(L2 - L1);
    v5Fact = L3/(L3 - L1);
    
    Gap2 = 2*(L1/L3*Ubar3 + Ubar1);
    Gap4 = Ubar2*(1 - L1/L2);
    Gap6 = Ubar3*(1 - L1/L3);
    
    // initialize other variables
    this->revertToStart();
    
    return 0;
}


int TripleFrictionPendulum::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    int errCode = 0;
    
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;
    
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    
    if (displayMode >= 0)  {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i=0; i<3; i++)  {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v3(i) = end2Crd(i) + end2Disp(i)*fact;
        }
        v2(0) = end1Crd(0) + (end2Disp(0) + xp(1)*end2Disp(5) - xp(2)*end2Disp(4))*fact;
        v2(1) = end1Crd(1) + (end2Disp(1) - xp(0)*end2Disp(5) + xp(2)*end2Disp(3))*fact;
        v2(2) = end1Crd(2) + (end2Disp(2) + xp(0)*end2Disp(4) - xp(1)*end2Disp(3))*fact;
    } else  {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        
        if (eigen1.noCols() >= mode)  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v3(i) = end2Crd(i) + eigen2(i,mode-1)*fact;
            }
            v2(0) = end1Crd(0) + (eigen2(0,mode-1) + xp(1)*eigen2(5,mode-1) - xp(2)*eigen2(4,mode-1))*fact;
            v2(1) = end1Crd(1) + (eigen2(1,mode-1) - xp(0)*eigen2(5,mode-1) + xp(2)*eigen2(3,mode-1))*fact;
            v2(2) = end1Crd(2) + (eigen2(2,mode-1) + xp(0)*eigen2(4,mode-1) - xp(1)*eigen2(3,mode-1))*fact;
        } else  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i);
                v2(i) = end1Crd(i);
                v3(i) = end2Crd(i);
            }
        }
    }
    
    errCode += theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag(), 0);
    errCode += theViewer.drawLine (v2, v3, 1.0, 1.0, this->getTag(), 0);
    
    return errCode;
}


void TripleFrictionPendulum::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag(); 
        s << "  type: TripleFrictionPendulum, iNode: " << externalNodes(0);
        s << ", jNode: " << externalNodes(1) << endln;
        s << "  FrictionModels: " << theFrnMdls[0]->getTag() << ", ";
        s << theFrnMdls[1]->getTag() << ", " << theFrnMdls[2]->getTag() << endln;
        s << "  Materials: " << theMaterials[0]->getTag() << ", ";
        s << theMaterials[1]->getTag() << ", " << theMaterials[2]->getTag();
        s << ", " << theMaterials[3]->getTag() << endln;
        s << "  L1: " << L1 << ", L2: " << L2 << ", L3: " << L3 << endln;
        s << "  d1: " << Ubar1 << ", d2: " << Ubar2 << ", d3: " << Ubar3 << endln;
        s << "  uy: " << Uy << ", kvt: " << Kvt << ",  minFv: " << MinFv << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"TripleFrictionPendulum\", ";
        s << "\"nodes\": [" << externalNodes(0) << ", " << externalNodes(1) << "], ";
        s << "\"frictionModels\": [\"";
        s << theFrnMdls[0]->getTag() << "\", \"";
        s << theFrnMdls[1]->getTag() << "\", \"";
        s << theFrnMdls[2]->getTag() << "\"], ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\", \"";
        s << theMaterials[2]->getTag() << "\", \"";
        s << theMaterials[3]->getTag() << "\"], ";
        s << "\"L1\": " << L1 << ", ";
        s << "\"L2\": " << L2 << ", ";
        s << "\"L3\": " << L3 << ", ";
        s << "\"d1\": " << Ubar1 << ", ";
        s << "\"d2\": " << Ubar2 << ", ";
        s << "\"d3\": " << Ubar3 << ", ";
        s << "\"uy\": " << Uy << ", ";
        s << "\"kvt\": " << Kvt << ", ";
        s << "\"minFv\": " << MinFv << "}";
    }
}


Response* TripleFrictionPendulum::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","TripleFrictionPendulum");
    output.attr("eleTag",this->getTag());
    output.attr("node1",externalNodes[0]);
    output.attr("node2",externalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 ||
        strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0)
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
        
        theResponse = new ElementResponse(this, 1, eleR);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
    {
        output.tag("ResponseType","N_1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Vy_2");
        output.tag("ResponseType","Vz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 2, Vector(12));
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 ||
        strcmp(argv[0],"basicForces") == 0)
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
        
        theResponse = new ElementResponse(this, 4, Vector(12));
    }
    // basic displacements
    else if (strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || 
        strcmp(argv[0],"basicDeformation") == 0 ||
        strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 ||
        strcmp(argv[0],"basicDisplacements") == 0)
    {
        output.tag("ResponseType","ub1");
        output.tag("ResponseType","ub2");
        output.tag("ResponseType","ub3");
        output.tag("ResponseType","ub4");
        output.tag("ResponseType","ub5");
        output.tag("ResponseType","ub6");
        
        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // displacement components
    else if (strcmp(argv[0],"compDeformation") == 0 ||
        strcmp(argv[0],"compDeformations") == 0 ||
        strcmp(argv[0],"compDisplacement") == 0 ||
        strcmp(argv[0],"compDisplacements") == 0)
    {
        output.tag("ResponseType","d1x");
        output.tag("ResponseType","d1y");
        output.tag("ResponseType","d3x");
        output.tag("ResponseType","d3y");
        output.tag("ResponseType","d5x");
        output.tag("ResponseType","d5y");
        
        theResponse = new ElementResponse(this, 6, Vector(6));
    }
    // friction model output
    else if (strcmp(argv[0],"frictionModel") == 0 || strcmp(argv[0],"frnMdl") == 0 ||
        strcmp(argv[0],"frictionMdl") == 0 || strcmp(argv[0],"frnModel") == 0)  {
        if (argc > 2)  {
            int frnNum = atoi(argv[1]);
            if (frnNum >= 1 && frnNum <= 3)
                theResponse = theFrnMdls[frnNum-1]->setResponse(&argv[2], argc-2, output);
        }
    }
    // material output
    else if (strcmp(argv[0],"material") == 0)  {
        if (argc > 2)  {
            int matNum = atoi(argv[1]);
            if (matNum >= 1 && matNum <= 4)
                theResponse = theMaterials[matNum-1]->setResponse(&argv[2], argc-2, output);
        }
    }
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int TripleFrictionPendulum::getResponse(int responseID, Information &eleInfo)
{
    Vector locForce(12), locDisp(12), basForce(6), basDisp(6), cmpDisp(6);
    switch (responseID)  {
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
        cmpDisp(0) = d1(0);
        cmpDisp(1) = d1(1);
        cmpDisp(2) = d3(0);
        cmpDisp(3) = d3(1);
        cmpDisp(4) = d5(0);
        cmpDisp(5) = d5(1);
        return eleInfo.setVector(cmpDisp);
        
    default:
        return -1;
    }
}


void TripleFrictionPendulum::CircularElasticGap(Matrix &kj, Vector &fj, double Ej, double Gapj, Vector di)
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


void TripleFrictionPendulum::BidirectionalPlastic(Matrix &ki, Vector &fi, Vector &epitmp, Vector &qitmp, double Fyi, double Ei, double Hi, Vector epi, Vector qi, Vector di)
{
    // This part was adapted to the source code of Bidirectional section in OpenSees
    Vector xsi;
    Vector ntmp(2);
    double normxsi;
    double fn;
    fi = Ei*(di - epi);  //trial stress
    xsi = fi - qi;
    normxsi = xsi.Norm();
    fn = normxsi - Fyi;  //yield function
    
    // Elastic step
    if (fn <= 0)
    {
        ki(0,0) = ki(1,1) = Ei;
        ki(1,0) = ki(0,1) = 0.0;
        epitmp = epi;
        qitmp = qi;
    }
    else {
        double dlam = fn/(Ei+Hi);
        double n1 = xsi(0)/normxsi;
        double n2 = xsi(1)/normxsi;
        double A = Ei*Ei/(Ei+Hi);
        double B = Ei*Ei*dlam/normxsi;
        double EB = Ei-B;
        double BA = B-A;
        ki(0,0) = EB + BA*n1*n1;
        ki(1,1) = EB + BA*n2*n2;
        ki(0,1) = ki(1,0) = BA*n1*n2;
        
        n1 = n1*dlam;
        n2 = n2*dlam;
        fi(0) -= Ei*n1;
        fi(1) -= Ei*n2;	
        ntmp(0) = n1;
        ntmp(1) = n2;
        epitmp = epi + ntmp;
        qitmp = qi + ntmp*Hi;
    }
}


void TripleFrictionPendulum::Segment(Vector &epitmp, Vector &qitmp, bool &conv, Matrix &kij, Vector &di, Vector epi, Vector qi, Vector f, Vector df, double Fyi, double Ei, double Hi, double Ej, double Gapj, double Tol, int Niter)
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
    
    while (((dd.Norm() > 0.01*Tol) && (iter <= Niter)) || (!enterloop))
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


void TripleFrictionPendulum::TFPElement(bool &Conv, Vector &ep1tmp, Vector &ep3tmp, Vector &ep5tmp, Vector &q1tmp, Vector &q3tmp, Vector &q5tmp, Matrix &K, Vector &f, Matrix &k12, Matrix &k34, Matrix &k56, Vector &d1, Vector &d3, Vector &d5, Vector ep1, Vector ep3, Vector ep5, Vector q1, Vector q3, Vector q5, Vector u, Vector dusub, double Fy1, double Fy3, double Fy5, double E1, double E3, double E5, double H1, double H3, double H5, double E2, double E4, double E6, double Gap2, double Gap4, double Gap6, double Tol, int Niter)
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


void TripleFrictionPendulum::StiffnessForm(Matrix &K, Matrix k12, Matrix k34, Matrix k56)
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
