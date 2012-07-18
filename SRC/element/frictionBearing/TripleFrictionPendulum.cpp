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
// TripleFrictionPendulum element
// Written by Nhan Dao, nhan.unr@gmail.com

#include "TripleFrictionPendulum.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix TripleFrictionPendulum::eleK(12,12);
Matrix TripleFrictionPendulum::eleKinit(12,12);
Matrix TripleFrictionPendulum::eleM(12,12);
Matrix TripleFrictionPendulum::eleD(12,12);
Vector TripleFrictionPendulum::eleR(12);

static int numTripleFrictionPendulum = 0;

void *
OPS_TripleFrictionPendulum()
{
  if (numTripleFrictionPendulum == 0) {
    numTripleFrictionPendulum++;
    opserr << "TripleFrictionPendulum element v1.0.1 - Written by Nhan@unr\n";
  }

  // get the id and end nodes 
  int iData[3];
  double dData[34];
  int numData;

  numData = 3;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data";
    return 0;
  }

  int eleTag = iData[0];

  numData = 34;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element" << eleTag << endln;
    return 0;
  }

  // create the element and add it to the Domain
  Element *theTripleFrictionPendulum = new TripleFrictionPendulum(eleTag, iData[1], iData[2], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20], dData[21], dData[22], dData[23], dData[24], dData[25], dData[26], dData[27], dData[28], dData[29], dData[30], dData[31], dData[32],dData[33]);


  if (theTripleFrictionPendulum == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }
  //opserr << "Finish initializing" << endln;
  return theTripleFrictionPendulum;
}

// typical constructor
TripleFrictionPendulum::TripleFrictionPendulum(int tag,
					       int Nd1, int Nd2,
					       double l1,
					       double l2,
					       double l3,
					       double ubar1,
					       double ubar2,
					       double ubar3,
					       double mu1slow,
					       double mu1fast,
					       double n1slow,
					       double n1fast,
					       double alpha10,
					       double alpha11,
					       double alpha12,
					       double mu2slow,
					       double mu2fast,
					       double n2slow,
					       double n2fast,
					       double alpha20,
					       double alpha21,
					       double alpha22,
					       double mu3slow,
					       double mu3fast,
					       double n3slow,
					       double n3fast,
					       double alpha30,
					       double alpha31,
					       double alpha32,
					       double w,
					       double uy,
					       double kvc,
					       double kvt,
					       double minFv,
					       double maxMuFac,
					       double tol)
 :Element(tag,ELE_TAG_TripleFrictionPendulum),
  externalNodes(2), L1(l1), L2(l2), L3(l3), Ubar1(ubar1), Ubar2(ubar2), Ubar3(ubar3),
  Mu1slow(mu1slow), Mu1fast(mu1fast), N1slow(n1slow), N1fast(n1fast), Alpha10(alpha10), Alpha11(alpha11), Alpha12(alpha12),
  Mu2slow(mu2slow), Mu2fast(mu2fast), N2slow(n2slow), N2fast(n2fast), Alpha20(alpha20), Alpha21(alpha21), Alpha22(alpha22),
  Mu3slow(mu3slow), Mu3fast(mu3fast), N3slow(n3slow), N3fast(n3fast), Alpha30(alpha30), Alpha31(alpha31), Alpha32(alpha32),
  W(w), Uy(uy), Kvc(kvc), Kvt(kvt), MinFv(minFv), MaxMuFac(maxMuFac), TOL(tol),
  K(2,2), Kpr(2,2), f(2), fpr(2),
  k12(2,2), k12pr(2,2), k34(2,2), k34pr(2,2), k56(2,2), k56pr(2,2),
  d1(2), d1pr(2), d3(2), d3pr(2), d5(2), d5pr(2),
  ep1(2), ep1pr(2), ep3(2), ep3pr(2), ep5(2), ep5pr(2),
  q1(2), q1pr(2), q3(2), q3pr(2), q5(2), q5pr(2),
  ep1tmp(2), ep3tmp(2), ep5tmp(2)
{
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL TripleFrictionPendulum::TripleFrictionPendulum() - out of memory, could not create an ID of size 2\n";
    exit(-1);
  }
  Niter = 20;
  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;
  theNodes[0] = 0; 
  theNodes[1] = 0;
  
  A1slow = Mu1slow*pow(W,(1-N1slow));
  A1fast = Mu1fast*pow(W,(1-N1fast));
  av1 = Alpha10 + Alpha11*W + Alpha12*W*W;
  A2slow = Mu2slow*pow(W,(1-N2slow));
  A2fast = Mu2fast*pow(W,(1-N2fast));
  av2 = Alpha20 + Alpha21*W + Alpha22*W*W;
  A3slow = Mu3slow*pow(W,(1-N3slow));
  A3fast = Mu3fast*pow(W,(1-N3fast));
  av3 = Alpha30 + Alpha31*W + Alpha32*W*W;
  avbar1 = av1/2;
  avbar2 = L2*av2/(L2-L1);
  avbar3 = L3*av3/(L3-L1);
  Vel1Avg = Vel3Avg = Vel5Avg =0.;
  Fy1pr = 0.; Fy3pr = 0.; Fy5pr = 0.;

  Wpr = W;
  Wcr = W;
  Wavg = W;
  Fy1 = A1fast*pow(Wavg,(N1fast-1))*(1-exp(-avbar1*Vel1Avg))+A1slow*pow(Wavg,(N1slow-1))*exp(-avbar1*Vel1Avg);
  Fy3 = A2fast*pow(Wavg,(N2fast-1))*(1-exp(-avbar2*Vel3Avg))+A2slow*pow(Wavg,(N2slow-1))*exp(-avbar2*Vel3Avg);
  Fy5 = A3fast*pow(Wavg,(N3fast-1))*(1-exp(-avbar3*Vel5Avg))+A3slow*pow(Wavg,(N3slow-1))*exp(-avbar3*Vel5Avg);
  
  E1 = E2 = 3*Fy1/Uy;
  E3 = E4 = 3*Fy1/Uy;
  E5 = E6 = 3*Fy1/Uy;
  
  double E1p = 1./2./L1;
  double E3p = 1./(L2 - L1);
  double E5p = 1./(L3 - L1);
  H1 = E1*E1p/(E1-E1p);
  H3 = E3*E3p/(E3-E3p);
  H5 = E5*E5p/(E5-E5p);
  Gap2 = 2*(L1/L3*Ubar3 + Ubar1);
  Gap4 = Ubar2*(1-L1/L2);
  Gap6 = Ubar3*(1-L1/L3);	
  
  Vector tmp1(2), tmp2(2), tmp3(2);	
  
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
  
  //opserr << "Finish Constructor" << endln;
}

// constructor which should be invoked by an FE_ObjectBroker only
TripleFrictionPendulum::TripleFrictionPendulum()
 :Element(0,ELE_TAG_TripleFrictionPendulum),
  externalNodes(2), L1(0.0), L2(0.0), L3(0.0), Ubar1(0.0), Ubar2(0.0), Ubar3(0.0),
  Mu1slow(0.0), Mu1fast(0.0), N1slow(0.0), N1fast(0.0), Alpha10(0.0), Alpha11(0.0), Alpha12(0.0),
  Mu2slow(0.0), Mu2fast(0.0), N2slow(0.0), N2fast(0.0), Alpha20(0.0), Alpha21(0.0), Alpha22(0.0),
  Mu3slow(0.0), Mu3fast(0.0), N3slow(0.0), N3fast(0.0), Alpha30(0.0), Alpha31(0.0), Alpha32(0.0),
  W(0.0), Uy(0.0), Kvc(0.0), Kvt(0.0), MaxMuFac(0.0), MinFv(0.0)
{
  //opserr << "Acessing FE_OjbectBroker" << endln;
  theNodes[0] = 0; 
  theNodes[1] = 0;
  //opserr << "Finish FE_OjbectBroker" << endln;
}

//  destructor - provided to clean up any memory
TripleFrictionPendulum::~TripleFrictionPendulum()
{
  //opserr << "Acessing Destructor" << endln;
  
  //opserr << "Finish Destructor" << endln;
}

int
TripleFrictionPendulum::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
TripleFrictionPendulum::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
TripleFrictionPendulum::getNodePtrs(void) 
{
  return theNodes;
}

int
TripleFrictionPendulum::getNumDOF(void) {
  return 12;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
TripleFrictionPendulum::setDomain(Domain *theDomain)
{
  //opserr << "Acessing setDomain" << endln;	
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    opserr << "Domain doesnot exist" << endln;	
    exit(0);
  }
  
  // first ensure nodes exist in Domain and set the node pointers
  Node *end1Ptr, *end2Ptr;    
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  end1Ptr=(*theDomain).getNode(Nd1);
  end2Ptr=(*theDomain).getNode(Nd2);
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
  int dofNd1 = (*end1Ptr).getNumberDOF();
  int dofNd2 = (*end2Ptr).getNumberDOF();
  
  if ((dofNd1 != 6) || (dofNd2 != 6)) {
    opserr << "TripleFrictionPendulum::setDomain(): 6 dof required at nodes\n";
    return;
  }	
  //opserr << "Finish setDomain" << endln;	
}


int
TripleFrictionPendulum::commitState()
{
  // opserr << "Access commitState" << endln;
  this->Element::commitState();
  Wpr = Wcr;
  Fy1pr = Fy1; Fy3pr = Fy3; Fy5pr = Fy5;
  Kpr = K;
  fpr = f;
  k12pr = k12; k34pr = k34; k56pr = k56;
  d1pr = d1; d3pr = d3; d5pr = d5;
  ep1pr = ep1tmp; ep3pr = ep3tmp; ep5pr = ep5tmp;
  q1pr = q1tmp; q3pr = q3tmp; q5pr = q5tmp;
  //opserr << "Finish commitState" << endln;
  return 0;
}

int
TripleFrictionPendulum::revertToLastCommit()
{
  //opserr << "Access revertToLastCommit" << endln;
  //opserr << "Finish revertToLastCommit" << endln;
  return 0;
}

int
TripleFrictionPendulum::revertToStart()
{
  //opserr << "Acess revertToStart" << endln;	
  Vector tmp1(2), tmp2(2), tmp3(2);
  Vel1Avg = Vel3Avg = Vel5Avg = 0.;
  Wpr = W; Wavg = W; Wcr = W;
  av1 = Alpha10 + Alpha11*Wavg + Alpha12*Wavg*Wavg;
  av2 = Alpha20 + Alpha21*Wavg + Alpha22*Wavg*Wavg;
  av3 = Alpha30 + Alpha31*Wavg + Alpha32*Wavg*Wavg;
  avbar1 = av1/2;
  avbar2 = L2*av2/(L2-L1);
  avbar3 = L3*av3/(L3-L1);
  Fy1pr = 0.; Fy3pr = 0.; Fy5pr = 0.;
  Fy1 = A1fast*pow(Wavg,(N1fast-1))*(1-exp(-avbar1*Vel1Avg))+A1slow*pow(Wavg,(N1slow-1))*exp(-avbar1*Vel1Avg);
  Fy3 = A2fast*pow(Wavg,(N2fast-1))*(1-exp(-avbar2*Vel3Avg))+A2slow*pow(Wavg,(N2slow-1))*exp(-avbar2*Vel3Avg);
  Fy5 = A3fast*pow(Wavg,(N3fast-1))*(1-exp(-avbar3*Vel5Avg))+A3slow*pow(Wavg,(N3slow-1))*exp(-avbar3*Vel5Avg);
  
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
  //opserr << "Finish revertToStart" << endln;	
  return 0;
}

int
TripleFrictionPendulum::update()
{
  //opserr << "Access update" << endln;
  const Vector &duNd1 = (*theNodes[0]).getIncrDisp();
  const Vector &duNd2 = (*theNodes[1]).getIncrDisp();
  const Vector &utrialNd1 = (*theNodes[0]).getTrialDisp();
  const Vector &utrialNd2 = (*theNodes[1]).getTrialDisp();
  const Vector &uNd1 = (*theNodes[0]).getDisp();
  const Vector &uNd2 = (*theNodes[1]).getDisp();
  const Vector &end1Crd = (*theNodes[0]).getCrds();
  const Vector &end2Crd = (*theNodes[1]).getCrds();
  
  Vector u(2);
  u(0) = uNd2(0) - uNd1(0);	// converged displacement from previous step.
  u(1) = uNd2(1) - uNd1(1);
  Vector utrial(2);
  Dx = utrial(0) = utrialNd2(0) - utrialNd1(0); // trial displacement (target displacement)
  Dy = utrial(1) = utrialNd2(1) - utrialNd1(1);
  
  Vector dusub(2);
  dusub(0) = duNd2(0) - duNd1(0); // incremental displacement
  dusub(1) = duNd2(1) - duNd1(1);
  
  double uvert = utrialNd2(2) - utrialNd1(2); // vertical displacement
  // vertical force and stiffness
  if (uvert > 0) {
    Fvert = uvert*Kvt;
    Kvert = Kvt;
  } else {
    Fvert = uvert*Kvc;
    Kvert = Kvc;
  }		
  // vertical force for computing friction
  if (Fvert < -MinFv) {
    Wcr = -Fvert;
  } else {
    Wcr = MinFv;
  }
  // Isolator height
  Hisolator = end2Crd(2) - end1Crd(2);
  
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
  
  // compute velocity
  const Vector &vNd1 = (*theNodes[0]).getVel();
  const Vector &vNd2 = (*theNodes[1]).getVel();
  Vector Vel(2);
  Vel(0) = vNd2(0) - vNd1(0);	// converged velocity.
  Vel(1) = vNd2(1) - vNd1(1);
  
  // converged velocity of element groups from previous time step
  Vector TotalDisp(2);
  Vector Vel1pr(2), Vel3pr(2), Vel5pr(2);
  TotalDisp = d1pr + d3pr + d5pr;
  if (TotalDisp(0) != 0) {
    Vel1pr(0) = d1pr(0)/TotalDisp(0)*Vel(0);
    Vel3pr(0) = d3pr(0)/TotalDisp(0)*Vel(0);
    Vel5pr(0) = d5pr(0)/TotalDisp(0)*Vel(0);
  } else {
    Vel1pr(0) = Vel3pr(0) = Vel5pr(0) = 0.;
  }
  if (TotalDisp(1) != 0) {
    Vel1pr(1) = d1pr(1)/TotalDisp(1)*Vel(1);
    Vel3pr(1) = d3pr(1)/TotalDisp(1)*Vel(1);
    Vel5pr(1) = d5pr(1)/TotalDisp(1)*Vel(1);
  } else {
    Vel1pr(1) = Vel3pr(1) = Vel5pr(1) = 0.;
  }
  
  Vel1Avg = Vel1pr.Norm();
  Vel3Avg = Vel3pr.Norm();
  Vel5Avg = Vel5pr.Norm();
  
  // averaging vertical force
  Wavg = (Wpr + Wcr)/2;
  
  av1 = Alpha10 + Alpha11*Wavg + Alpha12*Wavg*Wavg;
  av2 = Alpha20 + Alpha21*Wavg + Alpha22*Wavg*Wavg;
  av3 = Alpha30 + Alpha31*Wavg + Alpha32*Wavg*Wavg;
  avbar1 = av1/2;
  avbar2 = L2*av2/(L2-L1);
  avbar3 = L3*av3/(L3-L1);
  double Fy1cr, Fy3cr, Fy5cr, dFy1, dFy3, dFy5;
  Fy1cr = A1fast*pow(Wavg,(N1fast-1))*(1-exp(-avbar1*Vel1Avg))+A1slow*pow(Wavg,(N1slow-1))*exp(-avbar1*Vel1Avg);
  Fy3cr = A2fast*pow(Wavg,(N2fast-1))*(1-exp(-avbar2*Vel3Avg))+A2slow*pow(Wavg,(N2slow-1))*exp(-avbar2*Vel3Avg);
  Fy5cr = A3fast*pow(Wavg,(N3fast-1))*(1-exp(-avbar3*Vel5Avg))+A3slow*pow(Wavg,(N3slow-1))*exp(-avbar3*Vel5Avg);
  
  if (Fy1cr > MaxMuFac*Mu1fast) {
    Fy1cr = MaxMuFac*Mu1fast;
  }
  if (Fy3cr > MaxMuFac*Mu2fast) {
    Fy3cr = MaxMuFac*Mu2fast;
  }
  if (Fy5cr > MaxMuFac*Mu3fast) {
    Fy5cr = MaxMuFac*Mu3fast;
  }
  dFy1 = Fy1cr - Fy1pr; dFy3 = Fy3cr - Fy3pr;dFy5 = Fy5cr - Fy5pr;
  Fy1 = Fy1pr; Fy3 = Fy3pr; Fy5 = Fy5pr;
  
  int nDiv = 0; int nWhileIter = 0;
  double TolOriginal = Tol;
  while ((nDiv < 10) && (ErrDisp.Norm() > TolOriginal)) {
    Fy1 = Fy1 + dFy1; Fy3 = Fy3 + dFy3; Fy5 = Fy5 + dFy5;
    TFPElement(Conv, ep1tmp, ep3tmp, ep5tmp, q1tmp, q3tmp, q5tmp, K, f, k12, k34, k56, d1, d3, d5, ep1, ep3, ep5, q1, q3, q5, u, dusub, Fy1, Fy3, Fy5, E1, E3, E5, H1, H3, H5, E2, E4, E6, Gap2, Gap4, Gap6, Tol, Niter);
    if ((!Conv) && (nDiv < 7)){
      u(0) = uNd2(0) - uNd1(0);
      u(1) = uNd2(1) - uNd1(1);
      dFy1 /= 2; dFy3 /= 2;dFy5 /= 2;
      Fy1 = Fy1pr; Fy3 = Fy3pr; Fy5 = Fy5pr;
      K = Kpr;
      f = fpr;
      k12 = k12pr; k34 = k34pr; k56 = k56pr;
      d1 = d1pr; d3 = d3pr; d5 = d5pr;
      ep1 = ep1pr; ep3 = ep3pr; ep5 = ep5pr;
      q1 = q1pr; q3 = q3pr; q5 = q5pr;
      dusub /= 2;
      nDiv++;
      nWhileIter = 0;
    } else {
      if (nWhileIter >= pow(2.,nDiv)) {
	break;
      }
      ep1 = ep1tmp; ep3 = ep3tmp; ep5 = ep5tmp;
      q1 = q1tmp; q3 = q3tmp; q5 = q5tmp;
      u += dusub;
      ErrDisp(0) = utrial(0) - u(0);
      ErrDisp(1) = utrial(1) - u(1);
      nWhileIter++;
    }
  }
  if (nDiv == 10) {
    if ((!Conv) || (Tol < ErrDisp.Norm())) {
      opserr <<"Warning: isolator " << this->getTag() << " has not converged, ErrDisp = "<< ErrDisp << endln;
    }
  }
  return 0;
}

const Matrix &
TripleFrictionPendulum::getTangentStiff(void)
{
	//opserr << "Access getTangentStiff" << endln;
	Matrix a(2,12);
	Matrix aT(12,2);
	
	a.Zero();
	aT.Zero();	
	a(0,0) = a(1,1) = -1;
	a(0,6) = a(1,7) = 1;
	aT(0,0)=aT(1,1)=-1;
	aT(6,0)=aT(7,1)=1;
	eleK= aT*K*a;
	if (Fvert < -MinFv) {
		eleK *= -Fvert;
	} else {
		eleK *= MinFv;
	}
	if (Fvert > 0.0) {
		eleK(2,2) = eleK(8,8) = Kvt;
	} else {
		eleK(2,2) = eleK(8,8) = Kvc;
	}
	//opserr << "Finish getTangentStiff" << endln;
	return eleK;
}

const Matrix &
TripleFrictionPendulum::getInitialStiff(void)
{
	//opserr << "Access getInitialStiff" << endln;
	Matrix a(2,12);
	Matrix aT(12,2);
	Matrix Kinit(2,2);
	
	Kinit.Zero();
	Kinit(0,0) = Kinit(1,1) = E1/3;
	a.Zero();
	aT.Zero();	
	a(0,0) = a(1,1) = -1;
	a(0,6) = a(1,7) = 1;
	aT(0,0)=aT(1,1)=-1;
	aT(6,0)=aT(7,1)=1;
	eleKinit = aT*Kinit*a;
	eleKinit *= W;
	eleKinit(2,2) = eleKinit(8,8) = Kvc;
	//opserr << "Finish getInitialStiff" << endln;
	return eleKinit;
}

const Matrix &
TripleFrictionPendulum::getDamp(void)
{
	//opserr << "Access getDamp" << endln;
	eleD.Zero();
	//opserr << "Finish getDamp" << endln;
  	return eleD;  
}


const Matrix &
TripleFrictionPendulum::getMass(void)
{ 
	//opserr << "Access getMass" << endln;
  eleM.Zero();
  //opserr << "Finish getMass" << endln;
  return eleM;
}

const Vector &
TripleFrictionPendulum::getResistingForce()
{
	//opserr << "Access getResistingForce" << endln;
	Matrix aT(12,2);
	aT.Zero();
	aT(0,0) = aT(1,1) = -1;
	aT(6,0) = aT(7,1) = 1;
	eleR = aT*f;	
	double Mx;
	double My;
	double Mz;
	if (Fvert < -MinFv) {
		eleR *= -Fvert;
		Mx = -Fvert*Dy + eleR(7)*Hisolator;
		My = Fvert*Dx - eleR(6)*Hisolator;
		Mz = eleR(6)*Dy - eleR(7)*Dx;
		eleR(3) = eleR(9) = Mx/2;
		eleR(4) = eleR(10)= My/2;
		eleR(5) = eleR(11)= Mz/2;
	} else {
		eleR *= MinFv;
		Mx = eleR(7)*Hisolator;
		My = - eleR(6)*Hisolator;
		Mz = eleR(6)*Dy - eleR(7)*Dx;
		eleR(3) = eleR(9) = Mx/2;
		eleR(4) = eleR(10)= My/2;
		eleR(5) = eleR(11)= Mz/2;
	}
	eleR(2) = -Fvert;
	eleR(8) = Fvert;
	//opserr << "Finish getResistingForce, f(0)=" << f(0) << endln;
	return eleR;
}

Element *
TripleFrictionPendulum::getCopy(void)
{
	//opserr << "Access getCopy" << endln;
    TripleFrictionPendulum *theCopy = new TripleFrictionPendulum(this->getTag(), externalNodes(0), externalNodes(1), L1, L2, L3, Ubar1, Ubar2, Ubar3, Mu1slow, Mu1fast, N1slow, N1fast, Alpha10, Alpha11, Alpha12, Mu2slow, Mu2fast, N2slow, N2fast, Alpha20, Alpha21, Alpha22, Mu3slow, Mu3fast, N3slow, N3fast, Alpha30, Alpha31, Alpha32, W, Uy, Kvc, Kvt, MinFv, MaxMuFac, TOL);
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
    //opserr << "Finish getCopy" << endln;
    return theCopy;
}

int
TripleFrictionPendulum::sendSelf(int commitTag, Channel &theChannel)
{
	//opserr << "Access sendSelf" << endln;
    int res;
    int dataTag = this->getDbTag();
    Vector data(39);
    data(0) = this->getTag();
    data(1) = L1;
    data(2) = L2;
    data(3) = L3;
    data(7) = Ubar1;
    data(8) = Ubar2;
    data(9) = Ubar3;
    data(10) = Mu1slow;
	data(11) = Mu1fast;
	data(12) = N1slow;
	data(13) = N1fast;
	data(14) = Alpha10;
	data(15) = Alpha11;
	data(16) = Alpha12;
	data(17) = Mu2slow;
	data(18) = Mu2fast;
	data(19) = N2slow;
	data(20) = N2fast;
	data(21) = Alpha20;
	data(22) = Alpha21;
	data(23) = Alpha22;
	data(24) = Mu3slow;
	data(25) = Mu3fast;
	data(26) = N3slow;
	data(27) = N3fast;
	data(28) = Alpha30;
	data(29) = Alpha31;
	data(30) = Alpha32;
    data(31)= W;
    data(32)= Uy;
    data(33)= Kvc;
    data(34)= Kvt;
	data(35)= MinFv;
	data(36) = MaxMuFac;
	data(37)= TOL;
    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send Vector\n";
      return -1;
    }	      

    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING TripleFrictionPendulum::sendSelf() - failed to send ID\n";
      return -2;
    }
	//opserr << "Finish sendSelf" << endln;
    return 0;
}

int
TripleFrictionPendulum::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	//opserr << "Access recvSelf" << endln;
    int res;
    int dataTag = this->getDbTag();
    Vector data(17);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive Vector\n";
      return -1;
    }	      

    this->setTag((int)data(0));
    L1 = data(1);
    L2 = data(2);
    L3 = data(3);
    Ubar1 = data(7);
    Ubar2 = data(8);
    Ubar3 = data(9);
    Mu1slow = data(10);
	Mu1fast = data(11);
	N1slow = data(12);
	N1fast = data(13);
	Alpha10 = data(14);
	Alpha11 = data(15);
	Alpha12 = data(16);
	Mu2slow = data(17);
	Mu2fast = data(18);
	N2slow = data(19);
	N2fast = data(20);
	Alpha20 = data(21);
	Alpha21 = data(22);
	Alpha22 = data(23);
	Mu3slow = data(24);
	Mu3fast = data(25);
	N3slow = data(26);
	N3fast = data(27);
	Alpha30 = data(28);
	Alpha31 = data(29);
	Alpha32 = data(30);
    W = data(31);
    Uy = data(32);
    Kvc = data(33);
    Kvt = data(34);
	MinFv = data(35);
	MaxMuFac = data(36);
	TOL = data(37);
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING TripleFrictionPendulum::recvSelf() - failed to receive ID\n";
      return -2;
    }
	//opserr << "Finish recvSelf" << endln;
    return 0;
}

void
TripleFrictionPendulum::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << ". Type: TripleFrictionPendulum,  iNode: " << externalNodes(0);
  s << ", jNode: " << externalNodes(1);
}



Response *
TripleFrictionPendulum::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	//opserr << "ACCESSED setResponse\n";
	return 0;
}



int 
TripleFrictionPendulum::getResponse(int responseID, Information &eleInfo)
{
	return 0;
}

//////////////////////////////////////////////
void
TripleFrictionPendulum::CircularElasticGap(Matrix &kj, Vector &fj, double Ej,double Gapj,Vector di)
{
	double r = di.Norm();
	if (r==0)
	{
		kj.Zero();
		fj.Zero();
	}
	else
	{
		double sn = di(1)/r;
		double cs = di(0)/r;
		if (r <= Gapj)
		{
			kj.Zero();
			fj.Zero();
		}
		else
		{
			kj(0,0) = Ej*(1 - Gapj/r*sn*sn);
			kj(0,1) = kj(1,0) = Ej*Gapj/r*sn*cs;
			kj(1,1) = Ej*(1 - Gapj/r*cs*cs);
			fj(0) = Ej*(r-Gapj)*cs;
			fj(1) = Ej*(r-Gapj)*sn;
		}
	}	
}

void
TripleFrictionPendulum::BidirectionalPlastic(Matrix &ki, Vector &fi, Vector &epitmp, Vector &qitmp, double Fyi, double Ei, double Hi, Vector epi, Vector qi, Vector di)
{
	// This part was adapted to the source code of Bidirectional section in OpenSees
	Vector xsi;
	Vector ntmp(2);
	double normxsi;
	double fn;
	fi = Ei*(di - epi); 	//trial stress
	xsi = fi - qi;
	normxsi=xsi.Norm();
	fn=normxsi-Fyi; 	//yield function
	
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

void
TripleFrictionPendulum::Segment(Vector &epitmp, Vector &qitmp, bool &conv, Matrix &kij, Vector &di, Vector epi, Vector qi, Vector f, Vector df, double Fyi, double Ei, double Hi, double Ej, double Gapj, double Tol, int Niter)
{
	Vector dftmp = df;
	Vector dd;
	Matrix ki(2,2);
	Matrix kj(2,2);
	Vector fi(2);
	Vector fj(2);
	Vector fprime(2);
	Matrix invkij(2,2);
	
	MInverse(kij,invkij,2);	
	dd=invkij*dftmp;
	register int iter = 1;
	epitmp = epi;
	qitmp = qi;
	
	bool enterloop = false;
	
	while (((dd.Norm() > 0.01*Tol) && (iter <= Niter)) || (!enterloop))
	{
		enterloop = true;
		iter++;
		di = di + dd;
		BidirectionalPlastic(ki, fi, epitmp, qitmp, Fyi, Ei, Hi, epi, qi, di);
		CircularElasticGap(kj, fj, Ej,Gapj,di);
		kij = ki + kj;
		fprime = fi + fj;
		dftmp = f + df - fprime;
		MInverse(kij,invkij,2);
		dd = invkij*dftmp;
	}
	if (iter > Niter)
	{
		conv = false;
	}
	else
	{
		conv = true;
	}
}

void
TripleFrictionPendulum::TFPElement(bool &Conv, Vector &ep1tmp, Vector &ep3tmp, Vector &ep5tmp, Vector &q1tmp, Vector &q3tmp, Vector &q5tmp, Matrix &K, Vector &f, Matrix &k12, Matrix &k34, Matrix &k56, Vector &d1, Vector &d3, Vector &d5, Vector ep1, Vector ep3, Vector ep5, Vector q1, Vector q3, Vector q5, Vector u, Vector dusub, double Fy1, double Fy3, double Fy5, double E1, double E3, double E5, double H1, double H3, double H5, double E2, double E4, double E6, double Gap2, double Gap4, double Gap6, double Tol, int Niter)
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
		f=f + df;
		//opserr << "finside=" << f(0) <<endln;
		uprime(0) = d1(0) + d3(0) + d5(0);
		uprime(1) = d1(1) + d3(1) + d5(1);
		du(0) = u(0) + dusub(0) - uprime(0);
		du(1) = u(1) + dusub(1) - uprime(1);
		StiffnessForm(K,k12,k34,k56);
	}
	if (iter > Niter)
	{
		Conv = false;
	}
}

void
TripleFrictionPendulum::StiffnessForm(Matrix &K, Matrix k12, Matrix k34, Matrix k56)
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
	K88(4,4) = k12(0,0)+k34(0,0);
	K88(4,5) = K88(5,4) = k12(0,1)+k34(0,1);
	K88(4,6) = K88(6,4) = -k34(0,0);
	K88(4,7) = K88(7,4) = -k34(0,1);
	K88(5,5) = k12(1,1)+k34(1,1);
	K88(5,6) = K88(6,5) = -k34(0,1);
	K88(5,7) = K88(7,5) = -k34(1,1);
	K88(6,6) = k34(0,0)+k56(0,0);
	K88(6,7) = K88(7,6) = k34(0,1)+k56(0,1);
	K88(7,7) = k34(1,1)+k56(1,1);

	
	for (int i=0; i < 4; i++) {
		for (int j=0; j < 4; j++) {
			ktt(i,j)=K88(i+4,j+4);
			kot(i,j) = kto(j,i) = K88(i+4,j);			
			Ktmp1(i,j)=K88(i,j);
		}
	}
	invktt.Zero();
	MInverse(ktt,invktt,4);
	Ktmp2 = Ktmp1 - kot*invktt*kto;
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++) {
			K(i,j) = Ktmp2(i+2,j+2);
		}
	}
}

void
TripleFrictionPendulum::MInverse(Matrix MInput,Matrix &Mout,int ncol)
{
	// This part was adapted to the current code in OpenSees
	double tmp;
	int i,ii,j;	
	Matrix MIn = MInput;
	Mout.Zero();
	for (i=0;i<ncol;i++) {
		Mout(i,i)=1.0;
	}
	// Go forward to find U
	for (i=0;i<ncol-1;i++) {
		if (MIn(i,i) == 0.) {MIn(i,i)=1.e-12;}
		tmp = MIn(i,i);
		for (j=i;j<ncol;j++) {MIn(i,j)=MIn(i,j)/tmp;}
		for (j=0;j<ncol;j++) {Mout(i,j)=Mout(i,j)/tmp;}
		for (ii=i+1;ii<ncol;ii++) {
			tmp = MIn(ii,i);
			for (j=i;j<ncol;j++) {MIn(ii,j)=MIn(ii,j)-tmp*MIn(i,j);}
			for (j=0;j<ncol;j++) {Mout(ii,j)=Mout(ii,j)-tmp*Mout(i,j);}
		}
	}
	if (MIn(ncol-1,ncol-1)==0.) {MIn(ncol-1,ncol-1)=1.e-12;}
	tmp=MIn(ncol-1,ncol-1);
	MIn(ncol-1,ncol-1)=1.;
	for (j=0;j<ncol;j++) {Mout(ncol-1,j)=Mout(ncol-1,j)/tmp;}
	//Go backward to find L
	for (i=ncol-1;i>0;i--) {
		for (ii=i-1;ii>=0;ii--) {
			for (j=0;j<ncol;j++) {Mout(ii,j)=Mout(ii,j)-MIn(ii,i)*Mout(i,j);}
		}
	}
	
	// Check results
	Matrix I = Mout;
	I=MInput*Mout;
	for (i=0;i<ncol;i++) {
		if ((I(i,i)<0.99999) || (I(i,i)>1.00001)) {
			opserr << "Invert Matrix failed!" << I(i,i);
			exit(1);
		}
	}
}
