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

// Created: C.McGann, UW, 12.2011
//          adopted from Kathy Petek's BeamContact3D element
//
// Description: This file contains the implementation for the BeamContact3Dp class.

#include "BeamContact3Dp.h"
#include <Information.h>
#include <ElementResponse.h>
#include <CrdTransf.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <Parameter.h>
#include <ContactMaterial3D.h>
#include <elementAPI.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define OPS_Export

static int num_BeamContact3Dp = 0;

OPS_Export void *
OPS_BeamContact3Dp(void)
{
  if (num_BeamContact3Dp == 0) {
    num_BeamContact3Dp++;
    //OPS_Error("BeamContact3Dp element - Written: K.Petek, C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "BeamContact3Dp element - Written: K.Petek, C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args,  want: element BeamContact3Dp eleTag?  iNode? jNode? secondaryNode? radius? crdTransf? matTag? penalty? <cSwitch>?\n";
    return 0;
  }
   
  int    iData[6];
  double dData[2];
  int icSwitch = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0] << endln;
    return 0;  
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement" << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0] << endln;
    return 0;  
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 8;
  while (numRemainingInputArgs >= 1) {
	  numData = 1;
	  if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
		  opserr << "WARNING invalid initial contact flag: element BeamContact3Dp " << iData[0] << endln;
	  	  return 0;
      }
	  numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact3Dp(iData[0], iData[1], iData[2], iData[3], dData[0], *theTransf, *theMaterial,
                                  dData[1], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact3DpElement\n";
    return 0;
  }

  return theElement;
}

// constructors:
BeamContact3Dp::BeamContact3Dp(int tag, int Nd1, int Nd2, int NdS, 
			       double rad, 
			       CrdTransf &coordTransf,
                               NDMaterial &theMat, 
			       double pen, 
			       int cSwitch)
 :Element(tag,ELE_TAG_BeamContact3Dp),    
  crdTransf(0),
  theMaterial(0),
  externalNodes(BC3Dp_NUM_NODE),
  mTangentStiffness(BC3Dp_NUM_DOF, BC3Dp_NUM_DOF),
  mInternalForces(BC3Dp_NUM_DOF),
  meye1(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mg1(BC3Dp_NUM_NDM),
  mg2(BC3Dp_NUM_NDM),
  mg_metric(2,2),
  mn(BC3Dp_NUM_NDM),
  mH(4),
  mIcrd_a(BC3Dp_NUM_NDM),
  mIcrd_b(BC3Dp_NUM_NDM),
  mIcrd_s(BC3Dp_NUM_NDM),
  mDcrd_a(BC3Dp_NUM_NDM),
  mDcrd_b(BC3Dp_NUM_NDM),
  mDcrd_s(BC3Dp_NUM_NDM),
  mDisp_a_n(6),
  mDisp_b_n(6),
  mDisp_s_n(3),
  mQa(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mQb(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mQc(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mc1(BC3Dp_NUM_NDM),
  mBn(BC3Dp_NUM_DOF),
  mBs(BC3Dp_NUM_DOF,2),
  mBphi(3,12),
  mSlip(2)
{
    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;
    externalNodes(2) = NdS;

    mRadius = rad;
    mPenalty = pen;
    mIniContact = cSwitch;
    
    if (mIniContact == 0) {
      inContact          = true;
      was_inContact      = true;
      in_bounds          = true;
    } else {
      inContact          = false;
      was_inContact      = false;
      in_bounds          = true;
    }
    
    // initialize penetration function and contact force
    mGap    = 0.0;
    mLambda = 0.0;
    
    // get copy of the transformation & material object  
    crdTransf = coordTransf.getCopy3d();
    NDMaterial *theMatCopy = theMat.getCopy("ContactMaterial3D");
    if (theMatCopy != 0) {
      theMaterial = (ContactMaterial3D *)theMatCopy;
    } else {
      opserr << "BeamContact3Dp::BeamContact3Dp - material needs to be of type Contact3D for ele: " << this->getTag() << endln;
    }      
    
    // check them:         
    if (!crdTransf) {
      opserr << "Error: BeamContact3d::BeamContact3d: could not create copy of coordinate transformation object" << endln;
      exit(-1);
    } 
    if (theMaterial == 0) {
      opserr << "BeamContact3Dp::BeamContact3Dp - failed allocate material model pointer\n";
      exit(-1);
    }

	// set initialization to true for setDomain function
	mInitialize = true;
}


BeamContact3Dp::BeamContact3Dp()
 :Element(0,ELE_TAG_BeamContact3Dp),    
  crdTransf(0),
  theMaterial(0),
  externalNodes(BC3Dp_NUM_NODE),
  mTangentStiffness(BC3Dp_NUM_DOF, BC3Dp_NUM_DOF),
  mInternalForces(BC3Dp_NUM_DOF),
  meye1(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mg1(BC3Dp_NUM_NDM),
  mg2(BC3Dp_NUM_NDM),
  mg_metric(2,2),
  mn(BC3Dp_NUM_NDM),
  mH(4),
  mIcrd_a(BC3Dp_NUM_NDM),
  mIcrd_b(BC3Dp_NUM_NDM),
  mIcrd_s(BC3Dp_NUM_NDM),
  mDcrd_a(BC3Dp_NUM_NDM),
  mDcrd_b(BC3Dp_NUM_NDM),
  mDcrd_s(BC3Dp_NUM_NDM),
  mDisp_a_n(6),
  mDisp_b_n(6),
  mDisp_s_n(3),
  mQa(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mQb(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mQc(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM),
  mc1(BC3Dp_NUM_NDM),
  mBn(BC3Dp_NUM_DOF),
  mBs(BC3Dp_NUM_DOF,2),
  mBphi(3,12),
  mSlip(2)
{
	mRadius     = 0.0;
    mPenalty    = 0.0;
    mIniContact = 0;
	mGap        = 0.0;
    mLambda     = 0.0;

    inContact     = true;
    was_inContact = true;
    in_bounds     = true;

	mInitialize = false;
}

//  destructor:
BeamContact3Dp::~BeamContact3Dp()
{
}

int
BeamContact3Dp::getNumExternalNodes(void) const
{
  return BC3Dp_NUM_NODE;
}

const ID &
BeamContact3Dp::getExternalNodes(void)
{
  return externalNodes;
}


Node **
BeamContact3Dp::getNodePtrs(void)
{
    return theNodes;                        
}

int
BeamContact3Dp::getNumDOF(void)
{
    return BC3Dp_NUM_DOF;
}

#include <ElementIter.h>

void
BeamContact3Dp::setDomain(Domain *theDomain)
{
  meye1.Zero();
  meye1(0,0) = 1.0;
  meye1(1,1) = 1.0;
  meye1(2,2) = 1.0;
  
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  int NdS = externalNodes(2);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);
  theNodes[2] = theDomain->getNode(NdS);
  
  for (int i = 0; i < 3; i++) {
    if (theNodes[i] == 0) {
      opserr << "BeamContact3Dp::setDomain() - no node with tag: " << theNodes[i] << endln;
      return;  // don't go any further - otherwise segmentation fault
    }
  }  

  // only perform these steps during initial creation of element
  if (mInitialize) {
    //  initialize coordinate vectors and set to initial state
    mIcrd_a = theNodes[0]->getCrds();
    mIcrd_b = theNodes[1]->getCrds();
    mIcrd_s = theNodes[2]->getCrds();    
    mDcrd_a = mIcrd_a;
    mDcrd_b = mIcrd_b;
    mDcrd_s = mIcrd_s;
    mDisp_a_n.Zero();
    mDisp_b_n.Zero();
    mDisp_s_n.Zero();
    
    // Fill the coordinate transformation matrices mQa, mQb
    // initialize the transformation
    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
      opserr << "BeamContact3Dp::setDomain(): Error initializing coordinate transformation";  
      exit(0);
    }
    // Note: the following initialization is already performed in crdTransf
    //       but it doesn't return the values to beamContact element
    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++) {
      mQa(i,0) = initXAxis(i);
      mQa(i,1) = initYAxis(i);
      mQa(i,2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    // if mQa = mQb -> mchi = 0
    mQb = mQa;
    mchi = 0;  
  
    // length of primary segment L
    mL = (mDcrd_b - mDcrd_a).Norm();  
    
    // perform projection to update local coordinate along centerline of primary segment
    mxi = ((mDcrd_b - mDcrd_s)^(mDcrd_b - mDcrd_a)) / ((mDcrd_b - mDcrd_a)^(mDcrd_b - mDcrd_a)) ;  // initial approx

	// adjust cohesion force
    //double area = 1.0*sqrt(mg_metric(0,0)*mg_metric(1,1)-mg_metric(0,1)*mg_metric(1,0));
    double area = 1.0;
  
    theMaterial->ScaleCohesion(area);
    theMaterial->ScaleTensileStrength(area);

    // initial basis function values for use in projection
    mxi = project(mxi);
  
    // check contact state based on projection
    in_bounds = ( (mxi>0.000) && (mxi<1.0000) );
    inContact = ( was_inContact && in_bounds );
  
    // update base vectors g1, g2 for determination of mg_metric
    UpdateBase(mxi);
  }
    // compute mBn, mBs
    ComputeB();

  // call the base class method
  this->DomainComponent::setDomain(theDomain);
}
  
int
BeamContact3Dp::commitState()
{
    // update projection (includes update of Qa, Qb & recalc of mn, mc1)
    mxi = project(mxi);

    // update tangent plane & metric tensor (includes update of Qc, rho2, rho3, g1, g2)
    UpdateBase(mxi);

    // update Bn, Bs for next step
    ComputeB();

    // update Boolean Variables for contact condition
    double tol = 0.000001*mRadius;
    was_inContact  = ( mGap < tol );   
    in_bounds      = ( (mxi > 0.000) && (mxi < 1.000) );
    inContact      = ( was_inContact && in_bounds );

    int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "BeamContact3Dp::commitState () - failed in base class";
    }    
    retVal = theMaterial->commitState();

    return retVal;
}

int
BeamContact3Dp::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
BeamContact3Dp::revertToStart()
{
  if (mIniContact == 0) {
    inContact          = true;
    was_inContact      = true;
    in_bounds          = true;
  } else {
    inContact          = false;
    was_inContact      = false;
    in_bounds          = true;
  }
  
  // reset applicable member variables
  mDcrd_a = mIcrd_a;
  mDcrd_b = mIcrd_b;
  mDcrd_s = mIcrd_s;
  mDisp_a_n.Zero();
  mDisp_b_n.Zero();
  mDisp_s_n.Zero();
  
  mL = (mDcrd_b - mDcrd_a).Norm();  
  
  mxi = ((mDcrd_b - mDcrd_s)^(mDcrd_b - mDcrd_a))/((mDcrd_b - mDcrd_a)^(mDcrd_b - mDcrd_a));
  mxi = project(mxi);
  
  in_bounds      = ( (mxi>0.000) && (mxi<1.0000) );
  inContact      = ( was_inContact && in_bounds );
  
  UpdateBase(mxi);
  ComputeB();
  
  return theMaterial->revertToStart();
}

int
BeamContact3Dp::update(void)
{
  if (mInitialize) {
    double tensileStrength;
    Vector a1(BC3Dp_NUM_NDM);
    Vector b1(BC3Dp_NUM_NDM);
    Vector a1_n(BC3Dp_NUM_NDM);
    Vector b1_n(BC3Dp_NUM_NDM);
    Vector disp_a(6);
    Vector disp_b(6);
    Vector disp_s(3);
    Vector disp_L(BC3Dp_NUM_NDM);
    Vector rot_a(BC3Dp_NUM_NDM);
    Vector rot_b(BC3Dp_NUM_NDM);
    Vector x_c(BC3Dp_NUM_NDM);
    Vector d(BC3Dp_NUM_NDM);
    
    // update secondary node coordinate
    mDcrd_s = mIcrd_s + theNodes[2]->getTrialDisp();
    
    // update nodes a, b  coordinates & rotations
    disp_a = theNodes[0]->getTrialDisp();
    disp_b = theNodes[1]->getTrialDisp();
    disp_s = theNodes[2]->getTrialDisp();
  
    for (int i = 0; i < 3; i++) {
      // compute updated nodal coordinates a, b
      mDcrd_a(i) = mIcrd_a(i) +  disp_a(i);
      mDcrd_b(i) = mIcrd_b(i) +  disp_b(i);
      // compute incremental rotations from step n to n+1
      rot_a(i) = disp_a(i+3) - mDisp_a_n(i+3);
      rot_b(i) = disp_b(i+3) - mDisp_b_n(i+3);
    }
    
    // Get tangent vectors from previous converged step
    a1_n = Geta1();
    b1_n = Getb1();
    
    // Perform linear update of tangent vectors such that
    //  a1_{n+1} = a1_n + theta_{n->n+1}
    //  theta_{n->n+1} = rot_a x a1_n
    a1 = a1_n + CrossProduct(rot_a, a1_n);
    b1 = b1_n + CrossProduct(rot_b, b1_n);
  
    // update centerline projection coordinate
    // based upon coords & tangents at step n+1,  xi value @ step n
    x_c = mDcrd_a*mH(0) + a1*mH(1) + mDcrd_b*mH(2) + b1*mH(3);
  
    // update the penetration function
    d = mDcrd_s - x_c;
    mGap = (mn^d) - mRadius;
    double tol = 0.000001*mRadius;
    if (mGap < tol && in_bounds) {
      inContact = true;
    } else {
      mGap = 0.0;
      inContact = false;
    }
    
    // update normal contact force
    if (was_inContact) {
      mLambda = mPenalty*mGap;
    } else {
      mLambda = 0.0;
    }
    
    // get tensile strength from material
    tensileStrength = theMaterial->getTensileStrength();
    
    // determine trial strain vector based on contact state
    if (inContact) {        
      Vector strain(4);
      Vector slip(2);
      Vector phi_c(3);
      Vector c2n(3);
      Vector c3n(3);
      Vector c2n1(3);
      Vector c3n1(3);
      Vector incDisp_ab(12);
      Vector incDisp_s(3);
      Vector dstar(3);
      
      for (int i = 0; i<3; i++) {
        c2n(i) = mQc(i,1);
        c3n(i) = mQc(i,2);
        incDisp_ab(i)   = disp_a(i)   - mDisp_a_n(i);
        incDisp_ab(i+3) = rot_a(i);
        incDisp_ab(i+6) = disp_b(i)   - mDisp_b_n(i);
        incDisp_ab(i+9) = rot_b(i);
      }
      
      incDisp_s = disp_s - mDisp_s_n;
      
      phi_c = mBphi * incDisp_ab;
      
      c2n1 = c2n + CrossProduct(phi_c, c2n);
      c3n1 = c3n + CrossProduct(phi_c, c3n);
      
      dstar = mDcrd_s - x_c - mrho2*c2n1 - mrho3*c3n1;
      
      slip(0) = mg1 ^ dstar;
      slip(1) = mg2 ^ dstar;
      
      // slip is treated in local (curvilinear) coordinates
      strain(0) = mGap;
      strain(1) = slip(0);
      strain(2) = slip(1);
      strain(3) = -mLambda;      
      
      theMaterial->setTrialStrain(strain);
      mSlip = slip;
    } else {
      Vector strain(4);
      
      strain(0) = mGap;
      strain(1) = 0.0;
      strain(2) = 0.0;
      strain(3) = -mLambda;
  
      theMaterial->setTrialStrain(strain);
      mSlip.Zero();
    }
  }
  mInitialize = true;
  
  return 0;
}

double
BeamContact3Dp::project(double xi)
{
    double xi_P;                            // local value of xi
    double H1;                              // local value of Hermitian Function, H1
    double H2;                              // local value of Hermitian Function, H2
    double H3;                              // local value of Hermitian Function, H3
    double H4;                              // local value of Hermitian Function, H4
    double R;
    double DR;
    double dxi;                             // change in xi
    double xi_P_squared;                   
    double xi_P_cubed;                             
    Vector a1(BC3Dp_NUM_NDM);                // tangent at end a
    Vector b1(BC3Dp_NUM_NDM);                // tangent at end b
    Vector x_c_P(BC3Dp_NUM_NDM);             // current centerline porjection coordinate
    Vector d(BC3Dp_NUM_NDM);                 // distance from secondary node to centerline coord
    Vector tc(BC3Dp_NUM_NDM);                // tangent at projection point = 1st deriv of x_c
    Vector ddx_c(BC3Dp_NUM_NDM);             // 2nd derivative of x_c

    // initialize xi_P to previous value of xi
    xi_P = xi;

    // Perform exponential update of coordinate transforms Qa, Qb
    UpdateTransforms();
    // set tangent vectors from Qa, Qb for calc of x_c_P
    a1 = Geta1();
    b1 = Getb1();

    // calculate current projection location
    xi_P_squared = xi_P*xi_P;
    xi_P_cubed = xi_P*xi_P_squared;
    H1 = 1.0-3.0*xi_P_squared + 2.0*xi_P_cubed;
    H2 = (xi_P - 2.0*xi_P_squared + xi_P_cubed) * mL;
    H3 = 1.0 - H1;
    H4 = (-xi_P_squared + xi_P_cubed) * mL;
    x_c_P = mDcrd_a*H1 + a1*H2 + mDcrd_b*H3 + b1*H4;

    d = mDcrd_s - x_c_P;
    tc = Getdx_c(xi_P);
    R = (d^tc);

    // Iterate to determine a new value of xi
    //    such that:  normal ^ tangent = d ^ tc = 0
    int Gapcount = 0;
    while (fabs(R/mL) > 1.0e-10 && Gapcount < 50) {
        ddx_c = Getddx_c(xi_P);
        DR = (d^ddx_c) - (tc^tc);
        dxi = -R / DR;
        xi_P = xi_P + dxi;

        xi_P_squared = xi_P*xi_P;
        xi_P_cubed = xi_P*xi_P_squared;
        H1 = 1.0-3.0*xi_P_squared + 2.0*xi_P_cubed;
        H2 = (xi_P - 2.0*xi_P_squared + xi_P_cubed) * mL;
        H3 = 1.0 - H1;
        H4 = (-xi_P_squared + xi_P_cubed) * mL;
               
        x_c_P = mDcrd_a*H1 + a1*H2 + mDcrd_b*H3 + b1*H4;
        d = mDcrd_s - x_c_P;
        tc = Getdx_c(xi_P);
        R = (d^tc);
        Gapcount +=1;
    }

    // update norm, n, for current projection
    mn = (mDcrd_s - x_c_P) / ( (mDcrd_s - x_c_P).Norm() );

    // set value of c1 for current projection for use in ComputeQc
    Setc1( (tc/tc.Norm()) );

    // update Hermitian Basis functions for use in ComputeB
    //    and for use in next time step in function update()
    mH(0) = H1;
    mH(1) = H2;
    mH(2) = H3;
    mH(3) = H4;

    return xi_P;
}


int
BeamContact3Dp::UpdateBase(double xi)
// this function calculates g1, g2, mg_metric, and sends value of metric tensor to the material
{
    Vector c1(BC3Dp_NUM_NDM);       
    Vector c2(BC3Dp_NUM_NDM);        // c1, c2, c3 are coord transf. vectors at projected point
    Vector c3(BC3Dp_NUM_NDM);
    Vector tc(BC3Dp_NUM_NDM);    // tangent at projected point = 1st deriv of x_c
    Vector ddx_c(BC3Dp_NUM_NDM); // 2nd deriv of x_c w.r.t. xi
    Vector d_c1(BC3Dp_NUM_NDM);  // deriv of c1 w.r.t xi
    Vector d_c2(BC3Dp_NUM_NDM);  // deriv of c2 w.r.t xi
    Vector d_c3(BC3Dp_NUM_NDM);  // deriv of c3 w.r.t xi
    Matrix Qc(BC3Dp_NUM_NDM, BC3Dp_NUM_NDM);  // coord transf at projected point
    Vector g1(BC3Dp_NUM_NDM);    // temporary vector for g1 (co-variant)
    Vector g2(BC3Dp_NUM_NDM);    // temporary vector for g2 (co-variant)
    Vector test(3);

    // Update interpolated coordinate transform, Qc
    ComputeQc(xi);  // sets new value of mQc
    Qc = mQc;

    // obtain c1, c2, c3 vectors
    c1 = Getc1();
    for (int i = 0; i<3; i++) {
        test(i) = c1(i) - Qc(i,0);
        c2(i) = Qc(i,1);
        c3(i) = Qc(i,2);
    }

    // calculate local coordinates rho2, rho3 for radial vector, r_vec
    //  where: r_vec = rho2*c2 + rho3*c3
    mrho2 = mRadius * (mn^c2); //  = cos(psi)*mRadius
    mrho3 = mRadius * (mn^c3); //  = sin(psi)*mRadius

    // calculate all derivatives
    tc = Getdx_c(xi);
    ddx_c = Getddx_c(xi);
    d_c1 = ddx_c - ((c1^ddx_c)*c1);
    d_c1 = d_c1/(tc.Norm());
    d_c2 = -(d_c1 ^ c2)*c1 + mchi*c3;
    d_c3 = -(d_c1 ^ c3)*c1 - mchi*c2;

    // calculate tangent plane vectors, g1 & g2
    g1 = tc + mrho2*d_c2 + mrho3*d_c3;
    g2 = -mrho3 * c2 + mrho2 * c3;

    // fill metric tensor (covariant)
    mg_metric(0,0) = g1^g1;
    mg_metric(0,1) = g1^g2;
    mg_metric(1,0) = mg_metric(0,1);
    mg_metric(1,1) = g2^g2;

    theMaterial->setMetricTensor(mg_metric);

    Matrix G_metric(2,2);  // contravariant
    double det = (mg_metric(0,0) * mg_metric(1,1)) - (mg_metric(0,1) * mg_metric(1,0)) ;
    G_metric(0,0) =  mg_metric(1,1);
    G_metric(1,0) = -mg_metric(1,0);
    G_metric(0,1) = -mg_metric(0,1);
    G_metric(1,1) =  mg_metric(0,0);
    G_metric = G_metric / det;

    // transform mg1, mg2 to contravariant form
    mg1 = G_metric(0,0)*g1 + G_metric(0,1)*g2;
    mg2 = G_metric(1,0)*g1 + G_metric(1,1)*g2;

    return 0;
}

void
BeamContact3Dp::ComputeB(void)
{
    mBn.Zero();
    mBs.Zero();

    int i;

    Vector a1(3);
    Vector b1(3);
    Vector a1xn(3);
    Vector b1xn(3);
    Vector r(3);
    Vector rxg1(3);
    Vector rxg2(3);

    Matrix Bx(3,12);

    a1 = Geta1();
    b1 = Getb1();
    a1xn = CrossProduct(a1, mn);
    b1xn = CrossProduct(b1, mn);  

    r(0) = mrho2*mQc(0,1) + mrho3*mQc(0,2);
    r(1) = mrho2*mQc(1,1) + mrho3*mQc(1,2);
    r(2) = mrho2*mQc(2,1) + mrho3*mQc(2,2);

    rxg1 = CrossProduct(r,mg1);
    rxg2 = CrossProduct(r,mg2);

    mBn(0)  = mn(0)*mH(0);
    mBn(1)  = mn(1)*mH(0);
    mBn(2)  = mn(2)*mH(0);

    mBn(3)  = a1xn(0)*mH(1);
    mBn(4)  = a1xn(1)*mH(1);
    mBn(5)  = a1xn(2)*mH(1);

    mBn(6)  = mn(0)*mH(2);
    mBn(7)  = mn(1)*mH(2);
    mBn(8)  = mn(2)*mH(2);

    mBn(9)  = b1xn(0)*mH(3);
    mBn(10) = b1xn(1)*mH(3);
    mBn(11) = b1xn(2)*mH(3);

    mBn(12) =  -mn(0);
    mBn(13) =  -mn(1);
    mBn(14) =  -mn(2);

    // Derivatives of H
    double dH1 =      - 6.0*mxi + 6.0*mxi*mxi;
    double dH2 =  1.0 - 4.0*mxi + 3.0*mxi*mxi;
    double dH3 =        6.0*mxi - 6.0*mxi*mxi;
    double dH4 =      - 2.0*mxi + 3.0*mxi*mxi;

    Matrix At(3,3);
    Matrix Bt(3,3);
    Matrix Ct(3,3);
    Matrix Dt(3,3);
    Matrix Et(3,3);
    Matrix Ft(3,3);
    Matrix Gt(3,3);
    Matrix Ht(3,3);
    Matrix Qct(3,3);
    Matrix temp1(3,3);
    Matrix temp2(3,3);
    Matrix temp3(3,3);

    At.Zero();
    At(0,0) = mH(0);
    At(1,1) = mH(0);
    At(2,2) = mH(0);

    Bt.Zero();
    temp1.Zero();
    temp1 = ComputeSkew(a1);
    Bt = 1*mH(1) * temp1;  // sign should be +:  -(a1_skew^T) = -(-a1_skew) = + a1_skew

    Ct.Zero();
    Ct(0,0) = mH(2);
    Ct(1,1) = mH(2);
    Ct(2,2) = mH(2);

    Dt.Zero();
    temp1.Zero();
    temp1 = ComputeSkew(b1);
    Dt = 1*mH(3) * temp1;  // sign should be +:  -(b1_skew^T) = -(-b1_skew) = + b1_skew

    Qct = Transpose(3, 3, mQc);

    Et.Zero();
    temp1.Zero();
    temp1(0,1) = -mQa(0,2);
    temp1(0,2) =  mQa(0,1);
    temp1(1,1) = -mQa(1,2);
    temp1(1,2) =  mQa(1,1);
    temp1(2,1) = -mQa(2,2);
    temp1(2,2) =  mQa(2,1);
    Et = (dH1/mL) * temp1 * Qct;

    temp1.Zero();
    temp1(0,0) = (1-mxi)* mQa(0,0);
    temp1(0,1) =    dH2 * mQa(0,1);
    temp1(0,2) =    dH2 * mQa(0,2);
    temp1(1,0) = (1-mxi)* mQa(1,0);
    temp1(1,1) =    dH2 * mQa(1,1);
    temp1(1,2) =    dH2 * mQa(1,2);
    temp1(2,0) = (1-mxi)* mQa(2,0);
    temp1(2,1) =    dH2 * mQa(2,1);
    temp1(2,2) =    dH2 * mQa(2,2);
    Ft = temp1 * Qct;

    Gt.Zero();
    temp1.Zero();
    temp1(0,1) = -mQb(0,2);
    temp1(0,2) =  mQb(0,1);
    temp1(1,1) = -mQb(1,2);
    temp1(1,2) =  mQb(1,1);
    temp1(2,1) = -mQb(2,2);
    temp1(2,2) =  mQb(2,1);
    Gt = (dH3/mL) * temp1 * Qct;

    Ht.Zero();
    temp1(0,0) = mxi * mQb(0,0);
    temp1(0,1) = dH4 * mQb(0,1);
    temp1(0,2) = dH4 * mQb(0,2);
    temp1(1,0) = mxi * mQb(1,0);
    temp1(1,1) = dH4 * mQb(1,1);
    temp1(1,2) = dH4 * mQb(1,2);
    temp1(2,0) = mxi * mQb(2,0);
    temp1(2,1) = dH4 * mQb(2,1);
    temp1(2,2) = dH4 * mQb(2,2);
    Ht = temp1 * Qct;

    // Fill Bs(0:2,0), Bs(0:2,1)
    Vector col1(3);
    Vector col2(3);
    col1.Zero();
    col2.Zero();
    col1 = At*mg1 + Et*rxg1;
    col2 = At*mg2 + Et*rxg2;
    for (i = 0; i < 3; i++) {
        mBs(i,0) = -col1(i);
        mBs(i,1) = -col2(i);
    }

    // Fill Bs(3:5,0), Bs(3:5,1)
    col1 = Bt*mg1 + Ft*rxg1;
    col2 = Bt*mg2 + Ft*rxg2;
    for (i = 0; i < 3; i++) {
        mBs(i+3,0) = -col1(i);
        mBs(i+3,1) = -col2(i);
    }

    // Fill Bs(6:8,0), Bs(6:8,1)
    col1 = Ct*mg1 + Gt*rxg1;
    col2 = Ct*mg2 + Gt*rxg2;
    for (i = 0; i < 3; i++) {
        mBs(i+6,0) = -col1(i);
        mBs(i+6,1) = -col2(i);
    }

    // Fill Bs(9:11,0), Bs(9:11,1)
    col1 = Dt*mg1 + Ht*rxg1;
    col2 = Dt*mg2 + Ht*rxg2;
    for (i = 0; i < 3; i++) {
        mBs(i+9,0) = -col1(i);
        mBs(i+9,1) = -col2(i);
    }

    // Fill Bs(12:14,0), Bs(12:14,1)
    for (i = 0; i < 3; i++) {
        mBs(i+12,0) = mg1(i);
        mBs(i+12,1) = mg2(i);
    }

    mBphi = ComputeBphi();

    return;
}


Matrix
BeamContact3Dp::ComputeBphi(void)
{
    int i, j;
    Matrix dummy1(3,3);
    Matrix dummy2(3,3);
    Matrix dummy3(3,3);
    Matrix Bphi(3,12);

    Bphi.Zero();

    // Derivatives of H
    double dH1 =      - 6.0*mxi + 6.0*mxi*mxi;
    double dH2 =  1.0 - 4.0*mxi + 3.0*mxi*mxi;
    double dH3 =        6.0*mxi - 6.0*mxi*mxi;
    double dH4 =      - 2.0*mxi + 3.0*mxi*mxi;

    // Compute Bphi(0:2, 3:5)
    dummy1.Zero();
    dummy2.Zero();
    dummy3.Zero();
    dummy1(0,0) = (1-mxi)*mQc(0,0);  // N1 * mQc(0:2,0)
    dummy1(1,0) = (1-mxi)*mQc(1,0);
    dummy1(2,0) = (1-mxi)*mQc(2,0);

    dummy1(0,1) = dH2*mQc(0,1);  // dH2 * mQc(0:2,1:2)
    dummy1(1,1) = dH2*mQc(1,1);
    dummy1(2,1) = dH2*mQc(2,1);
    dummy1(0,2) = dH2*mQc(0,2);
    dummy1(1,2) = dH2*mQc(1,2);
    dummy1(2,2) = dH2*mQc(2,2);
    dummy2 = Transpose(BC3Dp_NUM_NDM, BC3Dp_NUM_NDM, mQa);
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Bphi(i,3+j) = dummy3(i,j);
        }
    }

    // Reuse dummy2 and Compute Bphi(0:2, 0:2)
    dummy1.Zero();
    dummy3.Zero();
    dummy1(0,1) =  mQc(0,2);
    dummy1(0,2) = -mQc(0,1);
    dummy1(1,1) =  mQc(1,2);
    dummy1(1,2) = -mQc(1,1);
    dummy1(2,1) =  mQc(2,2);
    dummy1(2,2) = -mQc(2,1);
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Bphi(i,j) = dH1/mL * dummy3(i,j);
        }
    }

    // Reuse dummy1 and Compute Bphi(0:2, 6:8)
    dummy2.Zero();
    dummy3.Zero();
    dummy2 = Transpose(BC3Dp_NUM_NDM, BC3Dp_NUM_NDM, mQb);
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Bphi(i,6+j) = dH3/mL * dummy3(i,j);
        }
    }

    // Ruse dummy2 and Compute Bphi(0:2, 9:11)
    dummy1.Zero();
    dummy3.Zero();
    dummy1(0,0) = mxi*mQc(0,0);  // N2 * mQc(0:2,0)
    dummy1(1,0) = mxi*mQc(1,0);
    dummy1(2,0) = mxi*mQc(2,0);

    dummy1(0,1) = dH4*mQc(0,1);     // dH4 * mQc(0:2,1:2)
    dummy1(1,1) = dH4*mQc(1,1);
    dummy1(2,1) = dH4*mQc(2,1);
    dummy1(0,2) = dH4*mQc(0,2);
    dummy1(1,2) = dH4*mQc(1,2);
    dummy1(2,2) = dH4*mQc(2,2);
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Bphi(i,9+j) = dummy3(i,j);
        }
    }

    return Bphi;
}

void
BeamContact3Dp::UpdateTransforms(void)
{
    Vector disp_a(6);           // trial disp/rot vector at a(total disp/rot)
    Vector disp_b(6);           // trial disp/rot vector at a(total disp/rot)  
    Vector rot_a(BC3Dp_NUM_NDM);  // incr. rot vector at a (from n->n+1)
    Vector rot_b(BC3Dp_NUM_NDM);  // incr. rot vector at a (from n->n+1)
    Matrix Omega(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);  // Matrix used for Exponential Map

    // Recalculate incremental rotations from step n to n+1
    disp_a = theNodes[0]->getTrialDisp();
    disp_b = theNodes[1]->getTrialDisp();
    int i;
    //for (i=0; i<3; i++) {
    for (i = 3; i<6; i++) {
        rot_a(i-3) = disp_a(i) - mDisp_a_n(i);
        rot_b(i-3) = disp_b(i) - mDisp_b_n(i);
    }
               
    // Perform exponential update of Qa
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_a);
    //   calculate new Qa
    mQa = Omega*mQa;

    // Perform exponential update of Qb
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_b);
    //   calculate new Qb
    mQb = Omega*mQb;
               
    // Reset total disp & rotation vectors for calculation of
    // incremental values in subsequent step
    for (i=0; i<6; i++) {
        mDisp_a_n(i) = disp_a(i);
        mDisp_b_n(i) = disp_b(i);
    }

    mDisp_s_n = theNodes[2]->getTrialDisp();

    return;
}

void
BeamContact3Dp::ComputeQc(double xi)
{      
    Vector c1(BC3Dp_NUM_NDM);        // tangent vector at projection point, c
    Vector a1(BC3Dp_NUM_NDM);        // tangent vector at a
    Vector b1(BC3Dp_NUM_NDM);        // tangent vector at b
    Vector temp(BC3Dp_NUM_NDM);          // dummy vector for use in calcs
    Matrix Qc_df(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);   // Drill free transformation matrix for c
    Matrix Qc_chi(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);  // Twist transformation matrix for c
    Matrix Qb_df(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);   // Drill free transf. matrix from a to b

    // Fill tangent vectors
    a1 = Geta1();
    b1 = Getb1();
    c1 = Getc1();
    temp.Zero();

    // Calculate the drill free transformation from a to c, Qc_df
    temp = CrossProduct(a1,c1);
    Qc_df = ExponentialMap(temp);

    // Calculate the drill free transformation from a to b, Qb_df
    // for determination of twist angle mchi
    temp = CrossProduct(a1,b1);
    Qb_df = ExponentialMap(temp);
    Qb_df = Qb_df * mQa;

    // mchi = arcsin( b3 dot b2_df) = arcsin(mQb(:,2) dot Qb_df(:,1))
    // WATCH SIGN!!!!
    mchi = mQb(0,2)*Qb_df(0,1) + mQb(1,2)*Qb_df(1,1) + mQb(2,2)*Qb_df(2,1);
    mchi = -asin(mchi);

    // Calculate twist transformation from a to c, Qc_df
    // based upon linear scaling of twist angle: mxi * mchi * c1
    temp = mxi*mchi*c1;
    Qc_chi = ExponentialMap(temp);
   
    mQc = (Qc_chi * Qc_df) * mQa;

    return;
}


Matrix
BeamContact3Dp::ExponentialMap(Vector th)
{
// This function transforms a vector th to a skew symmetric matrix
// and performs an exponential map to transform the skew symmetric
// matrix into an orthogonal matrix.

//  exp[sk_theta] = cos(theta)*1 + sin(theta)/theta * sk_theta +
//                   (1-cos(theta))/(theta^2) * theta_vec^T*theta_vec
//
//                = sf1*1 + sf2*k_theta + sf3*theta_vec^T*theta_vec
//
//  approximations are used for small theta when necessary

    double theta;                                  // vector norm
    Vector theta_vec(BC3Dp_NUM_NDM);                // input vector
    Matrix sk_theta(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);    // skew of vector
    Matrix theta_theta(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM); // dyadic product of vector
    Matrix Q(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);           // Exonential Map Vector  

    double sf1;  // local variables for calculation of exp. map terms
    double sf2;
    double sf3;

    Q.Zero();
    sk_theta.Zero();
    theta_theta.Zero();
       
    theta_vec = th;
    theta = theta_vec.Norm();
    sk_theta = ComputeSkew(theta_vec);
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            theta_theta(i,j) = theta_vec(i) * theta_vec(j);
        }
    }

    sf1 = cos(theta);

    if (theta > 5.0e-3) {
        sf2 = (sin(theta))/theta;
    } else {
        // small theta approximation
        sf2 = 1.0 - theta*theta/6.0 + pow(theta,4.0)/120.0;
    }

    if (theta > 0.1) {
        sf3 = (1.0-cos(theta))/(theta*theta);
    } else {
        // small theta approximation
        sf3 = 0.5 - theta*theta/24.0 + pow(theta, 4.0)/720.0 - pow(theta,6.0)/40320.0 + pow(theta, 8.0)/3628800.0;
    }


    Q = sf1*meye1 + sf2*sk_theta + sf3*theta_theta;

    return Q;
}

Matrix
BeamContact3Dp::ComputeSkew(Vector th)
{
    Matrix skew_th(BC3Dp_NUM_NDM,BC3Dp_NUM_NDM);

    skew_th(0,0) =  0.0;
    skew_th(0,1) = -th(2);
    skew_th(0,2) =  th(1);
    skew_th(1,0) =  th(2);
    skew_th(1,1) =  0.0;
    skew_th(1,2) = -th(0);
    skew_th(2,0) = -th(1);
    skew_th(2,1) =  th(0);
    skew_th(2,2) =  0.0;

    return skew_th;
}

Vector
BeamContact3Dp::CrossProduct(Vector &V1, Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Matrix  
BeamContact3Dp::Transpose( int dim1, int dim2, const Matrix &M )
{
    Matrix Mtran(dim2, dim1);

    for ( int i = 0; i < dim1; i++ ) {
        for ( int j = 0; j < dim2; j++ ) {
            Mtran(j,i) = M(i,j);
		}
    }

    return Mtran;
}

Vector
BeamContact3Dp::Geta1(void)
{
    Vector a1(BC3Dp_NUM_NDM);
    int i;
    for (i=0; i<3; i++) {
        a1(i) = mQa(i,0);
    }

    return a1;
}

Vector
BeamContact3Dp::Getb1(void)
{
    Vector b1(BC3Dp_NUM_NDM);
    int i;
    for (i=0; i<3; i++) {
        b1(i) = mQb(i,0);
    }

    return b1;
}

void
BeamContact3Dp::Setc1(Vector c1_vec)
{
    mc1 = c1_vec;

    return;
}

Vector
BeamContact3Dp::Getc1(void)
{
    return mc1;
}

Vector
BeamContact3Dp::Getdx_c(double xi)
{
    double xi_squared;
    Vector a1(BC3Dp_NUM_NDM);
    Vector b1(BC3Dp_NUM_NDM);
    Vector deriv1(BC3Dp_NUM_NDM);

    a1 = Geta1();
    b1 = Getb1();
    xi_squared = xi*xi;

    deriv1 = ( - 6.0*xi + 6.0*xi_squared)*mDcrd_a
               + (1.0 - 4.0*xi + 3.0*xi_squared)*mL*a1
               + (      6.0*xi - 6.0*xi_squared)*mDcrd_b
               + (    - 2.0*xi + 3.0*xi_squared)*mL*b1;

    return deriv1;
}

Vector
BeamContact3Dp::Getddx_c(double xi)
{
    Vector a1(BC3Dp_NUM_NDM);
    Vector b1(BC3Dp_NUM_NDM);
    Vector deriv2(BC3Dp_NUM_NDM);

    a1 = Geta1();
    b1 = Getb1();

    deriv2 = ( - 6.0 + 12.0*xi)*mDcrd_a
                   + ( - 4.0 +  6.0*xi)*mL*a1
                   + (   6.0 - 12.0*xi)*mDcrd_b
                   + ( - 2.0 +  6.0*xi)*mL*b1;

    return deriv2;
}


const Matrix &
BeamContact3Dp::getTangentStiff(void)
{
    mTangentStiffness.Zero();

    if (inContact) { 
        Matrix Cmat = theMaterial->getTangent();
           
		// constitutive terms
        double Cn1  = Cmat(1,3);
        double Cn2  = Cmat(2,3);        
        double Cs11 = Cmat(1,1);
        double Cs12 = Cmat(1,2);
        double Cs21 = Cmat(2,1);
        double Cs22 = Cmat(2,2);

		// build the tangent
		for (int i = 0; i < BC3Dp_NUM_DOF; i++) {
			for (int j = 0; j < BC3Dp_NUM_DOF; j++) {

				mTangentStiffness(i,j) = mPenalty*mBn(i)*mBn(j) + (mBs(i,0)*Cs11 + mBs(i,1)*Cs21)*mBs(j,0)
				                         + (mBs(i,0)*Cs21 + mBs(i,1)*Cs22)*mBs(j,1) 
										 + mPenalty*(mBs(i,0)*Cn1 + mBs(i,1)*Cn2)*mBn(j);
            }
		}
	}

    return mTangentStiffness;
}

const Matrix &
BeamContact3Dp::getInitialStiff(void)
{
    return getTangentStiff();
}
   
void
BeamContact3Dp::zeroLoad(void)
{
    return;
}

int
BeamContact3Dp::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    return 0;
}

int
BeamContact3Dp::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector &
BeamContact3Dp::getResistingForce()
{
    mInternalForces.Zero();

    // get contact "stress":  stress(0) = lambda, stress(1) = t_s(0), stress(2) = t_s(1)
    Vector stress = theMaterial->getStress();

    if (inContact) {
        for (int i = 0; i < BC3Dp_NUM_DOF; i++) {
            mInternalForces(i) = -mLambda*mBn(i) + stress(1)*mBs(i,0) + stress(2)*mBs(i,1);
        }
    }

    return mInternalForces;
}


const Vector &
BeamContact3Dp::getResistingForceIncInertia()
{      
    return getResistingForce();
}

int
BeamContact3Dp::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // BeamContact3Dp packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments
  
  static Vector data(88);
  data(0) = this->getTag();
  data(1) = mRadius;
  data(2) = mPenalty;
  data(3) = mIniContact;
  data(4) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;
  data(6) = crdTransf->getClassTag();
  int crdDbTag = crdTransf->getDbTag();
  
  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (crdDbTag == 0) {
    crdDbTag = theChannel.getDbTag();
    if (crdDbTag != 0)
      crdTransf->setDbTag(crdDbTag);
  }
  
  data(7) = crdDbTag;
  data(8) = inContact;
  data(9) = was_inContact;
  data(10) = in_bounds;
  data(11) = mGap;
  data(12) = mLambda;
  data(13) = mxi;
  data(14) = mchi;
  data(15) = mL;
  data(16) = mrho2;
  data(17) = mrho3;

  data(18) = mH(0);
  data(19) = mH(1);
  data(20) = mH(2);
  data(21) = mH(3);

  data(22) = mn(0);
  data(23) = mn(1);
  data(24) = mn(2);

  data(25) = mg1(0);
  data(26) = mg1(1);
  data(27) = mg1(2);
  data(28) = mg2(0);
  data(29) = mg2(1);
  data(30) = mg2(2);
  
  data(31) = mQa(0,0);
  data(32) = mQa(0,1);
  data(33) = mQa(0,2);
  data(34) = mQa(1,0);
  data(35) = mQa(1,1);
  data(36) = mQa(1,2);
  data(37) = mQa(2,0);
  data(38) = mQa(2,1);
  data(39) = mQa(2,2);
  
  data(40) = mQb(0,0);
  data(41) = mQb(0,1);
  data(42) = mQb(0,2);
  data(43) = mQb(1,0);
  data(44) = mQb(1,1);
  data(45) = mQb(1,2);
  data(46) = mQb(2,0);
  data(47) = mQb(2,1);
  data(48) = mQb(2,2);

  data(49) = mQc(0,0);
  data(50) = mQc(0,1);
  data(51) = mQc(0,2);
  data(52) = mQc(1,0);
  data(53) = mQc(1,1);
  data(54) = mQc(1,2);
  data(55) = mQc(2,0);
  data(56) = mQc(2,1);
  data(57) = mQc(2,2);

  data(58) = mIcrd_a(0);
  data(59) = mIcrd_a(1);
  data(60) = mIcrd_a(2);
  data(61) = mIcrd_b(0);
  data(62) = mIcrd_b(1);
  data(63) = mIcrd_b(2);
  data(64) = mIcrd_s(0);
  data(65) = mIcrd_s(1);
  data(66) = mIcrd_s(2);

  data(67) = mDcrd_a(0);
  data(68) = mDcrd_a(1);
  data(69) = mDcrd_a(2);
  data(70) = mDcrd_b(0);
  data(71) = mDcrd_b(1);
  data(72) = mDcrd_b(2);

  data(73) = mDisp_a_n(0);
  data(74) = mDisp_a_n(1);
  data(75) = mDisp_a_n(2);
  data(76) = mDisp_a_n(3);
  data(77) = mDisp_a_n(4);
  data(78) = mDisp_a_n(5);
  data(79) = mDisp_b_n(0);
  data(80) = mDisp_b_n(1);
  data(81) = mDisp_b_n(2);
  data(82) = mDisp_b_n(3);
  data(83) = mDisp_b_n(4);
  data(84) = mDisp_b_n(5);
  data(85) = mDisp_s_n(0);
  data(86) = mDisp_s_n(1);
  data(87) = mDisp_s_n(2);

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }          

  // BeamContact3Dp then sends the tags of it's four nodes
  res = theChannel.sendID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }


  // finally BeamContact3Dp asks it's material object to send itself
  res = crdTransf->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::sendSelf() - " << this->getTag() << " failed to send its crdTransf\n";
    return -4;
  }

  // finally BeamContact3Dp asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }
  
  return 0;
}

int
BeamContact3Dp::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();
  
  // BeamContact3Dp creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(88);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::recvSelf() - failed to receive Vector\n";
    return -1;
  }          
  
  this->setTag((int)data(0));
  mRadius       = data(1);
  mPenalty      = data(2);
  mIniContact   = (int)data(3);
  inContact     = (int)data(8);
  was_inContact = (int)data(9);
  in_bounds     = (int)data(10);
  mGap          = data(11);
  mLambda       = data(12);
  mxi           = data(13);
  mchi          = data(14);
  mL            = data(15);
  mrho2         = data(16);
  mrho3         = data(17);

  mH(0)         = data(18);
  mH(1)         = data(19);
  mH(2)         = data(20);
  mH(3)         = data(21);

  mn(0)         = data(22);
  mn(1)         = data(23);
  mn(2)         = data(24);

  mg1(0)        = data(25);
  mg1(1)        = data(26);
  mg1(2)        = data(27);
  mg2(0)        = data(28);
  mg2(1)        = data(29);
  mg2(2)        = data(30);

  mQa(0,0)      = data(31);
  mQa(0,1)      = data(32);
  mQa(0,2)      = data(33);
  mQa(1,0)      = data(34);
  mQa(1,1)      = data(35);
  mQa(1,2)      = data(36);
  mQa(2,0)      = data(37);
  mQa(2,1)      = data(38);
  mQa(2,2)      = data(39);

  mQb(0,0)      = data(40);
  mQb(0,1)      = data(41);
  mQb(0,2)      = data(42);
  mQb(1,0)      = data(43);
  mQb(1,1)      = data(44);
  mQb(1,2)      = data(45);
  mQb(2,0)      = data(46);
  mQb(2,1)      = data(47);
  mQb(2,2)      = data(48);

  mQc(0,0)      = data(49);
  mQc(0,1)      = data(50);
  mQc(0,2)      = data(51);
  mQc(1,0)      = data(52);
  mQc(1,1)      = data(53);
  mQc(1,2)      = data(54);
  mQc(2,0)      = data(55);
  mQc(2,1)      = data(56);
  mQc(2,2)      = data(57);

  mIcrd_a(0)    = data(58);
  mIcrd_a(1)    = data(59);
  mIcrd_a(2)    = data(60);
  mIcrd_b(0)    = data(61);
  mIcrd_b(1)    = data(62);
  mIcrd_b(2)    = data(63);
  mIcrd_s(0)    = data(64);
  mIcrd_s(1)    = data(65);
  mIcrd_s(2)    = data(66);

  mDcrd_a(0)    = data(67);
  mDcrd_a(1)    = data(68);
  mDcrd_a(2)    = data(69);
  mDcrd_b(0)    = data(70);
  mDcrd_b(1)    = data(71);
  mDcrd_b(2)    = data(72);

  mDisp_a_n(0)  = data(73);
  mDisp_a_n(1)  = data(74);
  mDisp_a_n(2)  = data(75);
  mDisp_a_n(3)  = data(76);
  mDisp_a_n(4)  = data(77);
  mDisp_a_n(5)  = data(78);
  mDisp_b_n(0)  = data(79);
  mDisp_b_n(1)  = data(80);
  mDisp_b_n(2)  = data(81);
  mDisp_b_n(3)  = data(82);
  mDisp_b_n(4)  = data(83);
  mDisp_b_n(5)  = data(84);
  mDisp_s_n(0)  = data(85);
  mDisp_s_n(1)  = data(86);
  mDisp_s_n(2)  = data(87);

  // BeamContact3Dp now receives the tags of its external nodes
  res = theChannel.recvID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }
  
  int crdClass = (int)data(6);
  int crdDb = (int)data(7);
  
  // check if we have a material object already & if we do if of right type
  if ((crdTransf == 0) || (crdTransf->getClassTag() != crdClass)) {
    
    // if old one .. delete it
    if (crdTransf != 0)
      delete crdTransf;
    
    // create a new material object
    crdTransf = theBroker.getNewCrdTransf(crdClass);
    
    if (crdTransf == 0) {
      opserr <<"WARNING BeamContact3Dp::recvSelf() - " << this->getTag()
	     << " failed to get a blank CrdTransf of type " << crdClass << endln;
      return -3;
    }
  }
  
  crdTransf->setDbTag(crdDb); // note: we set the dbTag before we receive the material
  
  res = crdTransf->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }
  
  // finally BeamContact3Dp creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.
  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {
    
    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;
    
    // create a new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (ContactMaterial3D *)theMatCopy;
    
    if (theMaterial == 0) {
      opserr <<"WARNING BeamContact3Dp::recvSelf() - " << this->getTag()
	     << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }
  
  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING BeamContact3Dp::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }
  
  return 0;
}


int
BeamContact3Dp::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}


void
BeamContact3Dp::Print(OPS_Stream &s, int flag)
{
    opserr << "BeamContact3Dp, element id:  " << this->getTag() << endln;
    opserr << "   Connected external nodes:  " << externalNodes;
    opserr << "   Transformation: ";
    if (crdTransf != 0)
      crdTransf->Print(s,flag);
    opserr << "\n    Material: ";
    if (theMaterial != 0) 
      theMaterial->Print(s,flag);
    opserr << "\n";

    return;
}

Response*
BeamContact3Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
        return new ElementResponse(this, 1, Vector(3));

    } else if (strcmp(argv[0],"frictionforce") == 0 || strcmp(argv[0],"frictionforces") == 0) {
        return new ElementResponse(this, 2, Vector(3));
   
    } else if (strcmp(argv[0],"forcescalar") == 0 || strcmp(argv[0],"forcescalars") == 0) {
        return new ElementResponse(this, 3, Vector(3));

    } else if (strcmp(argv[0],"masterforce") == 0 || strcmp(argv[0],"masterforces") == 0 ||
	       strcmp(argv[0],"primaryforce") == 0 || strcmp(argv[0],"primaryforces") == 0) {
        return new ElementResponse(this, 4, Vector(6));

    } else if (strcmp(argv[0],"mastermoment") == 0 || strcmp(argv[0],"mastermoments") == 0 ||
	       strcmp(argv[0],"primarymoment") == 0 || strcmp(argv[0],"primarymoments") == 0) {
        return new ElementResponse(this, 5, Vector(6));

    } else if (strcmp(argv[0],"masterreaction") == 0 || strcmp(argv[0],"masterreactions") == 0 ||
	       strcmp(argv[0],"primaryreaction") == 0 || strcmp(argv[0],"primaryreactions") == 0) {
        return new ElementResponse(this, 6, Vector(12));

	} else if (strcmp(argv[0],"slip") == 0) {
	    return new ElementResponse(this, 7, Vector(2));

	} else {
        // otherwise response quantity is unknown for the BeamContact3Dp class
        opserr << "BeamContact3Dp::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): " << argv[0] << " unknown request" << endln;
        return 0;
	}
}

int
BeamContact3Dp::getResponse(int responseID, Information &eleInfo)
{
    Vector force(3);
    Vector sForce(3);
    Vector mForces(6);
    Vector mMoments(6);
    Vector mReactions(12);
	Vector theSlip(2);
     
    // get contact stresses/forces
    Vector stress = theMaterial->getStress();

    if (responseID == 1) {
     
        // forces on secondary node
        for (int ii=0; ii<3; ii++) {
            sForce(ii)   = -mInternalForces(BC3Dp_NUM_DOF - 6 + ii);
        }
        return eleInfo.setVector(sForce);
 
    } else if (responseID == 2) {
         
        force = stress(1)*mg1 + stress(2)*mg2;
        return eleInfo.setVector(force);

    } else if (responseID == 3) {

        force(0) = stress(0);
        force(1) = stress(1);
        force(2) = stress(2);
        return eleInfo.setVector(force);
 
    } else if (responseID == 4) {

        // forces on primary nodes
        for (int ii=0; ii<3; ii++) {
            mForces(ii)   = -mInternalForces(ii);
            mForces(ii+3) = -mInternalForces(ii+6);
        }
        return eleInfo.setVector(mForces);
 
    } else if (responseID == 5) {

        // moments on primary nodes
        for (int ii=0; ii<3; ii++) {
            mMoments(ii)   = -mInternalForces(ii+3);
            mMoments(ii+3) = -mInternalForces(ii+9);
        }
        return eleInfo.setVector(mMoments);
 
    } else if (responseID == 6) {

        // full reactions on primary nodes
        for (int ii=0; ii<6; ii++) {
            mReactions(ii)   = -mInternalForces(ii);
            mReactions(ii+6) = -mInternalForces(ii+6);
        }
        return eleInfo.setVector(mReactions);
  
    } else if (responseID == 7) {

	    // slip vector on tangent plane
	    return eleInfo.setVector(mSlip);
 
    } else {
        // responsibility is unknown for the BeamContact3Dp class
        opserr << "BeamContact3Dp::getResponse(int responseID=" << responseID << ", Information &eleInfo): " << " unknown request" << endln;
        return -1;
	}
}

int
BeamContact3Dp::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"friction") == 0) {
		return param.addObject(1, this);
	}
	
	return -1;
}

int
BeamContact3Dp::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes =  theMaterial->updateParameter(parameterID, info);
	if (matRes != -1) {
		res = matRes;
	}
	return res;
}
