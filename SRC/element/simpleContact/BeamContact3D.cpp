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

// Written: Kathryn A. Petek
// Created: 06/05
//
// Revisions
//    06/05 created
//    10/07 Peter Mackenzie-Helnwein: adding recorder features
//    11/10 F.Mckenna and C.McGann: changes for incorporation into main source code
//    02/11 C.McGann: added initial contact switch (default is inContact)
//
// Description: This file contains the implementation for the BeamContact3D class.

#include "BeamContact3D.h"
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#include <elementAPI.h>
#define OPS_Export

static int num_BeamContact3D = 0;

OPS_Export void *
OPS_BeamContact3D(void)
{
  if (num_BeamContact3D == 0) {
    num_BeamContact3D++;
    OPS_Error("BeamContact3D element - Written: K.Petek, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 10) {
    opserr << "Invalid #args,  want: element BeamContact3D eleTag?  iNode? jNode? slaveNode? lambdaNode? radius? crdTransf? matTag? tolGap? tolF? <cSwitch>?\n";
    return 0;
  }
   
  int    iData[7];
  double dData[3];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0] << endln;
    return 0;  
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement" << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0] << endln;
    return 0;  
  }

  int transfTag = iData[5];
  CrdTransf *theTransf = OPS_GetCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[6];
  NDMaterial *theMaterial = OPS_GetNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 10;
  while (numRemainingInputArgs >= 1) {
	  numData = 1;
	  if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
		  opserr << "WARNING invalid initial contact flag: element BeamContact3D " << iData[0] << endln;
	  	  return 0;
      }
	  numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact3D(iData[0], iData[1], iData[2], iData[3], iData[4],
                                 dData[0], *theTransf, *theMaterial,
                                 dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact3DElement\n";
    return 0;
  }

  return theElement;
}



// constructors:
BeamContact3D::BeamContact3D(int tag, int Nd1, int Nd2,
                             int NdS, int NdL, double rad, CrdTransf &coordTransf,
                             NDMaterial &theMat, double tolG, double tolF, int cSwitch)
 :Element(tag,ELE_TAG_BeamContact3D),    
   crdTransf(0),
   theMaterial(0),
   externalNodes(BC3D_NUM_NODE),
   mTangentStiffness(BC3D_NUM_DOF, BC3D_NUM_DOF),
   mInternalForces(BC3D_NUM_DOF),
   meye1(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mg1(BC3D_NUM_NDM),
   mg2(BC3D_NUM_NDM),
   mg_metric(2,2),
//   mG_metric(2,2),
   mn(BC3D_NUM_NDM),
   mH(4),
   mIcrd_a(BC3D_NUM_NDM),
   mIcrd_b(BC3D_NUM_NDM),
   mIcrd_s(BC3D_NUM_NDM),
   mDcrd_a(BC3D_NUM_NDM),
   mDcrd_b(BC3D_NUM_NDM),
   mDcrd_s(BC3D_NUM_NDM),
   //mRot_a_n(BC3D_NUM_NDM),
  // mRot_b_n(BC3D_NUM_NDM),
   mDisp_a_n(6),
   mDisp_b_n(6),
   mDisp_s_n(3),
   mQa(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mQb(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mQc(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mc1(BC3D_NUM_NDM),
   mBn(BC3D_NUM_DOF-3),
   mBs(BC3D_NUM_DOF-3,2),
   mBphi(3,12),
   mIniContact(cSwitch)
{
#ifdef DEBUG
        opserr << "BeamContact3D::BeamContact3D(): " << MyTag << endln;
#endif
        externalNodes(0) = Nd1;
        externalNodes(1) = Nd2;
        externalNodes(2) = NdS;
        externalNodes(3) = NdL;

        mRadius = rad;
        mTolGap = tolG;
        mTolForce = tolF;
		mIniContact = cSwitch;

		if (mIniContact == 0) {
			inContact          = true;
			was_inContact      = true;
			to_be_released     = false;
			should_be_released = false;
			in_bounds          = true;
		} else {
			inContact          = false;
			was_inContact      = false;
			to_be_released     = false;
			should_be_released = false;
			in_bounds          = true;
		}

        mGap    = 0.0;
        mLambda = 0.0;
       
        // get copy of the transformation & material object  
        crdTransf = coordTransf.getCopy3d();
        NDMaterial *theMatCopy = theMat.getCopy("ContactMaterial3D");
        if (theMatCopy != 0) {
          theMaterial = (ContactMaterial3D *)theMatCopy;
        } else {
          opserr << "BeamContact3D::BeamContact3D - material needs to be of type Contact3D for ele: " << this->getTag() << endln;
        }      
       
        // check them:         
        if (!crdTransf)
          {
            opserr << "Error: BeamContact3d::BeamContact3d: could not create copy of coordinate transformation object" << endln;
            exit(-1);
          }
       
        if (theMaterial == 0) {
          opserr << "BeamContact3D::BeamContact3D - failed allocate material model pointer\n";
          exit(-1);
        }
       
        MyTag = tag;
}

BeamContact3D::BeamContact3D()
 :Element(0,ELE_TAG_BeamContact3D),    
   crdTransf(0),
   theMaterial(0),
   externalNodes(BC3D_NUM_NODE),
   mTangentStiffness(BC3D_NUM_DOF, BC3D_NUM_DOF),
   mInternalForces(BC3D_NUM_DOF),
   meye1(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mg1(BC3D_NUM_NDM),
   mg2(BC3D_NUM_NDM),
   mg_metric(2,2),
//   mG_metric(2,2),
   mn(BC3D_NUM_NDM),
   mH(4),
   mIcrd_a(BC3D_NUM_NDM),
   mIcrd_b(BC3D_NUM_NDM),
   mIcrd_s(BC3D_NUM_NDM),
   mDcrd_a(BC3D_NUM_NDM),
   mDcrd_b(BC3D_NUM_NDM),
   mDcrd_s(BC3D_NUM_NDM),
   //mRot_a_n(BC3D_NUM_NDM),
  // mRot_b_n(BC3D_NUM_NDM),
   mDisp_a_n(6),
   mDisp_b_n(6),
   mDisp_s_n(3),
   mQa(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mQb(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mQc(BC3D_NUM_NDM,BC3D_NUM_NDM),
   mc1(BC3D_NUM_NDM),
   mBn(BC3D_NUM_DOF-3),
   mBs(BC3D_NUM_DOF-3,2),
   mBphi(3,12)
{
#ifdef DEBUG
        opserr << "BeamContact3D::BeamContact3D(): " << MyTag << endln;
#endif
}


//  destructor:
BeamContact3D::~BeamContact3D()
{
#ifdef DEBUG
        opserr << "BeamContact3D::~BeamContact3D(): " << MyTag << endln;
#endif
        if (theMaterial != 0)
          delete theMaterial;
       
        delete crdTransf;
}


int
BeamContact3D::getNumExternalNodes(void) const
{
#ifdef DEBUG
        opserr << "BeamContact3D::getNumExternalNodes(): " << MyTag << endln;
#endif
    return BC3D_NUM_NODE;
}

const ID &
BeamContact3D::getExternalNodes(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getExternalNodes(): " << MyTag << endln;
#endif
    return externalNodes;
}


Node **
BeamContact3D::getNodePtrs(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getNodePtrs(): " << MyTag << endln;
#endif
        return theNodes;                        
}

int
BeamContact3D::getNumDOF(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getNumDOF(): " << MyTag << endln;
#endif
    return BC3D_NUM_DOF;
}


void
BeamContact3D::setDomain(Domain *theDomain)
{
#ifdef DEBUG
        opserr << "BeamContact3D::setDomain(Domain *theDomain): " << MyTag << endln;
#endif

                Vector a1(BC3D_NUM_NDM);       // tangent vector at end a
                Vector b1(BC3D_NUM_NDM);       // tangent vector at end b
                Vector x_c(BC3D_NUM_NDM);      // coordinate at projected point on centerline

                meye1.Zero();
                meye1(0,0) = 1.0;
                meye1(1,1) = 1.0;
                meye1(2,2) = 1.0;

        int Nd1 = externalNodes(0);
        int Nd2 = externalNodes(1);
        int NdS = externalNodes(2);
        int NdL = externalNodes(3);

        theNodes[0] = theDomain->getNode(Nd1);
        theNodes[1] = theDomain->getNode(Nd2);
        theNodes[2] = theDomain->getNode(NdS);
        theNodes[3] = theDomain->getNode(NdL);
               
                int i;
        for (i = 0; i < 4; i++) {
           if (theNodes[i] == 0)
                return;  // don't go any further - otherwise segmentation fault
        }  

                //  initialize coordinate vectors and set to initial state
                mIcrd_a = theNodes[0]->getCrds();
        mIcrd_b = theNodes[1]->getCrds();
                mIcrd_s = theNodes[2]->getCrds();    
                mDcrd_a = mIcrd_a;
                mDcrd_b = mIcrd_b;
                mDcrd_s = mIcrd_s;
                //mRot_a_n.Zero();
                //mRot_b_n.Zero();
                mDisp_a_n.Zero();
                mDisp_b_n.Zero();
                mDisp_s_n.Zero();

                // Fill the coordinate transformation matrices mQa, mQb
                // initialize the transformation
                if (crdTransf->initialize(theNodes[0], theNodes[1]))
                {
                        opserr << "BeamContact3D::setDomain(): Error initializing coordinate transformation";  
                        exit(0);
                }
                // Note: the following intialization is already performed in crdTransf
                //       but it doesnt return the values to beamContact element
                Vector initXAxis(3);
                Vector initYAxis(3);
                Vector initZAxis(3);
                crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
                // fill mQa
                for (i = 0; i < 3; i++) {
                        mQa(i,0) = initXAxis(i);
                        mQa(i,1) = initYAxis(i);
                        mQa(i,2) = initZAxis(i);
                }
                // set mQb = mQa : beam column element requires zero initial twist
                // if mQa = mQb -> mchi = 0
                mQb = mQa;
                mchi = 0;  
                // fill tangent vectors = first column mQa, mQb = initXAxis
                a1 = initXAxis;
                b1 = a1;

                // length of master segment L
                mL = (mDcrd_b - mDcrd_a).Norm();  
       
                // perform projection to update local coordinate along centerline
                //  of master segment.  projection function also sets mn, mc1
                mxi = ((mDcrd_b - mDcrd_s)^(mDcrd_b - mDcrd_a))
                          / ((mDcrd_b - mDcrd_a)^(mDcrd_b - mDcrd_a)) ;  // initial approx
                // initial basis function values for use in projection
        mxi = project(mxi);

                // initialize contact state based on projection
                // in_bounds      = ( (mxi>=0.000) && (mxi<=1.0000) );
                in_bounds      = ( (mxi>0.000) && (mxi<1.0000) );
        inContact      = ( was_inContact && in_bounds );

                // compute centerline projection coordinate
                x_c = mDcrd_a*mH(0) + a1*mH(1) + mDcrd_b*mH(2) + b1*mH(3);

                // update base vectors g1, g2 for determination of mg_metric
                UpdateBase(mxi);
       
        // adjust cohesion force
        //double area = 1.0*sqrt(mg_metric(0,0)*mg_metric(1,1)-mg_metric(0,1)*mg_metric(1,0));
        double area = 1.0;
        theMaterial->ScaleCohesion(area);
        theMaterial->ScaleTensileStrength(area);
       
        // compute mBn, mBs
        ComputeB();

        // call the base class method
        this->DomainComponent::setDomain(theDomain);

}        


int
BeamContact3D::commitState()
{
#ifdef DEBUG
        opserr << "BeamContact3D::commitState(): " << MyTag << endln;
#endif


                // update projection
                // (includes update of Qa, Qb & recalc of mn, mc1)
                mxi = project(mxi);

                // update tangent plane & metric tensor
                // (includes update of Qc, rho2, rho3, g1, g2)
                UpdateBase(mxi);

                // update Bn, Bs for next step
                ComputeB();

                // update Boolean Variables for contact condition
        was_inContact  = ( mGap < mTolGap );   
                // in_bounds      = ( (mxi>=-0.0001) && (mxi<=1.0001) );
                in_bounds      = ( (mxi>0.000) && (mxi<1.000) );
                to_be_released = ( should_be_released || !in_bounds );
        inContact      = ( was_inContact && !to_be_released && in_bounds );

/*              double tol = 1e-10;
                if ( (mQa(0,0)-mQc(0,0)>tol) || (mQa(0,1)-mQc(0,1)>tol) ||(mQa(0,2)-mQc(0,2)>tol) ||
                         (mQa(1,0)-mQc(1,0)>tol) || (mQa(1,1)-mQc(1,1)>tol) ||(mQa(1,2)-mQc(1,2)>tol) ||
                         (mQa(2,0)-mQc(2,0)>tol) || (mQa(2,1)-mQc(2,1)>tol) ||(mQa(2,2)-mQc(2,2)>tol) ) {
                        opserr << "ELE: " << MyTag << endln;
                        opserr << "*** mQa neq. mQc *** "<< endln;
                        opserr << "mQa = "  << mQa << endln;
                        opserr << "mQc = " << mQc << endln;
                }

                if ( (mQa(0,0)-mQb(0,0)>tol) || (mQa(0,1)-mQb(0,1)>tol) ||(mQa(0,2)-mQb(0,2)>tol) ||
                         (mQa(1,0)-mQb(1,0)>tol) || (mQa(1,1)-mQb(1,1)>tol) ||(mQa(1,2)-mQb(1,2)>tol) ||
                         (mQa(2,0)-mQb(2,0)>tol) || (mQa(2,1)-mQb(2,1)>tol) ||(mQa(2,2)-mQb(2,2)>tol) ) {
                        opserr << "ELE: " << MyTag << endln;
                        opserr << "*** mQa neq. mQb *** "<< endln;
                        opserr << "mQa = "  << mQa << endln;
                        opserr << "mQb = " << mQb << endln;
                }
                */

/*              Vector a1(3);
                Vector b1(3);
                Vector c2(3);
                Vector c3(3);
                Vector x_c_G(3);
                double check1;
                double check2;
               
                a1 = Geta1();
                b1 = Getb1();
                for (int i=0; i<3; i++) {
                                c2(i) = mQc(i,1);
                                c3(i) = mQc(i,2);
                }

                x_c_G = mDcrd_a*mH(0) + a1*mH(1) + mDcrd_b*mH(2) + b1*mH(3) + mrho2*c2 + mrho3*c3;

                check1 = mg1 ^ (mDcrd_s - x_c_G);
                check2 = mg2 ^ (mDcrd_s - x_c_G);

                opserr << "ELE = " << MyTag << endln;
                opserr << "mDcrd_s = [ " << mDcrd_s << "];" << endln;
                opserr << "x_c    = [ " << (mDcrd_a*mH(0) + a1*mH(1) + mDcrd_b*mH(2) + b1*mH(3)) << "];" << endln;
                opserr << "r = [ " << (mrho2*c2 + mrho3*c3) << "];" << endln;
                opserr << "x_c_G  = [ " << x_c_G << "];" << endln;
        opserr << "(x_s - x_c_G) = [" << (mDcrd_s - x_c_G) << " ];" << endln;
                opserr << "check 1 = "  << check1 << ";      check 2 = " << check2 << ";" << endln;
                opserr << endln; */


        int retVal = 0;
        // call element commitState to do any base class stuff
        if ((retVal = this->Element::commitState()) != 0) {
                opserr << "BeamContact3D::commitState () - failed in base class";
                }    
        retVal = theMaterial->commitState();

                return retVal;

}


int
BeamContact3D::revertToLastCommit()
{
#ifdef DEBUG
        opserr << "BeamContact3D::revertToLastCommit(): " << MyTag << endln;
#endif
        return theMaterial->revertToLastCommit();
}

int
BeamContact3D::revertToStart()
{
	if (mIniContact == 0) {
		inContact          = true;
		was_inContact      = true;
		to_be_released     = false;
		should_be_released = false;
		in_bounds          = true;
	} else {
		inContact          = false;
		was_inContact      = false;
		to_be_released     = false;
		should_be_released = false;
		in_bounds          = true;
	}

    return theMaterial->revertToStart();
}

int
BeamContact3D::update(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::update(): " << MyTag << endln;
#endif

                double tensileStrength;
                Vector a1(BC3D_NUM_NDM);
                Vector b1(BC3D_NUM_NDM);
                Vector a1_n(BC3D_NUM_NDM);
                Vector b1_n(BC3D_NUM_NDM);
                Vector disp_a(6);
                Vector disp_b(6);
                Vector disp_s(3);
                Vector disp_L(BC3D_NUM_NDM);
                Vector rot_a(BC3D_NUM_NDM);
                Vector rot_b(BC3D_NUM_NDM);
                Vector x_c(BC3D_NUM_NDM);
                Vector d(BC3D_NUM_NDM);

                // update slave node coordinate
                mDcrd_s = mIcrd_s + theNodes[2]->getTrialDisp();

                // update Lagrange Multiplier Value
                disp_L = theNodes[3]->getTrialDisp();
                mLambda = disp_L(0);

                // update nodes a, b  coordinates & rotations
                disp_a = theNodes[0]->getTrialDisp();
                disp_b = theNodes[1]->getTrialDisp();
                disp_s = theNodes[2]->getTrialDisp();

                int i;
                for (i=0; i<3; i++) {
                        // compute updated nodal coordinates a, b
                        mDcrd_a(i) = mIcrd_a(i) +  disp_a(i);
                        mDcrd_b(i) = mIcrd_b(i) +  disp_b(i);
                        // compute incremental rotations from step n to n+1
                        //rot_a(i) = disp_a(i+3) - mRot_a_n(i);
                        //rot_b(i) = disp_b(i+3) - mRot_b_n(i);
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

                // update value of the gap
                d = mDcrd_s - x_c;
                mGap = (mn^d) - mRadius;

                // get tensile strength from material
                tensileStrength = theMaterial->getTensileStrength();

                // set boolean variable for release condition
                //should_be_released = ( mLambda <= -(tensileStrength + mTolForce ) );
                should_be_released = ( mLambda <= - mTolForce );

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
            opserr << "   CONTACT:            " << inContact << endln;
            opserr << "   should be released: " << should_be_released << endln;
            }
#endif

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

                        for (i=0; i<3; i++) {
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


Vector dummySlip(2);
dummySlip.Zero();
for (i = 0; i < 12; i++) {
        dummySlip(0) += incDisp_ab(i) * mBs(i,0);
        dummySlip(1) += incDisp_ab(i) * mBs(i,1);
}

for (i = 0; i < 3; i++) {
        dummySlip(0) += incDisp_s(i) * mBs(12+i,0);
        dummySlip(1) += incDisp_s(i) * mBs(12+i,1);
}



if ( (slip - dummySlip).Norm() > 1e-4 ) {
#ifdef DEBUG
        //opserr << "      slip = " << slip;
        //opserr << " dummySlip = " << dummySlip;
        //opserr << " ELE: " << MyTag << "      Norm = " << (slip - dummySlip).Norm()
        //               << "      DIFF = " << (slip - dummySlip)  << endln;
#endif
}

//slip(0) = -dummySlip(0);
//slip(1) = -dummySlip(1);
 
/*
opserr << "ELE = " << MyTag << endln;
opserr << "mxi = " << mxi << endln;
opserr << "mL = "  << mL  << ";" << endln;
opserr << "mQa = [ " << mQa << " ];" << endln;
opserr << "mQb = [ " << mQb << " ];" << endln;
opserr << "mQc = [ " << mQc << " ];" << endln;
opserr << "mDcrd_s = [ " << mDcrd_s << " ]';" << endln;
opserr << "x_c  = [ " << x_c  << " ]';" endln;
opserr << "disp_a = [ " << disp_a << " ]';" << endln;
opserr << "disp_b = [ " << disp_b << " ]';" << endln;
opserr << "mDisp_a_n = [ " << mDisp_a_n << " ]';" << endln;
opserr << "mDisp_b_n = [ " << mDisp_b_n << " ]';" << endln;
opserr << "incDisp_ab = [ "  << incDisp_ab << " ]'" << endln;
opserr << "mrho2 = " << mrho2 << ";" << endln;
opserr << "mrho3 = " << mrho3 << ";" << endln;
opserr << "mn = [ " << mn << "]';" << endln;
opserr << "mg1 = [ " << mg1 << " ]';" << endln;
opserr << "mg2 = [ " << mg2 << " ]';" << endln;
opserr << "phi_c = [ " << phi_c << " ]'" << endln;
opserr << "c2n1 = [ " << c2n1 << "]'" << endln;
opserr << "c3n1 = [ " << c3n1 << "]'" << endln;
opserr << "dstar= [ " << dstar << " ]'" << endln;
opserr << "mBphi = [ " << mBphi << " ]" << endln;
opserr << "mBs = [ " << mBs << " ]" << endln;
//opserr << "mBn = [ " << mBn << " ]" << endln;
opserr << "slip = [ " << slip << " ]' " << endln;
opserr << "dummySlip = [ " << dummySlip << " ]' " <<endln;
opserr << endln; */


            // slip is treated in local (curvilinear) coordinates
            strain(0) = mGap;
            strain(1) = slip(0);
            strain(2) = slip(1);
            strain(3) = mLambda;      
            theMaterial->setTrialStrain(strain);

        }

        else if (to_be_released) {
                // prevents sliding & stabilizes behavior in lag step
                        Vector strain(4);

            strain(0) = mGap;
            strain(1) = 0.0;
            strain(2) = 0.0;
            strain(3) = mLambda;    
            theMaterial->setTrialStrain(strain);
        }



  return 0;
}

double
BeamContact3D::project(double xi)
{
#ifdef DEBUG
        opserr << "BeamContact3D::project(): " << MyTag << endln;
#endif

                double xi_P;                                    // local value of xi
                double H1;                                              // local value of Hermitian Function, H1
                double H2;                                              // local value of Hermitian Function, H2
                double H3;                                              // local value of Hermitian Function, H3
                double H4;                                              // local value of Hermitian Function, H4
                double R;
                double DR;
                double dxi;                                             // change in xi
                double xi_P_squared;                   
                double xi_P_cubed;                             
                Vector a1(BC3D_NUM_NDM);                // tangent at end a
                Vector b1(BC3D_NUM_NDM);                // tangent at end b
                Vector x_c_P(BC3D_NUM_NDM);             // current centerline porjection coordinate
                Vector d(BC3D_NUM_NDM);                 // distance from slave node to centerline coord
                Vector tc(BC3D_NUM_NDM);                // tangent at projection point = 1st deriv of x_c
                Vector ddx_c(BC3D_NUM_NDM);             // 2nd derivative of x_c

                // initialize xi_P to previous value of xi
                xi_P = xi;

                // Perform exponential update of coordinate tranforms Qa, Qb
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
                while (fabs(R/mL) > mTolGap && Gapcount < 50) {
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

                // set value of c1 for current proejction for use in ComputeQc
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
BeamContact3D::UpdateBase(double xi)
// this function calculates g1, g2, mg_metric, and sends value of metric tensor to the material
{
#ifdef DEBUG
        opserr << "BeamContact3D::UpdateBase(): " << MyTag << endln;
#endif

                //double rho2;                          // angular coord relating c2 to radial vector
                //double rho3;                          // angular coord relating c3 to radial vector
                Vector c1(BC3D_NUM_NDM);       
                Vector c2(BC3D_NUM_NDM);        // c1, c2, c3 are coord transf. vectors at projected point
                Vector c3(BC3D_NUM_NDM);
                Vector tc(BC3D_NUM_NDM);    // tangent at projected point = 1st deriv of x_c
                Vector ddx_c(BC3D_NUM_NDM); // 2nd deriv of x_c w.r.t. xi
                Vector d_c1(BC3D_NUM_NDM);  // deriv of c1 w.r.t xi
                Vector d_c2(BC3D_NUM_NDM);  // deriv of c2 w.r.t xi
                Vector d_c3(BC3D_NUM_NDM);  // deriv of c3 w.r.t xi
                //Vector g1(BC3D_NUM_NDM);      // tangent plane basis vector, g_xi
                //Vector g2(BC3D_NUM_NDM);    // tangent plane basis vector, g_psi
                Matrix Qc(BC3D_NUM_NDM, BC3D_NUM_NDM);  // coord transf at projected point
                Vector g1(BC3D_NUM_NDM);    // temporary vector for g1 (co-variant)
                Vector g2(BC3D_NUM_NDM);    // temporary vector for g2 (co-variant)
                Vector test(3);

                // Update interpolated coordinate transform, Qc
                ComputeQc(xi);  // sets new value of mQc
                Qc = mQc;

                // obtain c1, c2, c3 vectors
                c1 = Getc1();
                int i;
                for (i=0; i<3; i++) {
                        // c1(i) = Qc(i,0);  if //: c1 = tc/norm(tc);
                        test(i) = c1(i) - Qc(i,0);
                        c2(i) = Qc(i,1);
                        c3(i) = Qc(i,2);
                }

//opserr << "ELE : " << MyTag << "     test.norm() = " << test.Norm() << ";   test = [" << test << "];" << endln;

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
                // g1 = d(x_c + r_vec)/dxi
                // g2 = d(x_c + r_vec)/dpsi
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

/*              opserr << "ELE = " << MyTag << endln;
                opserr << "g1 = [ " << g1 << "]" << endln;
                opserr << "g2 = [ " << g2 << "]" << endln;
                opserr << "mg1 = [ " << mg1 << "]" << endln;
                opserr << "mg2 = [ " << mg2 << "]" << endln; */

                return 0;
}


// Bn and Bs are computed at CommitState for the next step
// this allows to only be calculated once
/*void
BeamContact3D::ComputeB(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::ComputeB(): " << MyTag << endln;
#endif

        // vector n used below should be of same step

                // initialize Bn, Bs;
                mBn.Zero();
                mBs.Zero();

                int i, j;

                Vector a1(3);
                Vector b1(3);
                Vector a1xn(3);
                Vector b1xn(3);
                Vector r(3);
                Vector rxg1(3);
                Vector rxg2(3);

                Matrix Bx(3,12);

                a1= Geta1();
                b1 = Getb1();
                a1xn = CrossProduct(a1, mn);
                b1xn = CrossProduct(b1, mn);  

                r(0) = mrho2*mQc(0,1) + mrho3*mQc(0,2);
                r(1) = mrho2*mQc(1,1) + mrho3*mQc(1,2);
                r(2) = mrho2*mQc(2,1) + mrho3*mQc(2,2);

        // dgap = Bn' * dq
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
               
                // Compute Bs in components
                // for al = 1, 2
                // Bs = [ -(g_n^al'*Bx + (r_n x g_n^al)'*Bphi)       g_n^al'   0'  ]
                Bx.Zero();
                Bx(0,0) =  mH(0);
                Bx(1,1) =  mH(0);
                Bx(2,2) =  mH(0);
                Bx(1,3) = -mH(1)*mL*a1(2);
                Bx(2,3) =  mH(1)*mL*a1(1);
                Bx(0,4) =  mH(1)*mL*a1(2);
                Bx(2,4) = -mH(1)*mL*a1(0);
                Bx(0,5) = -mH(1)*mL*a1(1);
                Bx(1,5) =  mH(1)*mL*a1(0);
                Bx(0,6) =  mH(2);
                Bx(1,7) =  mH(2);
                Bx(2,8) =  mH(2);
                Bx(1,9)  = -mH(3)*mL*b1(2);
                Bx(2,9)  =  mH(3)*mL*b1(1);
                Bx(0,10) =  mH(3)*mL*b1(2);
                Bx(2,10) = -mH(3)*mL*b1(0);
                Bx(0,11) = -mH(3)*mL*b1(1);
                Bx(1,11) =  mH(3)*mL*b1(0);


                // Compute Bphi using function
                mBphi = ComputeBphi();
               
                // Assemble mBs
                rxg1 = CrossProduct(r,mg1);
                rxg2 = CrossProduct(r,mg2);


                for (i = 0; i < 12; i++) {
                        for (j = 0; j<3; j++) {
                                mBs(i,0) -= mg1(j) * Bx(j,i) + rxg1(j) * mBphi(j,i);
                                mBs(i,1) -= mg2(j) * Bx(j,i) + rxg2(j) * mBphi(j,i);
                        }
                }

                for (j = 0; j < 3; j++) {
                        mBs(12+j,0) = mg1(j);
                        mBs(12+j,1) = mg2(j);
                }

        return;
} */



void
BeamContact3D::ComputeB(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::ComputeB(): " << MyTag << endln;
#endif

        // vector n used below should be of same step

                // initialize Bn, Bs;
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
                //rxg1 = -1*rxg1;
                //rxg2 = -1*rxg2;

        // dgap = Bn' * dq
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
BeamContact3D::ComputeBphi(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::ComputeB(): " << MyTag << endln;
#endif
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
//              dummy1(0,0) = mxi*mQc(0,0);  // N1 * mQc(0:2,0)
//              dummy1(1,0) = mxi*mQc(1,0);
//              dummy1(2,0) = mxi*mQc(2,0);

                dummy1(0,1) = dH2*mQc(0,1);  // dH2 * mQc(0:2,1:2)
                dummy1(1,1) = dH2*mQc(1,1);
                dummy1(2,1) = dH2*mQc(2,1);
                dummy1(0,2) = dH2*mQc(0,2);
                dummy1(1,2) = dH2*mQc(1,2);
                dummy1(2,2) = dH2*mQc(2,2);
                dummy2 = Transpose(BC3D_NUM_NDM, BC3D_NUM_NDM, mQa);
                dummy3 = dummy1*dummy2;
                for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                                Bphi(i,3+j) = dummy3(i,j);
                                //Bphi(i,3+j) = dummy3(j,i);
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
                                //Bphi(i,j) = dH1/mL * dummy3(j,i);
                        }
                }

                // Reuse dummy1 and Compute Bphi(0:2, 6:8)
                dummy2.Zero();
                dummy3.Zero();
                dummy2 = Transpose(BC3D_NUM_NDM, BC3D_NUM_NDM, mQb);
                dummy3 = dummy1*dummy2;
                for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                                Bphi(i,6+j) = dH3/mL * dummy3(i,j);
                                //Bphi(i,6+j) = dH3/mL * dummy3(j,i);
                        }
                }

                // Ruse dummy2 and Compute Bphi(0:2, 9:11)
                dummy1.Zero();
                dummy3.Zero();
                dummy1(0,0) = mxi*mQc(0,0);  // N2 * mQc(0:2,0)
                dummy1(1,0) = mxi*mQc(1,0);
                dummy1(2,0) = mxi*mQc(2,0);
//              dummy1(0,0) = (1-mxi)*mQc(0,0);  // N2 * mQc(0:2,0)
//              dummy1(1,0) = (1-mxi)*mQc(1,0);
//              dummy1(2,0) = (1-mxi)*mQc(2,0);

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
                                //Bphi(i,9+j) = dummy3(i,j);
                        }
                }

                return Bphi;
}

void
BeamContact3D::UpdateTransforms(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::UpdateTransforms(): " << MyTag << endln;
#endif

                Vector disp_a(6);           // trial disp/rot vector at a(total disp/rot)
                Vector disp_b(6);           // trial disp/rot vector at a(total disp/rot)  
                Vector rot_a(BC3D_NUM_NDM);  // incr. rot vector at a (from n->n+1)
                Vector rot_b(BC3D_NUM_NDM);  // incr. rot vector at a (from n->n+1)
                Matrix Omega(BC3D_NUM_NDM,BC3D_NUM_NDM);  // Matrix used for Exponential Map

                // Recalculate incremental rotations from step n to n+1
                disp_a = theNodes[0]->getTrialDisp();
                disp_b = theNodes[1]->getTrialDisp();
                int i;
                //for (i=0; i<3; i++) {
                for (i = 3; i<6; i++) {
                        //rot_a(i) = disp_a(i+3) - mRot_a_n(i);
                        //rot_b(i) = disp_b(i+3) - mRot_b_n(i);
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
               
                /* ************************************************************
                // TRY LINEAR UPDATE INSTEAD!  
                Omega = ComputeSkew(rot_a);
                Omega = Omega + meye1;
                mQa = Omega*mQa;

                Omega = ComputeSkew(rot_b);
                Omega = Omega + meye1;
                mQb = Omega*mQb;
                ***************************************************************************************** */

                // Reset total disp & rotation vectors for calculation of
                // incremental values in subsequent step
                //for (i=0; i<3; i++) {
                for (i=0; i<6; i++) {
                        //mRot_a_n(i) = disp_a(i+3);
                        //mRot_b_n(i) = disp_b(i+3);
                        mDisp_a_n(i) = disp_a(i);
                        mDisp_b_n(i) = disp_b(i);
                }

                mDisp_s_n = theNodes[2]->getTrialDisp();

        return;
}

void
BeamContact3D::ComputeQc(double xi)
{      
#ifdef DEBUG
        opserr << "BeamContact3D::ComputeQc(): " << MyTag << endln;
#endif
                Vector c1(BC3D_NUM_NDM);        // tangent vector at projection point, c
                Vector a1(BC3D_NUM_NDM);        // tangent vector at a
                Vector b1(BC3D_NUM_NDM);        // tangent vector at b
                Vector temp(BC3D_NUM_NDM);          // dummy vector for use in calcs
                Matrix Qc_df(BC3D_NUM_NDM,BC3D_NUM_NDM);   // Drill free transformation matrix for c
                Matrix Qc_chi(BC3D_NUM_NDM,BC3D_NUM_NDM);  // Twist transformation matrix for c
                Matrix Qb_df(BC3D_NUM_NDM,BC3D_NUM_NDM);   // Drill free transf. matrix from a to b
                //Matrix Qc(BC3D_NUM_NDM,BC3D_NUM_NDM);      // Total transformation matrix at projection point, c

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

/*opserr << "ELE: " << MyTag << endln;
opserr << "xi = " << mxi << "     mchi = " << mchi << "     testchi = " << testchi << endln;
opserr << "c1 = "  << c1;
opserr << "temp = " << temp;
opserr << "mQa = " << mQa;
opserr << "mQb = " << mQb;
opserr << "Qb_df = " << Qb_df;
opserr << "Qc_df = " << Qc_df;
opserr << "Qc_chi = " << Qc_chi;
opserr << "mQc = " << mQc; */
                return;
}


Matrix
BeamContact3D::ExponentialMap(Vector th)
{
//#ifdef DEBUG
//        opserr << "BeamContact3D::ExponentialMap(): " << MyTag << endln;
//#endif

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
        Vector theta_vec(BC3D_NUM_NDM);                // input vector
        Matrix sk_theta(BC3D_NUM_NDM,BC3D_NUM_NDM);    // skew of vector
        Matrix theta_theta(BC3D_NUM_NDM,BC3D_NUM_NDM); // dyadic product of vector
        Matrix Q(BC3D_NUM_NDM,BC3D_NUM_NDM);           // Exonential Map Vector  

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
                sf3 = 0.5 - theta*theta/24.0 + pow(theta, 4.0)/720.0
                            -pow(theta,6.0)/40320.0 + pow(theta, 8.0)/3628800.0;
        }


        Q = sf1*meye1 + sf2*sk_theta + sf3*theta_theta;

        return Q;
}



Matrix
BeamContact3D::ComputeSkew(Vector th)
{
//#ifdef DEBUG
//        opserr << "BeamContact3D::ComputeSkew(): " << MyTag << endln;
//#endif

        Matrix skew_th(BC3D_NUM_NDM,BC3D_NUM_NDM);

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
BeamContact3D::CrossProduct(Vector &V1, Vector &V2)
{
//#ifdef DEBUG
//        opserr << "BeamContact3D::CrossProduct(): " << MyTag << endln;
//#endif

        Vector V3(3);

        V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
        V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
        V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

        return V3;

}

Matrix  
BeamContact3D::Transpose( int dim1, int dim2, const Matrix &M )
{
// copied from transpose function in Brick.cpp

  Matrix Mtran( dim2, dim1 ) ;

  for ( int i = 0; i < dim1; i++ ) {
     for ( int j = 0; j < dim2; j++ )
         Mtran(j,i) = M(i,j) ;
  } // end for i

  return Mtran ;
}



Vector
BeamContact3D::Geta1(void){
//#ifdef DEBUG
//        opserr << "BeamContact3D::Geta1(): " << MyTag << endln;
//#endif

    Vector a1(BC3D_NUM_NDM);
        int i;

        for (i=0; i<3; i++) {
                a1(i) = mQa(i,0);
        }

        return a1;
}

Vector
BeamContact3D::Getb1(void){
//#ifdef DEBUG
//        opserr << "BeamContact3D::Getb1(): " << MyTag << endln;
//#endif

    Vector b1(BC3D_NUM_NDM);
        int i;

        for (i=0; i<3; i++) {
                b1(i) = mQb(i,0);
        }

        return b1;
}

void
BeamContact3D::Setc1(Vector c1_vec){
//#ifdef DEBUG
//        opserr << "BeamContact3D::Setc1(): " << MyTag << endln;
//#endif

        mc1 = c1_vec;

        return;
}

Vector
BeamContact3D::Getc1(void){
//#ifdef DEBUG
//        opserr << "BeamContact3D::Getc1(): " << MyTag << endln;
//#endif

        return mc1;
}



Vector
BeamContact3D::Getdx_c(double xi){
//#ifdef DEBUG
//        opserr << "BeamContact3D::Getdx_c(): " << MyTag << endln;
//#endif
// ? Is this function necessary?  Calculate directly at place used instead?    
        double xi_squared;
        Vector a1(BC3D_NUM_NDM);
        Vector b1(BC3D_NUM_NDM);
        Vector deriv1(BC3D_NUM_NDM);

        a1 = Geta1();
        b1 = Getb1();
        xi_squared = xi*xi;

        deriv1 = (    - 6.0*xi + 6.0*xi_squared)*mDcrd_a
                   + (1.0 - 4.0*xi + 3.0*xi_squared)*mL*a1
                   + (      6.0*xi - 6.0*xi_squared)*mDcrd_b
                   + (    - 2.0*xi + 3.0*xi_squared)*mL*b1;


        return deriv1;
}

Vector
BeamContact3D::Getddx_c(double xi){
//#ifdef DEBUG
//       opserr << "BeamContact3D::Getddx_c(): " << MyTag << endln;
//#endif
// ? Is this function necessary?  Calculate directly at place used instead?            
        Vector a1(BC3D_NUM_NDM);
        Vector b1(BC3D_NUM_NDM);
        Vector deriv2(BC3D_NUM_NDM);

        a1 = Geta1();
        b1 = Getb1();

        deriv2 = ( - 6.0 + 12.0*xi)*mDcrd_a
                   + ( - 4.0 +  6.0*xi)*mL*a1
                   + (   6.0 - 12.0*xi)*mDcrd_b
                   + ( - 2.0 +  6.0*xi)*mL*b1;

        return deriv2;
}


const Matrix &
BeamContact3D::getTangentStiff(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getTangentStiff(): " << MyTag << endln;
#endif
        int i, j;
                // initialize Kt
        mTangentStiffness.Zero();

        if (inContact) {                // in contact
            Matrix Cmat = theMaterial->getTangent();
           
            double Cnl;
            Vector Csl(2);
            Matrix Css(2,2);
           
            Cnl      = Cmat(0,3);
            Csl(0)   = Cmat(1,3);
            Csl(1)   = Cmat(2,3);        
            Css(0,0) = Cmat(1,1);
            Css(0,1) = Cmat(1,2);
            Css(1,0) = Cmat(2,1);
            Css(1,1) = Cmat(2,2);

            // frictionless contact part
            if (Cnl != 0.0) {
#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                opserr << "   ** tangent: normal" << endln;
                }
#endif
                // assume Cnl == 1 (!)
                for (i = 0; i < BC3D_NUM_DOF - 3; i++) {
                    mTangentStiffness(i,BC3D_NUM_DOF-3) = mBn(i);
                                        mTangentStiffness(BC3D_NUM_DOF-3,i) = mBn(i);
                                }

                //mTangentStiffness(BC3D_NUM_DOF-3,BC3D_NUM_DOF-3) = -1.0e-16;
                                mTangentStiffness(BC3D_NUM_DOF-3,BC3D_NUM_DOF-3) = 0.0;
                                mTangentStiffness(BC3D_NUM_DOF-2,BC3D_NUM_DOF-2) =  1.0;
                                mTangentStiffness(BC3D_NUM_DOF-1,BC3D_NUM_DOF-1) =  1.0;
             }

            // contribution from friction
                       
                        //  sliding: Csl = mu * r   &    Css != 0
                        // sticking: Csl = 0        &    Css != 0
           
            bool inSlip;
            inSlip = theMaterial->getContactState();
            if (inSlip) {              
            //if (Csl(0) != 0.0 || Csl(1) != 0.0 ) {

            // sliding
#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                opserr << "   ** tangent: sliding" << endln; }
#endif
//opserr << "   ** tangent: sliding" << endln;

                for (i = 0; i < BC3D_NUM_DOF-3; i++) {              
                    // assume first is the row
                    mTangentStiffness(i,BC3D_NUM_DOF-3) += mBs(i,0)*Csl(0)
                                                        +  mBs(i,1)*Csl(1);
                    }
                        // end sliding
            } else {
            // sticking
#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                opserr << "   ** tangent: sticking" << endln;  }
#endif
         }
                        // in 3D: both sticking & sliding cases have Css coefficients
                        // in 2D: Css = 0 for sliding
            for (i = 0; i < BC3D_NUM_DOF-3; i++) {
                for (j = 0; j < BC3D_NUM_DOF-3; j++) {
                    mTangentStiffness(i,j) += mBs(i,0)*mBs(j,0)*Css(0,0)
                                           +  mBs(i,1)*mBs(j,0)*Css(1,0)
                                           +  mBs(i,0)*mBs(j,1)*Css(0,1)
                                           +  mBs(i,1)*mBs(j,1)*Css(1,1);
                                }
                        }

         } else {    
             // not in contact
                mTangentStiffness(BC3D_NUM_DOF-3,BC3D_NUM_DOF-3) = 1.0;
                mTangentStiffness(BC3D_NUM_DOF-2,BC3D_NUM_DOF-2) = 1.0;
                mTangentStiffness(BC3D_NUM_DOF-1,BC3D_NUM_DOF-1) = 1.0;
        }

//              opserr << "   K_T: " << mTangentStiffness;

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                        opserr << "   inContact:   " << inContact << endln;
            opserr << "   K_T: " << mTangentStiffness;
            }
#endif


        return mTangentStiffness;
}

const Matrix &
BeamContact3D::getInitialStiff(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getInitialStiff(): " << MyTag << endln;
#endif
  return getTangentStiff();
}
   
void
BeamContact3D::zeroLoad(void)
{
#ifdef DEBUG
        opserr << "BeamContact3D::zeroLoad(): " << MyTag << endln;
#endif
  return;
}

int
BeamContact3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
#ifdef DEBUG
        opserr << "BeamContact3D::addLoad(ElementalLoad *theLoad, double loadFactor): " << MyTag << endln;
#endif
  return 0;
}

int
BeamContact3D::addInertiaLoadToUnbalance(const Vector &accel)
{
#ifdef DEBUG
        opserr << "BeamContact3D::addInertiaLoadToUnbalance(const Vector &accel): " << MyTag << endln;
#endif
  return 0;
}

const Vector &
BeamContact3D::getResistingForce()
{
#ifdef DEBUG
        opserr << "BeamContact3D::getResistingForce(): " << MyTag << endln;
#endif

      int     i;                
      // initialize F
      mInternalForces.Zero();

         // get contact stresses/forces:  stress(0) = lambda, stress(1) = t_s(0), stress(2) = t_s(1)
     Vector stress = theMaterial->getStress();

      if (inContact) {
                  // in contact
          for (i=0; i < BC3D_NUM_DOF - 3; i++) {
                          mInternalForces(i) = mLambda*mBn(i) + stress(1)*mBs(i,0) + stress(2)*mBs(i,1);
                  }
                  mInternalForces(BC3D_NUM_DOF - 3) = -mGap;
         
          } else {
                 
                  mInternalForces(BC3D_NUM_DOF - 3) = mLambda;
       
          }

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                        opserr << "inContact:  " << inContact << endln;
            opserr << "   F_int  = " << mInternalForces;
            }
#endif

        return mInternalForces;
}


const Vector &
BeamContact3D::getResistingForceIncInertia()
{      
#ifdef DEBUG
        opserr << "BeamContact3D::getResistingForceIncInertia(): " << MyTag << endln;
#endif
  //return theVector;
                return mInternalForces;
}




int
BeamContact3D::sendSelf(int commitTag, Channel &theChannel)
{
#ifdef DEBUG
        opserr << "BeamContact3D::sendSelf(int commitTag, Channel &theChannel)" << endln;
#endif
       
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // BeamContact3D packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(5);
  data(0) = this->getTag();
  data(1) = BC3D_NUM_DOF;
  data(2) = mTolGap;
  data(3) = theMaterial->getClassTag();

  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(4) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }          

  // BeamContact3D then sends the tags of it's four nodes

  res = theChannel.sendID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally BeamContact3D asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
BeamContact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
#ifdef DEBUG
        opserr << "BeamContact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)" << endln;
#endif

  int res;
  int dataTag = this->getDbTag();

  // BeamContact3D creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector

  static Vector data(5);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::recvSelf() - failed to receive Vector\n";
    return -1;
  }          

  this->setTag((int)data(0));
  // BC3D_NUM_DOF = (int)data(1);       // this must not be since BC3D_NUM_DOF is used to initialize variables and thus must not change
  mTolGap = data(2);
 
  // BeamContact3D now receives the tags of it's four external nodes
  res = theChannel.recvID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally BeamContact3D creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(3);
  int matDb = (int)data(4);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (ContactMaterial3D *)theMatCopy;

    if (theMaterial == 0) {
      opserr <<"WARNING BeamContact3D::recvSelf() - " << this->getTag()
        << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING BeamContact3D::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
BeamContact3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
#ifdef DEBUG
        opserr << "BeamContact3D::displaySelf(Renderer &theViewer, int displayMode, float fact)" << endln;
#endif
  return 0;
}


void
BeamContact3D::Print(OPS_Stream &s, int flag)
{
#ifdef DEBUG
        opserr << "BeamContact3D::Print(OPS_Stream &s, int flag)" << endln;
#endif
        opserr << "BeamContact3D, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  " ;
        for (int i = 0; i<BC3D_NUM_NODE; i++)
        {
                opserr << externalNodes(i)<< " ";
        }
        return;
}


Response*
BeamContact3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
#ifdef DEBUG
        opserr << "BeamContact3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): " << MyTag << endln;
#endif
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
      return new ElementResponse(this, 1, Vector(3));

    else if (strcmp(argv[0],"frictionforce") == 0 || strcmp(argv[0],"frictionforces") == 0)
      return new ElementResponse(this, 2, Vector(3));
   
    else if (strcmp(argv[0],"forcescalar") == 0 || strcmp(argv[0],"forcescalars") == 0)
      return new ElementResponse(this, 3, Vector(3));

    else if (strcmp(argv[0],"masterforce") == 0 || strcmp(argv[0],"masterforces") == 0)
      return new ElementResponse(this, 4, Vector(6));

    else if (strcmp(argv[0],"mastermoment") == 0 || strcmp(argv[0],"mastermoments") == 0)
      return new ElementResponse(this, 5, Vector(6));

    else if (strcmp(argv[0],"masterreaction") == 0 || strcmp(argv[0],"masterreactions") == 0)
      return new ElementResponse(this, 6, Vector(12));

    // otherwise response quantity is unknown for the BeamContact3D class
    else
        opserr << "BeamContact3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): " << argv[0] << " unknown request" << endln;
      return 0;
}


int
BeamContact3D::getResponse(int responseID, Information &eleInfo)
{
#ifdef DEBUG
        opserr << "BeamContact3D::getResponse(int responseID =" << responseID << ", Information &eleInfo): " << MyTag << endln;
#endif
          Vector force(3);
          Vector sForce(3);
          Vector mForces(6);
          Vector mMoments(6);
          Vector mReactions(12);
     
          // get contact stresse/forces
      Vector stress = theMaterial->getStress();

  if (responseID == 1) {
     
        // force = stress(0)*mn + stress(1)*mg1 + stress(2)*mg2;

        // forces on slave node
          for (int ii=0; ii<3; ii++) {
                  sForce(ii)   = -mInternalForces(BC3D_NUM_DOF - 6 + ii);
/*
                  if ( sForce(ii) != force(ii) ) {
                        opserr << "** problem in recorder: " << sForce(ii) << " != " << force(ii) << endln ;
                        opserr << "   stress = " << stress;
            opserr << "   Fn     = " << stress(0)*mn;
            opserr << "   Fs1    = " << stress(1)*mg1;
            opserr << "   Fs2    = " << stress(2)*mg2;
            opserr << "   force  = " << force;
            opserr << "   F_int  = " << mInternalForces;
                        }
*/
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

        // forces on master nodes
          for (int ii=0; ii<3; ii++) {
                  mForces(ii)   = -mInternalForces(ii);
                  mForces(ii+3) = -mInternalForces(ii+6);
                }

          return eleInfo.setVector(mForces);
 
  } else if (responseID == 5) {

        // moments on master nodes
          for (int ii=0; ii<3; ii++) {
                  mMoments(ii)   = -mInternalForces(ii+3);
                  mMoments(ii+3) = -mInternalForces(ii+9);
                }

          return eleInfo.setVector(mMoments);
 
  } else if (responseID == 6) {

        // full reactions on master nodes
          for (int ii=0; ii<6; ii++) {
                  mReactions(ii)   = -mInternalForces(ii);
                  mReactions(ii+6) = -mInternalForces(ii+6);
                }

          return eleInfo.setVector(mReactions);
 
  } else

        opserr << "BeamContact3D::getResponse(int responseID=" << responseID << ", Information &eleInfo): " << " unknown request" << endln;
    return -1;
}

int
BeamContact3D::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"friction") == 0) {
		return param.addObject(1, this);
	}
	
	return -1;
}

int
BeamContact3D::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes =  theMaterial->updateParameter(parameterID, info);
	if (matRes != -1) {
		res = matRes;
	}
	return res;
}
