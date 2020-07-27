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
                                                                        
// $Revision: 1.1
// $Date: 
// $Source: /OpenSees/SRC/element/SimpleContact/SimpleContact3D.cpp,v $
                                                                        
// Written: Kathryn A. Petek
// Created: 02/04
//
// Revisions
//    02/04 created
//    11/10 F.Mckenna and C.McGann: changes for incorporation into main source code
//
// Description: This file contains the implementation for the SimpleContact3D class.
//

#include "SimpleContact3D.h"
#include <Information.h>
#include <ElementResponse.h>

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

static int num_SimpleContact3D = 0;

OPS_Export void *
OPS_SimpleContact3D(void)
{
  if (num_SimpleContact3D == 0) {
    num_SimpleContact3D++;
    //OPS_Error("SimpleContact3D element - Written: K.Petek, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr<<"SimpleContact3D element - Written: K.Petek, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 10) {
    opserr << "Invalid #args,  want: element SimpleContact3D eleTag? iNode? jNode? kNode? lNode? secondaryNode? lambdaNode? matTag? tolGap? tolForce?\n";
    return 0;
  }
    
  int    iData[8];
  double dData[2];

  int numData = 8;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact3DElement" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact3D " << iData[0] << endln;
    return 0;	
  }

  int matID = iData[7];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new SimpleContact3D(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6],
				   *theMaterial,
				   dData[0], dData[1]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SimpleContact3DElement\n";
    return 0;
  }

  return theElement;
}



// constructors:
SimpleContact3D::SimpleContact3D(int tag, int Nd1, int Nd2, int Nd3, int Nd4, 
                int NdS, int NdL, NDMaterial &theMat, double tolG, double tolF)
 :Element(tag,ELE_TAG_SimpleContact3D),     
   externalNodes(SC3D_NUM_NODE),
   tangentStiffness(SC3D_NUM_DOF, SC3D_NUM_DOF),
   internalForces(SC3D_NUM_DOF),
   d(SC3D_NUM_NDF), 
   x(SC3D_NUM_NDF,5),
   g_metric(2,2),
   G_metric(2,2),
   XiEta_n(2),
   XiEta_nplus1(2),
   x_XiCrd(SC3D_NUM_NDF),
   slip(2),
   g1(SC3D_NUM_NDF), 
   g2(SC3D_NUM_NDF),
   n(SC3D_NUM_NDF), 
   Kinv(2,2),
   KinvLin(2,2),
   Bn(SC3D_NUM_DDOF), Bs(SC3D_NUM_DDOF,2),
   dcrd1(SC3D_NUM_NDF),
   dcrd2(SC3D_NUM_NDF),
   dcrd3(SC3D_NUM_NDF),
   dcrd4(SC3D_NUM_NDF),
   dcrdS(SC3D_NUM_NDF),
   dispL(SC3D_NUM_NDF)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::SimpleContact3D(): " << MyTag << endln;
#endif
        externalNodes(0) = Nd1;
        externalNodes(1) = Nd2;
        externalNodes(2) = Nd3;
        externalNodes(3) = Nd4;
        externalNodes(4) = NdS;
        externalNodes(5) = NdL;

		MyTag = tag;

        tolGap = tolG;
        tolForce = tolF;

        inContact          = false;
        was_inContact      = false;
        should_be_released = false;
        to_be_released     = false;
		in_bounds          = false;

        gap    = 0.0;
        lambda = 0.0;
        slip.Zero();
        
	NDMaterial *theMatCopy = theMat.getCopy("ContactMaterial3D");
	if (theMatCopy != 0) {
	  theMaterial = (ContactMaterial3D *)theMatCopy;
	} else {
	  opserr << "SimpleContact3D::SimpleContact3D - material needs to be of type Contact3D for ele: " << this->getTag() << endln;
	}

        if (theMaterial == 0) {
	  opserr << "SimpleContact3D::SimpleContact3D - failed allocate material model pointer\n";
	  exit(-1);
	}

}

SimpleContact3D::SimpleContact3D()
 :Element(0,ELE_TAG_SimpleContact3D),     
   externalNodes(SC3D_NUM_NODE),
   tangentStiffness(SC3D_NUM_DOF, SC3D_NUM_DOF),
   internalForces(SC3D_NUM_DOF),
   d(SC3D_NUM_NDF), 
   x(SC3D_NUM_NDF,5),
   g_metric(2,2),
   G_metric(2,2),
   XiEta_n(2),
   XiEta_nplus1(2),
   x_XiCrd(SC3D_NUM_NDF),
   slip(2),
   g1(SC3D_NUM_NDF), 
   g2(SC3D_NUM_NDF),
   n(SC3D_NUM_NDF), 
   Kinv(2,2),
   KinvLin(2,2),
   Bn(SC3D_NUM_DDOF), Bs(SC3D_NUM_DDOF,2),
   dcrd1(SC3D_NUM_NDF),
   dcrd2(SC3D_NUM_NDF),
   dcrd3(SC3D_NUM_NDF),
   dcrd4(SC3D_NUM_NDF),
   dcrdS(SC3D_NUM_NDF),
   dispL(SC3D_NUM_NDF)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::SimpleContact3D(): " << MyTag << endln;
#endif
}


//  destructor:
SimpleContact3D::~SimpleContact3D()
{
#ifdef DEBUG
        opserr << "SimpleContact3D::~SimpleContact3D(): " << MyTag << endln;
#endif
        delete theMaterial;
}


int
SimpleContact3D::getNumExternalNodes(void) const
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getNumExternalNodes(): " << MyTag << endln;
#endif
    return SC3D_NUM_NODE;
}

const ID &
SimpleContact3D::getExternalNodes(void) 
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getExternalNodes(): " << MyTag << endln;
#endif
    return externalNodes;
}


Node **
SimpleContact3D::getNodePtrs(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getNodePtrs(): " << MyTag << endln;
#endif
        return theNodes;                        
}

int
SimpleContact3D::getNumDOF(void) 
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getNumDOF(): " << MyTag << endln;
#endif
    return SC3D_NUM_DOF;
}


void
SimpleContact3D::setDomain(Domain *theDomain)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::setDomain(Domain *theDomain): " << MyTag << endln;
#endif

        int Nd1 = externalNodes(0);
        int Nd2 = externalNodes(1);
        int Nd3 = externalNodes(2);
        int Nd4 = externalNodes(3);
        int NdS = externalNodes(4);
        int NdL = externalNodes(5);

        theNodes[0] = theDomain->getNode(Nd1);
        theNodes[1] = theDomain->getNode(Nd2);
        theNodes[2] = theDomain->getNode(Nd3);
        theNodes[3] = theDomain->getNode(Nd4);
        theNodes[4] = theDomain->getNode(NdS);
        theNodes[5] = theDomain->getNode(NdL);

        int i;
        for (i = 0; i < 6; i++) {
           if (theNodes[i] == 0)
                return;  // don't go any further - otherwise segmentation fault
        }

        dcrd1 = theNodes[0]->getCrds();
        dcrd2 = theNodes[1]->getCrds();
        dcrd3 = theNodes[2]->getCrds();
        dcrd4 = theNodes[3]->getCrds();
        dcrdS = theNodes[4]->getCrds();
        dispL.Zero();
        
        // length of primary segments
        Vector L1       = (dcrd2 - dcrd1);
        // double L1primary = L1.Norm();
        Vector L2       = (dcrd3 - dcrd2);
        // double L2primary = L2.Norm();
        Vector L3       = (dcrd4 - dcrd3);
        // double L3primary      = L3.Norm();
        Vector L4       = (dcrd1 - dcrd4);
        // double L4primary = L4.Norm();

        // error check that nodes are not in same location
        if (fabs(L1(0)) < tolGap && fabs(L1(1)) < tolGap && fabs(L1(2)) < tolGap ) { 
                opserr << "SimpleContact3D::SimpleContact3D - node 1 and node 2 share same coordinates\n";
                opserr << "Program Terminated\n";
                exit(-1);
                }
        if (fabs(L2(0)) < tolGap && fabs(L2(1)) < tolGap && fabs(L2(2)) < tolGap ) { 
                opserr << "SimpleContact3D::SimpleContact3D - node 2 and node 3 share same coordinates\n";
                opserr << "Program Terminated\n";
                exit(-1);
                }
        if (fabs(L3(0)) < tolGap && fabs(L3(1)) < tolGap && fabs(L3(2)) < tolGap ) { 
                opserr << "SimpleContact3D::SimpleContact3D - node 3 and node 4 share same coordinates\n";
                opserr << "Program Terminated\n";
                exit(-1);
                }
        if (fabs(L4(0)) < tolGap && fabs(L4(1)) < tolGap && fabs(L4(2)) < tolGap ) { 
	  opserr << "SimpleContact3D::SimpleContact3D - node 1 and node 4 share same coordinates\n";
	  opserr << "Program Terminated\n";
	  exit(-1);
	}

        XiEta_n.Zero();  // initially xi and eta = 0

        XiEta_n = this->project(XiEta_n);
        // this->UpdateBase(XiEta_n);

        // Set initial values of g_metric, n, Bn, and Bs 
        // (calculated based on XiEta_nplus1)
        XiEta_nplus1 = XiEta_n;

        // calculation of metric tensor
        g_metric(0,0) =   g1 ^ g1;
        g_metric(0,1) =   g1 ^ g2;
        g_metric(1,0) =   g_metric(0,1);
        g_metric(1,1) =   g2 ^ g2;

                // calculation of Kinvlin (contra-variant metric tensor)
                // for linear projection in update()
                double detKLin = g_metric(0,0)*g_metric(1,1) - g_metric(1,0)*g_metric(0,1);
                KinvLin(0,0) =  g_metric(1,1)/detKLin;
                KinvLin(1,0) = -g_metric(1,0)/detKLin;
                KinvLin(0,1) = -g_metric(0,1)/detKLin;
                KinvLin(1,1) =  g_metric(0,0)/detKLin;

        // normal vector to primary surface
        n(0) = g1(1)*g2(2) - g1(2)*g2(1);
        n(1) = g1(2)*g2(0) - g1(0)*g2(2);
        n(2) = g1(0)*g2(1) - g1(1)*g2(0);
        n = n / n.Norm();

        theMaterial->setMetricTensor(g_metric); 
        
        // adjust cohesion force
        double area = 4.0*sqrt(g_metric(0,0)*g_metric(1,1)-g_metric(0,1)*g_metric(1,0));
        theMaterial->ScaleCohesion(area);
		theMaterial->ScaleTensileStrength(area);
        
        // compute Bn, Bs
        this->ComputeB();

        // call the base class method
        this->DomainComponent::setDomain(theDomain);

#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
                opserr << "   g1: " << g1 ;
                opserr << "   g2: " << g2 ;
                opserr << "   n: " << n ;
                }
#endif

}        


int
SimpleContact3D::commitState()
{
#ifdef DEBUG
        opserr << "SimpleContact3D::commitState(): " << MyTag << endln;
#endif
        
        was_inContact  = ( gap < tolGap );


// --- START --- for geometry update between increments
            
        XiEta_nplus1 = this->project(XiEta_nplus1);
        // this->UpdateBase(XiEta_nplus1);

        // calculate g_metric, KinvLin n, Bn, and Bs for next step
        // calculation of metric tensor
        g_metric(0,0) =   g1 ^ g1;
        g_metric(0,1) =   g1 ^ g2;
        g_metric(1,0) =   g_metric(0,1);
        g_metric(1,1) =   g2 ^ g2;

        // calculation of Kinvlin (contra-variant metric tensor)
        // for linear projection in update()
        double detKLin = g_metric(0,0)*g_metric(1,1) - g_metric(1,0)*g_metric(0,1);
        KinvLin(0,0) =  g_metric(1,1)/detKLin;
        KinvLin(1,0) = -g_metric(1,0)/detKLin;
        KinvLin(0,1) = -g_metric(0,1)/detKLin;
        KinvLin(1,1) =  g_metric(0,0)/detKLin;


        // normal vector to primary surface
        n(0) = g1(1)*g2(2) - g1(2)*g2(1);
        n(1) = g1(2)*g2(0) - g1(0)*g2(2);
        n(2) = g1(0)*g2(1) - g1(1)*g2(0);
        n = n / n.Norm();

        theMaterial->setMetricTensor(g_metric); 

        this->ComputeB();

// --- END --- for geometry update between increments

        XiEta_n = XiEta_nplus1; 


		in_bounds      = ( (fabs(XiEta_n(0)) <= 1) && 
			               (fabs(XiEta_n(1)) <= 1) );

		to_be_released = ( should_be_released || !in_bounds );
        inContact = ( was_inContact && !to_be_released && in_bounds );


        
        int retVal = 0;
        // call element commitState to do any base class stuff
        if ((retVal = this->Element::commitState()) != 0) {
                opserr << "SimpleContact3D::commitState () - failed in base class";
                }    
        retVal = theMaterial->commitState();
        return retVal; 
}


int
SimpleContact3D::revertToLastCommit()
{
#ifdef DEBUG
        opserr << "SimpleContact3D::revertToLastCommit(): " << MyTag << endln;
#endif
        return theMaterial->revertToLastCommit();
}

int
SimpleContact3D::revertToStart()
{
#ifdef DEBUG
        opserr << "SimpleContact3D::revertToStart(): " << MyTag << endln;
#endif
          inContact = false;
          was_inContact = false;
          should_be_released = false;
          to_be_released = false;
		  in_bounds = false;
          
          return theMaterial->revertToStart();
}

int
SimpleContact3D::update(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::update(): " << MyTag << endln;
#endif
//        double cohesion;
		double tensileStrength;
                        
        Vector deltaXi(2);          // vector of delta xi and eta
        Vector strain(4);           // vector to send "strains" to material

        dcrd1 = theNodes[0]->getCrds() +  theNodes[0]->getTrialDisp();
        dcrd2 = theNodes[1]->getCrds() +  theNodes[1]->getTrialDisp();
        dcrd3 = theNodes[2]->getCrds() +  theNodes[2]->getTrialDisp();
        dcrd4 = theNodes[3]->getCrds() +  theNodes[3]->getTrialDisp();
        dcrdS = theNodes[4]->getCrds() +  theNodes[4]->getTrialDisp();
        dispL = theNodes[5]->getTrialDisp();
        
        // calculate initial g1, g2, x_Xi @ step n+1
        x_XiCrd = GetPoint(XiEta_n);
        // this->UpdateBase(XiEta_n); // do not update during iteration -> commitState()

        d = dcrdS - x_XiCrd;
        gap = n ^ d;

        Vector R(2);
        R(0) = d ^ g1;  // delta xi  (co-variant)
        R(1) = d ^ g2;  // delta eta (co-variant)

                /*
        //Matrix K(2,2);
        double detK = g_metric(0,0)*g_metric(1,1) - g_metric(1,0)*g_metric(0,1);
      
        // calculation of Kinv (contra-variant metric tensor)
        Kinv(0,0) =  g_metric(1,1)/detK;
        Kinv(1,0) = -g_metric(1,0)/detK;
        Kinv(0,1) = -g_metric(0,1)/detK;
        Kinv(1,1) =  g_metric(0,0)/detK;

                deltaXi = Kinv * R;     // conta-variant components of slip
                */

        deltaXi = KinvLin * R;     // conta-variant components of slip

        XiEta_nplus1 = XiEta_n + deltaXi;

        
        // n, N1-N4 use undeformed coords, rest is in deformed coords
        
        // note: if need N1-N4 here, take out of ComputeB and put into CommitState
        //gap = n ^ (endSCrd + dispS - N1*(end1Crd + disp1) - N2*(end2Crd + disp2)
        //                         - N3*(end3Crd + disp3) - N4*(end4Crd + disp4) );
         
        lambda = dispL(0);
        
		tensileStrength = theMaterial->getTensileStrength();
		

        // default: not in contact ( inContact = false)
        slip.Zero();                            
        
        // define state of the contact element
        should_be_released = ( lambda <= -(tensileStrength + tolForce ) );


#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
            opserr << "   CONTACT:            " << inContact << endln;
            opserr << "   should be released: " << should_be_released << endln;
            }
#endif

        if (inContact) {        
                  
            slip = XiEta_nplus1 - XiEta_n;

            // slip is treated in local (curvilinear) coordinates
            strain(0) = gap; 
            strain(1) = slip(0);
            strain(2) = slip(1);
            strain(3) = lambda;
       
            theMaterial->setTrialStrain(strain);
        }

        else if (to_be_released) {
		// prevents sliding & stabilizes behavior in lag step

            strain(0) = gap; 
            strain(1) = 0.0;
            strain(2) = 0.0;
            strain(3) = lambda;
       
            theMaterial->setTrialStrain(strain);
        }

#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
            opserr << "end of update\n";
            }
#endif

  return 0;
}

Vector
SimpleContact3D::project(Vector XiEta0)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::project(): " << MyTag << endln;
#endif
        Vector XiEta_P(2);
        
        XiEta_P = XiEta0;

        Vector x_XiCrd(SC3D_NUM_NDF);
        Vector deltaXi(2);          // vector of delta xi and eta

        // calculate initial g1, g2, x_Xi @ step n+1
        x_XiCrd = GetPoint(XiEta_P);  
        this->UpdateBase(XiEta_P);  

        d = dcrdS - x_XiCrd;
        gap = n ^ d;

        Vector R(2);
        R(0) = d ^ g1;
        R(1) = d ^ g2;

        // loop to calculate g1, g2, deltaXi, x_Xi @ step n+1
        while (R.Norm() > tolGap) {
                //Matrix K(2,2);
                double K2dbl;
                double detK;
              
            // calculation of Kinv
                Kinv(0,0) =   g1 ^ g1;
                Kinv(0,1) =   g1 ^ g2;
                Kinv(1,0) =   Kinv(0,1);
                Kinv(1,1) =   g2 ^ g2;
        
                K2dbl  = -( d ^ ( dcrd1 - dcrd2 + dcrd3 - dcrd4 )) * 0.25;
                Kinv(0,1) += K2dbl;
                Kinv(1,0) += K2dbl;

                detK = Kinv(0,0)*Kinv(1,1) - Kinv(1,0)*Kinv(0,1);
                K2dbl = Kinv(0,0);
                Kinv(0,0) =  Kinv(1,1)/detK;
                Kinv(1,0) = -Kinv(1,0)/detK;
                Kinv(0,1) = -Kinv(0,1)/detK;
                Kinv(1,1) =  K2dbl/detK;

                deltaXi = Kinv * R;
    
                XiEta_P = XiEta_P + deltaXi;

            // calculates g1, g2, x_XiCrd
                x_XiCrd = GetPoint(XiEta_P);
                this->UpdateBase(XiEta_P);

                d = dcrdS - x_XiCrd;

                R(0) = d ^ g1;
                R(1) = d ^ g2;

                }

  return XiEta_P;
}

Vector
SimpleContact3D::GetPoint(Vector XiEta)
// this function calculates x_Xi 
{
#ifdef DEBUG
        opserr << "SimpleContact3D::GetPoint(): " << MyTag << endln;
#endif

        Vector x_Xi(SC3D_NUM_NDF);

        double oneMinusEta = 1 - XiEta(1);
        double onePlusEta  = 1 + XiEta(1);
        double oneMinusXi  = 1 - XiEta(0);
        double onePlusXi   = 1 + XiEta(0);

        // x_Xi = N1*x1 + N2*x2 + N3*x3 + N4*x4
        x_Xi  = 0.25 * (oneMinusXi * oneMinusEta * dcrd1
                      + onePlusXi  * oneMinusEta * dcrd2
                      + onePlusXi  * onePlusEta  * dcrd3
                      + oneMinusXi * onePlusEta  * dcrd4);    
            
        return x_Xi;
}


int
SimpleContact3D::UpdateBase(Vector XiEta)
// this function calculates g1, g2
{
#ifdef DEBUG
        opserr << "SimpleContact3D::UpdateBase(): " << MyTag << endln;
#endif

        double oneMinusEta = 1 - XiEta(1);
        double onePlusEta  = 1 + XiEta(1);
        double oneMinusXi  = 1 - XiEta(0);
        double onePlusXi   = 1 + XiEta(0);

        // calculate vectors g1 and g2
        // g1 = d(x_Xi)/dXi, g2 = d(x_Xi)/dEta
        g1 = (oneMinusEta * (dcrd2 - dcrd1) + onePlusEta  * (dcrd3 - dcrd4)) * 0.25;
        g2 = (onePlusXi   * (dcrd3 - dcrd2) + oneMinusXi  * (dcrd4 - dcrd1)) * 0.25;

        return 0;
}


// Bn and Bs are computed at CommitState for the next step
// this allows to only be calculated once
void
SimpleContact3D::ComputeB(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::ComputeB(): " << MyTag << endln;
#endif

        // Xi_Eta_nplus1, g1, g2, Kinv all from converged state (n+1)
        // vector n used below should be of same step (calculated in CommitState)
        
        N1 = 0.25*(1 - XiEta_nplus1(0))*(1 - XiEta_nplus1(1));
        N2 = 0.25*(1 + XiEta_nplus1(0))*(1 - XiEta_nplus1(1));
        N3 = 0.25*(1 + XiEta_nplus1(0))*(1 + XiEta_nplus1(1));
        N4 = 0.25*(1 - XiEta_nplus1(0))*(1 + XiEta_nplus1(1));

#ifdef DEBUG
        if (DEBUG_LEVEL > 2) {
              opserr << "\tN1+N2+N3+N4 = " << N1+N2+N3+N4 << endln;
              }
#endif
        // dgap = Bn' * dq
        Bn(0)  = -n(0)*N1;
        Bn(1)  = -n(1)*N1;
        Bn(2)  = -n(2)*N1;

        Bn(3)  = -n(0)*N2;
        Bn(4)  = -n(1)*N2;
        Bn(5)  = -n(2)*N2;

        Bn(6)  = -n(0)*N3;
        Bn(7)  = -n(1)*N3;
        Bn(8)  = -n(2)*N3;

        Bn(9)  = -n(0)*N4;
        Bn(10) = -n(1)*N4;
        Bn(11) = -n(2)*N4;

        Bn(12) =  n(0);
        Bn(13) =  n(1);
        Bn(14) =  n(2);


        // dslip = Bs' * dq
        //   Bs' = (-dR/dXi)^-1 * (dR/dq)
        //   Bs' =   Kinv    * (dR/dq)
        //   Bs  =  (dR/dq)' * Kinv'
        Matrix dR_dqT(18,2);
        
        dR_dqT(0,0)  = -g1(0)*N1;
        dR_dqT(1,0)  = -g1(1)*N1;
        dR_dqT(2,0)  = -g1(2)*N1;
        dR_dqT(3,0)  = -g1(0)*N2;
        dR_dqT(4,0)  = -g1(1)*N2;
        dR_dqT(5,0)  = -g1(2)*N2;
        dR_dqT(6,0)  = -g1(0)*N3;
        dR_dqT(7,0)  = -g1(1)*N3;
        dR_dqT(8,0)  = -g1(2)*N3;
        dR_dqT(9,0)  = -g1(0)*N4;
        dR_dqT(10,0) = -g1(1)*N4;
        dR_dqT(11,0) = -g1(2)*N4;
        dR_dqT(12,0) =  g1(0);
        dR_dqT(13,0) =  g1(1);
        dR_dqT(14,0) =  g1(2);
        
        dR_dqT(0,1)  = -g2(0)*N1;
        dR_dqT(1,1)  = -g2(1)*N1;
        dR_dqT(2,1)  = -g2(2)*N1;
        dR_dqT(3,1)  = -g2(0)*N2;
        dR_dqT(4,1)  = -g2(1)*N2;
        dR_dqT(5,1)  = -g2(2)*N2;
        dR_dqT(6,1)  = -g2(0)*N3;
        dR_dqT(7,1)  = -g2(1)*N3;
        dR_dqT(8,1)  = -g2(2)*N3;
        dR_dqT(9,1)  = -g2(0)*N4;
        dR_dqT(10,1) = -g2(1)*N4;
        dR_dqT(11,1) = -g2(2)*N4;
        dR_dqT(12,1) =  g2(0);
        dR_dqT(13,1) =  g2(1);
        dR_dqT(14,1) =  g2(2);
        
        int i;
        double detK = g_metric(0,0)*g_metric(1,1) - g_metric(1,0)*g_metric(0,1);
        for (i = 0; i < SC3D_NUM_DDOF; i++) {
            // Bs(i,0) = dR_dqT(i,0)*Kinv(0,0) + dR_dqT(i,1)*Kinv(0,1);
            // Bs(i,1) = dR_dqT(i,0)*Kinv(1,0) + dR_dqT(i,1)*Kinv(1,1);
            Bs(i,0) = ( dR_dqT(i,0)*g_metric(1,1) - dR_dqT(i,1)*g_metric(0,1))/detK;
            Bs(i,1) = (-dR_dqT(i,0)*g_metric(1,0) + dR_dqT(i,1)*g_metric(0,0))/detK;
        }
        
        return;
}

const Matrix &
SimpleContact3D::getTangentStiff(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getTangentStiff(): " << MyTag << endln;

        if (DEBUG_LEVEL > 1) {
            opserr << "   inContact:   " << inContact << endln;
            }
#endif

        // initialize Kt
        tangentStiffness.Zero();

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

#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
            opserr << "   Cmat = " << Cmat
                   << "   C_nl = " << Cnl
                   << endln
                   << "   C_ss = " << Css
                   << "   C_sl = " << Csl
                   << endln;
            }
#endif

            int i, j;

            // frictionless contact part

            if (Cnl != 0.0) {
#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
                opserr << "   ** tangent: normal" << endln;
                }
#endif
                // assume Cnl == 1 (!)

                for (i = 0; i < SC3D_NUM_DDOF; i++) {

                    tangentStiffness(i,SC3D_NUM_DDOF) -= Bn(i);
                    tangentStiffness(SC3D_NUM_DDOF,i) -= Bn(i);
                    }
                tangentStiffness(SC3D_NUM_DOF-3,SC3D_NUM_DOF-3) = -1.0e-10;
                tangentStiffness(SC3D_NUM_DOF-2,SC3D_NUM_DOF-2) =  1.0;
                tangentStiffness(SC3D_NUM_DOF-1,SC3D_NUM_DOF-1) =  1.0;
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
                opserr << "   ** tangent: sliding" << endln;
                }
#endif
opserr << "   ** tangent: sliding" << endln;

                for (i = 0; i < SC3D_NUM_DDOF; i++) {               
                    // assume first is the row
                    tangentStiffness(i,SC3D_NUM_DDOF) += Bs(i,0)*Csl(0) 
                                                      +  Bs(i,1)*Csl(1);
                    }
         }       // end sliding

         else {
            // sticking
#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
                opserr << "   ** tangent: sticking" << endln;
                }
#endif
         }

         
         // in 3D: both sticking & sliding cases have Css coefficients
         // in 2D: Css = 0 for sliding
         for (i = 0; i < SC3D_NUM_DDOF; i++) {
            for (j = 0; j < SC3D_NUM_DDOF; j++) {
               tangentStiffness(i,j) += Bs(i,0)*Bs(j,0)*Css(0,0)
                                     +  Bs(i,1)*Bs(j,0)*Css(1,0)
                                     +  Bs(i,0)*Bs(j,1)*Css(0,1)
                                     +  Bs(i,1)*Bs(j,1)*Css(1,1);
            }
         } 
 

        } else {     
             // not in contact
                tangentStiffness(SC3D_NUM_DOF-3,SC3D_NUM_DOF-3) = 1.0;
                tangentStiffness(SC3D_NUM_DOF-2,SC3D_NUM_DOF-2) = 1.0;
                tangentStiffness(SC3D_NUM_DOF-1,SC3D_NUM_DOF-1) = 1.0;
        }

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
            opserr << "   K_T: " << tangentStiffness;
            }
#endif
        return tangentStiffness;
}

const Matrix &
SimpleContact3D::getInitialStiff(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getInitialStiff(): " << MyTag << endln;
#endif
  return getTangentStiff();
}
    
void 
SimpleContact3D::zeroLoad(void)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::zeroLoad(): " << MyTag << endln;
#endif
  return;
}

int 
SimpleContact3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::addLoad(ElementalLoad *theLoad, double loadFactor): " << MyTag << endln;
#endif
  return 0;
}

int 
SimpleContact3D::addInertiaLoadToUnbalance(const Vector &accel)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::addInertiaLoadToUnbalance(const Vector &accel): " << MyTag << endln;
#endif
  return 0;
}

const Vector &
SimpleContact3D::getResistingForce()
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getResistingForce(): " << MyTag << endln;
#endif
        double  t_n;            // normal force
        Vector  t_s(2);         // shear force in xi direction

        int     i;      // loop counter
                
        // initialize F
        internalForces.Zero();

        // get contact stresse/forces
        Vector stress = theMaterial->getStress();

        t_n     = stress(0);

        if (inContact) {        // in contact

            t_s(0)  = stress(1);
            t_s(1)  = stress(2);

            for (i=0; i<SC3D_NUM_DDOF; i++) {
                internalForces(i) = -t_n*Bn(i) + t_s(0)*Bs(i,0) + t_s(1)*Bs(i,1);
            }
            internalForces(SC3D_NUM_DDOF) = -gap;

        } else {
            
            // internalForces(SC3D_NUM_DDOF) = lambda;
            internalForces(SC3D_NUM_DDOF) = t_n;
        }

#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
            opserr << "   Gap    = " << gap << endln;
            opserr << "   Slip   = " << slip ;
            Vector checkDisp = theNodes[5]->getTrialDisp();
            opserr << "   dispL(0) = " << checkDisp(0) << endln;
            opserr << "   lambda   = " << lambda << endln;
            opserr << "   T_n      = " << t_n << endln;
            opserr << "   T_s      = " << t_s;
            }
        if (DEBUG_LEVEL > 2) {
            opserr << "   F_int  = " << internalForces;
            }
#endif

        return internalForces;
}


const Vector &
SimpleContact3D::getResistingForceIncInertia()
{       
#ifdef DEBUG
        opserr << "SimpleContact3D::getResistingForceIncInertia(): " << MyTag << endln;
#endif
  return getResistingForce();
}




int
SimpleContact3D::sendSelf(int commitTag, Channel &theChannel)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::sendSelf(int commitTag, Channel &theChannel)" << endln;
#endif
        
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // SimpleContact3D packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(5);
  data(0) = this->getTag();
  data(1) = SC3D_NUM_DOF;
  data(2) = tolGap;
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
    opserr <<"WARNING SimpleContact3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }           

  // SimpleContact3D then sends the tags of it's four nodes

  res = theChannel.sendID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING SimpleContact3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally SimpleContact3D asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING SimpleContact3D::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
SimpleContact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)" << endln;
#endif

  int res;
  int dataTag = this->getDbTag();

  // SimpleContact3D creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(5);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING SimpleContact3D::recvSelf() - failed to receive Vector\n";
    return -1;
  }           

  this->setTag((int)data(0));
  // SC3D_NUM_DOF = (int)data(1);       // this must not be since SC3D_NUM_DOF is used to initialize variables and thus must not change
  tolGap = data(2);
  
  // SimpleContact3D now receives the tags of it's four external nodes
  res = theChannel.recvID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING SimpleContact3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally SimpleContact3D creates a material object of the correct type,
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
      opserr <<"WARNING SimpleContact3D::recvSelf() - " << this->getTag() 
        << " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING SimpleContact3D::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
SimpleContact3D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::displaySelf(Renderer &theViewer, int displayMode, float fact)" << endln;
#endif
  return 0;
}


void
SimpleContact3D::Print(OPS_Stream &s, int flag)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::Print(OPS_Stream &s, int flag)" << endln;
#endif
        opserr << "SimpleContact3D, element id:  " << this->getTag() << endln;
        opserr << "   Connected external nodes:  " ; 
        for (int i = 0; i<SC3D_NUM_NODE; i++)
        {
                opserr << externalNodes(i)<< " ";
        }
        return;
}


Response*
SimpleContact3D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::setResponse(const char **argv, int argc, Information &eleInfo): " << MyTag << endln;
#endif
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
      return new ElementResponse(this, 1, Vector(3));

    if (strcmp(argv[0],"frictionforce") == 0 || strcmp(argv[0],"frictionforces") == 0)
      return new ElementResponse(this, 2, Vector(3));
    
    else if (strcmp(argv[0],"forcescalar") == 0 || strcmp(argv[0],"forcescalars") == 0)
      return new ElementResponse(this, 3, Vector(3));

    // otherwise response quantity is unknown for the SimpleContact3D class
    else
      return 0;
}


int 
SimpleContact3D::getResponse(int responseID, Information &eleInfo)
{
#ifdef DEBUG
        opserr << "SimpleContact3D::getResponse(int responseID, Information &eleInfo): " << MyTag << endln;
#endif

	  Vector force(3);
      Vector G1(3);
      Vector G2(3);
      
	  // get contact stresse/forces
      Vector stress = theMaterial->getStress();
	  
	  // Calculate contravariant metric tensor G = inv(g)
      double det = g_metric(0,0)*g_metric(1,1) - g_metric(0,1)*g_metric(1,0);

      G_metric(0,0) =  g_metric(1,1);
      G_metric(1,0) = -g_metric(1,0);
      G_metric(0,1) = -g_metric(0,1);
      G_metric(1,1) =  g_metric(0,0);
      G_metric = G_metric / det;

	  // Calculate contravariant base vectors, G1 & G2 from g1 & g2
	  // Gj = G_metric_ji*gi
	  G1 = G_metric(0,0) * g1 + G_metric(0,1) * g2;
	  G2 = G_metric(1,0) * g1 + G_metric(1,1) * g2;

  if (responseID == 1) {
      
	  force = stress(0)*n + stress(1)*G1 + stress(2)*G2;

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
			opserr << "Total force vector = " << force ;
			opserr << "stress vector = " << stress;
		}
#endif

    return eleInfo.setVector(force);

  } else if (responseID == 2) {
	  
	  force = stress(1)*G1 + stress(2)*G2;

#ifdef DEBUG
        if (DEBUG_LEVEL > 0) {
			opserr << "Frictional force vector = " << force ;
			opserr << "stress vector = " << stress;
		}
#endif
		return eleInfo.setVector(force);

  } else if (responseID == 3) {

	  force(0) = stress(0);
	  force(1) = stress(1);
	  force(2) = stress(2);

	  return eleInfo.setVector(force);

  } else

    return -1;
}

int
SimpleContact3D::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"friction") == 0) {
		return param.addObject(1, this);
	}
	
	return -1;
}

int
SimpleContact3D::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes =  theMaterial->updateParameter(parameterID, info);
	if (matRes != -1) {
		res = matRes;
	}
	return res;
}
