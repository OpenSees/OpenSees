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

// Written by: Manish Kumar (mkumar2@buffalo.edu)
// Credits: This element extends the formulation of elastomericBearing element written by Andreas Schellenberg 
// Created: 02/29/2012
//
// Description: This file contains the implementation of the
// ElastomericX class.
//
// What: "@(#) ElastomericX.cpp, revA"

#include "ElastomericX.h"

#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

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
Matrix ElastomericX::theMatrix(12,12);
Vector ElastomericX::theVector(12);
Vector ElastomericX::theLoad(12);

static int numMyBearing = 0;
void *OPS_ElastomericX(void)
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyBearing == 0) {
    opserr << "ElastomericX element - Written by Manish Kumar, University at Buffalo, 2012\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs == 0) { // parallel processing
    theEle = new ElastomericX();
    return theEle;
  }
  
  if (numArgs !=12 && numArgs !=18 && numArgs !=19 && numArgs !=20
	  && numArgs !=21 && numArgs !=22 && numArgs !=23 && numArgs !=24 && numArgs !=25) {
    opserr << "ERROR - ElastomericX incorrect # args provided";
    return theEle;
  }
  
  // get the id and end nodes
  int iData[3];
  double dData[22];
  int numData;					// specify the number of arguments to be read from command line
  // every time an argument is read through OPS_Get.... OPS_GetNumRemainingInputArgs() is increased by one
  
  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  
  int eleTag = iData[0];
  
  numData=9;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  
  // Get the orientation vector
  Vector x(0);
  Vector y(3); y(0)=-1.0; y(1)=0.0; y(2)=0.0;

  // The default values of the parameters
  double kl=10.0;				// Cavitation parameter
  double phi=0.5;				// Damage index
  double al=1.0;				// Strength reduction parameter
  double sDratio=0.5;			// Shear distance ratio
  double m=0.0;					// Mass of the bearing
  double cd1=0;					// Viscous damping parameter
  double tc=0;					// Cover thickness

  if(numArgs>=18) {
    double value;
    x.resize(3);
    numData=1;
    for(int i=0; i<3; i++) {
      if(OPS_GetDoubleInput(&numData, &value) != 0) {
	opserr << "WARNING invalid orientation value for element" << eleTag << endln;
	return 0;
      } else {
	x(i)=value;
      }
    }
    for(int i=0; i<3; i++){
      if (OPS_GetDoubleInput(&numData, &value) != 0) {
	opserr << "WARNING invalid orientation value for element" << eleTag << endln;
	return 0;
      } else {
	y(i)=value;
      }
    }
    if(numArgs>=19){
      numData=1;
      if (OPS_GetDoubleInput(&numData, &kl) != 0) {
	opserr << "WARNING error reading element property cavitation parameter for element" << eleTag << endln;
	return 0;
      }
      if(numArgs>=20){
	numData=1;
	if (OPS_GetDoubleInput(&numData, &phi) != 0) {
	  opserr << "WARNING error reading element property damage index for element" << eleTag << endln;
	  return 0;
	}
	if(numArgs>=21){
	  numData=1;
	  if (OPS_GetDoubleInput(&numData, &al) != 0) {
	    opserr << "WARNING error reading element property strength degradation parameter for element" << eleTag << endln;
	    return 0;
	  }
	  if(numArgs>=22){
	    numData=1;
	    if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
	      opserr << "WARNING error reading element property shear distance ratio for element" << eleTag << endln;
	      return 0;
	    }
	    if(numArgs>=23){
	      numData=1;
	      if (OPS_GetDoubleInput(&numData, &m) != 0) {
		opserr << "WARNING error reading element property mass for element" << eleTag << endln;
		return 0;
	      }
	      if(numArgs>=24){
		numData=1;
		if (OPS_GetDoubleInput(&numData, &cd1) != 0) {
		  opserr << "WARNING error reading element property viscous damping parameter for element" << eleTag << endln;
		  return 0;
		}
		if(numArgs==25){
		  numData=1;
		  if (OPS_GetDoubleInput(&numData, &tc) != 0) {
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
  }
  
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm == 3) {
    
    // check space frame problem has 6 dof per node
    if (ndf != 6)  {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - ElastomericX \n"; 
    }
    
    theEle = new ElastomericX(iData[0], iData[1], iData[2], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], y, x, kl, phi, al, sDratio, m, cd1, tc);
    
  } 
  
  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }
  
  return theEle;
}

  
  ElastomericX::ElastomericX(int tag, int Nd1, int Nd2, double qRubber, double uy, double Gr, double Kbulk, double Di, double Do, 
			     double ts, double tr, int n, const Vector _y, const Vector _x, double kl, double PhiMax, double al, double sDratio, double m, double cd1, double tc)
    :Element(tag, ELE_TAG_ElastomericX), connectedExternalNodes(2), G(Gr), x(_x), y(_y), PhiM(PhiMax), ac(al), shearDistI(sDratio), mass(m), cd(cd1), 
	  L(0.0), D1(Di), D2(Do), ub(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(12), Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6) 
{
	// Ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "ElastomericX::setUp() - element: "
            << this->getTag() << " failed to create an ID of size 2\n";
    }
   
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    // Set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;

	// Vertical motion
	A=(PI/4)*((D2+tc)*(D2+tc)-D1*D1);
	S=(D2-D1)/(4*tr);
	Tr=n*tr;
	h=Tr + (n-1)*ts;
	double F;
	if(D1<DBL_EPSILON) {
		 F=1.0;
	} else {
		double r=D2/D1;											// Outer to inner diameter ratio
		F=(r*r+1)/((r-1)*(r-1)) + (1+r)/((1-r)*log(r));	// Dimension modification factor
	}
	Ec=1.0/((1/(6*G*S*S*F))+(4.0/3.0)*(1/Kbulk));      // Compressive modulus of elastomeric bearing
	double E=3*G;											// Elastic modulus
	double I=(PI/64)*(pow((D2+tc),4)-pow(D1,4));            // Moment of inertia of bearing
	rg=sqrt(I/A);                                           // Radius of gyration 
	Kv0=A*Ec/Tr;											// Pre-cavitation stiffness at zero lateral displacement 
	Kv=Kv0;													// Pre-cavitation stiffness
	Fc=3*G*A;                                               // Cavitation force
	double Er=Ec/3;											// Rotation modulus of bearing
	double As=A*h/Tr;										// Adjusted shear area of bearing
	double Is=I*h/Tr;										// Adjusted moment of intertia of bearing
	double Pe=PI*PI*Er*Is/(h*h);							// Euler buckling load of bearing
	Fcr=-sqrt(Pe*G*As);										// Critical buckling load in compression
	Fcrn=Fcr;												// Initial value of critical buckling load
	Fcrmin=Fcr;											// Initial value of critical buckling load during loading
	ucr=Fcrn/Kv0;											// Initial value of critical buckling deformation
	uc=Fc/Kv0;												// Initial cavitation deformation
	Fcn=Fc;													// Initial value of cavitation deformation
	umax=uc;												// Initial value of maximum tensile deformation

	if (kl<DBL_EPSILON) {
		kc=0.0001;
	} else {
		kc=kl;
	}
	// Horizontal motion
	qYield=qRubber;
	ke=G*A/Tr;
	k0=qRubber/uy;

	// Rotation
	Kr= Er*Is/Tr;

	// Torsion
	Kt=G*(2*Is)/Tr;
	 
    // Initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = k0+ke;
    kbInit(2,2) = k0+ke;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;

    // Initialize variables
    this->revertToStart();

	//cout<<"Fcr: "<< Fcr << "ucr: "<<ucr<<endln;

}

ElastomericX::ElastomericX()
    : Element(0, 0),
    connectedExternalNodes(2),
    k0(0.0), qYield(0.0), ke(0.0), x(0), y(0), shearDistI(0.5),	Kv0(0.0), Kv(0.0), Fc(0.0), Fcr(0.0), Tr(0.0), kc(0.0), 
	PhiM(0.0), ac(0.0), Fcn(0.0), umax(0.0), D1(0.0), D2(0.0), rg(0.0), Ar(0.0), mass(0.0), cd(0.0), L(0.0), 
	tCurrent(0.0), tCommit(0.0), ub(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(6), Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6)
{      
    // ensure the connectedExternalNode ID is of correct size & set values
        if (connectedExternalNodes.Size() != 2)  {
                opserr << "ElastomericX::ElastomericX() - "
                        <<  "failed to create an ID of size 2\n";
                exit(-1);
    }
   
    // set node pointers to NULL
        for (int i=0; i<2; i++)
                theNodes[i] = 0;
}

ElastomericX::~ElastomericX()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
}

int ElastomericX::getNumExternalNodes() const
{
    return 2;
}


const ID& ElastomericX::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** ElastomericX::getNodePtrs()
{
        return theNodes;
}


int ElastomericX::getNumDOF()
{
    return 12;
}


void ElastomericX::setDomain(Domain *theDomain)
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
                        opserr << "WARNING ElastomericX::setDomain() - Nd1: "
                                << connectedExternalNodes(0) << " does not exist in the model for ";
                } else  {
                        opserr << "WARNING ElastomericX::setDomain() - Nd2: "
                                << connectedExternalNodes(1) << " does not exist in the model for ";
                }
                opserr << "ElastomericX ele: " << this->getTag() << endln;
               
                return;
    }
       
        // now determine the number of dof and the dimension    
        int dofNd1 = theNodes[0]->getNumberDOF();
        int dofNd2 = theNodes[1]->getNumberDOF();      
       
        // if differing dof at the ends - print a warning message
    if (dofNd1 != 6)  {
                opserr << "ElastomericX::setDomain() - node 1: "
                        << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
                return;
    }
    if (dofNd2 != 6)  {
                opserr << "ElastomericX::setDomain() - node 2: "
                        << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
                return;
    }
       
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
   
    // set up the transformation matrix for orientation
    this->setUp();
}        


int ElastomericX::commitState()
{
	int errCode = 0;

	double uh=sqrt(ub(1)*ub(1)+ub(2)*ub(2));
	//cout<<"uh: "<<uh<<"rg: "<<rg<<"PI: "<<PI<<endln;
	// Vertical motion
	Kv=Kv0*(1.0/(1.0+(3.0/(PI*PI))*(uh/rg)*(uh/rg)));
	if (uh>DBL_EPSILON) uc=Fc/Kv;
	
	// Tension
	if(ub(0)>umax) {
		umax =ub(0);
		Fcn=Fc*(1-PhiM*(1-exp(-ac*(ub(0)-uc)/uc)));
	}
	// Compression
	double Delta = 2*acos(uh/D2);
	Ar=(D2*D2/4)*(Delta-sin(Delta));
	if(Ar/A>0.2) {
		Fcrn=Fcr*Ar/A;
	} else {
		Fcrn=0.2*Fcr;
	}
	ucr = Fcrn/Kv;
	Fcrmin=max(Fcrmin,Fcrn);

	// Horizontal motion
	//ke=(G*A/Tr)*(1-pow(qb(0)/Fcrn,2));
	//if(ke<0) opserr<<"Negative horizontal stiffness\n";
	// commit trial history variables for horizontal direction
	ubC = ub;
    zC = z;
       
    return errCode;
}


int ElastomericX::revertToLastCommit()
{
	int errCode = 0;
    return errCode;
}


int ElastomericX::revertToStart()
{  
	int errCode=0;

    // reset trial history variables
    ub.Zero();
    z.Zero();
    qb.Zero();
   
    // reset committed history variables
	ubC.Zero();
    zC.Zero();

	// reset tangent of hysteretic evolution parameters
	double A=1.0;
    dzdu(0,0) = dzdu(1,1) = A*k0/qYield;
    dzdu(1,0) = dzdu(0,1) = 0.0;
   
    // reset stiffness matrix in basic system
    kb = kbInit;
   
    return errCode;
}


int ElastomericX::update()
{
	// get global trial displacements and velocities
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
   
    static Vector ug(12), ugdot(12), uldot(12), ubdot(6);
    for (int i=0; i<6; i++)  {
        ug(i)   = dsp1(i);  ugdot(i)   = vel1(i);
        ug(i+6) = dsp2(i);  ugdot(i+6) = vel2(i);
    }
   
    // Transform response from the global to the local system
    ul = Tgl*ug;
    uldot = Tgl*ugdot;
   
    // Transform response from the local to the basic system
    ub = Tlb*ul;
    ubdot = Tlb*uldot;

	//get displacement increments (trial-commited)
	Vector delta_ub = ub - ubC;

	// 1) Get axial force and stiffness in basic x-direction

	ucr = Fcrn/Kv;
	double ucn=Fcn/Kv;
	double Fmax=Fc*(1+(1.0/(Tr*kc))*(1-exp(-kc*(umax-uc))));
		
   if (ub(0)<=ucr) {
		kb(0,0)=Kv/1000;
		qb(0)=Fcrn+kb(0,0)*(ub(0)-ucr);
		//qb(0)=qb(0)+kb(0,0)*delta_ub(0);
		//qb(0)=Fcrn;
			//opserr<<"Elastomer failed in buckling\n";
			//opserr<<"ucr: "<<ucr<<"trialStrain: "<<ub(0)<<"Kv: "<<Kv<<"Fcrn: "<<Fcrn<<endln;
			//exit(-1);
   }
   else if (ub(0)<=ucn) {
		qb(0)=Kv*ub(0);
		kb(0,0)=Kv;
	}
   else if(ub(0)<=umax) {
		qb(0)=Fcn+((Fmax-Fcn)/(umax-ucn))*(ub(0)-ucn);
		kb(0,0)=((Fmax-Fcn)/(umax-ucn));
	}
   else {
		qb(0)=Fc*(1+(1.0/(Tr*kc))*(1-exp(-kc*(ub(0)-uc))));
		kb(0,0)=((Fc/Tr)*exp(-kc*(ub(0)-uc)));
	}
   

	//2) calculate shear forces and stiffnesses in basic y- and z-direction
	
    if (sqrt(pow(delta_ub(1),2)+pow(delta_ub(2),2)) > 0.0)  {
        // get yield displacement
        double uy = qYield/k0;

        // calculate hysteretic evolution parameter z using Newton-Raphson
        int iter = 0;
		int maxIter = 25;
		double tol = 1E-12;
		double A=1;
		double beta=0.1;
		double gamma=0.9;
        double zNrm, tmp1, tmp2, tmp3;
        Vector f(2), delta_z(2);
        Matrix Df(2,2);
        do  {
            zNrm = z.Norm();
            if (zNrm == 0.0)  // check because of negative exponents
                zNrm = DBL_EPSILON;
            tmp1 = z(0)*delta_ub(1) + z(1)*delta_ub(2);
            tmp2 = beta + gamma*sgn(tmp1);
            tmp3 = tmp1*tmp2;
            //cout<<"iter: "<<iter<<" z(0): "<<z(0)<<" z(1): "<<z(1)<<" zC(0): "<<zC(0)<<" zC(1): "<<zC(1)<<endln;
            // function and derivative
            f(0) = z(0) - zC(0) - 1.0/uy*(A*delta_ub(1) - z(0)*tmp3);
            f(1) = z(1) - zC(1) - 1.0/uy*(A*delta_ub(2) - z(1)*tmp3);
            
            Df(0,0) = 1.0 + (tmp2/uy)*(2*z(0)*delta_ub(1)+z(1)*delta_ub(2));
            Df(1,0) = (tmp2/uy)*z(1)*delta_ub(1);
            Df(0,1) = (tmp2/uy)*z(0)*delta_ub(2);
            Df(1,1) = 1.0 + (tmp2/uy)*(2*z(1)*delta_ub(2)+z(0)*delta_ub(1));
            
            // issue warning if diagonal of derivative Df is zero
            if ((fabs(Df(0,0)) <= DBL_EPSILON) || (fabs(Df(1,1)) <= DBL_EPSILON))  {
                opserr << "WARNING: ElastomericBearingBoucWen3d::update() - "
                    << "zero Jacobian in Newton-Raphson scheme for hysteretic "
                    << "evolution parameter z.\n";
                return -1;
            }
            
            // advance one step
			// delta_z = f/Df; either write a function to do matrix devision or use the solution below
            delta_z(0) = (f(0)*Df(1,1)-f(1)*Df(0,1))/(Df(0,0)*Df(1,1)-Df(0,1)*Df(1,0));
			delta_z(1) = (f(0)*Df(1,0)-f(1)*Df(0,0))/(Df(0,1)*Df(1,0)-Df(0,0)*Df(1,1));
            z -= delta_z;
            iter++;
        } while ((delta_z.Norm() >= tol) && (iter < maxIter));
        
        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter)   {
            opserr << "WARNING: ElastomericBearingBoucWen3d::update() - "
                << "did not find the hysteretic evolution parameters z after "
                << iter << " iterations and norm: " << delta_z.Norm() << endln;
            return -2;
        }
        
        // get derivative of hysteretic evolution parameter
        delta_z = z-zC;
        if (fabs(delta_ub(1)) > DBL_EPSILON)  {
            dzdu(0,0) = delta_z(0)/delta_ub(1);
            dzdu(1,0) = delta_z(1)/delta_ub(1);
        }
        if (fabs(delta_ub(2)) > DBL_EPSILON)  {
            dzdu(0,1) = delta_z(0)/delta_ub(2);
            dzdu(1,1) = delta_z(1)/delta_ub(2);
        }

        tCurrent=(this->getDomain())->getCurrentTime();
		if(tCurrent<tCommit) {
			tCommit=0.0;
		}
		double dT=tCurrent-tCommit;

		// set shear force
        qb(1) = cd*ubdot(1) + qYield*z(0) + ke*ub(1);
        qb(2) = cd*ubdot(2) + qYield*z(1) + ke*ub(2);
        // set tangent stiffness
        kb(1,1) = cd/dT+qYield*dzdu(0,0) + ke;
        kb(1,2) = qYield*dzdu(0,1);
        kb(2,1) = qYield*dzdu(1,0);
        kb(2,2) = cd/dT+qYield*dzdu(1,1) + ke;
		//cout<<" L "<< L<<" kc: "<< kc<<" phiM "<<PhiM<<" cd:"<<cd<<endl;
    }
   
    // 3) get moment and stiffness in basic x-direction
    qb(3) = Kt*ub(3);
    kb(3,3) = Kt;
   
    // 4) get moment and stiffness in basic y-direction
    qb(4) = Kr*ub(4);
    kb(4,4) = Kr;
   
    // 5) get moment and stiffness in basic z-direction
    qb(5) = Kt*ub(5);
    kb(5,5) = Kr;
    return 0;
}


const Matrix& ElastomericX::getTangentStiff()
{

	// zero the matrix
    theMatrix.Zero();
   
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
   
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
   
    return theMatrix;

}


const Matrix& ElastomericX::getInitialStiff()
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


const Matrix& ElastomericX::getMass()
{
        // zero the matrix
    theMatrix.Zero();
   
        // check for quick return
        if (mass == 0.0)  {
                return theMatrix;
        }    
   
        double m = 0.5*mass;
        for (int i = 0; i < 3; i++)  {
                theMatrix(i,i)     = m;
                theMatrix(i+3,i+3) = m;
        }
       
    return theMatrix;
}


void ElastomericX::zeroLoad()
{
    theLoad.Zero();
}


int ElastomericX::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
        opserr <<"ElastomericX::addLoad() - "
                << "load type unknown for element: "
                << this->getTag() << endln;
   
        return -1;
}


int ElastomericX::addInertiaLoadToUnbalance(const Vector &accel)
{
	// check for quick return
        if (mass == 0.0)  {
                return 0;
        }
   
        // get R * accel from the nodes
        const Vector &Raccel1 = theNodes[0]->getRV(accel);
        const Vector &Raccel2 = theNodes[1]->getRV(accel);
       
        if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
                opserr << "ElastomericX::addInertiaLoadToUnbalance() - "
                        << "matrix and vector sizes are incompatible\n";
                return -1;
        }
   
        // want to add ( - fact * M R * accel ) to unbalance
        // take advantage of lumped mass matrix
        double m = 0.5*mass;
    for (int i = 0; i < 3; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+3) -= m * Raccel2(i);
    }
   
        return 0;
}


const Vector& ElastomericX::getResistingForce()
{
	// zero the residual
    theVector.Zero();
   
    // determine resisting forces in local system
    static Vector ql(12);
    ql = Tlb^qb;

    // determine resisting forces in global system
    theVector = Tgl^ql;
   
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    return theVector;
}


const Vector& ElastomericX::getResistingForceIncInertia()
{      
	theVector = this->getResistingForce();
       
        // add the damping forces if rayleigh damping
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
                theVector += this->getRayleighDampingForces();
   
        // now include the mass portion
        if (mass != 0.0)  {
                const Vector &accel1 = theNodes[0]->getTrialAccel();
                const Vector &accel2 = theNodes[1]->getTrialAccel();    
               
                double m = 0.5*mass;
                for (int i = 0; i < 3; i++)  {
                        theVector(i)   += m * accel1(i);
                        theVector(i+3) += m * accel2(i);
                }
        }
       
        return theVector;
}


int ElastomericX::sendSelf(int commitTag, Channel &sChannel)
{
	// send element parameters
    static Vector data(20);
    data(0) = this->getTag();
    data(1) = k0;
    data(2) = qYield;
    data(3) = ke;
	data(4) = Kv0;
	data(5) = kc;
	data(6) = PhiM;
	data(7) = Fcr;
	data(8) = Fc;
	data(9) = Kt;
	data(10) = Kr;
	data(11) = shearDistI;
    data(12) = mass;
    data(13) = x.Size();
    data(14) = y.Size();
	data(15) = Tr;
	data(16) = D1;
	data(17) = D2;
	data(18) = L;
	data(19) = rg;

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


int ElastomericX::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
   
    // receive element parameters
    static Vector data(23);
    rChannel.recvVector(0, commitTag, data);    
    this->setTag((int)data(0));
    k0 = data(1);
    qYield = data(2);
    ke = data(3);
	Kv0 = data(4);
	kc = data(5);
	PhiM = data(6);
	Fcr = data(7);
	Fc = data(8);
	Kt = data(9);
	Kr = data(10);
	shearDistI = data(11);
    mass = data(12);
	Tr = data(15);
	D1 = data(16);
	D2 = data(17);
	L = data(18);
	rg = data(19);
   
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
   
   
    // receive remaining data
    if ((int)data(14) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(15) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }


    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = ke + k0;
    kbInit(2,2) = ke + k0;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
   
    // initialize variables
    this->revertToStart();
   
    return 0;
}


int ElastomericX::displaySelf(Renderer &theViewer,
    int displayMode, float fact)
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

    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}


void ElastomericX::Print(OPS_Stream &s, int flag)
{
	if (flag == 0)  {
        // print everything
        s << "************************************************************" << endln;
        s << "Element: " << this->getTag();
        s << "  type: ElastomericX  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
		s << "************************************************************" << endln;
		s << "GEOMETRIC PROPERTIES" << endln;
		s << "D1: "<< D1 << " D2: "<< D2 <<" L: "<< L<<" Tr: "<< Tr <<" S: "<< S <<" A: "<< A <<endln;
		s << "MATERIAL PROPERTIES" << endln;
		s << "kc: "<< kc <<" ac: "<<ac<<" PhiM: "<< PhiM <<"  shearDistI: " << shearDistI << "  mass: " << mass<<endln;
		s << "MECHANICAL PROPERTIES: HORIZONTAL MOTION\n"<<endln;
		s << "k0: " << k0 << "  ke: " << ke << " qYield: " << qYield  << " Fcrmin: "<< Fcrmin <<endln;
		s << "MECHANICAL PROPERTIES: VERTICAL MOTION"<<endln;
		s << "Ec: "<<Ec<<" Kv0: "<< Kv0<<" Kv: "<< Kv <<" uc: "<< uc <<" Fcr: "<< Fcr <<" Fcrn: "<< Fcrn <<" ucr: "<< ucr <<" umax: "<< umax <<endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
		s << "************************************************************" << endln;
    } else if (flag == 1)  {
                // does nothing
    }
}


Response* ElastomericX::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
   
    output.tag("ElementOutput");
    output.attr("eleType","ElastomericX");
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
    output.endTag(); // ElementOutput
    return theResponse;
}


int ElastomericX::getResponse(int responseID, Information &eleInfo)
{
	double kGeo1, MpDelta1, MpDelta2, MpDelta3, MpDelta4, MpDelta5, MpDelta6;
   
    switch (responseID)  {
        case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
       
        case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector = Tlb^qb;
        return eleInfo.setVector(theVector);
       
        case 3:  // basic forces
        return eleInfo.setVector(qb);
       
        case 4:  // local displacements
        return eleInfo.setVector(ul);
       
        case 5:  // basic displacements
        return eleInfo.setVector(ub);
       
    default:
                return -1;
        }
}


// set up the transformation matrix for orientation
void ElastomericX::setUp()
{
	const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();    
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();

    if (L > DBL_EPSILON)  {
                if (x.Size() == 0)  {
                    x.resize(3);
                    x = xp;
        } else  {
            opserr << "WARNING ElastomericX::setUp() - "
                << "element: " << this->getTag() << endln
                << "ignoring nodes and using specified "
                << "local x vector to determine orientation\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "ElastomericX::setUp() - "
            << "element: " << this->getTag() << endln
            << "incorrect dimension of orientation vectors\n";
        exit(-1);
    }

	// establish orientation of element for the tranformation matrix
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
        opserr << "ElastomericX::setUp() - "
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


double ElastomericX::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
