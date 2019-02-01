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

// Written: Manish Kumar (mkumar3@buffalo.edu)
// Created: 2014/12/22
//
// Description: This file contains the implementation of the
// FPBearingPTV class. This class is similar to the existing
// singleFPBearing class. The key difference lies in 
// consideration of dependence of coefficient of friction on
// instantaneous values of sliding velocity, axial pressure on 
// the bearing, and temperature at the sliding surface.



#include "FPBearingPTV.h"
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

// initialize the class wide variables
Matrix FPBearingPTV::theMatrix(12,12);
Vector FPBearingPTV::theVector(12);
Vector FPBearingPTV::theLoad(12);

static int numMyBearing = 0;
void *OPS_FPBearingPTV()
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyBearing == 0) {
    opserr << "FPBearingPTV element - Written by Manish Kumar, University at Buffalo Copyright 2013\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new FPBearingPTV();
    return theEle;
  }

  if (numRemainingArgs < 30) {
    opserr << "ERROR - FPBearingPTV incorrect # args provided";
    return theEle;
  } 
  
  int iData[13];
  double dData[17];
  int numData;  
    
  // get the id and end nodes
  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];
  int iNode = iData[1];
  int jNode = iData[2];
  
  //MuRef
  double aDouble;
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;	
  }
  dData[0]=aDouble;
  double MuRef = dData[0];

  //IsPressureDependent
  numData = 1;
  int aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[3]=aInt;
  int IsPressureDependent = iData[3];

  //Reference Pressure
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  dData[1]=aDouble;
  double pRef = dData[1];

  //IsTemperatureDependent
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[4]=aInt;
  int IsTemperatureDependent = iData[4];

  //Thermal diffusivity of Steel
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  dData[2]=aDouble;
  double Diffusivity = dData[2];

  //Thermal conductivity of Steel
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  dData[3]=aDouble;
  double Conductivity = dData[3];

  //IsVelocityDependent
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[5]=aInt;
  int IsVelocityDependent = iData[5];

  //rateParameter
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[4]=aDouble;
  double rateParameter = dData[4];

  //Effective radius of curvature of sliding surface
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[5]=aDouble;
  double ReffectiveFP = dData[5];

  //Radius of contact area
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[6]=aDouble;
  double Radius_Contact = dData[6];

  //kInitial
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[7]=aDouble;
  double kInitial = dData[7];

  //Material tags
  numData = 1;
  int aNumber;
  for(int i=0; i<4; i++){
	if (OPS_GetIntInput(&numData, &aNumber) != 0) {
		opserr << "WARNING invalid material information\n";
		return 0;
	} else {
	  iData[6+i]=aNumber;
	}
  }
   
  int matIDP = iData[6];
  int matIDT = iData[7];
  int matIDMy = iData[8];
  int matIDMz = iData[9]; 

  UniaxialMaterial *theMaterialA = OPS_GetUniaxialMaterial(matIDP);
  UniaxialMaterial *theMaterialB = OPS_GetUniaxialMaterial(matIDT);
  UniaxialMaterial *theMaterialC = OPS_GetUniaxialMaterial(matIDMy);
  UniaxialMaterial *theMaterialD = OPS_GetUniaxialMaterial(matIDMz);

  //Orientation vector
  Vector x(3); 
  Vector y(3); 
  //x(0)=0.0; x(1)=0.0; x(2)=1.0;
  //y(0)=1.0; y(1)=0.0; y(2)=0.0;
  numData = 1;  
  for(int i=0; i<6; i++){
	  if (OPS_GetDoubleInput(&numData, &aDouble) !=0) {
		  opserr << "WARNING invalid element data\n";
		  return 0;
	  } else {
		  dData[8+i]=aDouble;
	  }
  }
  x(0)=dData[8]; x(1)=dData[9]; x(2)=dData[10];
  y(0)=dData[11]; y(1)=dData[12]; y(2)=dData[13];

  //Shear distance
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING invalid element data" << eleTag << endln;
    return 0;
  }
  dData[14]=aDouble;
  double shearDist = dData[14];

  //doRayleigh
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[10]=aInt;
  int doRayleigh = iData[10];

  //mass
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  dData[15]=aDouble;
  double mass = dData[15];

  //NumberOfIterations
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[11]=aInt;
  int iter = iData[11];

  //tol
  numData = 1;  
  if (OPS_GetDoubleInput(&numData, &aDouble) != 0) {
    opserr << "WARNING error reading element properties for element" << eleTag << endln;
    return 0;
  }
  dData[16]=aDouble;
  double tol = dData[16];

  //Units 1: N,m,s,C; 2: kN,m,s,C; 3: N,mm,s,C; 4: kN,mm,s,C; 5: lb,in,s,C; 6: kip,in,s,C; 7: lb,ft,s,C; 8: kip,ft,s,C
  numData = 1;
  aInt = 0;
  if (OPS_GetIntInput(&numData, &aInt) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
  iData[12]=aInt;
  int unit = iData[12];

  
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();
  if (ndm == 3) {

	  // check space frame problem has 6 dof per node
        if (ndf != 6)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for space problem need 6 - FPBearingPTV \n"; 
		}
			 
    theEle = new FPBearingPTV(eleTag, iNode, jNode, MuRef, IsPressureDependent, pRef, IsTemperatureDependent, Diffusivity, Conductivity, IsVelocityDependent, rateParameter, ReffectiveFP, Radius_Contact, kInitial, *theMaterialA, *theMaterialB, *theMaterialC, *theMaterialD, x, y, shearDist, doRayleigh, mass, iter, tol, unit); 
	
  } 
  

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
	delete theMaterialA;
	delete theMaterialB;
	delete theMaterialC;
	delete theMaterialD;
    return 0;
  }
 
  return theEle;
    
}

//Typical constructor

FPBearingPTV::FPBearingPTV(int tag, int Nd1, int Nd2, double MuReference, int IsPDependent, 
	double refP, int IsTDependent, double Diffusivity_Steel, double Conductivity_Steel, 
	int IsVDependent, double rate_v_mu, double Reff, double r_Contact, double kInit, 
	UniaxialMaterial &theMatA, UniaxialMaterial &theMatB,
	UniaxialMaterial &theMatC, UniaxialMaterial &theMatD, 
	const Vector _x, const Vector _y, double sdI, int addRay,
    double m, int maxiter, double _tol, int _unit)
    : Element(tag, ELE_TAG_FPBearingPTV),
    connectedExternalNodes(2), muRef(MuReference), 
	kpFactor(IsPDependent), refPressure(refP), 
	kTFactor(IsTDependent), diffuse(Diffusivity_Steel),
	conduct(Conductivity_Steel), kvFactor(IsVDependent),
	rateParam(rate_v_mu), Reffective(Reff), rContact(r_Contact),
	k0(kInit), 
	x(_x), y(_y),
    shearDistI(sdI), addRayleigh(addRay),
    mass(m), maxIter(maxiter), tol(_tol), unit(_unit),
    L(0.0), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6), DomainTime(1), DomainHeatFlux(1), 
	DomainTimeTemp(0), DomainHeatFluxTemp(0), iCountTime(0), kpFTemp(1), kTFTemp(1), kvFTemp(1),
	DomainDisp(2,3), TemperatureCenter(1), MuFactors(3), MuAdjusted(1), HeatFluxCenter(1)
{
	// get a copy of the material object for our own use
	theMaterials[0] = theMatA.getCopy();
	theMaterials[1] = theMatB.getCopy();
	theMaterials[2] = theMatC.getCopy();
	theMaterials[3] = theMatD.getCopy();	

	// check material input
    if (theMaterials[0] == 0)  {
        opserr << "FPBearingPTV::FPBearingPTV() - "
            << "null material array passed.\n";
        exit(-1);
    }
	if (theMaterials[1] == 0)  {
        opserr << "FPBearingPTV::FPBearingPTV() - "
            << "null material array passed.\n";
        exit(-1);
    }
	if (theMaterials[2] == 0)  {
        opserr << "FPBearingPTV::FPBearingPTV() - "
            << "null material array passed.\n";
        exit(-1);
    }
	if (theMaterials[3] == 0)  {
        opserr << "FPBearingPTV::FPBearingPTV() - "
            << "null material array passed.\n";
        exit(-1);
    }	

	// ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "FPBearingPTV::FPBearingPTV() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;	
	
	// set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
       
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = k0;
    kbInit(2,2) = k0;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
		  	
	//Initialize DomainTime
	DomainTime(0) = 0.0;
	DomainHeatFlux(0) = 0.0;	
	int iCountTime = 0;

    // initialize other variables
    this->revertToStart();	
}

// constructor which should be invoked by an FE_ObjectBroker only

FPBearingPTV::FPBearingPTV()
    : Element(0, ELE_TAG_FPBearingPTV),
    connectedExternalNodes(2), muRef(0.0), 
	kpFactor(0), refPressure(0.0), 
	kTFactor(0), diffuse(0.0),
	conduct(0.0), kvFactor(0),
	rateParam(0.0), Reffective(0.0), rContact(0.0),
	k0(0.0), 
	x(0), y(0),
    shearDistI(0.0), addRayleigh(0),
    mass(0.0), maxIter(25), tol(1E-12), unit(0),
    L(0.0), ub(6), ubPlastic(2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubPlasticC(2), kbInit(6,6), DomainTime(1), DomainHeatFlux(1), 
	DomainTimeTemp(0), DomainHeatFluxTemp(0), iCountTime(0), kpFTemp(1), kTFTemp(1), kvFTemp(1),
	DomainDisp(2,3), TemperatureCenter(1), MuFactors(3), MuAdjusted(1), HeatFluxCenter(1)
{    

	// ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)  {
        opserr << "FPBearingPTV::FPBearingPTV() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }
    
    // set node pointers to NULL
    for (int i=0; i<2; i++)
        theNodes[i] = 0;
	
    
    // set material pointers to NULL    
	theMaterials[0] = 0;
	theMaterials[1] = 0;
	theMaterials[2] = 0;
	theMaterials[3] = 0;	

	DomainTime(0) = 0.0;
	DomainHeatFlux(0) = 0.0;	
	int iCountTime = 0;
}


FPBearingPTV::~FPBearingPTV()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to    
    
    if (theMaterials[0] != 0)
       delete theMaterials[0];
	if (theMaterials[1] != 0)
       delete theMaterials[1];
	if (theMaterials[2] != 0)
       delete theMaterials[2];
	if (theMaterials[3] != 0)
       delete theMaterials[3];		
}


int FPBearingPTV::getNumExternalNodes() const
{
    return 2;	
}


const ID& FPBearingPTV::getExternalNodes() 
{
    return connectedExternalNodes;
}


Node** FPBearingPTV::getNodePtrs() 
{
    return theNodes;		
}



int FPBearingPTV::getNumDOF() 
{    
	return 12;	
}


void FPBearingPTV::setDomain(Domain *theDomain)
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
            opserr << "WARNING FPBearingPTV::setDomain() - Nd1: " 
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING FPBearingPTV::setDomain() - Nd2: " 
                << connectedExternalNodes(1)
                << " does not exist in the model for";
        }
        opserr << " element: " << this->getTag() << ".\n";        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != 6)  {
        opserr << "FPBearingPTV::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "FPBearingPTV::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
	    
    // set up the transformation matrix for orientation
    this->setUp();
}



int FPBearingPTV::commitState()
{	
    int errCode = 0;
    
    // commit trial history variables
    ubPlasticC = ubPlastic;
    
    // commit material models    
	errCode += theMaterials[0]->commitState();
	errCode += theMaterials[1]->commitState();
	errCode += theMaterials[2]->commitState();
	errCode += theMaterials[3]->commitState();	
    
    // commit the base class
    errCode += this->Element::commitState();  

    return errCode;
}


int FPBearingPTV::revertToLastCommit()
{
    int errCode = 0;    
    
    // revert material models    
	errCode += theMaterials[0]->revertToLastCommit();
	errCode += theMaterials[1]->revertToLastCommit();
	errCode += theMaterials[2]->revertToLastCommit();
	errCode += theMaterials[3]->revertToLastCommit();    
    return errCode;
}



int FPBearingPTV::revertToStart()
{
    int errCode=0;
    
    // reset trial history variables
    ub.Zero();
    ubPlastic.Zero();
    qb.Zero();
    
    // reset committed history variables
    ubPlasticC.Zero();
    
    // reset stiffness matrix in basic system
    kb = kbInit;    
    
    // revert material models    
	errCode += theMaterials[0]->revertToStart();
	errCode += theMaterials[1]->revertToStart();
	errCode += theMaterials[2]->revertToStart();
	errCode += theMaterials[3]->revertToStart();	
    
    return errCode;
}


int FPBearingPTV::update()
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
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
	    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);

	// Compute resultant displacement in horizontal direction
	double resultantU =  sqrt(ub(1)*ub(1)+ub(2)*ub(2));
    
	// get radii in basic y- and z- direction
	double Ry = sqrt(pow(Reffective,2)-pow(ub(1),2));
	double Rz = sqrt(pow(Reffective,2)-pow(ub(2),2));		
	
    // get absolute velocity
    //double ubdotAbs = sqrt(pow(ubdot(1),2) + pow(ubdot(2),2));		
	double ubdotAbs = sqrt(pow(ubdot(1)/Ry*ub(1) + ubdot(2)/Rz*ub(2),2) + pow(ubdot(1),2) + pow(ubdot(2),2));
    
    // 1) get axial force and stiffness in basic x-direction	
    double ub0Old = theMaterials[0]->getStrain();    
    theMaterials[0]->setTrialStrain(ub(0), ubdot(0));    
    qb(0) = theMaterials[0]->getStress();
    kb(0,0) = theMaterials[0]->getTangent();

	//Check for uplift	
	double IsUplift = 0.0;
	double kFactUplift = 1.0E-6;
	if (qb(0) >= 0.0) {
		qb(0) = -0.0001*qb(0);
		IsUplift = 1.0;
	}        
    if (qb(0) >= 0.0)  {
        //Do nothing		
        if (qb(0) > 0.0)  {
            //Do nothing		
        }		 			
        return 0;		
    }		
	
	//Factors to account for pressure, temperature and velocity dependences of friction
	double kpF; double kTF; double kvF;
	double Mu_Adj;	
	if (iCountTime == 0) {
		kpF=1.0;	kTF=1.0;	kvF=1.0;
		Mu_Adj=muRef*kpF*kTF*kvF;
	} else {
		kpF=kpFTemp(0); kTF=kTFTemp(0); kvF=kvFTemp(0);
		Mu_Adj=muRef*kpF*kTF*kvF;
	}	
			
	//Unit conversion for pressure to be used in the pressure factor computation
	double p_Unit_Convert; //To convert to MPa
	if (unit == 1) { p_Unit_Convert=0.000001; }
	if (unit == 2) { p_Unit_Convert=0.001; }
	if (unit == 3) { p_Unit_Convert=1.0; }
	if (unit == 4) { p_Unit_Convert=1000.0; }
	if (unit == 5) { p_Unit_Convert=0.006894; }
	if (unit == 6) { p_Unit_Convert=6.894; }
	if (unit == 7) { p_Unit_Convert=0.00004788; }
	if (unit == 8) { p_Unit_Convert=0.04788; }

	// get normal forces
	double N = -qb(0);
    	
	//Instantaneous pressure
	double Inst_Pressure=fabs(N)/(22.0*rContact*rContact/7.0);	

	//Current domain time
	double tCurrent=(this->getDomain())->getCurrentTime();
	double Disp_CurrentStep;
	double Vel_CurrentStep;	

	//Compute factors for pressure, temperature and velocity dependencies	
	if (tCurrent > DomainTime(iCountTime-1)) {
		//Initiate the DomainTimeVector
		if (iCountTime == 0) {
				DomainTime.resize(2000);								
				for (int iTemp=1; iTemp <=2000; iTemp++) {
					DomainTime(iTemp-1)=0.0;										
				}
			}	
		DomainTime(iCountTime)=tCurrent;			
		//Resize and repopulate the DomainTime Vector
		for (int iResize=1; iResize<=1000; iResize++) { //Vectors can be as long as 2000*1000
			if (iCountTime == (iResize*2000-3)) {
				DomainTimeTemp.resize(iResize*2000-2);									
				//Populate DomainTimeTemp 
				for (int iTemp=1; iTemp<=iResize*2000-2; iTemp++) {
					DomainTimeTemp(iTemp-1)=DomainTime(iTemp-1);											
				}
				DomainTime.resize((iResize+1)*2000);									
				for (int iTemp=1; iTemp<=(iResize+1)*2000; iTemp++) {
					DomainTime(iTemp-1)=0.0;										
				}
				//Populate DomainTime 
				for (int iTemp=1; iTemp<=iResize*2000-2; iTemp++) {
					DomainTime(iTemp-1)=DomainTimeTemp(iTemp-1);											
				}
				DomainTimeTemp.resize(0);									
				break;
			}
		}
		//Displacement in the previous step
		if (iCountTime == 0) {
			DomainDisp(0,0)=DomainDisp(0,1)=DomainDisp(0,2)=0.0;
			DomainDisp(1,0)=ub(0); DomainDisp(1,1)=ub(1); DomainDisp(1,2)=ub(2);
		} else {
			DomainDisp(0,0)=DomainDisp(1,0); 
			DomainDisp(0,1)=DomainDisp(1,1);
			DomainDisp(0,2)=DomainDisp(1,2);
			DomainDisp(1,0)=ub(0); DomainDisp(1,1)=ub(1); DomainDisp(1,2)=ub(2);
		}		
		Disp_CurrentStep = sqrt(pow((DomainDisp(1,1)-DomainDisp(0,1)),2)+pow((DomainDisp(1,2)-DomainDisp(0,2)),2));
		//Sliding velocity in the previous step (ubDotAbs will be 0 in static analysis)		
		if (iCountTime == 0) {
			Vel_CurrentStep=Disp_CurrentStep/tCurrent;
		} else {
			Vel_CurrentStep=Disp_CurrentStep/(DomainTime(iCountTime)-DomainTime(iCountTime-1));					
		}			
		//set Vel_CurrentStep=ubdotAbs for dynamic analysis
		if (ubdotAbs > 0.0) { Vel_CurrentStep=ubdotAbs; }			
		//Temperature factor
		if (fabs(kTFactor-1.0) <=10.00001) {
			//Initialize DomainHeatFlux
			if (iCountTime == 0) {				
				DomainHeatFlux.resize(2000);				
				for (int iTemp=1; iTemp <=2000; iTemp++) {					
					DomainHeatFlux(iTemp-1)=0.0;					
				}
			}									
			//Resize and repopulate the DomainHeatFlux Vector
			for (int iResize=1; iResize<=1000; iResize++) { //Vectors can be as long as 2000*1000
				if (iCountTime == (iResize*2000-3)) {					
					DomainHeatFluxTemp.resize(iResize*2000-2);					
					//Populate DomainHeatFluxTemp
					for (int iTemp=1; iTemp<=iResize*2000-2; iTemp++) {						
						DomainHeatFluxTemp(iTemp-1)=DomainHeatFlux(iTemp-1);						
					}					
					DomainHeatFlux.resize((iResize+1)*2000);					
					for (int iTemp=1; iTemp<=(iResize+1)*2000; iTemp++) {						
						DomainHeatFlux(iTemp-1)=0.0;						
					}
					//Populate DomainHeatFlux
					for (int iTemp=1; iTemp<=iResize*2000-2; iTemp++) {						
						DomainHeatFlux(iTemp-1)=DomainHeatFluxTemp(iTemp-1);						
					}					
					DomainHeatFluxTemp.resize(0);					
					break;
				}
			}			
			//Heat flux			
			if (resultantU < (rContact*sqrt(22.0/(7.0*4.0)))) {
				DomainHeatFlux(iCountTime)=Mu_Adj*(fabs(N)/(22.0*rContact*rContact/7.0))*Vel_CurrentStep;
			} else DomainHeatFlux(iCountTime)=0.0;	
			HeatFluxCenter(0) = DomainHeatFlux(iCountTime);
			//Temperature at the sliding surface	
			double Temperature_Surface;
			double Temperature_Change = 0.0;	

			if (iCountTime > 1 && resultantU > 0.0) { //When there is some lateral displacement
				double DtAnalysis=0.000001;
				double tau=0.000001;
				for (int jTemp=0; jTemp<=iCountTime; jTemp++) {										
					if (iCountTime == 0) {
						DtAnalysis = DomainTime(iCountTime);
						tau = DtAnalysis;					
					} else {
						DtAnalysis = DomainTime(iCountTime) - DomainTime(iCountTime-1);					
						tau = DomainTime(jTemp);						
					}						
					Temperature_Change=Temperature_Change+sqrt(diffuse*7.0/22.0)*DomainHeatFlux(iCountTime-jTemp)*DtAnalysis/(conduct*sqrt(tau));					
				}				
			}
			Temperature_Surface = 20.0+Temperature_Change;	
			TemperatureCenter(0) = Temperature_Surface;
			if (fabs(kTFactor-1.0) <= 0.00001) { kTF = 0.789*(pow(0.7, (Temperature_Surface/50.0))+0.4); }			
		} 	
					
		//Pressure factor
		if (fabs(kpFactor-1.0) <= 0.00001) { kpF=pow(0.7, ((Inst_Pressure-refPressure)*p_Unit_Convert/50.0)); } 
		//Velocity factor
		if (fabs(kvFactor-1.0) <= 0.00001) { kvF=1-0.5*exp(-rateParam * Vel_CurrentStep);	}
		//Coefficients for all the stages of an iteration during a time-step
		kpFTemp(0)=kpF; kTFTemp(0)=kTF; kvFTemp(0)=kvF; 
		//Update the counter
		iCountTime=iCountTime+1;
	}
	//When domain time is reset
	if (tCurrent < DomainTime(iCountTime-1)) {
		iCountTime=0;		
		DomainTime(iCountTime)=tCurrent;
		iCountTime=iCountTime+1;
		kpF=1.0; kTF=1.0; kvF=1.0;
		Mu_Adj=muRef*kpF*kTF*kvF;
		kpFTemp(0)=kpF; kTFTemp(0)=kTF; kvFTemp(0)=kvF;			
	}
	//Factors used during iterations to reach the next time-step
	if (fabs(tCurrent-DomainTime(iCountTime-1)) < 0.00000001) {
		kpF=kpFTemp(0); kTF=kTFTemp(0); kvF=kvFTemp(0);
	}
	Mu_Adj=muRef*kpF*kTF*kvF;
	MuFactors(0) = kpF;	MuFactors(1) = kTF;	MuFactors(2) = kvF;
	MuAdjusted(0) = Mu_Adj;	
	
    // 2) calculate shear forces and stiffnesses in basic y- and z-direction
    int iter = 0;
    Vector qbOld(2);
    do  {
        // save old shear forces
        qbOld(0) = qb(1);
        qbOld(1) = qb(2);
        
        // get normal forces			
		N = -qb(0);
	
		//Yield force 		
		double qYield = Mu_Adj*N;
		
        // get stiffnesses of elastic components
		double k2y = N/Ry;
		double k2z = N/Rz;		

		// get initial stiffnesses of hysteretic components
		double k0y = k0 - k2y;
		double k0z = k0 - k2z;   

		//Account for uplift
		kFactUplift = 1.0E-6;
		if (fabs(IsUplift - 1.0) <= 0.00001) {
			k0y = kFactUplift*k0y;
			k0z = kFactUplift*k0z;
		}
				
		// get trial shear forces of hysteretic component
        Vector qTrial(2);
        qTrial(0) = k0y*(ub(1) - ubPlasticC(0));
        qTrial(1) = k0z*(ub(2) - ubPlasticC(1));
        
        // compute yield criterion of hysteretic component
        double qTrialNorm = qTrial.Norm();
        double Y = qTrialNorm - qYield;
        
        // elastic step -> no updates required
        if (Y <= 0.0)  {
            // set shear forces
            qb(1) = qTrial(0) - N*ul(5) + k2y*ub(1);
            qb(2) = qTrial(1) + N*ul(4) + k2z*ub(2);
            // set tangent stiffnesses
            kb(1,1) = kb(2,2) = k0;
            kb(1,2) = kb(2,1) = 0.0;	
			//Account for uplift
			if (fabs(IsUplift - 1.0) <= 0.00001) {
				kb(1,1) = kFactUplift*k0;
				kb(2,2) = kFactUplift*k0;
			}
        }
        // plastic step -> return mapping
        else  {
            // compute consistency parameters
            double dGammaY = Y/k0y;
			double dGammaZ = Y/k0z;
            // update plastic displacements
            ubPlastic(0) = ubPlasticC(0) + dGammaY*qTrial(0)/qTrialNorm;
            ubPlastic(1) = ubPlasticC(1) + dGammaZ*qTrial(1)/qTrialNorm;
            // set shear forces
            qb(1) = qYield*qTrial(0)/qTrialNorm - N*ul(5) + k2y*ub(1);
            qb(2) = qYield*qTrial(1)/qTrialNorm + N*ul(4) + k2z*ub(2);			
            // set tangent stiffnesses
            //Account for uplift
			double k0Plastic = k0;
			if (fabs(IsUplift - 1.0) <= 0.00001) {
				k0Plastic = kFactUplift*k0;				
			}
            double D = pow(qTrialNorm,3);
            kb(1,1) =  k0Plastic*(1.0-(qYield*qTrial(1)*qTrial(1)/D)) + k2y;
            kb(1,2) = -qYield*k0Plastic*qTrial(0)*qTrial(1)/D;
            kb(2,1) =  kb(1,2);
            kb(2,2) =  k0Plastic*(1.0-(qYield*qTrial(0)*qTrial(0)/D)) + k2z;				
        }
        iter++;
    } while ((sqrt(pow(qb(1)-qbOld(0),2)+pow(qb(2)-qbOld(1),2)) >= tol) && (iter <= maxIter));
		
	

    // issue warning if iteration did not converge
    if (iter >= maxIter)   {
        opserr << "WARNING: FPBearingPTV::update() - element: "
            << this->getTag() << " - did not find the shear force after "
            << iter << " iterations and norm: "
            << sqrt(pow(qb(1)-qbOld(0),2)+pow(qb(2)-qbOld(1),2)) << ".\n";
        return -1;
    }
    
    // 3) get moment and stiffness in basic x-direction
    theMaterials[1]->setTrialStrain(ub(3),ubdot(3));
    qb(3) = theMaterials[1]->getStress();
    kb(3,3) = theMaterials[1]->getTangent();
    
    // 4) get moment and stiffness in basic y-direction
    theMaterials[2]->setTrialStrain(ub(4),ubdot(4));
    qb(4) = theMaterials[2]->getStress();
    kb(4,4) = theMaterials[2]->getTangent();
    
    // 5) get moment and stiffness in basic z-direction	
	if (vel1(0) == 0.0 && vel1(1) == 0.0 && vel1(2) == 0.0 && vel2(0) == 0.0 && vel2(1) == 0.0 && vel2(2) == 0.0) { //Under static loading
		ub(5)=0.0; ubdot(5)=0.0; 
	}	
    theMaterials[3]->setTrialStrain(ub(5),ubdot(5));
    qb(5) = theMaterials[3]->getStress();
    kb(5,5) = theMaterials[3]->getTangent();    	
    return 0;
}


const Matrix& FPBearingPTV::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    double Ls = (1.0 - shearDistI)*L;
    // add P-Delta moment stiffness terms
    kl(5,1)   -= qb(0);
    kl(5,7)   += qb(0);
    kl(5,11)  -= qb(0)*Ls;
    kl(11,11) += qb(0)*Ls;
    kl(4,2)   += qb(0);
    kl(4,8)   -= qb(0);
    kl(4,10)  -= qb(0)*Ls;
    kl(10,10) += qb(0)*Ls;
    // add V-Delta torsion stiffness terms
    kl(3,1)   += qb(2);
    kl(3,2)   -= qb(1);
    kl(3,7)   -= qb(2);
    kl(3,8)   += qb(1);
    kl(3,10)  += qb(1)*Ls;
    kl(3,11)  += qb(2)*Ls;
    kl(9,10)  -= qb(1)*Ls;
    kl(9,11)  -= qb(2)*Ls;
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& FPBearingPTV::getInitialStiff()
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


const Matrix& FPBearingPTV::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    // call base class to setup Rayleigh damping
    double factThis = 0.0;
    if (addRayleigh == 1)  {
        theMatrix = this->Element::getDamp();
        factThis = 1.0;
    }
    
    // now add damping tangent from materials
    static Matrix cb(6,6);
    cb.Zero();
    cb(0,0) = theMaterials[0]->getDampTangent();
    cb(3,3) = theMaterials[1]->getDampTangent();
    cb(4,4) = theMaterials[2]->getDampTangent();
    cb(5,5) = theMaterials[3]->getDampTangent();
    
    // transform from basic to local system
    static Matrix cl(12,12);
    cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);
    
    // transform from local to global system and add to cg
    theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    
    return theMatrix;
}


const Matrix& FPBearingPTV::getMass()
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


void FPBearingPTV::zeroLoad()
{
    theLoad.Zero();
}


int FPBearingPTV::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"FPBearingPTV::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << ".\n";
    
    return -1;
}


int FPBearingPTV::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
        opserr << "FPBearingPTV::addInertiaLoadToUnbalance() - "
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


const Vector& FPBearingPTV::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(12);
	ql.Zero();	
    ql = Tlb^qb;
    
    // add P-Delta moments to local forces	
    double MpDelta1 = qb(0)*(ul(7)-ul(1));
    ql(5)  += MpDelta1;	
    double MpDelta2 = qb(0)*(1.0 - shearDistI)*L*ul(11);
    ql(5)  -= MpDelta2;
    ql(11) += MpDelta2;	
    double MpDelta3 = qb(0)*(ul(8)-ul(2));
    ql(4)  -= MpDelta3;
    double MpDelta4 = qb(0)*(1.0 - shearDistI)*L*ul(10);
    ql(4)  -= MpDelta4;
    ql(10) += MpDelta4;
    
    // add V-Delta torsion to local forces
    double Vdelta1 = qb(1)*(ul(8)-ul(2)) - qb(2)*(ul(7)-ul(1));
    ql(3) += Vdelta1;
    double Vdelta2 = (1.0 - shearDistI)*L*(qb(1)*ul(10) + qb(2)*ul(11));
    ql(3) += Vdelta2;
    ql(9) -= Vdelta2;
    
    // determine resisting forces in global system
    theVector = Tgl^ql;
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    return theVector;
}


const Vector& FPBearingPTV::getResistingForceIncInertia()
{
    // this already includes damping forces from materials
    theVector = this->getResistingForce();
    
    // add the damping forces from rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector += this->getRayleighDampingForces();
    }
    
    // add inertia forces from element mass
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


int FPBearingPTV::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(20);
    data(0) = this->getTag();
	data(1)=muRef;
	data(2)=kpFactor;
	data(3)=refPressure;
	data(4)=kTFactor;
	data(5)=diffuse;
	data(6)=conduct;
	data(7)=kvFactor;
	data(8)=rateParam;
	data(9)=Reffective;
	data(10)=rContact;
	data(11)=k0;
	data(12)=x.Size();
	data(13)=y.Size();
	data(14)=shearDistI;
	data(15)=addRayleigh;
	data(16)=mass;
	data(17)=maxIter;
	data(18)=tol;
	data(19)=unit;    
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);    
    
    // send the material class tags
    ID matClassTags(4);    
	matClassTags(0) = theMaterials[0]->getClassTag();
	matClassTags(1) = theMaterials[1]->getClassTag();
	matClassTags(2) = theMaterials[2]->getClassTag();
	matClassTags(3) = theMaterials[3]->getClassTag();
    sChannel.sendID(0, commitTag, matClassTags);
    
    // send the material models   
	theMaterials[0]->sendSelf(commitTag, sChannel);
	theMaterials[1]->sendSelf(commitTag, sChannel);
	theMaterials[2]->sendSelf(commitTag, sChannel);
	theMaterials[3]->sendSelf(commitTag, sChannel);	
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    
    return 0;
}


int FPBearingPTV::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    // delete material memory    
	if (theMaterials[0] != 0)
		delete theMaterials[0];
	if (theMaterials[1] != 0)
		delete theMaterials[1];
	if (theMaterials[2] != 0)
		delete theMaterials[2];
	if (theMaterials[3] != 0)
		delete theMaterials[3];	
    
    // receive element parameters
    static Vector data(20);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    muRef=data(1);
	kpFactor=data(2);
	refPressure=data(3);
	kTFactor=data(4);
	diffuse=data(5);
	conduct=data(6);
	kvFactor=data(7);
	rateParam=data(8);
	Reffective=data(9);
	rContact=data(10);
	k0=data(11);	
	shearDistI=data(14);
	addRayleigh=data(15);
	mass=data(16);
	maxIter=data(17);
	tol=data(18);
	unit=data(19);  
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);    
    
    // receive the material class tags
    ID matClassTags(4);
    rChannel.recvID(0, commitTag, matClassTags);
    
    // receive the material models
	theMaterials[0] = theBroker.getNewUniaxialMaterial(matClassTags(0));
        if (theMaterials[0] == 0) {
            opserr << "FPBearingPTV::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
    theMaterials[0]->recvSelf(commitTag, rChannel, theBroker);

	theMaterials[1] = theBroker.getNewUniaxialMaterial(matClassTags(1));
        if (theMaterials[1] == 0) {
            opserr << "FPBearingPTV::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
    theMaterials[1]->recvSelf(commitTag, rChannel, theBroker);

	theMaterials[2] = theBroker.getNewUniaxialMaterial(matClassTags(2));
        if (theMaterials[2] == 0) {
            opserr << "FPBearingPTV::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
    theMaterials[2]->recvSelf(commitTag, rChannel, theBroker);

	theMaterials[3] = theBroker.getNewUniaxialMaterial(matClassTags(3));
        if (theMaterials[3] == 0) {
            opserr << "FPBearingPTV::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
    theMaterials[3]->recvSelf(commitTag, rChannel, theBroker);	
    
    // receive remaining data
    if ((int)data(12) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(13) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = theMaterials[0]->getInitialTangent();
    kbInit(1,1) = k0;
    kbInit(2,2) = k0;
    kbInit(3,3) = theMaterials[1]->getInitialTangent();
    kbInit(4,4) = theMaterials[2]->getInitialTangent();
    kbInit(5,5) = theMaterials[3]->getInitialTangent();
    
    // initialize variables
    this->revertToStart();
    
    return 0;
}


int FPBearingPTV::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
    int errCode = 0;
    
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();     
    
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    
    if (displayMode >= 0)  {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i=0; i<3; i++)  {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v2(i) = end1Crd(i) + (end1Disp(i) + end2Disp(i))*fact;
            v3(i) = end2Crd(i) + end2Disp(i)*fact;    
        }
    } else  {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        
        if (eigen1.noCols() >= mode)  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v2(i) = end1Crd(i) + (eigen1(i,mode-1) + eigen2(i,mode-1))*fact;
                v3(i) = end2Crd(i) + eigen2(i,mode-1)*fact;
            }
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


void FPBearingPTV::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag(); 
        s << "  type: FPBearingPTV  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        //s << "  FrictionModel: " << theFrnMdl->getTag() << endln;
        s << "  Reff: " << Reffective << "  kInit: " << k0 << endln;
        s << "  Material ux: " << theMaterials[0]->getTag() << endln;
        s << "  Material rx: " << theMaterials[1]->getTag() << endln;
        s << "  Material ry: " << theMaterials[2]->getTag() << endln;
        s << "  Material rz: " << theMaterials[3]->getTag() << endln;		
        s << "  shearDistI: " << shearDistI << "  addRayleigh: "
            << addRayleigh << "  mass: " << mass << endln;
        s << "  maxIter: " << maxIter << "  tol: " << tol << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FPBearingPTV\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        //s << "\"frictionModel\": \"" << theFrnMdl->getTag() << "\", ";
        s << "\"Reff\": " << Reffective << ", ";
        s << "\"kInit\": " << k0 << ", ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\", \"";
        s << theMaterials[2]->getTag() << "\", \"";
        s << theMaterials[3]->getTag() << "\"], ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"addRayleigh\": " << addRayleigh << ", ";
        s << "\"mass\": " << mass << ", ";
        s << "\"maxIter\": " << maxIter << ", ";
        s << "\"tol\": " << tol << "}";
    }
}


Response* FPBearingPTV::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;    
    output.tag("ElementOutput");
    output.attr("eleType","FPBearingPTV");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
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
        
        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 ||
        strcmp(argv[0],"localForces") == 0)
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
        
        theResponse = new ElementResponse(this, 4, theVector);
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
	// Temperature at the center
    else if (strcmp(argv[0],"temperature") == 0 ||
        strcmp(argv[0],"Temperature") == 0 || 
        strcmp(argv[0],"temp") == 0 ||
        strcmp(argv[0],"Temp") == 0 )
    {
        output.tag("ResponseType","TemperatureCenter");        
        
        theResponse = new ElementResponse(this, 6, Vector(1));
    }
	// Pressure, temperature, velocity factors
    else if (strcmp(argv[0],"MuFactors") == 0 ||
        strcmp(argv[0],"mufactors") == 0 || 
        strcmp(argv[0],"mufactor") == 0 ||
        strcmp(argv[0],"FrictionFactors") == 0 )
    {
        output.tag("ResponseType","kPressure");  
		output.tag("ResponseType","kTemperature");
		output.tag("ResponseType","kVelocity");
        
        theResponse = new ElementResponse(this, 7, Vector(3));
    }
	// Adjusted coefficient of friction
    else if (strcmp(argv[0],"MuAdj") == 0 ||
        strcmp(argv[0],"muadj") == 0 || 
        strcmp(argv[0],"MuAdjusted") == 0 ||
        strcmp(argv[0],"muadjusted") == 0 )
    {
        output.tag("ResponseType","MuAdjusted");		
        
        theResponse = new ElementResponse(this, 8, Vector(1));
    }
	// Heat flux at center of sliding surface
    else if (strcmp(argv[0],"HeatFlux") == 0 ||
        strcmp(argv[0],"heatflux") == 0 || 
        strcmp(argv[0],"heatFlux") == 0 ||
        strcmp(argv[0],"Heatflux") == 0 )
    {
        output.tag("ResponseType","HeatFluxCenter");		
        
        theResponse = new ElementResponse(this, 9, Vector(1));
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


int FPBearingPTV::getResponse(int responseID, Information &eleInfo)
{
    double MpDelta1, MpDelta2, MpDelta3, MpDelta4, Vdelta1, Vdelta2;
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector = Tlb^qb;
        // add P-Delta moments
        MpDelta1 = qb(0)*(ul(7)-ul(1));
        theVector(5)  += MpDelta1;
        MpDelta2 = qb(0)*(1.0 - shearDistI)*L*ul(11);
        theVector(5)  -= MpDelta2;
        theVector(11) += MpDelta2;
        MpDelta3 = qb(0)*(ul(8)-ul(2));
        theVector(4)  -= MpDelta3;
        MpDelta4 = qb(0)*(1.0 - shearDistI)*L*ul(10);
        theVector(4)  -= MpDelta4;
        theVector(10) += MpDelta4;
        // add V-Delta torsion
        Vdelta1 = qb(1)*(ul(8)-ul(2)) - qb(2)*(ul(7)-ul(1));
        theVector(3)  += Vdelta1;
        Vdelta2 = (1.0 - shearDistI)*L*(qb(1)*ul(10) + qb(2)*ul(11));
        theVector(3)  += Vdelta2;
        theVector(9)  -= Vdelta2;
        return eleInfo.setVector(theVector);
        
    case 3:  // basic forces
        return eleInfo.setVector(qb);
        
    case 4:  // local displacements
        return eleInfo.setVector(ul);
        
    case 5:  // basic displacements
        return eleInfo.setVector(ub);

	case 6:  // Temperature at the center of sliding surface
		return eleInfo.setVector(TemperatureCenter);

	case 7:  // Pressure, temperature and velocity factors
		return eleInfo.setVector(MuFactors);

	case 8:  // Adjusted coefficient of friction
		return eleInfo.setVector(MuAdjusted);

	case 9:  // Heat flux at the center of the sliding surface
		return eleInfo.setVector(HeatFluxCenter);
        
    default:
        return -1;
    }
}

// Code for using updateParameter command in the tcl file
int
FPBearingPTV::setParameter(const char **argv, int argc, Parameter &param)
{
  int result = -1;  
  int numMaterials1d = 4;
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "material") == 0) {
      if (argc > 2) {
        int matNum = atoi(argv[1]);
        if (matNum >= 1 && matNum <= numMaterials1d)    
          return theMaterials[matNum-1]->setParameter(&argv[2], argc-2, param);
      } else {
        return -1;
      }
  }

  for (int i=0; i<numMaterials1d; i++) {
    int res = theMaterials[i]->setParameter(argv, argc, param);
    if (res != -1) {
      result = res;
    }
  }  
  return result;
}

// Establish the external nodes and set up the transformation matrix for orientation
void FPBearingPTV::setUp()
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
            opserr << "WARNING FPBearingPTV::setUp() - " 
                << "element: " << this->getTag()
                << " - ignoring nodes and using specified "
                << "local x vector to determine orientation.\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "FPBearingPTV::setUp() - "
            << "element: " << this->getTag()
            << " - incorrect dimension of orientation vectors.\n";
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
        opserr << "FPBearingPTV::setUp() - "
            << "element: " << this->getTag()
            << " - invalid orientation vectors.\n";
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


double FPBearingPTV::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}

