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
                                                                        
//
// $Date: 2016/05/29 12:19
// 
// Written: Erkan Bicici 
//
// Description: Axial deflection model with respect ot lateral deflection 
//				02/11/2018 Applicaiton of New features. 
//				Eliminate zero stiffness.

#include <elementAPI.h>
#include "Axial01.h"

#include <stdlib.h>
#include <G3Globals.h>
#include <math.h>
#include <ElementResponse.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <DomainComponent.h>
#include <ArrayOfTaggedObjects.h>
#include <Information.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Parameter.h>


#include <DummyStream.h>
#define PI (3.141592653589793)
#include <fstream>
#include <iomanip>
#include <iostream>


using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;




static int numAxial01 = 0;


void *
OPS_Axial01()
{
  // print out some KUDO's
  if (numAxial01 == 0) {
    opserr << "Axial01 unaxial material\n";
    numAxial01 =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  int    iNodeData[2];
  double dData[3];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, iNodeData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;	
  }
  
  Domain *theDomain = 0;
  theDomain = OPS_GetDomain();
  if (theDomain == 0) {
	  opserr << "WARNING Axial01 -- Pointer to Domain was not returned\n" << endln;
	  return 0;
  }
  Node *theNodeT = 0;
  theNodeT = theDomain->getNode(iNodeData[0]);
  if (theNodeT == 0) {
	  opserr << "WARNING Axial01 -- nodeT with tag \n" << endln;
	  return 0;
  }
  Node *theNodeB = 0;
  theNodeB = theDomain->getNode(iNodeData[1]);
  if (theNodeB == 0) {
	  opserr << "WARNING Axial01 -- nodeB with tag \n" << endln;
	  return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
	  opserr << "WARNING invalid E & ep\n";
	  return 0;
  }

  
  theMaterial = new Axial01(	iData[0], 
							iNodeData[0], iNodeData[1], 
							dData[0], dData[1], dData[2], 
	  theDomain, theNodeT, theNodeB);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Axial01\n";
    return 0;
  }

  // return the material
  return theMaterial;
}




Axial01::Axial01(int tag, int nodeT, int nodeB, double kax, double dsh, double da,  Domain *theDom, Node *theNdT, Node *theNdB)
	:UniaxialMaterial(tag, 0),
	nodeTop(nodeT), nodeBot(nodeB), Kax(kax), Dsh(dsh), Da(da),
	theDomain(theDom), theNodeT(theNdT), theNodeB(theNdB),

trialStrain(0.0), trialStress(0.0), trialTangent(9.9e9), dof(1),
commitStrain(0.0), commitStress(0.0), commitTangent(9.9e9), failure(0.0)
{
	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}


Axial01::~Axial01()
{
  // does nothing
}


int 
Axial01::setTrialStrain(double strain, double strainRate)
{
	Node *theNodeT = theDomain->getNode(nodeTop);
	Node *theNodeB = theDomain->getNode(nodeBot);

	double lt = getLateralDrift();
	
	//Definition of maximum and minimum

	Tlt = lt;
	trialStrain = strain;

	// Same results

	if (trialStrain == commitStrain) { trialTangent = commitTangent; }


	// Degradation on the stiffness
	double DegK = Kax / (Da - Dsh);


	if (lt > C_max_lt) { T_max_lt = lt; }
	if (lt < C_min_lt) { T_min_lt = lt; }


	if ((C_max_lt < Dsh) && (C_min_lt > -Dsh)) {
		trialTangent = 9.9e9;

	}

	else {

		if (lt > C_max_lt) {

						if (commitTangent > Kax){
						trialTangent = Kax - DegK*(lt - Dsh);
						//MaxCheck();
						}
						else {
						trialTangent = commitTangent - (commitTangent / (Da - abs(Clt))) *	(abs(lt) - abs(Clt));
						MaxCheck();
						}
		}
		else if (lt < C_min_lt) {

						if (commitTangent > Kax) {
						trialTangent = Kax - DegK* (abs(lt) - Dsh);
						//MaxCheck();
						}
						else {
						trialTangent = commitTangent - (commitTangent / (Da - abs(Clt))) *	(abs(lt) - abs(Clt));
						MaxCheck();
						}	
		}
		else {
		trialTangent = commitTangent;
		}
	}
	if (trialTangent < CMinTangent) { TMinTangent = trialTangent; }

	trialStress = trialStrain * trialTangent;




	// opserr << "nodeT = " << theNodeT << " " << endln;
	// opserr << "nodeB = " << theNodeB << " " << endln;
	//opserr << "		di = " << dispI(0) << " " << endln;
	// opserr << "		di = " << T_max_lt << " " << endln;
	//opserr << "		di = " << lt << " " << endln;
	//opserr << "		di = " << Tlt << " " << endln;
	//opserr << "		di = " << Clt << " " << endln;

	//opserr << "		tangent = " << trialTangent << " " << endln;
	//opserr << "		tangent = " << commitTangent << " " << endln;
	//opserr << "		failure = " << failure << " " << endln;
	//opserr << "		di = " << lta << " " << endln;


    return 0;
}



double
Axial01::getLateralDrift(void) {

	// get displacements
	const Vector &dispI = theNodeT->getTrialDisp();
	const Vector &dispJ = theNodeB->getTrialDisp();

	double dx = dispI(0) - dispJ(0);

	return dx;  
}

double
Axial01::getVerticalDrift(void) {

	// get displacements
	const Vector &dispI = theNodeT->getTrialDisp();
	const Vector &dispJ = theNodeB->getTrialDisp();

	double dy = dispI(1) - dispJ(1);

	return dy;
}




double 
Axial01::getStrain(void)
{
  return trialStrain;
}

double 
Axial01::getStress(void)
{
  return trialStress;
}


double 
Axial01::getTangent(void)
{
  return trialTangent;
}

int 
Axial01::commitState(void)
{
    commitStrain = trialStrain;
    commitTangent=trialTangent;
    commitStress = trialStress;
	C_max_lt = T_max_lt;
	C_min_lt = T_min_lt;
	CMinTangent = TMinTangent;
	Clt = Tlt;

    return 0;
}	


int 
Axial01::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;
  T_max_lt = C_max_lt;
  T_min_lt = C_min_lt;
  TMinTangent = CMinTangent;
  Tlt = Clt;

  return 0;
}


int 
Axial01::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialTangent = commitTangent = 9.9e9;
  trialStress = commitStress = 0.0;
  C_max_lt = 0.0;
  C_min_lt = 0.0;
  Clt = 0.0;
  CMinTangent = 9.9e9;

  return 0;
}


UniaxialMaterial *
Axial01::getCopy(void)
{
  Axial01 *theCopy =
    new Axial01(this->getTag(), 
					nodeTop, nodeBot,
					Kax, Dsh, Da,
		theDomain, theNodeT, theNodeB);

  theCopy->nodeTop = nodeTop;
  theCopy->nodeBot = nodeBot;
  theCopy->Kax = Kax;
  theCopy->Dsh = Dsh;
  theCopy->Da = Da;
  

  theCopy->commitStrain = commitStrain;
  theCopy->commitStress = commitStress;
  theCopy->commitTangent = commitTangent;


  return theCopy;
}


int 
Axial01::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = nodeTop;
  data(2) = nodeBot;
  data(3) = Kax;
  data(4) = Dsh;
  data(5) = Da;
  data(6) = commitStrain;
  data(7) = commitStress;
  data(8) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Axial01::sendSelf() - failed to send data\n";

  return res;
}

int 
Axial01::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Axial01::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
	nodeTop = data(1);
	nodeBot = data(2);
    Kax       = data(3);
	Dsh		= data(4);
	Da		= data(5);
    commitStrain=data(6);
    commitStress=data(7);
    commitTangent=data(8);
    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
  }

  return res;
}

void 
Axial01::Print(OPS_Stream &s, int flag)
{

}

void 
Axial01::MaxCheck(void) {

	//if (trialTangent > CMinTangent) {  trialTangent = CMinTangent; }
	if (trialTangent < 1) { trialTangent = 1; }
}