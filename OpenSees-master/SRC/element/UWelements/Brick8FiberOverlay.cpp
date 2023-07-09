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
                                                                       
// Created: M. Chiaramonte,  P. Arduino,  P.Mackenzie-Helnwein, UW, 03.29.2011
// Modified: Alborz Ghofrani, UW, 2011
//
// Description: This file contains the implementation of the Brick8FiberOverlay class

#include <Brick8FiberOverlay.h>

#include <elementAPI.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Channel.h>
#include <G3Globals.h>
#include <ErrorHandler.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

double Brick8FiberOverlay::pts[2][3];
double Brick8FiberOverlay::wts[2];
Matrix Brick8FiberOverlay::FiberK(24,24);
Vector Brick8FiberOverlay::P(24);

static int num_Brick8FiberOverlay = 0;
	 

void *
OPS_Brick8FiberOverlay(void)  
{
if (num_Brick8FiberOverlay == 0) {
	num_Brick8FiberOverlay++;
	opserr<<"Brick8FiberOverlay element - Written: M.Chiaramonte, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 15) {
    opserr << "Want: Brick8FiberOverlay tag? nd1? nd2? nd3? nd4? nd5? nd6? nd7? nd8? matTag? AreaFiber? B1? B2? B3? B4?\n";
    return 0;
  }
    
  int		iData[9];
  double	dData[5];
  int		matTag = 0;
  int		eleTag = 0;
  int		numData = 9;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Brick8FiberOverlay" << endln;
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Brick8FiberOverlay: invalid matTag for element: " << eleTag << "\n";
    return 0;
  }
  numData = 5;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Brick8FiberOverlay " << eleTag << endln;
    return 0;	
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matTag);
  
  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element " << eleTag << endln;
    return 0;

  }

  // Parsing was successful, allocate the material
  theElement = new Brick8FiberOverlay(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8], *theMaterial, dData[0], dData[1], dData[2],  dData[3], dData[4]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Brick8FiberOverlay\n";
    return 0;
  }

  return theElement;
}



Brick8FiberOverlay::Brick8FiberOverlay(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8,
										UniaxialMaterial &m, double AreaFiber,double B1, double B2, double B3, double B4)
					:Element(tag, ELE_TAG_Brick8FiberOverlay), 
					externalNodes(SL_NUM_NODE), 
					g1(SL_NUM_NDF), dualg1(SL_NUM_NDF), 
					g2(SL_NUM_NDF), dualg2(SL_NUM_NDF), 
					g3(SL_NUM_NDF), dualg3(SL_NUM_NDF),
					beta1(B1), beta2(B2), beta3(B3), beta4(B4), 
					dNidxAlphai(SL_NUM_NODE,SL_NUM_NDF), 
					Q1(SL_NUM_NDF), Q2(SL_NUM_NDF), Q3(SL_NUM_NDF), Q4(SL_NUM_NDF),
					Q5(SL_NUM_NDF), Q6(SL_NUM_NDF), Q7(SL_NUM_NDF), Q8(SL_NUM_NDF),
					Qfi(SL_NUM_NDF), Qfj(SL_NUM_NDF), nFi(SL_NUM_NDF), nFj(SL_NUM_NDF),
					Vf(SL_NUM_NDF),
					A(SL_NUM_NDF), AA(SL_NUM_NDF*2), Bb(SL_NUM_DOF),
					Af(AreaFiber), u(SL_NUM_DOF)
{	
	//determine int pts location based on fiber orientation

	nFi(0) = (beta1 - 0.5) * 2.0;
	nFi(1) = -1.0;
	nFi(2) = (beta2 - 0.5) * 2.0;

	nFj(0) = (beta3 - 0.5) * 2.0;
	nFj(1) = 1.0;
	nFj(2) = (beta4 - 0.5) * 2.0;

	Vf = nFj - nFi; 
	Vf.Normalize();

	//set up integration parameters (2 intgr pts)
	pts[0][0] = nFi(0)+Vf(0)*(1-0.5773502691896258);   
	pts[0][1] = nFi(1)+Vf(1)*(1-0.5773502691896258);	
	pts[0][2] = nFi(2)+Vf(2)*(1-0.5773502691896258);		
	pts[1][0] = nFj(0)-Vf(0)*(1-0.5773502691896258);    
	pts[1][1] = nFj(1)-Vf(1)*(1-0.5773502691896258);	
	pts[1][2] = nFj(2)-Vf(2)*(1-0.5773502691896258);		
	wts[0] = 1.0;
	wts[1] = 1.0;

	externalNodes(0) = nd1;
	externalNodes(1) = nd2;
	externalNodes(2) = nd3;
	externalNodes(3) = nd4;
	externalNodes(4) = nd5;
	externalNodes(5) = nd6;
	externalNodes(6) = nd7;
	externalNodes(7) = nd8;

	for(int i = 0; i < 2; i++)
		theMaterial[i] = m.getCopy();
	
	for(int i = 0; i < SL_NUM_NODE; i++)
		theNodes[i] = 0;
}


Brick8FiberOverlay::Brick8FiberOverlay()
				: Element(0, ELE_TAG_Brick8FiberOverlay),
				externalNodes(0),
				g1(0), dualg1(0),
				g2(0), dualg2(0),
				g3(0), dualg3(0),
				beta1(0), beta2(0), beta3(0), beta4(0),
				dNidxAlphai(0,0),
				Q1(0), Q2(0), Q3(0), Q4(0),
				Q5(0), Q6(0), Q7(0), Q8(0),
				Qfi(0), Qfj(0),
				Vf(0), A(0), AA(0), 
				Bb(0), Af(0), u(0)
{	
	for(int i = 0; i < SL_NUM_NODE; i++)
		theNodes[i] = 0;
	
	for(int i = 0; i < 2; i++)
		theMaterial[i] = 0;
}	



Brick8FiberOverlay::~Brick8FiberOverlay()
{
	for(int i = 0; i < 2; i++)
		delete theMaterial[i];
}


void
Brick8FiberOverlay::setDomain(Domain *theDomain)	
{	

	 // Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		theNodes[4] = 0;
		theNodes[5] = 0;
		theNodes[6] = 0;
		theNodes[7] = 0;
		return;
     }

     int Nd1 = externalNodes(0);
     int Nd2 = externalNodes(1);
     int Nd3 = externalNodes(2);
     int Nd4 = externalNodes(3);
	 int Nd5 = externalNodes(4);
     int Nd6 = externalNodes(5);
     int Nd7 = externalNodes(6);
     int Nd8 = externalNodes(7);
	
     theNodes[0] = theDomain->getNode(Nd1);
     theNodes[1] = theDomain->getNode(Nd2);
     theNodes[2] = theDomain->getNode(Nd3);
     theNodes[3] = theDomain->getNode(Nd4);
     theNodes[4] = theDomain->getNode(Nd5);
     theNodes[5] = theDomain->getNode(Nd6);
     theNodes[6] = theDomain->getNode(Nd7);
     theNodes[7] = theDomain->getNode(Nd8);

	 Q1 = theNodes[0]->getCrds();
	 Q2 = theNodes[1]->getCrds();
	 Q3 = theNodes[2]->getCrds();
	 Q4 = theNodes[3]->getCrds();
	 Q5 = theNodes[4]->getCrds();
	 Q6 = theNodes[5]->getCrds();
	 Q7 = theNodes[6]->getCrds();
	 Q8 = theNodes[7]->getCrds();

	 //UpdateBase(nFi(0), nFi(1), nFi(2));
	 //for(int i = 0; i < SL_NUM_NDF; i++) {
	 //	 Qfi(i) += g1(i)*(nFi(0)+1); Qfi(i) += g2(i)*(nFi(1)+1) ; Qfi(i) += g3(i)*(nFi(2)+1) ;
	 //}
	 //Qfi += Q1;
	 //UpdateBase(nFj(0), nFj(1), nFj(2));
	 //for(int i = 0; i < SL_NUM_NDF; i++) {
	 //	 Qfj(i) += g1(i)*(nFj(0)+1); Qfj(i) += g2(i)*(nFj(1)+1) ; Qfj(i) += g3(i)*(nFj(2)+1) ;
	 //}
	 //Qfj += Q1;
	 
	double xe1 = nFi(0), eta1 = nFi(1), zet1 = nFi(2);
	double xe2 = nFj(0), eta2 = nFj(1), zet2 = nFj(2);
	Matrix xAlphai(SL_NUM_NDF,SL_NUM_NODE);

	xAlphai(0,0) = -1; xAlphai(0,1) = +1; xAlphai(0,2) = +1; xAlphai(0,3) = -1; xAlphai(0,4) = -1; xAlphai(0,5) = +1; xAlphai(0,6) = +1; xAlphai(0,7) = -1; 
	xAlphai(1,0) = -1; xAlphai(1,1) = -1; xAlphai(1,2) = +1; xAlphai(1,3) = +1; xAlphai(1,4) = -1; xAlphai(1,5) = -1; xAlphai(1,6) = +1; xAlphai(1,7) = +1; 
	xAlphai(2,0) = -1; xAlphai(2,1) = -1; xAlphai(2,2) = -1; xAlphai(2,3) = -1; xAlphai(2,4) = +1; xAlphai(2,5) = +1; xAlphai(2,6) = +1; xAlphai(2,7) = +1; 

	Qfi = 0.125 * (((1+xAlphai(0,0)*xe1)*(1+xAlphai(1,0)*eta1)*(1+xAlphai(2,0)*zet1)) * Q1	+
				   ((1+xAlphai(0,1)*xe1)*(1+xAlphai(1,1)*eta1)*(1+xAlphai(2,1)*zet1)) * Q2	+
				   ((1+xAlphai(0,2)*xe1)*(1+xAlphai(1,2)*eta1)*(1+xAlphai(2,2)*zet1)) * Q3	+
				   ((1+xAlphai(0,3)*xe1)*(1+xAlphai(1,3)*eta1)*(1+xAlphai(2,3)*zet1)) * Q4	+
				   ((1+xAlphai(0,4)*xe1)*(1+xAlphai(1,4)*eta1)*(1+xAlphai(2,4)*zet1)) * Q5	+
				   ((1+xAlphai(0,5)*xe1)*(1+xAlphai(1,5)*eta1)*(1+xAlphai(2,5)*zet1)) * Q6	+
				   ((1+xAlphai(0,6)*xe1)*(1+xAlphai(1,6)*eta1)*(1+xAlphai(2,6)*zet1)) * Q7	+
				   ((1+xAlphai(0,7)*xe1)*(1+xAlphai(1,7)*eta1)*(1+xAlphai(2,7)*zet1)) * Q8);

	Qfj = 0.125 * (((1+xAlphai(0,0)*xe2)*(1+xAlphai(1,0)*eta2)*(1+xAlphai(2,0)*zet2)) * Q1	+
				   ((1+xAlphai(0,1)*xe2)*(1+xAlphai(1,1)*eta2)*(1+xAlphai(2,1)*zet2)) * Q2	+
				   ((1+xAlphai(0,2)*xe2)*(1+xAlphai(1,2)*eta2)*(1+xAlphai(2,2)*zet2)) * Q3	+
				   ((1+xAlphai(0,3)*xe2)*(1+xAlphai(1,3)*eta2)*(1+xAlphai(2,3)*zet2)) * Q4	+
				   ((1+xAlphai(0,4)*xe2)*(1+xAlphai(1,4)*eta2)*(1+xAlphai(2,4)*zet2)) * Q5	+
				   ((1+xAlphai(0,5)*xe2)*(1+xAlphai(1,5)*eta2)*(1+xAlphai(2,5)*zet2)) * Q6	+
				   ((1+xAlphai(0,6)*xe2)*(1+xAlphai(1,6)*eta2)*(1+xAlphai(2,6)*zet2)) * Q7	+
				   ((1+xAlphai(0,7)*xe2)*(1+xAlphai(1,7)*eta2)*(1+xAlphai(2,7)*zet2)) * Q8);


	 A = Qfj - Qfi; 
	 Lf = A.Norm();
	 A.Normalize();
	 AA(0) = A(0)*A(0); AA(1) = A(1)*A(1); AA(2) = A(2)*A(2);
	 AA(3) = A(0)*A(1); AA(4) = A(2)*A(1); AA(5) = A(0)*A(2);
	 this->DomainComponent::setDomain(theDomain);
 }


 int
Brick8FiberOverlay::getNumExternalNodes(void) const
 {
    return SL_NUM_NODE;
 }


const ID &
Brick8FiberOverlay::getExternalNodes(void)
 {
   return externalNodes;
 }
 

Node **
Brick8FiberOverlay::getNodePtrs(void)
 {
   return theNodes;
 }
 

 int
Brick8FiberOverlay::getNumDOF(void) {
     return SL_NUM_DOF;
 }


int
Brick8FiberOverlay::commitState()
{
	int success = 0 ;

	// call element commitState to do any base class stuff
	if ((success = this->Element::commitState()) != 0) {
		opserr << "Brick8FiberOverlay::commitState () - failed in base class";
	}  

	for(int i = 0; i < 2; i++)
		success += theMaterial[i]->commitState();

	return success;
}
 

int
Brick8FiberOverlay::revertToLastCommit()
{
	int success=0;

	for(int i = 0; i < 2; i++)
		success += theMaterial[i]->revertToLastCommit();

	return success;
}
 

int
Brick8FiberOverlay::revertToStart()
{
	int success=0;

	for(int i = 0; i < 2; i++)
		success += theMaterial[i]->revertToStart();

	return success;
}


int
Brick8FiberOverlay::update()
{  
	int success = 0;

	for (int i = 0; i < 2; i++) {
		// Compute strain at Gauss point
		strain = computeCurrentStrain(pts[i][0], pts[i][1], pts[i][2]);
		// Set strain in material 
		success += theMaterial[i]->setTrialStrain(strain);
	}

	return success;
}

double
Brick8FiberOverlay::computeCurrentStrain(double Xi, double Eta,double Zeta)
{  // determine the current strain given trial displacements at nodes      
	strain = 0;

	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	const Vector &disp5 = theNodes[4]->getTrialDisp();
	const Vector &disp6 = theNodes[5]->getTrialDisp();
	const Vector &disp7 = theNodes[6]->getTrialDisp();
	const Vector &disp8 = theNodes[7]->getTrialDisp();

	 u(0) = disp1(0);  u(1) = disp1(1);  u(2) = disp1(2);
	 u(3) = disp2(0);  u(4) = disp2(1);  u(5) = disp2(2);
	 u(6) = disp3(0);  u(7) = disp3(1);  u(8) = disp3(2);
	 u(9) = disp4(0); u(10) = disp4(1); u(11) = disp4(2);
	u(12) = disp5(0); u(13) = disp5(1); u(14) = disp5(2);
	u(15) = disp6(0); u(16) = disp6(1); u(17) = disp6(2);
	u(18) = disp7(0); u(19) = disp7(1); u(20) = disp7(2);
	u(21) = disp8(0); u(22) = disp8(1); u(23) = disp8(2);

	this->getEltBb(Xi, Eta, Zeta);
	for (int i = 0; i < SL_NUM_DDOF; i++) {
		strain += Bb(i)*u(i);
	}

	return strain;
}


Vector
Brick8FiberOverlay::cross(Vector v1, Vector v2){
	Vector v3(SL_NUM_NDF);
	v3(0) = v1(1)*v2(2)-v1(2)*v2(1);
	v3(1) = v2(0)*v1(2)-v1(0)*v2(2);
	v3(2) = v1(0)*v2(1)-v1(1)*v2(0);
	return v3;
}

int
Brick8FiberOverlay::UpdateBase(double Xi, double Eta,double Zeta)
{	
	Matrix xAlphai(SL_NUM_NDF,SL_NUM_NODE);
	g1.Zero(); g2.Zero(); g3.Zero();

	xAlphai(0,0) = -1; xAlphai(0,1) = +1; xAlphai(0,2) = +1; xAlphai(0,3) = -1; xAlphai(0,4) = -1; xAlphai(0,5) = +1; xAlphai(0,6) = +1; xAlphai(0,7) = -1; 
	xAlphai(1,0) = -1; xAlphai(1,1) = -1; xAlphai(1,2) = +1; xAlphai(1,3) = +1; xAlphai(1,4) = -1; xAlphai(1,5) = -1; xAlphai(1,6) = +1; xAlphai(1,7) = +1; 
	xAlphai(2,0) = -1; xAlphai(2,1) = -1; xAlphai(2,2) = -1; xAlphai(2,3) = -1; xAlphai(2,4) = +1; xAlphai(2,5) = +1; xAlphai(2,6) = +1; xAlphai(2,7) = +1; 
	
	Vector Qi;
	for(int i = 0; i < SL_NUM_NODE; i++) {
		dNidxAlphai(i,0) = 0.125 * (xAlphai(0,i) * (1+Eta*xAlphai(1,i)) * (1+Zeta*xAlphai(2,i)));
		dNidxAlphai(i,1) = 0.125 * (xAlphai(1,i) * (1+Xi*xAlphai(0,i))  * (1+Zeta*xAlphai(2,i)));
		dNidxAlphai(i,2) = 0.125 * (xAlphai(2,i) * (1+Xi*xAlphai(0,i))  * (1+Eta*xAlphai(1,i)) );
		Qi = theNodes[i]->getCrds();
		g1 += dNidxAlphai(i,0)*Qi;
		g2 += dNidxAlphai(i,1)*Qi;
		g3 += dNidxAlphai(i,2)*Qi;
	}
    return 0;
}

int
Brick8FiberOverlay::Dual()
{	
	double vol = g1^(cross(g2,g3));
	dualg1 = cross(g2,g3)/vol;
	dualg2 = cross(g3,g1)/vol;
	dualg3 = cross(g1,g2)/vol;
	return 0;
}

int
Brick8FiberOverlay::getEltBb(double Xi, double Eta, double Zeta)
{		
	Matrix B(SL_NUM_NDF,SL_NUM_NODE);

	this->UpdateBase(Xi,Eta,Zeta);
	this->Dual();

	for(int i = 0; i < SL_NUM_NODE; i++) 
		for(int j = 0; j < SL_NUM_NDF; j++)
			B(j,i) = dNidxAlphai(i,0)*dualg1(j) + dNidxAlphai(i,1)*dualg2(j) + dNidxAlphai(i,2)*dualg3(j); 
	
	for(int i = 0; i < SL_NUM_NODE; i++) {
		Bb[i*SL_NUM_NDF  ] = AA(0)*B(0,i)+ AA(3)*B(1,i)+ AA(5)*B(2,i);
		Bb[i*SL_NUM_NDF+1] = AA(1)*B(1,i)+ AA(3)*B(0,i)+ AA(4)*B(2,i);
		Bb[i*SL_NUM_NDF+2] = AA(2)*B(2,i)+ AA(4)*B(1,i)+ AA(5)*B(0,i);
	}
	return 0;
}

const Matrix&
Brick8FiberOverlay::getTangentStiff()
{		
	double Ef;
	FiberK.Zero();
	for (int ip = 0; ip < 2; ip++) { 
		Ef = theMaterial[ip]->getTangent();
		this->getEltBb(pts[ip][0],pts[ip][1],pts[ip][2]);
		for (int i = 0; i < SL_NUM_DOF; i++)
			for (int j = 0; j < SL_NUM_DOF; j++)
				FiberK(i,j) += 0.5*Lf*Af*Ef*wts[ip]*Bb(i)*Bb(j);
	}
	return FiberK;
}

const Matrix &
Brick8FiberOverlay::getInitialStiff(void)
{
	return getTangentStiff();
}


const Vector&
Brick8FiberOverlay::getResistingForce()
{
	P.Zero();
	for (int ip = 0; ip<2; ip++) { 
		this->getEltBb(pts[ip][0],pts[ip][1],pts[ip][2]);
		
		for (int i = 0; i < SL_NUM_DOF; i++)
			P(i) += 0.5 * Lf*Af*wts[ip]*Bb(i)*theMaterial[ip]->getStress();
	}
	
	return P;
}


void
Brick8FiberOverlay::Print(OPS_Stream &s, int flag)
{
	if (flag == 2) {
	    s << "Brick8FiberOverlay element\n";
		const int numNodes = 4;
		const int nstress = 1 ;
		for (int i = 0; i < numNodes; i++) {
			const Vector &nodeCrd = theNodes[i]->getCrds();
			s << "Node " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
		}
    }
}

Response *
Brick8FiberOverlay::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType","Brick8FiberOverlay");
  output.attr("eleTag",this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];
 
  for (int i=0; i<numNodes; i++) {
    sprintf(nodeData,"node%d",i+1);
    output.attr(nodeData,nodes(i));
  }
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
		strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
		const Vector &force = this->getResistingForce();
		int size = force.Size();
		for (int i=0; i<size; i++) {
			sprintf(nodeData,"P%d",i+1);
			output.tag("ResponseType",nodeData);
		}
		theResponse = new ElementResponse(this, 1, Vector(SL_NUM_DOF));
  } else if (strcmp(argv[0],"axialForce") ==0) {
		theResponse = new ElementResponse(this, 2, 0.0);
  }
   output.endTag();
  return theResponse;
}
 
 
int
Brick8FiberOverlay::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
  case -1:
		return -1;
  case 1: // global forces                                                                                               
		return eleInfo.setVector(this->getResistingForce());
  case 2:
	  {
		// Compute average axial force
		double avAxForce = 0;	
		for(int i = 0; i < 2; i++)
			avAxForce += 0.5 * Af * theMaterial[i]->getStress();
		return eleInfo.setDouble(avAxForce);
	  }
  default:
	  return 0;
  }
}

int
Brick8FiberOverlay::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{	
	int dimension = 3;
    static Vector v1(3);
    static Vector v2(3);

    if (displayMode == 1 || displayMode == 2) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[3]->getDisp();    
	 
	 for (int i=0; i<dimension; i++) {
	    v1(i) = Qfi(i)+end1Disp(i)*fact;
	    v2(i) = Qfj(i)+end2Disp(i)*fact;    
	}
	
	// compute the strain and axial force in the member
	double avStrain = 0, avAxForce = 0;
	if (Lf == 0.0) {
	    avStrain = 0.0;
	    avAxForce = 0.0;
	} else {
		for(int i = 0; i < 2; i++){
			strain = this->computeCurrentStrain(pts[i][0], pts[i][1], pts[i][2]);
			theMaterial[i]->setTrialStrain(strain);
			avStrain += strain;
			avAxForce += Af * theMaterial[i]->getStress();
		}
		avAxForce *= 0.5; 
		avStrain *= 0.5;
	}
    
	if (displayMode == 2) // use the strain as the drawing measure
	  return theViewer.drawLine(v1, v2, (float)avStrain, (float)avStrain);	
	else { // otherwise use the axial force as measure
	  return theViewer.drawLine(v1,v2, (float)avAxForce, (float)avAxForce);
	}
    } else if (displayMode < 0) {
      int mode = displayMode  *  -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = Qfi(i) + eigen1(i,mode-1)*fact;
	  v2(i) = Qfj(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = Qfi(i);
	  v2(i) = Qfj(i);
	}    
      }
      return theViewer.drawLine(v1, v2, 1.0, 1.0);	
    }
    return 0;
    return 0;
}

// // database & parallel ////////////////////////////////////////////////////////////////////////
int
Brick8FiberOverlay::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
return 0;
}
int
Brick8FiberOverlay::sendSelf(int commitTag, Channel &theChannel)
{
return 0;
}


