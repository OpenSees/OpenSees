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
                                                                              
// Written: Alborz Ghofrani, Pedro Arduino, U.Washington
// Created: Oct 2014
// Description: This file contains the class definition for QuadBeamEmbedContact.

#include <QuadBeamEmbedContact.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <cmath>

static int num_QuadBeamEmbedContact = 0;
Matrix QuadBeamEmbedContact::m_ContactStiffness(QBEC_NUM_DOF,QBEC_NUM_DOF);
Vector QuadBeamEmbedContact::m_ContactForces(QBEC_NUM_DOF);
const double QuadBeamEmbedContact::m_Pi = 3.14159265359;

void *
OPS_QuadBeamEmbedContact(void)  
{
	if (num_QuadBeamEmbedContact == 0) {
        num_QuadBeamEmbedContact++;
	opserr<<"QuadBeamEmbedContact element - Written: A.Ghofrani, P.Arduino, U.Washington\n";
	}
	
	Element *theElement = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 10) {
		opserr << "Want: QuadBeamEmbedContact tag? Qnd1? Qnd2? Qnd3? Qnd4? Bnd1? Bnd2? radius? fricCoeff? normalPenalty? <tangentialPenalty?> \n";
		return 0;
	}
    
	int iData[7];
	int eleTag = 0;
	int numData = 7;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid integer data: element QuadBeamEmbedContact" << endln;
		return 0;
	}

	numData = 3;
	double dData[3];
	double oData[1];

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
		return 0;
	}

	oData[0] = dData[2];
	numData = numArgs - 10;
	if (numData != 0)
		if (OPS_GetDouble(&numData, oData) != 0) {
		opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
		return 0;
		}
		
	eleTag = iData[0];
	
	theElement = new QuadBeamEmbedContact(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], dData[0], dData[1], dData[2], oData[0]);

	if (theElement == 0) {
		opserr << "WARNING could not create element of type QuadBeamEmbedContact\n";
		return 0;
	}

	return theElement;
}


QuadBeamEmbedContact::QuadBeamEmbedContact(int tag, int Qnd1, int Qnd2, int Qnd3, int Qnd4,
		 int Bnd1, int Bnd2, double r, double mu, double ep, double et): Element(tag, ELE_TAG_QuadBeamEmbedContact),
		 externalNodes(6),
		 m_Q1(2), m_Q2(2), m_Q3(2), m_Q4(2), m_B1(3), m_B2(3),
		 m_Q1_c(2), m_Q2_c(2), m_Q3_c(2), m_Q4_c(2), m_B1_c(3), m_B2_c(3),
		 m_Q1_disp(2), m_Q2_disp(2), m_Q3_disp(2), m_Q4_disp(2), 
		 m_B1_disp(3), m_B2_disp(3), m_Bn(14), m_Bs(14),
		 m_Ba(2), m_Bb(2), m_Ba_c(2), m_Bb_c(2), 
		 m_bShape(4), m_DbShape(4), m_sShape(4),
		 m_nFi(2), m_nFj(2), m_beamIntersect(2),
		 m_eta_n(0.0), m_y_n(2), m_y_n1(2), m_yc0(2), 
		 m_xi_n(2), m_x_n(2), m_x_n1(2),
		 m_t1(2), m_n(2), m_tau1(2), m_nu(2)
{
	m_friction	= mu;
	m_radius	= r;
	m_ep		= ep;
	m_et		= et;

	inContact	= true;
	normDepOnY	= true;
	isStuck		= true;
	wasStuck	= true;
	createElement = false;

	m_SigmaN = m_SigmaT = m_SigmaT_n = 0.0;
	m_Gap = m_Slip = 0.0;

	externalNodes(0) = Qnd1;
	externalNodes(1) = Qnd2;
	externalNodes(2) = Qnd3;
	externalNodes(3) = Qnd4;
	externalNodes(4) = Bnd1;
	externalNodes(5) = Bnd2;

	for(int i = 0; i < 6; i++) 
		theNodes[i] = 0;
}

QuadBeamEmbedContact::QuadBeamEmbedContact()
		: Element(0, ELE_TAG_QuadBeamEmbedContact),
		 externalNodes(6),
		 m_Q1(2), m_Q2(2), m_Q3(2), m_Q4(2), m_B1(3), m_B2(3),
		 m_Q1_c(2), m_Q2_c(2), m_Q3_c(2), m_Q4_c(2), m_B1_c(3), m_B2_c(3),
		 m_Q1_disp(2), m_Q2_disp(2), m_Q3_disp(2), m_Q4_disp(2), 
		 m_B1_disp(3), m_B2_disp(3), m_Bn(14), m_Bs(14),
		 m_Ba(2), m_Bb(2), m_Ba_c(2), m_Bb_c(2), 
		 m_bShape(4), m_DbShape(4), m_sShape(4),
		 m_nFi(2), m_nFj(2), m_beamIntersect(2), 
		 m_eta_n(0.0), m_y_n(2), m_y_n1(2), m_yc0(2), 
		 m_xi_n(2), m_x_n(2), m_x_n1(2),
		 m_t1(2), m_n(2), m_tau1(2), m_nu(2)
{
	for(int i = 0; i < 6; i++) 
		theNodes[i] = 0;
}

QuadBeamEmbedContact::~QuadBeamEmbedContact()
{

}

int 
QuadBeamEmbedContact::getNumExternalNodes(void) const
{
	return QBEC_NUM_NODE;
}

const ID&
QuadBeamEmbedContact::getExternalNodes(void)
{
	return externalNodes;
}
    
Node **
QuadBeamEmbedContact::getNodePtrs(void)
{
	return theNodes;
}

int 
QuadBeamEmbedContact::getNumDOF(void)
{
	return QBEC_NUM_DOF;
}
    
void 
QuadBeamEmbedContact::setDomain(Domain *theDomain)
{
	// assign nodes
 	theNodes[0] = theDomain->getNode(externalNodes(0));
	theNodes[1] = theDomain->getNode(externalNodes(1));
	theNodes[2] = theDomain->getNode(externalNodes(2));
	theNodes[3] = theDomain->getNode(externalNodes(3));
	theNodes[4] = theDomain->getNode(externalNodes(4));
	theNodes[5] = theDomain->getNode(externalNodes(5));

	// check if nodes exist in domain
	if ((theNodes[0] == 0) || (theNodes[1] == 0) || (theNodes[2] == 0) || 
		(theNodes[3] == 0) || (theNodes[4] == 0) || (theNodes[5] == 0)) 
	{
		opserr << "FATAL ERROR QuadBeamEmbedContact (tag: " << this->getTag() << ") : " 
			   << "Node not found in the domain." << endln;
		return;
	}

	// check number of dof for each node
	if ((theNodes[0]->getNumberDOF() != 2) || (theNodes[1]->getNumberDOF() != 2) || 
		(theNodes[2]->getNumberDOF() != 2) || (theNodes[3]->getNumberDOF() != 2) || 
		(theNodes[4]->getNumberDOF() != 3) || (theNodes[5]->getNumberDOF() != 3)) 
	{			
		opserr << "FATAL ERROR QuadBeamEmbedContact (tag: " << this->getTag() << ") : " 
			   << "Node DOF not consistent." << endln;
		return;
	}
		
	// initialize node coordinates
	m_Q1_c = m_Q1 = theNodes[0]->getCrds();
	m_Q2_c = m_Q2 = theNodes[1]->getCrds();
	m_Q3_c = m_Q3 = theNodes[2]->getCrds();
	m_Q4_c = m_Q4 = theNodes[3]->getCrds();
	
	m_B1_c = m_B1 = theNodes[4]->getCrds();
	m_B2_c = m_B2 = theNodes[5]->getCrds();

	// initialize nodal displacements
	m_Q1_disp.Zero();
	m_Q2_disp.Zero();
	m_Q3_disp.Zero();
	m_Q4_disp.Zero();

	m_B1_disp.Zero();
	m_B2_disp.Zero();
	
	// initialize beam's start and end tangent vectors and beam length
	m_Ba = m_B2 - m_B1;
	m_Length = m_Ba.Norm();
	if (m_Length == 0) {
		opserr << "FATAL ERROR QuadBeamEmbedContact (tag: " << this->getTag() << ") : " 
			   << "Beam element has zero length." << endln;
		return;
	} else {
		m_Ba_c = m_Ba /= m_Length;
	}
	m_Bb_c = m_Bb = m_Ba;

	// find the contact point
	int res = getContactPt(m_nFi, m_nFj, m_xi_n, m_eta_n);
	// update shape functions
	updateShapeFuncs(m_xi_n, m_eta_n);
	// update the coordinate system at contact point
	updateBase(m_eta_n);
	
	// update contact point coordinate
	m_x_n1 = m_x_n = m_sShape(0) * m_Q1 + m_sShape(1) * m_Q2 + m_sShape(2) * m_Q3 + m_sShape(3) * m_Q4;
	m_y_n1 = m_y_n = m_bShape(0) * m_B1 + m_bShape(1) * m_Length * m_Ba + m_bShape(2) * m_B2 + m_bShape(3) * m_Length * m_Bb;
	m_yc0  = m_y_n;

	// update Bn and Bs
	computeB();

	//opserr << "xi = " << m_xi_n << ", eta = " << m_eta_n << endln;
	//opserr << "x= " << m_x_n(0) << ", " << m_x_n(1) << endln;
	//opserr << "y= " << m_y_n(0) << ", " << m_y_n(1) << endln;

	this->DomainComponent::setDomain(theDomain);
}

int 
QuadBeamEmbedContact::commitState(void)
{
	int retVal = 0;

	// update all step n variables
	m_y_n = m_y_n1;
	m_x_n = m_x_n1;
	m_Ba = m_Ba_c;
	m_Bb = m_Bb_c;
	m_SigmaT_n = m_SigmaT;
	wasStuck = isStuck;
	
	// update the coordinate system at contact point
	updateBase(m_eta_n);

	//opserr << "old xi = " << m_xi_n << endln;
	//project(m_xi_n, m_x_n);
	//opserr << "new xi = " << m_xi_n << endln;
	//m_x_n1 = m_x_n;

	//opserr << "new y = " << m_y_n1 << endln;
	//opserr << "new x = " << m_x_n1 << endln;
	
	// update nodal displacements
	m_Q1_disp = theNodes[0]->getTrialDisp();
	m_Q2_disp = theNodes[1]->getTrialDisp();
	m_Q3_disp = theNodes[2]->getTrialDisp();
	m_Q4_disp = theNodes[3]->getTrialDisp();
	m_B1_disp = theNodes[4]->getTrialDisp();
	m_B2_disp = theNodes[5]->getTrialDisp();
	
	// update Bn and Bs
	computeB();

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
	    opserr << "QuadBeamEmbedContact::commitState() - failed in base class";
	}
	
	return retVal;
}

int 
QuadBeamEmbedContact::revertToLastCommit(void)
{
	return 0;
}

int 
QuadBeamEmbedContact::revertToStart(void)
{
	return 0;
}

int 
QuadBeamEmbedContact::update(void)
{
	// check if the element was created
	if (createElement) {

		double rot_a;
		double rot_b;
		double sigmaTrial;
		Vector x_c(QBEC_NUM_DIM);
		Matrix mEyeS(2,2);

		mEyeS(0,1) = -1.0; mEyeS(1,0) = 1.0;
		
		// update trial displacements
		m_Q1_c = m_Q1 + theNodes[0]->getTrialDisp();
		m_Q2_c = m_Q2 + theNodes[1]->getTrialDisp();
		m_Q3_c = m_Q3 + theNodes[2]->getTrialDisp();
		m_Q4_c = m_Q4 + theNodes[3]->getTrialDisp();

		m_B1_c = m_B1 + theNodes[4]->getTrialDisp();
		m_B2_c = m_B2 + theNodes[5]->getTrialDisp();

		// get the incremental rotation of the beam nodes
		rot_a = theNodes[4]->getTrialDisp()(2) - m_B1_disp(2);
		rot_b = theNodes[5]->getTrialDisp()(2) - m_B2_disp(2);
		
		// update tangent vectors (linear update)
		m_Ba_c = m_Ba + rot_a*mEyeS*m_Ba;
		m_Bb_c = m_Bb + rot_b*mEyeS*m_Bb;

		// update contact point coordinates
		m_y_n1 = m_bShape(0)*m_B1_c + m_bShape(1)*m_Length*m_Ba_c + m_bShape(2)*m_B2_c + m_bShape(3)*m_Length*m_Bb_c;
		m_x_n1 = m_sShape(0)*m_Q1_c + m_sShape(1)*m_Q2_c + m_sShape(2)*m_Q3_c + m_sShape(3)*m_Q4_c;
		
		// update gap
		m_Gap = m_nu^(m_x_n1 - m_y_n1);

		// check if should be in contact
		if (m_Gap <= 0.0) { // in contact

			m_SigmaN = m_Gap * m_ep;
			inContact = true;

		} //%%<< if (m_Gap <= 0.0) >>%%
		else { // not in contact

			m_SigmaN = 0.0;
			m_SigmaT = 0.0;
			m_SignSigmaT = 0.0;
			inContact = false;

		} //%%<< if (m_Gap <= 0.0) >>%%

		if (inContact) { // in contact

			// update incremental slip and friction stress
			m_Slip	= m_tau1 ^ ((m_x_n1 - m_x_n) - (m_y_n1  - m_y_n));
			sigmaTrial = m_SigmaT_n + m_et * m_Slip;
			m_SignSigmaT = (m_SigmaT_n > 0.0) - (m_SigmaT_n < 0.0);
			m_Phi	= fabs(sigmaTrial) - fabs(m_friction * m_SigmaN);
			
			if (m_Phi <= 0.0) { // stick
				m_SigmaT = sigmaTrial;
				isStuck  = true;

			}  //%%<< if (m_Phi <= 0.0) >>%%
			else { // slip

				//m_SigmaT = sigmaTrial - m_Phi * m_SignSigmaT;
				//if (wasStuck)
					//m_SignSigmaT = (m_SigmaT_n > 0.0) - (m_SigmaT_n < 0.0);

				m_SigmaT = - m_friction * m_SigmaN * m_SignSigmaT;
				isStuck = false;
			}

		} //%%<< if (inContact) >>%%
		else { // not in contact
			

		} //%%<< if (inContact) >>%%


		return 0;

	}  //%%<< if (createElement) >>%%
	else { // element was not created

		inContact = false;
		return 0;

	} //%%<< if (createElement) >>%%
}

const Matrix&
QuadBeamEmbedContact::getTangentStiff(void)
{
	m_ContactStiffness.Zero();
	if (inContact) {
		double jacobian = getIntJacobian();

		if (isStuck)
			for (int i = 0; i < QBEC_NUM_DOF; i++) 
				for (int j = 0; j < QBEC_NUM_DOF; j++)
					m_ContactStiffness(i,j) = (m_ep * m_Bn(i) * m_Bn(j) -
							m_et * m_Bs(i) * m_Bs(j)) * 
							jacobian * m_Pi * m_radius * 2.0;
		else
			for (int i = 0; i < QBEC_NUM_DOF; i++) 
				for (int j = 0; j < QBEC_NUM_DOF; j++)
					m_ContactStiffness(i,j) = (m_ep * m_Bn(i) * m_Bn(j) +
							m_ep * m_SignSigmaT * m_friction * m_Bs(i) * m_Bn(j)) * 
							jacobian * m_Pi * m_radius * 2.0;
	}
	return m_ContactStiffness;
}

const Matrix&
QuadBeamEmbedContact::getInitialStiff(void)
{
	return this->getTangentStiff();
}

const Vector&
QuadBeamEmbedContact::getResistingForce(void)
{
	m_ContactForces.Zero();
	
	if(inContact) {
		double jacobian = getIntJacobian();
		
		for (int i = 0; i< QBEC_NUM_DOF; i++) 
			m_ContactForces(i) = (m_SigmaN * m_Bn(i) - m_SigmaT * m_Bs(i)) 
								* jacobian * m_Pi * m_radius* 2.0;
	}
	return m_ContactForces;
}

int 
QuadBeamEmbedContact::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int 
QuadBeamEmbedContact::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker)
{
	return 0;
}

int 
QuadBeamEmbedContact::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}
 
void 
QuadBeamEmbedContact::Print(OPS_Stream &s, int flag)
{
	return;
}

Response*
QuadBeamEmbedContact::setResponse(const char **argv, int argc, 
			  OPS_Stream &s)
{
	if (strcmp(argv[0],"beamForce") == 0 || strcmp(argv[0],"beamForces") == 0) {
		return new ElementResponse(this, 1, Vector(6));
		
	} else if (strcmp(argv[0],"solidForce") == 0 || strcmp(argv[0],"solidForces") == 0) {
		return new ElementResponse(this, 2, Vector(8));
		
	} else if (strcmp(argv[0],"Force") == 0 || strcmp(argv[0],"Forces") == 0) {
		return new ElementResponse(this, 3, Vector(2));

	} else if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
		return new ElementResponse(this, 4, Vector(3));

	} else {
		opserr << "QuadBeamEmbedContact Recorder, " << argv[0] << "is an unknown recorder request" 
			   << "  Element tag : " << this->getTag() << endln;
	return 0;
	}
}

int 
QuadBeamEmbedContact::getResponse(int responseID, Information &eleInformation)
{	
	if (responseID == 1) { // beam forces
		
		Vector tempForce(14);
		Vector bForce(6);
		tempForce = getResistingForce();
		for (int i = 0; i < 6; i++) 
			bForce(i) = tempForce(i+8);

		return eleInformation.setVector(bForce);
	
	} else if (responseID == 2) { // solid forces

		Vector tempForce(14);
		Vector sForce(8);
		tempForce = getResistingForce();
		for (int i = 0; i < 8; i++) 
			sForce(i) = tempForce(i);

		return eleInformation.setVector(sForce);
	
	} else if (responseID == 3) { // scalar contact forces on the beam

		double jacobian;
		Vector force(2);
		jacobian = getIntJacobian();
		
		force(0)	= m_SigmaN * jacobian * 2.0 * m_Pi * m_radius;
		force(1)	= m_SigmaT * jacobian * 2.0 * m_Pi * m_radius;

		return eleInformation.setVector(force);
	
	} else if (responseID == 4) { // scalar contact stresses on the beam 

		Vector stress(3);

		stress(0)	= m_SigmaN;
		stress(1)	= m_SigmaT;
		stress(2)	= isStuck;

		return eleInformation.setVector(stress);
	
	} else {
		opserr << "QuadBeamEmbedContact, tag = " << this->getTag()
		  << " -- unknown request" << endln;
		return -1;
	}
}

int 
QuadBeamEmbedContact::setParameter(const char **argv, int argc, Parameter &param)
{
	return 0;
}
    
int 
QuadBeamEmbedContact::updateParameter(int parameterID, Information &info)
{
	return 0;
}

int
QuadBeamEmbedContact::intersection(Vector& solidIntersect, Vector& beamIntersect)
{
	// output:
	// -2: No Intersection
	// -1: 1 Intersection
	//  0: 2 Intersections

	int		result = 0;
	double	tol = 1.0e-10;
	double	B_x[2], B_y[2];
	double	temp;

	solidIntersect(0) = solidIntersect(1) = -1.0;
	for(int i = 0; i < 2; i++) {
		B_x[i] = theNodes[i+4]->getCrds()[0];
		B_y[i] = theNodes[i+4]->getCrds()[1];
	}

	for(int i = 0; i < 4; i++){
		double Qxi = theNodes[i]->getCrds()[0];
		double Qyi = theNodes[i]->getCrds()[1];
		temp = std::fabs((Qyi-B_y[0])*(B_x[1]-B_x[0])-(B_y[1]-B_y[0])*(Qxi-B_x[0]));
		if (temp < tol) {
			solidIntersect(result) = (double)i;
			beamIntersect(result) = sqrt((pow(Qyi-B_y[0],2.0) + pow(Qxi-B_x[0],2.0))/
										 (pow(B_y[1]-B_y[0],2.0) + pow(B_x[1]-B_x[0],2.0)));
			result++;
			continue;
		} else {
			int j = i + 1;
			if (j == 4) j = 0;
			double Qxj = theNodes[j]->getCrds()[0];
			double Qyj = theNodes[j]->getCrds()[1];
			temp = ((B_x[1]-B_x[0])*(B_y[0]-Qyi)-(B_x[0]-Qxi)*(B_y[1]-B_y[0])) * 
				   ((B_x[1]-B_x[0])*(B_y[0]-Qyj)-(B_x[0]-Qxj)*(B_y[1]-B_y[0]));
			if (temp < -tol) {
				temp = ((B_x[1]-B_x[0])*(Qyi-B_y[0])-(B_y[1]-B_y[0])*(Qxi-B_x[0])) / 
					   ((B_y[1]-B_y[0])*(Qxj- Qxi)-(B_x[1]-B_x[0])*(Qyj- Qyi));
				solidIntersect(result) = (double)i + temp;
				double Qintersect[2];
				Qintersect[0] = Qxi + temp * (Qxj - Qxi);
				Qintersect[1] = Qyi + temp * (Qyj - Qyi);
				beamIntersect(result) = sqrt((pow(Qintersect[1]-B_y[0],2.0) + pow(Qintersect[0]-B_x[0],2.0))/
											 (pow(B_y[1]-B_y[0],2.0) + pow(B_x[1]-B_x[0],2.0)));
				result++;
				continue;
			}
		}
	}
	return (result-2);
}

int
QuadBeamEmbedContact::getContactPt(Vector& nFi, Vector& nFj, Vector& solidContactPt, double& beamContactPt)
{
	int result = 0;
	int iStartNode, jStartNode, iEndNode, jEndNode;
	Vector solidIntersect(2);

	result = intersection(solidIntersect, m_beamIntersect);
	if (result == -2){
		//opserr << "No intersection detected!" << endln;
		createElement = false;
	}
	else if(result == -1) {
		//opserr << "One intersection detected at node " << solidIntersect(0) << endln;
		//opserr << "Corresponding intersection " << m_beamIntersect(0) << endln;
		createElement = false; 
	}
	else {
		//opserr << "Two intersections detected at " << solidIntersect(0) << ", " << solidIntersect(1) << endln;
		//opserr << "Corresponding intersections " << m_beamIntersect(0) << ", " << m_beamIntersect(1) << endln;
		createElement = true;
	}
	
	if (createElement) {
		iStartNode = (int)(floor(solidIntersect(0)));
		iEndNode   = (int)(iStartNode + 1);
		jStartNode = (int)(floor(solidIntersect(1)));
		jEndNode   = (int)(jStartNode + 1);
	
		Matrix isoNodeCrds(2,5);
		
		isoNodeCrds(0,0) = -1.0;	isoNodeCrds(1,0) = -1.0;
		isoNodeCrds(0,1) =  1.0;	isoNodeCrds(1,1) = -1.0;
		isoNodeCrds(0,2) =  1.0;	isoNodeCrds(1,2) =  1.0;
		isoNodeCrds(0,3) = -1.0;	isoNodeCrds(1,3) =  1.0;
		isoNodeCrds(0,4) = -1.0;	isoNodeCrds(1,4) = -1.0;
	
		nFi(0) = isoNodeCrds(0,iStartNode) + (solidIntersect(0) - floor(solidIntersect(0))) * (isoNodeCrds(0,iEndNode) - isoNodeCrds(0,iStartNode));
		nFi(1) = isoNodeCrds(1,iStartNode) + (solidIntersect(0) - floor(solidIntersect(0))) * (isoNodeCrds(1,iEndNode) - isoNodeCrds(1,iStartNode));
		nFj(0) = isoNodeCrds(0,jStartNode) + (solidIntersect(1) - floor(solidIntersect(1))) * (isoNodeCrds(0,jEndNode) - isoNodeCrds(0,jStartNode));
		nFj(1) = isoNodeCrds(1,jStartNode) + (solidIntersect(1) - floor(solidIntersect(1))) * (isoNodeCrds(1,jEndNode) - isoNodeCrds(1,jStartNode));
		
		beamContactPt = 0.5 * (m_beamIntersect(0) + m_beamIntersect(1));
	
		solidContactPt(0) = 0.5*(nFi(0)+nFj(0));
		solidContactPt(1) = 0.5*(nFi(1)+nFj(1));
	
	}

	return result;
}

int
QuadBeamEmbedContact::updateShapeFuncs(Vector xi, double eta)
{
	int result = 0;

	if ((xi(0) > 1.0)||(xi(0) <-1.0)||(xi(1) > 1.0)
	  ||(xi(1) <-1.0)||(eta   > 1.0)||(eta   < 0.0)) 
	{
		opserr << "QuadBeamEmbedContact::Shape Function Parameter not in Range." << endln;
		opserr << "xi : " << xi << endln;
		opserr << "eta: " << eta << endln;
		result = -1;
	}

	double eta2 = eta*eta;
	double eta3 = eta2*eta;

	m_bShape(0) = 1 - 3.0*eta2 + 2.0*eta3;
	m_bShape(1) = eta - 2.0*eta2 + eta3;
	m_bShape(2) = 3.0*eta2 - 2.0*eta3;
	m_bShape(3) = -1.0*eta2 + eta3;

	m_DbShape(0) = -6.0*eta + 6.0*eta2;
	m_DbShape(1) = 1.0 - 4.0*eta + 3.0*eta2;
	m_DbShape(2) = 6.0*eta - 6.0*eta2;
	m_DbShape(3) = -2.0*eta + 3.0*eta2;

	m_sShape(0) = 0.25 * (1-xi(0)) * (1-xi(1));
	m_sShape(1) = 0.25 * (1+xi(0)) * (1-xi(1));
	m_sShape(2) = 0.25 * (1+xi(0)) * (1+xi(1));
	m_sShape(3) = 0.25 * (1-xi(0)) * (1+xi(1));

	return 0;
}

int
QuadBeamEmbedContact::updateBase(double eta)
{
	Matrix Projector(2,2);
	double normN;

	updateShapeFuncs(m_xi_n, eta);

	m_t1 = m_DbShape(0)*m_B1 + m_DbShape(1)*m_Length*m_Ba + m_DbShape(2)*m_B2 + m_DbShape(3)*m_Length*m_Bb;
	m_tau1 = m_t1 / m_t1.Norm();

	Projector(0,0) = 1-m_tau1(0)*m_tau1(0);
	Projector(0,1) = -m_tau1(0)*m_tau1(1);
	Projector(1,0) = -m_tau1(1)*m_tau1(0);
	Projector(1,1) = 1-m_tau1(1)*m_tau1(1);

	//if (normDepOnY) {
	//	m_n = Projector * (m_yc0 - m_y_n1);
	//	normN = m_n.Norm();
	//	if (normN > 1.0e-15) 
	//		m_nu = m_n / normN;
	//	else 
	//		normDepOnY = false;
	//}
	//if (!normDepOnY) {
	//	m_n = Projector * (m_x_n1 - m_yc0);
	//	normN = m_n.Norm();
	//	if (normN > 1.0e-15){
	//		m_nu = m_n / normN;
	//		normDepOnY = true;
	//	}
	//	else {
	//		m_nu.Zero();
	//		normDepOnY = true;
	//	}
	//}

	
	m_nu(0) = -m_tau1(1);
	m_nu(1) = m_tau1(0);

	//opserr << "tau = " << m_tau1;
	//opserr << " nu = " << m_nu;

	return 0;
}

double
QuadBeamEmbedContact::getIntJacobian(void)
{
	Vector iNode(2), jNode(2), intNodeD(2);
	Vector iNodeT(2), jNodeT(2);
	double length;
	double result = 0.0;

	updateShapeFuncs(m_xi_n, m_beamIntersect(0));
	iNode = m_bShape(0)*m_B1_c + m_bShape(1)*m_Length*m_Ba_c + m_bShape(2)*m_B2_c + m_bShape(3)*m_Length*m_Bb_c;
	iNodeT = m_DbShape(0)*m_B1_c + m_DbShape(1)*m_Length*m_Ba_c + m_DbShape(2)*m_B2_c + m_DbShape(3)*m_Length*m_Bb_c;
	
	updateShapeFuncs(m_xi_n, m_beamIntersect(1));
	jNode = m_bShape(0)*m_B1_c + m_bShape(1)*m_Length*m_Ba_c + m_bShape(2)*m_B2_c + m_bShape(3)*m_Length*m_Bb_c;
	jNodeT = m_DbShape(0)*m_B1_c + m_DbShape(1)*m_Length*m_Ba_c + m_DbShape(2)*m_B2_c + m_DbShape(3)*m_Length*m_Bb_c;

	length = (jNode - iNode).Norm();

	updateShapeFuncs(m_xi_n, 0.5);
	intNodeD = m_DbShape(0)*iNode + m_DbShape(1)*length*iNodeT + m_DbShape(2)*jNode + m_DbShape(3)*length*jNodeT;

	updateShapeFuncs(m_xi_n, m_eta_n);

	result = sqrt(intNodeD(0)*intNodeD(0)+intNodeD(1)*intNodeD(1));

	return result;
}

void
QuadBeamEmbedContact::computeB(void)
{
	Matrix mEyeS(2,2);
	mEyeS(0,1) = -1.0; mEyeS(1,0) = 1.0;
	
	m_Bn(0)		= m_sShape(0) * m_nu(0);
	m_Bn(1)		= m_sShape(0) * m_nu(1);
	m_Bn(2)		= m_sShape(1) * m_nu(0);
	m_Bn(3)		= m_sShape(1) * m_nu(1);
	m_Bn(4)		= m_sShape(2) * m_nu(0);
	m_Bn(5)		= m_sShape(2) * m_nu(1);
	m_Bn(6)		= m_sShape(3) * m_nu(0);
	m_Bn(7)		= m_sShape(3) * m_nu(1);
				   
	m_Bn(8)		= -m_bShape(0) * m_nu(0);
	m_Bn(9)		= -m_bShape(0) * m_nu(1);
	m_Bn(10)	= -m_bShape(1) * m_Length * (m_nu^(mEyeS * m_Ba));
	m_Bn(11)	= -m_bShape(2) * m_nu(0);
	m_Bn(12)	= -m_bShape(2) * m_nu(1);
	m_Bn(13)	= -m_bShape(3) * m_Length * (m_nu^(mEyeS * m_Bb));
	
	m_Bs(0)		= m_sShape(0) * m_tau1(0);
	m_Bs(1)		= m_sShape(0) * m_tau1(1);
	m_Bs(2)		= m_sShape(1) * m_tau1(0);
	m_Bs(3)		= m_sShape(1) * m_tau1(1);
	m_Bs(4)		= m_sShape(2) * m_tau1(0);
	m_Bs(5)		= m_sShape(2) * m_tau1(1);
	m_Bs(6)		= m_sShape(3) * m_tau1(0);
	m_Bs(7)		= m_sShape(3) * m_tau1(1);
		   
	m_Bs(8)		= -m_bShape(0) * m_tau1(0);
	m_Bs(9)		= -m_bShape(0) * m_tau1(1);
	m_Bs(10)	= -m_bShape(1) * m_Length * (m_tau1^(mEyeS * m_Ba));
	m_Bs(11)	= -m_bShape(2) * m_tau1(0);
	m_Bs(12)	= -m_bShape(2) * m_tau1(1);
	m_Bs(13)	= -m_bShape(3) * m_Length * (m_tau1^(mEyeS * m_Bb));

}

int
QuadBeamEmbedContact::project(Vector& xi, Vector& x_n1)
{
	double	tol		= 1.0e-10;
	int		maxIter = 50;
	Matrix	grad_inv(2,2);
	Vector	res(2);
	double	dx_dxi, dx_deta, dy_dxi, dy_deta, detG;

	updateShapeFuncs(xi, m_eta_n);

	x_n1 = m_sShape(0)*m_Q1_c + m_sShape(1)*m_Q2_c + m_sShape(2)*m_Q3_c + m_sShape(3)*m_Q4_c;
	res	 = x_n1 - m_y_n1;

	for (int i = 0; i < maxIter; i++) {
		dx_dxi	= 0.25 * ((xi(1)-1)*m_Q1_c(0)+(1-xi(1))*m_Q2_c(0)+(xi(1)+1)*m_Q3_c(0)-(xi(1)+1)*m_Q4_c(0));
		dy_dxi	= 0.25 * ((xi(1)-1)*m_Q1_c(1)+(1-xi(1))*m_Q2_c(1)+(xi(1)+1)*m_Q3_c(1)-(xi(1)+1)*m_Q4_c(1));
		dx_deta = 0.25 * ((xi(0)-1)*m_Q1_c(0)-(xi(0)+1)*m_Q2_c(0)+(xi(0)+1)*m_Q3_c(0)+(1-xi(0))*m_Q4_c(0));
		dy_deta	= 0.25 * ((xi(0)-1)*m_Q1_c(1)-(xi(0)+1)*m_Q2_c(1)+(xi(0)+1)*m_Q3_c(1)+(1-xi(0))*m_Q4_c(1));

		detG	= dx_dxi*dy_deta - dx_deta*dy_dxi;
		if (detG == 0) {
			opserr << "A problem here in Project()" << endln;
			return -1;
		}
		grad_inv(0,0) =  dy_deta;
		grad_inv(1,1) =  dx_dxi;
		grad_inv(0,1) = -dx_deta;
		grad_inv(1,0) = -dy_dxi;

		grad_inv /= detG;
		xi -= grad_inv * res;

		updateShapeFuncs(xi, m_eta_n);
		x_n1 = m_sShape(0)*m_Q1_c + m_sShape(1)*m_Q2_c + m_sShape(2)*m_Q3_c + m_sShape(3)*m_Q4_c;
		res	 = x_n1 - m_y_n1;
		
		if (i == 49) opserr << "maxIter reached!!!" << endln;

		if (res.Norm() < tol)
			break;
	}
	return 0;
}
