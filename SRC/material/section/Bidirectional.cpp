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
                                                                        
// $Revision: 1.6 $
// $Date: 2003-02-14 23:01:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/Bidirectional.cpp,v $

#include <Bidirectional.h>           
#include <Channel.h>

Vector Bidirectional::s(2);
Matrix Bidirectional::ks(2,2);
ID Bidirectional::code(2);

Bidirectional::Bidirectional
(int tag, double e, double sy, double Hi, double Hk) :
 SectionForceDeformation(tag, SEC_TAG_Bidirectional),
	 E(e), sigY(sy), Hiso(Hi), Hkin(Hk)
{
	for (int i = 0; i < 2; i++) {
		eP_n[i]  = 0.0;
		eP_n1[i] = 0.0;
		q_n[i]  = 0.0;
		q_n1[i] = 0.0;
	}

	alpha_n  = 0.0;
	alpha_n1 = 0.0;

	code(0) = SECTION_RESPONSE_VY;
	code(1) = SECTION_RESPONSE_P;
}

Bidirectional::Bidirectional():
 SectionForceDeformation(0, SEC_TAG_Bidirectional),
	 E(0.0), sigY(0.0), Hiso(0.0), Hkin(0.0)
{
	for (int i = 0; i < 2; i++) {
		eP_n[i]  = 0.0;
		eP_n1[i] = 0.0;
		q_n[i]  = 0.0;
		q_n1[i] = 0.0;
	}

	alpha_n  = 0.0;
	alpha_n1 = 0.0;

	code(0) = SECTION_RESPONSE_VY;
	code(1) = SECTION_RESPONSE_P;
}

Bidirectional::~Bidirectional()
{
	// Nothing to do here
}

int
Bidirectional::setTrialSectionDeformation(const Vector &e)
{
	e_n1[0] = e(0);
	e_n1[1] = e(1);

	return 0;
}

const Matrix&
Bidirectional::getSectionTangent(void)
{
	// Compute trial stress using elastic tangent
	s(0) = E*(e_n1[0]-eP_n[0]);
	s(1) = E*(e_n1[1]-eP_n[1]);

	static Vector xsi(2);

	// Predicted stress minus back stress
	xsi(0) = s(0) - q_n[0];
	xsi(1) = s(1) - q_n[1];

	double normxsi = xsi.Norm();

	// Current yield stress
	double sigY_n = sigY + alpha_n*Hiso;

	// Yield function
	double f_n1 = normxsi - sigY_n;

	// Elastic step
	if (f_n1 < 0.0) {
		ks(0,0) = ks(1,1) = E;
		ks(0,1) = ks(1,0) = 0.0;
	}

	// Plastic step
	else {
		// Consistency parameter
		double dlam = f_n1/(E+Hkin+Hiso);

		double n_n1[2];

		// Normal vector
		n_n1[0] = xsi(0)/normxsi;
		n_n1[1] = xsi(1)/normxsi;

		double A = E*(E/(Hiso+Hkin+E));
		double B = E*(E*dlam/normxsi);

		//ks(0,0) = E - A*n_n1[0]*n_n1[0] - B*(1.0 - n_n1[0]*n_n1[0]);
		//ks(1,1) = E - A*n_n1[1]*n_n1[1] - B*(1.0 - n_n1[1]*n_n1[1]);
		//ks(0,1) = -A*n_n1[0]*n_n1[1] - B*(-n_n1[0]*n_n1[1]);
		//ks(1,0) = ks(0,1);

		double EB = E-B;
		double BA = B-A;

		ks(0,0) = EB + BA*n_n1[0]*n_n1[0];
		ks(1,1) = EB + BA*n_n1[1]*n_n1[1];
		ks(0,1) = BA*n_n1[0]*n_n1[1];
		ks(1,0) = ks(0,1);

		//n_n1[0] *= dlam;
		//n_n1[1] *= dlam;

		// Update plastic strains
		//eP_n1[0] = eP_n[0] + n_n1[0];
		//eP_n1[1] = eP_n[1] + n_n1[1];

		// Update back stress
		//q_n1[0] = q_n[0] + Hkin*n_n1[0];
		//q_n1[1] = q_n[1] + Hkin*n_n1[1];

		// Update effective plastic strain
		//alpha_n1 = alpha_n + dlam;
	}

	return ks;
}

const Matrix&
Bidirectional::getInitialTangent(void)
{
  ks(0,0) = ks(1,1) = E;
  ks(0,1) = ks(1,0) = 0.0;

  return ks;
}

const Vector&
Bidirectional::getStressResultant(void)
{
	// Compute trial stress using elastic tangent
	s(0) = E*(e_n1[0]-eP_n[0]);
	s(1) = E*(e_n1[1]-eP_n[1]);

	static Vector xsi(2);

	// Predicted stress minus back stress
	xsi(0) = s(0) - q_n[0];
	xsi(1) = s(1) - q_n[1];

	double normxsi = xsi.Norm();

	// Current yield stress
	double sigY_n = sigY + alpha_n*Hiso;

	// Yield function
	double f_n1 = normxsi - sigY_n;

	// Elastic step
	if (f_n1 < 0.0) {
		// do nothing
	}

	// Plastic step
	else {
		// Consistency parameter
		double dlam = f_n1/(E+Hkin+Hiso);

		double n_n1[2];

		// Normal vector
		n_n1[0] = xsi(0)/normxsi;
		n_n1[1] = xsi(1)/normxsi;

		n_n1[0] *= dlam;
		n_n1[1] *= dlam;

		// Return stresses to yield surface
		s(0) -= E*n_n1[0];
		s(1) -= E*n_n1[1];

		// Update plastic strains
		eP_n1[0] = eP_n[0] + n_n1[0];
		eP_n1[1] = eP_n[1] + n_n1[1];

		// Update back stress
		q_n1[0] = q_n[0] + Hkin*n_n1[0];
		q_n1[1] = q_n[1] + Hkin*n_n1[1];

		// Update effective plastic strain
		alpha_n1 = alpha_n + dlam;
	}

	return s;
}

const Vector&
Bidirectional::getSectionDeformation(void)
{
	// Write to static variable for return
	s(0) = e_n1[0];
	s(1) = e_n1[1];

	return s;
}

int
Bidirectional::commitState(void)
{
	eP_n[0] = eP_n1[0];
	eP_n[1] = eP_n1[1];

	q_n[0] = q_n1[0];
	q_n[1] = q_n1[1];

	alpha_n = alpha_n1;

	return 0;
}

int
Bidirectional::revertToLastCommit(void)
{
	return 0;
}

int
Bidirectional::revertToStart(void)
{
	for (int i = 0; i < 2; i++) {
		eP_n[i]  = 0.0;
		eP_n1[i] = 0.0;
		q_n[i]  = 0.0;
		q_n1[i] = 0.0;
	}

	alpha_n  = 0.0;
	alpha_n1 = 0.0;

	return 0;
}

SectionForceDeformation*
Bidirectional::getCopy(void)
{
	Bidirectional *theCopy =
		new Bidirectional (this->getTag(), E, sigY, Hiso, Hkin);

	for (int i = 0; i < 2; i++) {
		theCopy->eP_n[i]  = eP_n[i];
		theCopy->eP_n1[i] = eP_n1[i];
		theCopy->q_n[i]  = q_n[i];
		theCopy->q_n1[i] = q_n1[i];
	}

	theCopy->alpha_n  = alpha_n;
	theCopy->alpha_n1 = alpha_n1;

	return theCopy;
}

const ID&
Bidirectional::getType(void)
{
	return code;
}

int
Bidirectional::getOrder(void) const
{
	return 2;
}

int 
Bidirectional::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigY;
  data(3) = Hiso;
  data(4) = Hkin;
  data(5) = eP_n[0];
  data(6) = eP_n[1];
  data(7) = q_n[0];
  data(8) = q_n[1];
  data(9) = alpha_n;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Bidirectional::sendSelf() - failed to send data\n";

  return res;
}

int 
Bidirectional::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(10);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "Bidirectional::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	E = data(1);
	sigY = data(2);
	Hiso = data(3);
	Hkin = data(4);
	eP_n[0] = data(5);
	eP_n[1] = data(6);
	q_n[0]  = data(7);
	q_n[1]  = data(8);
	alpha_n = data(9);

    // Set the trial state variables
    revertToLastCommit();
  }
    
  return res;
}

void
Bidirectional::Print(OPS_Stream &s, int flag)
{
	s << "Bidirectional, tag: " << this->getTag() << endln;
	s << "\tE:    " << E << endln;
	s << "\tsigY: " << sigY<< endln;
	s << "\tHiso: " << Hiso << endln;
	s << "\tHkin: " << Hkin << endln;
}
