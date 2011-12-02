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
                                                                        
// $Revision: 1.1 $
// $Date: 2000-12-18 09:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BidirectionalMaterial.cpp,v $

#include <BidirectionalMaterial.h>           
#include <Channel.h>
#include <Tensor.h>

Vector BidirectionalMaterial::sigma(2);
Matrix BidirectionalMaterial::D(2,2);

BidirectionalMaterial::BidirectionalMaterial
(int tag, double e, double s, double hi, double hk) :
 NDMaterial (tag, ND_TAG_Bidirectional),
	 E(e), sigY(s), Hi(hi), Hk(hk)
{
	for (int i = 0; i < 2; i++) {
		epsPn[i] = 0.0;
		epsPn1[i] = 0.0;
		backn[i] = 0.0;
		backn1[i] = 0.0;
	}

	effn = 0.0;
	effn1 = 0.0;
}

BidirectionalMaterial::BidirectionalMaterial():
 NDMaterial (0, ND_TAG_Bidirectional),
	 E(0.0), sigY(0.0), Hi(0.0), Hk(0.0)
{
	for (int i = 0; i < 2; i++) {
		epsPn[i] = 0.0;
		epsPn1[i] = 0.0;
		backn[i] = 0.0;
		backn1[i] = 0.0;
	}

	effn = 0.0;
	effn1 = 0.0;
}

BidirectionalMaterial::~BidirectionalMaterial()
{

}

int
BidirectionalMaterial::setTrialStrain(const Vector &v)
{
	epsn1[0] = v(0);
	epsn1[1] = v(1);

	return 0;
}

int
BidirectionalMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

int
BidirectionalMaterial::setTrialStrainIncr(const Vector &v)
{
	return 0;
}

int
BidirectionalMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
	return 0;
}

const Matrix&
BidirectionalMaterial::getTangent(void)
{
	// Compute trial stress using elastic tangent
	sigma(0) = E*(epsn1[0]-epsPn[0]);
	sigma(1) = E*(epsn1[1]-epsPn[1]);

	static Vector xsi(2);

	// Predicted stress minus back stress
	xsi(0) = sigma(0) - backn[0];
	xsi(1) = sigma(1) - backn[1];

	double normxsi = xsi.Norm();

	// Current yield stress
	double sigYn = sigY + effn*Hi;

	// Yield function
	double Fn1 = normxsi - sigYn;

	// Elastic step
	if (Fn1 < 0.0) {
		D(0,0) = D(1,1) = E;
		D(0,1) = D(1,0) = 0.0;

		return D;
	}

	// Plastic step
	else {
		// Consistency parameter
		double dlam = Fn1/(E+Hk+Hi);

		double nn1[2];

		// Normal vector
		nn1[0] = xsi(0)/normxsi;
		nn1[1] = xsi(1)/normxsi;

		double A = E*(E/(Hi+Hk+E));
		double B = E*(E*dlam/normxsi);

		//D(0,0) = E - A*nn1[0]*nn1[0] - B*(1.0 - nn1[0]*nn1[0]);
		//D(1,1) = E - A*nn1[1]*nn1[1] - B*(1.0 - nn1[1]*nn1[1]);
		//D(0,1) = -A*nn1[0]*nn1[1] - B*(-nn1[0]*nn1[1]);
		//D(1,0) = D(0,1);

		double EB = E-B;
		double BA = B-A;

		D(0,0) = EB + BA*nn1[0]*nn1[0];
		D(1,1) = EB + BA*nn1[1]*nn1[1];
		D(0,1) = BA*nn1[0]*nn1[1];
		D(1,0) = D(0,1);

		nn1[0] *= dlam;
		nn1[1] *= dlam;

		// Update plastic strains
		epsPn1[0] = epsPn[0] + nn1[0];
		epsPn1[1] = epsPn[1] + nn1[1];

		// Update back stress
		backn1[0] = backn[0] + Hk*nn1[0];
		backn1[1] = backn[1] + Hk*nn1[1];

		// Update effective plastic strain
		effn1 = effn + dlam;

		return D;
	}
}

const Vector&
BidirectionalMaterial::getStress(void)
{
	// Compute trial stress using elastic tangent
	sigma(0) = E*(epsn1[0]-epsPn[0]);
	sigma(1) = E*(epsn1[1]-epsPn[1]);

	static Vector xsi(2);

	// Predicted stress minus back stress
	xsi(0) = sigma(0) - backn[0];
	xsi(1) = sigma(1) - backn[1];

	double normxsi = xsi.Norm();

	// Current yield stress
	double sigYn = sigY + effn*Hi;

	// Yield function
	double Fn1 = normxsi - sigYn;

	// Elastic step
	if (Fn1 < 0.0)
		return sigma;

	// Plastic step
	else {
		// Consistency parameter
		double dlam = Fn1/(E+Hk+Hi);

		double nn1[2];

		// Normal vector
		nn1[0] = xsi(0)/normxsi;
		nn1[1] = xsi(1)/normxsi;

		nn1[0] *= dlam;
		nn1[1] *= dlam;

		// Return stresses to yield surface
		sigma(0) -= E*nn1[0];
		sigma(1) -= E*nn1[1];

		// Update plastic strains
		epsPn1[0] = epsPn[0] + nn1[0];
		epsPn1[1] = epsPn[1] + nn1[1];

		// Update back stress
		backn1[0] = backn[0] + Hk*nn1[0];
		backn1[1] = backn[1] + Hk*nn1[1];

		// Update effective plastic strain
		effn1 = effn + dlam;

		return sigma;
	}
}

const Vector&
BidirectionalMaterial::getStrain(void)
{
	// Write to static variable for return
	sigma(0) = epsn1[0];
	sigma(1) = epsn1[1];

	return sigma;
}

int
BidirectionalMaterial::commitState(void)
{
	for (int i = 0; i < 2; i++) {
		epsPn[i] = epsPn1[i];
		backn[i] = backn1[i];
	}

	effn = effn1;

	return 0;
}

int
BidirectionalMaterial::revertToLastCommit(void)
{
	return 0;
}

int
BidirectionalMaterial::revertToStart(void)
{
	for (int i = 0; i < 2; i++) {
		epsPn[i] = 0.0;
		epsPn1[i] = 0.0;
		backn[i] = 0.0;
		backn1[i] = 0.0;
	}

	effn = 0.0;
	effn1 = 0.0;

	return 0;
}

NDMaterial*
BidirectionalMaterial::getCopy(void)
{
	BidirectionalMaterial *theCopy =
		new BidirectionalMaterial (this->getTag(), E, sigY, Hi, Hk);

	for (int i = 0; i < 2; i++) {
		theCopy->epsPn[i] = epsPn[i];
		theCopy->backn[i] = backn[i];
	}

	theCopy->effn = effn;

	return theCopy;
}

NDMaterial*
BidirectionalMaterial::getCopy(const char *code)
{
    if (strcmp(code,this->getType()) == 0)
		return this->getCopy();

    else {
		g3ErrorHandler->fatal("%s -- failed to get copy of %s",
			"BidirectionalMaterial::getCopy", code);

		return 0;
    }
}

const char*
BidirectionalMaterial::getType(void) const
{
	return "ForceDeformation";
}

int
BidirectionalMaterial::getOrder(void) const
{
	return 2;
}

int 
BidirectionalMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigY;
  data(3) = Hi;
  data(4) = Hk;
  data(5) = epsPn[0];
  data(6) = epsPn[1];
  data(7) = backn[0];
  data(8) = backn[1];
  data(9) = effn;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "BidirectionalMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
BidirectionalMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(10);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      cerr << "BidirectionalMaterial::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	E = data(1);
	sigY = data(2);
	Hi = data(3);
	Hk = data(4);
	epsPn[0] = data(5);
	epsPn[1] = data(6);
	backn[0] = data(7);
	backn[1] = data(8);
	effn = data(9);

    // Set the trial state variables
    revertToLastCommit();
  }
    
  return res;
}
void
BidirectionalMaterial::Print(ostream &s, int flag)
{
	s << "BidirectionalMaterial, tag: " << this->getTag() << endl;
	s << "\tE:  " << E << endl;
	s << "\tsigY:  " << sigY<< endl;
	s << "\tHi:  " << Hi << endl;
	s << "\tHk:  " << Hk << endl;
}
