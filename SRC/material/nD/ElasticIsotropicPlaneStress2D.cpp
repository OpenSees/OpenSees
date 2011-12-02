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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.cpp,v $
                                                                        
                                                                        

#include <ElasticIsotropicPlaneStress2D.h>           
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicPlaneStress2D::sigma(3);

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D
(int tag, double E, double nu) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlaneStress2d, E, nu),
 D(3,3), Tepsilon(3), Cepsilon(3)
{
	this->update();
}

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlaneStress2d, 0.0, 0.0),
 D(3,3), Tepsilon(3), Cepsilon(3)
{

}

ElasticIsotropicPlaneStress2D::~ElasticIsotropicPlaneStress2D ()
{

}

int
ElasticIsotropicPlaneStress2D::setTrialStrain (const Vector &v)
{
	Tepsilon = v;

	return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrain (const Vector &v, const Vector &r)
{
	Tepsilon = v;

	return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrainIncr (const Vector &v)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, v, 1.0);

	return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	Tepsilon = Cepsilon;
	Tepsilon.addVector(1.0, v, 1.0);

	return 0;
}

const Matrix&
ElasticIsotropicPlaneStress2D::getTangent (void)
{
	return D;
}

const Vector&
ElasticIsotropicPlaneStress2D::getStress (void)
{
	//sigma = D*epsilon;
	sigma(0) = D(0,0)*Tepsilon(0) + D(0,1)*Tepsilon(1);
	sigma(1) = D(1,0)*Tepsilon(0) + D(1,1)*Tepsilon(1);
	sigma(2) = D(2,2)*Tepsilon(2);
	
	return sigma;
}

const Vector&
ElasticIsotropicPlaneStress2D::getStrain (void)
{
	return Tepsilon;
}

int
ElasticIsotropicPlaneStress2D::commitState (void)
{
	Cepsilon = Tepsilon;

	return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToLastCommit (void)
{
	Tepsilon = Cepsilon;

	return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToStart (void)
{
	Cepsilon.Zero();

	return 0;
}

NDMaterial*
ElasticIsotropicPlaneStress2D::getCopy (void)
{
	ElasticIsotropicPlaneStress2D *theCopy =
		new ElasticIsotropicPlaneStress2D (this->getTag(), E, v);

	theCopy->Cepsilon = Cepsilon;
	// D is created in the constructor call

	return theCopy;
}

const char*
ElasticIsotropicPlaneStress2D::getType (void) const
{
	return "ElasticIsotropicPlaneStress2D";
}

int
ElasticIsotropicPlaneStress2D::getOrder (void) const
{
	return 3;
}

int 
ElasticIsotropicPlaneStress2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(6);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
	data(3) = Cepsilon(0);
	data(4) = Cepsilon(1);
	data(5) = Cepsilon(2);

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropicPlaneStress2D::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropicPlaneStress2D::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(6);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropicPlaneStress2D::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    E = data(1);
	v = data(2);
	Cepsilon(0) = data(3);
	Cepsilon(1) = data(4);
	Cepsilon(2) = data(5);

	this->update();
	
	return res;
}

void 
ElasticIsotropicPlaneStress2D::update(void)
{
	// Set up the elastic constant matrix for plane stress
	D.Zero();
	D(0,0) = 1.0;
	D(0,1) = D(1,0) = v;
	D(1,1) = 1.0;
	D(2,2) = 0.5*(1.0-v);

	D = D * E/(1-v*v);
}
