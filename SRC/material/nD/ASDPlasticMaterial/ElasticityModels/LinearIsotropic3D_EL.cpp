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

// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial
//
// Fully general templated material class for plasticity modeling

#include "LinearIsotropic3D_EL.h"
#include "Vector.h"

VoigtMatrix LinearIsotropic3D_EL::Ee;


LinearIsotropic3D_EL::LinearIsotropic3D_EL(double E, double nu) : ElasticityBase<LinearIsotropic3D_EL>::ElasticityBase()  // Note the full-qualification of ElasticityBase through the scope resolution operator (::)
{
	lambda = ( nu * E ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
	mu = E / ( 2.0 * ( 1.0 + nu ) );
}


VoigtMatrix& LinearIsotropic3D_EL::operator()(const VoigtVector& stress) //See note on base class
{
	Ee.setZero(); //Zero it. It may have values from another instance with different parameters;

	const double mu2 = mu*mu; 

	Ee(0, 0) = Ee(1, 1) = Ee(2, 2) = mu2;
	Ee(0, 1) = Ee(1, 0) = Ee(0, 2) = Ee(2, 0) = Ee(1, 2) = Ee(2, 1) = lambda;
	Ee(3, 3) = mu;
	Ee(4, 4) = mu;
	Ee(5, 5) = mu;

	return Ee;
}

int LinearIsotropic3D_EL::sendSelf(int commitTag, Channel &theChannel)
{
	// static Vector data(2);
	// data(0) = lambda;
	// data(1) = mu;

	// if (theChannel.sendVector(0, commitTag, data) != 0)
	// {
	// 	cerr << "LinearIsotropic3D_EL::sendSelf() - Failed to send data. " << endl;
	// 	return -1;
	// }

	return 0;
}

int LinearIsotropic3D_EL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// static Vector data(2);

	// if (theChannel.receiveVector(0, commitTag, data) != 0)
	// {
	// 	cerr << "LinearIsotropic3D_EL::recvSelf() - Failed to receive data. " << endl;
	// 	return -1;
	// }

	// lambda = data(0);
	// mu = data(1);

	return 0;
}
