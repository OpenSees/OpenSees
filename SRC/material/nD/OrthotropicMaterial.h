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
                                                                        
// $Revision: 1.4 $
// $Date: 2020-04-19 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/OrthotropicMaterial.h,v $

// Davide Raino, Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#ifndef OrthotropicMaterial_h
#define OrthotropicMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class OrthotropicMaterial : public NDMaterial 
{
public:
	// life-cycle
	OrthotropicMaterial(
		int tag, 
		NDMaterial &theIsoMat,
		double Ex, double Ey, double Ez, double Gxy, double Gyz, double Gzx,
		double vxy, double vyz, double vzx,
		double Asigmaxx, double Asigmayy, double Asigmazz, double Asigmaxyxy, double Asigmayzyz, double Asigmaxzxz);
	OrthotropicMaterial();
	~OrthotropicMaterial();

	// info
	const char* getClassType(void) const { return "OrthotropicMaterial"; };

	// density
	double getRho(void);

	// set state
	int setTrialStrain(const Vector &strain);

	// get state
	const Vector &getStrain(void);
	const Vector &getStress(void);
	const Matrix &getTangent(void);
	const Matrix &getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	NDMaterial *getCopy(void);
	NDMaterial* getCopy(const char* code);
	const char *getType(void) const;
	int getOrder(void) const;
	void Print(OPS_Stream &s, int flag=0);

	// send/recv self
	virtual int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	Response* setResponse(const char** argv, int argc, OPS_Stream& s);

private:
	// the mapped isotropic material
	NDMaterial *theIsotropicMaterial = nullptr;
	// the strain in the real orthotropic space
	Vector epsilon = Vector(6);
	// strain tensor map
	Matrix Aepsilon = Matrix(6, 6);
	// inverse of stress tensor map (saved as vector, it's diagonal)
	Vector Asigma_inv = Vector(6);
};
#endif