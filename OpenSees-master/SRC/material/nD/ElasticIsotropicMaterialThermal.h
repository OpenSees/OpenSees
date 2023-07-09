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

// $Revision: 1.13 $
// $Date: 2006-09-05 21:21:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicMaterialThermal.h,v $


#ifndef ElasticIsotropicMaterialThermal_h
#define ElasticIsotropicMaterialThermal_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for ElasticIsotropicMaterialThermalModel.
// ElasticIsotropicMaterialThermalModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) ElasticIsotropicMaterialThermal.h, revA"

//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticIsotropicMaterialThermal : public NDMaterial
{
public:
	// Only called by subclasses to pass their tags to NDMaterialModel
	ElasticIsotropicMaterialThermal(int tag, int classTag, double E, double nu, double rho = 0.0, double alpha = 0, int softindex = 0);

	// Called by clients
	ElasticIsotropicMaterialThermal(int tag, double E, double nu, double rho = 0.0, double alpha = 0, int softindex=0);

	// For parallel processing
	ElasticIsotropicMaterialThermal(void);

	virtual ~ElasticIsotropicMaterialThermal(void);

	virtual const char *getClassType(void) const { return "ElasticIsotropicMaterialThermal"; };

	virtual double getRho();

	virtual int setTrialStrain(const Vector &v);
	virtual int setTrialStrain(const Vector &v, const Vector &r);
	virtual int setTrialStrainIncr(const Vector &v);
	virtual int setTrialStrainIncr(const Vector &v, const Vector &r);
	virtual const Matrix &getTangent(void);
	virtual const Matrix &getInitialTangent(void);
	virtual const Vector &getStress(void);
	virtual const Vector &getStrain(void);

	virtual int commitState(void);
	virtual int revertToLastCommit(void);
	virtual int revertToStart(void);

	// Create a copy of material parameters AND state variables
	// Called by GenericSectionXD
	virtual NDMaterial *getCopy(void);

	// Create a copy of just the material parameters
	// Called by the continuum elements
	virtual NDMaterial *getCopy(const char *type);

	// Return a string indicating the type of material model
	virtual const char *getType(void) const;

	virtual int getOrder(void) const;

	virtual int sendSelf(int commitTag, Channel &theChannel);
	virtual int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

	virtual int setParameter(const char **argv, int argc, Parameter &param);
	virtual int updateParameter(int parameterID, Information &info);



protected:
	double E;	// Elastic modulus
	double v;	// Poisson ratio
	double rho; //mass per unit 3D volume
	double Alpha;  // coefficient of thermal expansion
	int SoftIndex; // softening behaviour
private:
	// Uncomment when this material model is "sensitized"
	//int parameterID;
};


#endif