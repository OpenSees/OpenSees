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

// $Revision: 1.0 $
// $Date: 2025-01-03 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Alessia Casalucci - ASDEA Software, Italy
//
// todo...
//

#ifndef ASDSteel1DMaterial_h
#define ASDSteel1DMaterial_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>

class ASDSteel1DMaterial : public UniaxialMaterial
{

public:
	// life-cycle
	ASDSteel1DMaterial(
		int _tag,
		double _E,
		double _sy,
		double _H1,
		double _gamma1,
		double _H2,
		double _gamma2);
	ASDSteel1DMaterial();
	~ASDSteel1DMaterial();

	// info
	const char* getClassType(void) const { return "ASDSteel1DMaterial"; }

	// set strain
	int setTrialStrain(double v, double r = 0.0);

	// get state
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	UniaxialMaterial* getCopy(void);
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);
	double getEnergy(void);

private:
	int computeBaseSteel(double trial_strain);

 private:
	 // Young's modulus
	 double E = 0.0;
	 // Yield stress
	 double sy = 0.0;
	 // Chaboche kinematic hardening parameters
	 double H1 = 0.0;
	 double H2 = 0.0;
	 double gamma1 = 0.0;
	 double gamma2 = 0.0;
	 // state variables - backstresses
	 double alpha1 = 0.0;
	 double alpha1_commit = 0.0;
	 double alpha2 = 0.0;
	 double alpha2_commit = 0.0;
	 // state variables - plastic multiplier
	 double lambda = 0.0;
	 double lambda_commit = 0.0;
	 // strain, stress and tangent
	 double strain = 0.0;
	 double strain_commit = 0.0;
	 double stress = 0.0;
	 double stress_commit = 0.0;
	 double C = 0.0;
	 // other variables for output purposes
	 double energy = 0.0;
	 int niter = 0;
	 double Cnum = 0.0;
};

#endif
