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

/**
todo: global material
	flag computed_once -> if false -> call update
		revert to start -> call update (TO SET INITIAL DATA)

*/

class ASDSteel1DMaterialPIMPL;

class ASDSteel1DMaterial : public UniaxialMaterial
{
public:
	class InputParameters {
	public:
		// Young's modulus
		double E = 0.0;
		// Yield stress
		double sy = 0.0;
		// Chaboche kinematic hardening parameters
		double H1 = 0.0;
		double H2 = 0.0;
		double gamma1 = 0.0;
		double gamma2 = 0.0;
		// misc
		bool implex = false;
		// buckling
		double radius = 0.0;
		double length = 0.0;
		double lch_element = 0.0;

		//convergence
		double K_alpha = 0.0;
		double max_iter = 0.0;
		double tolU = 0.0;
		double tolR = 0.0;
		double p = 0.0;
		double n = 0.0;
		// counter
		static constexpr int NDATA = 16;
	};

public:
	// life-cycle
	ASDSteel1DMaterial(
		int _tag,
		const InputParameters& _params);
	ASDSteel1DMaterial();
	ASDSteel1DMaterial(const ASDSteel1DMaterial& other);
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
	int homogenize(bool do_implex);

public:
	void computeAlphaCr();

 private:	
	 // common input parameters
	 InputParameters params;
	 // state variables - implex
	 double dtime_n = 0.0;
	 double dtime_n_commit = 0.0;
	 bool commit_done = false;
	 // strain, stress and tangent (homogenized)
	 double strain = 0.0;
	 double strain_commit = 0.0;
	 double stress = 0.0;
	 double stress_commit = 0.0;
	 double C = 0.0;
	 double stress_rve = 0.0;
	 double stress_rve_commit = 0.0;
	 double C_rve = 0.0;
	 
	 // other variables for output purposes
	 double energy = 0.0;

	 //for regularization
	 bool use_regularization = true;
	 double alpha_cr = 1.0;
	 // private implementation
	 ASDSteel1DMaterialPIMPL* pdata = nullptr;
};

#endif
