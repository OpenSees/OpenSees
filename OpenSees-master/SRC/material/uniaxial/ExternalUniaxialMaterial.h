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
// $Date: 2019-01-26 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ExternalUniaxialMaterial.h,v $

// Written: M. Salehi
// Created: 26-1 2019
//


#ifndef ExternalUniaxialMaterial_h
#define ExternalUniaxialMaterial_h

#include <UniaxialMaterial.h>
typedef int(__stdcall *UMSetTrialStrain)(double strain, double strainRate);
typedef double(__stdcall *UMGetStress)(void);
typedef double(__stdcall *UMGetTangent)(void);
typedef double(__stdcall *UMGetInitialTangent)(void);
typedef double(__stdcall *UMGetStrain)(void);

typedef double(__stdcall *UMGetStrainRate)(void);
typedef double(__stdcall *UMGetDampTangent)(void);
typedef double(__stdcall *UMGetRho)(void);

typedef int(__stdcall *UMCommitState)();
typedef int(__stdcall *UMRevertToLastCommit)();
typedef int(__stdcall *UMRevertToStart)();
typedef UniaxialMaterial *(__stdcall *UMGetCopy)();
typedef void(__stdcall *UMPrint)(int);

class ExternalUniaxialMaterial : public UniaxialMaterial
{
public:
	ExternalUniaxialMaterial(int tag);


	void SetLinks(UMSetTrialStrain, UMGetStress, UMGetTangent,
		UMGetInitialTangent, UMGetDampTangent,
		UMGetStrain, UMGetStrainRate, UMGetRho, UMCommitState, UMRevertToLastCommit,
		UMRevertToStart, UMGetCopy, UMPrint);

	~ExternalUniaxialMaterial();

	const char *getClassType(void) const { return "ExternalUniaxialMaterial"; };

	int setTrialStrain(double strain, double strainRate = 0.0) {
		return _SetTrialStrain(strain, strainRate);
	};
	double getStrain(void) {
		return _GetStrain();
	};
	double getStress(void) {
		return _GetStress();
	};
	double getTangent(void) {
		return _GetTangent();
	};

	double getDampTangent(void) {
		return _GetDampTangent();
	};
	double getStrainRate(void) {
		return _GetStrainRate();
	};
	double getRho(void) {
		return _GetRho();
	};


	double getInitialTangent(void) {
		return _GetInitialTangent();
	};
	int commitState(void) {
		return _CommitState();
	};

	int revertToLastCommit(void) {
		return _RevertToLastCommit();
	};

	int revertToStart(void) {
		return _RevertToStart();
	};

	UniaxialMaterial *getCopy(void) {
		return _GetCopy();
	};

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0) {
		return _Print(flag);
	};

protected:

private:
	UMSetTrialStrain _SetTrialStrain;
	UMGetStress _GetStress;
	UMGetTangent _GetTangent;
	UMGetInitialTangent _GetInitialTangent;
	UMGetStrain _GetStrain;
	UMCommitState _CommitState;
	UMRevertToLastCommit _RevertToLastCommit;
	UMRevertToStart _RevertToStart;
	UMGetCopy _GetCopy;
	UMPrint _Print;
	UMGetDampTangent _GetDampTangent;
	UMGetStrainRate _GetStrainRate;
	UMGetRho _GetRho;
};

#endif

