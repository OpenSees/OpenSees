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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ExternalNDMaterial.h,v $

// Written: M. Salehi
// Created: 26-1 2019
//


#ifndef ExternalNDMaterial_h
#define ExternalNDMaterial_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

typedef NDMaterial *(__stdcall *NDM_GetCopy)();
typedef NDMaterial *(__stdcall *NDM_GetCopy_Type)(const char *type);
typedef void(__stdcall *NDM_Print)(int);
typedef double(__stdcall *NDM_GetRho)();
typedef int(__stdcall *NDM_SetTrialStrain_V)(const Vector &v);
typedef int(__stdcall *NDM_SetTrialStrain_VR)(const Vector &v, const Vector &r);
typedef int(__stdcall *NDM_SetTrialStrainIncr_V)(const Vector &v);
typedef int(__stdcall *NDM_SetTrialStrainIncr_VR)(const Vector &v, const Vector &r);
typedef Matrix * (__stdcall *NDM_GetTangent)(void);
typedef Matrix * (__stdcall *NDM_GetInitialTangent)(void);
typedef Vector * (__stdcall *NDM_GetStress)(void);
typedef Vector * (__stdcall *NDM_GetStrain)(void);
typedef int(__stdcall *NDM_RevertToStart)(void);
typedef int(__stdcall *NDM_CommitState)(void);
typedef int(__stdcall *NDM_RevertToLastCommit)(void);
typedef const char *(__stdcall *NDM_GetType) (void);
typedef int(__stdcall *NDM_GetOrder) (void);

class ExternalNDMaterial : public NDMaterial
{
public:

	// Constructor
	ExternalNDMaterial (int tag);		

	// Destructor
	~ExternalNDMaterial();				  

	double getRho(void);

	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector &v, const Vector &r);
	int setTrialStrainIncr(const Vector &v);
	int setTrialStrainIncr(const Vector &v, const Vector &r);
	const Matrix &getTangent(void);
	const Matrix &getInitialTangent(void);
	
	Response *setResponse(const char **argv, int argc, OPS_Stream &theOutputStream) { return 0; };
	int getResponse(int responseID, Information &matInformation) { return 0; };

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial *getCopy(void);
	NDMaterial *getCopy(const char *type);

	void Print(OPS_Stream &s, int flag = 0);
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	const char *getType(void) const;
	int getOrder(void) const;

	// Functions used for recorders
	const Vector &getStress(void);
	const Vector &getStrain(void);
	void SetLinks(NDM_GetCopy getCopy,
	NDM_GetCopy_Type getCopy_Type,
	NDM_Print print,
	NDM_GetRho getRho,
	NDM_SetTrialStrain_V setTrialStrain_V,
	NDM_SetTrialStrain_VR setTrialStrain_VR,
	NDM_SetTrialStrainIncr_V setTrialStrainIncr_V,
	NDM_SetTrialStrainIncr_VR setTrialStrainIncr_VR,
	NDM_GetTangent getTangent,
	NDM_GetInitialTangent getInitialTangent,
	NDM_GetStress getStress,
	NDM_GetStrain getStrain,
	NDM_RevertToStart revertToStart,
	NDM_CommitState commitState,
	NDM_RevertToLastCommit revertToLastCommit,
	NDM_GetType getType,
	NDM_GetOrder getOrder) {
		this->_NDMGetCopy = getCopy;
		this->_NDMGetCopy_Type = getCopy_Type;
		this->_NDMPrint = print;
		this->_NDMGetRho = getRho;
		this->_NDMSetTrialStrain_V = setTrialStrain_V;
		this->_NDMSetTrialStrain_VR = setTrialStrain_VR;
		this->_NDMSetTrialStrainIncr_V = setTrialStrainIncr_V;
		this->_NDMSetTrialStrainIncr_VR = setTrialStrainIncr_VR;
		this->_NDMGetTangent = getTangent;
		this->_NDMGetInitialTangent = getInitialTangent;
		this->_NDMGetStress = getStress;
		this->_NDMGetStrain = getStrain;
		this->_NDMRevertToStart = revertToStart;
		this->_NDMCommitState = commitState;
		this->_NDMRevertToLastCommit = revertToLastCommit;
		this->_NDMGetType = getType;
		this->_NDMGetOrder = getOrder;



	}
protected:

private:
	NDM_GetCopy _NDMGetCopy;
	NDM_GetCopy_Type _NDMGetCopy_Type;
	NDM_Print _NDMPrint;
	NDM_GetRho _NDMGetRho;
	NDM_SetTrialStrain_V _NDMSetTrialStrain_V;
	NDM_SetTrialStrain_VR _NDMSetTrialStrain_VR;
	NDM_SetTrialStrainIncr_V _NDMSetTrialStrainIncr_V;
	NDM_SetTrialStrainIncr_VR _NDMSetTrialStrainIncr_VR;
	NDM_GetTangent _NDMGetTangent;
	NDM_GetInitialTangent _NDMGetInitialTangent;
	NDM_GetStress _NDMGetStress;
	NDM_GetStrain _NDMGetStrain;
	NDM_RevertToStart _NDMRevertToStart;
	NDM_CommitState _NDMCommitState;
	NDM_RevertToLastCommit _NDMRevertToLastCommit;
	NDM_GetType _NDMGetType;
	NDM_GetOrder _NDMGetOrder;
};

#endif
