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

#ifndef PluginNDMaterial_h
#define PluginNDMaterial_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the PluginNDMaterial wrapper

#include <NDMaterial.h>
#include <PluginFrameworkAPI.h>

class PluginMaterialDescriptor;

/**
The PluginNDMaterial class is a wrapper for external materials
using the PluginFrameworkAPI
*/

class PluginNDMaterial : public NDMaterial
{
	friend void* OPS_PluginNDMaterial(void);

	// constructor and destructor
public:
	PluginNDMaterial();
	PluginNDMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d);
	~PluginNDMaterial();

	// from TaggedObject
public:
	virtual void Print(OPS_Stream& s, int flag = 0);

	// from MovableObject
public:
	const char* getClassType() const;
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	int activateParameter(int parameterID);
	int setVariable(const char* variable, Information& info);
	int getVariable(const char* variable, Information& info);

	// from NDMaterial
public:
	int setTrialStrain(const Vector& strain);
	int setTrialStrain(const Vector& strain, const Vector& strainRate);
	int setTrialStrainIncr(const Vector& strain);
	int setTrialStrainIncr(const Vector& strain, const Vector& strainRate);
	const Vector& getStrain();
	const Vector& getStress();
	const Matrix& getTangent();
	const Matrix& getInitialTangent();
	double getRho();
	double getThermalTangentAndElongation(double& TempT, double& ET, double& Elong);
	double setThermalTangentAndElongation(double& TempT, double& ET, double& Elong);
	const Vector& getTempAndElong();
	int commitState();
	int revertToLastCommit();
	int revertToStart();
	NDMaterial* getCopy();
	NDMaterial* getCopy(const char* code);
	const char* getType() const;
	int getOrder() const;
	Response* setResponse(const char** argv, int argc, OPS_Stream& theOutputStream);
	int getResponse(int responseID, Information& info);
	const Vector& getStressSensitivity(int gradIndex, bool conditional);
	const Vector& getStrainSensitivity(int gradIndex);
	const Matrix& getTangentSensitivity(int gradIndex);
	const Matrix& getInitialTangentSensitivity(int gradIndex);
	const Matrix& getDampTangentSensitivity(int gradIndex);
	double getRhoSensitivity(int gradIndex);
	int commitSensitivity(const Vector& strainGradient, int gradIndex, int numGrads);

private:
	PluginMaterialDescriptor* m_descriptor;
	PluginMaterialData* m_data;
	double m_lch;
	bool m_lch_calculated;
};

#endif // !PluginNDMaterial_h

