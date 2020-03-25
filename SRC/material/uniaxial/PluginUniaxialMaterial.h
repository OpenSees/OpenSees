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

#ifndef PluginUniaxialMaterial_h
#define PluginUniaxialMaterial_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the PluginUniaxialMaterial wrapper

#include <UniaxialMaterial.h>
#include <PluginFrameworkAPI.h>

class PluginMaterialDescriptor;

/**
The PluginUniaxialMaterial class is a wrapper for external materials
using the PluginFrameworkAPI
*/
class PluginUniaxialMaterial : public UniaxialMaterial
{
	friend void* OPS_PluginUniaxialMaterial(void);
	
	// constructor and destructor
public:
	PluginUniaxialMaterial();
	PluginUniaxialMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d);
	~PluginUniaxialMaterial();

	// from TaggedObject
public:
	virtual void Print(OPS_Stream& s, int flag = 0);

	// from MovableObject
public:
	virtual const char* getClassType() const;
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
	virtual int setParameter(const char** argv, int argc, Parameter& param);
	virtual int updateParameter(int parameterID, Information& info);
	virtual int activateParameter(int parameterID);
	virtual int setVariable(const char* variable, Information&);
	virtual int getVariable(const char* variable, Information&);

	// from UniaxialMaterial
public:
	virtual int setTrialStrain(double strain, double strainRate = 0);
	virtual int setTrialStrain(double strain, double temperature, double strainRate);
	virtual double getStrain();
	virtual double getStrainRate();
	virtual double getStress();
	virtual double getTangent();
	virtual double getInitialTangent();
	virtual double getDampTangent();
	virtual double getRho();
	virtual int commitState();
	virtual int revertToLastCommit();
	virtual int revertToStart();
	virtual UniaxialMaterial* getCopy();
	virtual Response* setResponse(const char** argv, int argc, OPS_Stream& theOutputStream);
	virtual int getResponse(int responseID, Information& matInformation);
	virtual int getResponseSensitivity(int responseID, int gradIndex, Information& info);
	virtual bool hasFailed();
	virtual double getStressSensitivity(int gradIndex, bool conditional);
	virtual double getStrainSensitivity(int gradIndex);
	virtual double getTangentSensitivity(int gradIndex);
	virtual double getInitialTangentSensitivity(int gradIndex);
	virtual double getDampTangentSensitivity(int gradIndex);
	virtual double getRhoSensitivity(int gradIndex);
	virtual int    commitSensitivity(double strainGradient, int gradIndex, int numGrads);
	virtual double getEnergy();

private:
	PluginMaterialDescriptor* m_descriptor;
	PluginMaterialData* m_data;
	double m_lch;
	bool m_lch_calculated;
};

#endif // PluginUniaxialMaterial_h
