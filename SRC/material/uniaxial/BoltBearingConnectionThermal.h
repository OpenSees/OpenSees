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

// $Revision: 1.7 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BoltBearingConnectionThermal.h,v $


#ifndef BoltBearingConnectionThermal_h
#define BoltBearingConnectionThermal_h

// Written: fmk 
// Created: 07/98
// Revision: A



#include <UniaxialMaterial.h>

class BoltBearingConnectionThermal : public UniaxialMaterial
{
public:
	BoltBearingConnectionThermal(int tag, double F, double Enot, double dbolt, double Abolt, double temp);
	BoltBearingConnectionThermal();
	~BoltBearingConnectionThermal();

	const char* getClassType(void) const { return "BoltBearingConnectionThermal"; };

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStress(void);
	double getStrain(void);
	double getTangent(void);
	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

protected:

private:
	double trialStrain;
	double trialStress;
	double trialTangent;
	double committedStrain;
	double committedStress;
	double committedTangent;
	double Fu0;
	double E0;
	double As;
	double db;
	double Temp;


};


#endif
