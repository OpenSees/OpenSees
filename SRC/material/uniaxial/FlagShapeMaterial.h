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

// $Revision: 0 $
// $Date: 2021-02-22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FlagShapeMaterial.h,v $

#ifndef FlagShapeMaterial_h
#define FlagShapeMaterial_h

// Written: CLLee
// Created: Feb 2021
// Reference: HardeningMaterial
//
// Description: This file contains the class definition for 
// FlagShapeMaterial.  FlagShapeMaterial provides the abstraction
// for a one-dimensional rate-independent plasticity model
// with kinematic hardening and a shift of backstress using beta variable

#include <UniaxialMaterial.h>
#include <Matrix.h>

class FlagShapeMaterial : public UniaxialMaterial
{
public:
    FlagShapeMaterial(int tag, double E, double fy,
        double Eh, double beta = 0.0);
    FlagShapeMaterial();

    const char* getClassType(void) const { return "FlagShapeMaterial"; };

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) { return E; };

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
    // Material parameters
    double E;	// Elastic modulus
    double fy;	// Yield stress
    double Eh;	// Post-yield modulus
    double beta;

    // Committed history variables
    double CplasticStrain;	// Committed plastic strain
    double CbackStress;		// Committed back stress

    // Trial history variables
    double TplasticStrain;	// Trial plastic strain
    double TbackStress;		// Trial internal hardening variable

    // Trial state variables
    double Tstrain;		// Trial strain
    double Tstress;		// Trial stress
    double Ttangent;	// Trial tangent

};

void* OPS_FlagShapeMaterial();
#endif

