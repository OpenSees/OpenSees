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

// $Revision: 1.3 $
// $Date: 2009-01-08 22:00:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HyperbolicGapMaterial.h,v $

// File: ~/material/HyperbolicGapMaterial.C
//
// Written: md
// Created: 04/2008
//
// Description: This file contains the class implementation for 
// HyperbolicGapMaterial.  This material is based on abutment stiffness
// models for bridge simulation proposed by Patrick Wilson and Ahmed Elgamal
// at UCSD.  The abutment stiffness models are based on large-scale abutment
// tests performed on the outdoor shaking table at UCSD.  The model is described
// for a 1.68 meter (5.5 ft) tall backwall height (typical size) and a 1 meter 
// wide section along the width of the abutment (to be scaled accordingly).
// The hyperbolic force-displacement model is based on work by Duncan and Mokwa
// (2001) and Shamsabadi et al. (2007) with calibrated parameters from UCSD
// abutment tests.  This model matches very well with test data up to 7.64 cm of 
// longitudinal displacement.
// Recommended values:
// Kmax = 20300 kN/m of abutment width
// Kur = Kmax for unloading/reloading stiffness
// Rf = 0.7
// Fult = 326 kN per meter of abutment width
// gap = 2.54 cm
// The model is implemented as a compression-only gap material, thus the values
// of Fult and gap should be negative.
//
#ifndef HyperbolicGapMaterial_h
#define HyperbolicGapMaterial_h

#include <UniaxialMaterial.h>

class HyperbolicGapMaterial : public UniaxialMaterial
{
public:
    HyperbolicGapMaterial(int tag, double Kmax, double Kur, double Rf, double Fult, double gap);
    HyperbolicGapMaterial();
    ~HyperbolicGapMaterial();

    const char *getClassType(void) const {return "HyperbolicGapMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

protected:

private:
    double Kmax;
    double Kur;
    double Rf;
    double Fult;
    double gap;
    double dStrain;
    double Tstress;
    double Tstrain;
    double Ttangent;
    double Cstrain;
    double Cstress;
    double TstrainMin;
    double CstrainMin;
    double TonsetOfUnloadingStrain;
    double ConsetOfUnloadingStrain;
    double TonsetOfUnloadingStress;
    double ConsetOfUnloadingStress;
    double TonsetOfReloadingStrain;
    double ConsetOfReloadingStrain;
    double TonsetOfReloadingStress;

    double negEnvStress(double strain);
    double negEnvTangent(double strain);

    void positiveIncrement(double dStrain);
    void negativeIncrement(double dStrain);
};

#endif
