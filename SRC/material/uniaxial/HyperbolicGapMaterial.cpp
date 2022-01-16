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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HyperbolicGapMaterial.cpp,v $

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
#include <stdlib.h>

#include <HyperbolicGapMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_HyperbolicGapMaterial)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 6) {
        opserr << "WARNING: Insufficient arguments\n";
        return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
        return 0;
    }

    double data[5];
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata,data)) {
        return 0;
    }

    UniaxialMaterial* mat = new HyperbolicGapMaterial(tag,data[0],data[1],data[2],data[3],data[4]);
    if (mat == 0) {
        opserr << "WARNING: failed to create Hyperbolicgapmaterial material\n";
        return 0;
    }

    return mat;
}

HyperbolicGapMaterial::HyperbolicGapMaterial(int tag, double kmax, double kur, double rf, double fult, double gap0)
    :UniaxialMaterial(tag,MAT_TAG_HyperbolicGapMaterial),
    Kmax(kmax), Kur(kur), Rf(rf), Fult(fult), gap(gap0)
{
    if (gap>=0) {
        opserr << "HyperbolicGapMaterial::HyperbolicGapMaterial -- Initial gap size must be negative for compression-only material, setting to negative\n";
        //exit(-1);
	gap = -gap;
    }
    if (Fult>0) {
        opserr << "HyperbolicGapMaterial::HyperbolicGapMaterial -- Fult must be negative for compression-only material, setting to negative\n";
        //exit(-1);
	Fult = -Fult;
    }
    if (Kmax == 0.0) {
        opserr << "HyperbolicGapMaterial::HyperbolicGapMaterial -- Kmax is zero, continuing with Kmax = Fult/0.002\n";
        if (Fult != 0.0)
            Kmax = fabs(Fult)/0.002;
        else {
            opserr << "HyperbolicGapMaterial::HyperbolicGapMaterial -- Kmax and Fult are zero\n";
            exit(-1);
        }
    }
    else

        // Initialize history variables
        this->revertToStart();
    this->revertToLastCommit();
}

HyperbolicGapMaterial::HyperbolicGapMaterial()
    :UniaxialMaterial(0,MAT_TAG_HyperbolicGapMaterial),
    Kmax(0.0), Kur(0.0), Rf(0.0), Fult(0.0), gap(0.0)
{
    // does nothing
}

HyperbolicGapMaterial::~HyperbolicGapMaterial()
{
    // does nothing
}

int 
HyperbolicGapMaterial::setTrialStrain(double strain, double strainRate)
{
    // set the trial strain
    Tstrain = strain;

    // set the strain increment
    dStrain = Tstrain - Cstrain;

    // loading on the envelope curve
    if (Tstrain <= CstrainMin) {
        TstrainMin = Tstrain;
        Ttangent = negEnvTangent(Tstrain);
        Tstress = negEnvStress(Tstrain);
    }
    else {
        // reloading
        if (dStrain < 0.0)
            negativeIncrement(dStrain);
        // unloading
        else if (dStrain > 0.0)
            positiveIncrement(dStrain);
    }
    return 0;
}

double 
HyperbolicGapMaterial::getStrain(void)
{
    return Tstrain;
}

double 
HyperbolicGapMaterial::getStress(void)
{
    return Tstress;
}

double 
HyperbolicGapMaterial::getTangent(void)
{
    return Ttangent;
}

double 
HyperbolicGapMaterial::getInitialTangent(void)
{
    return Kmax; 
}

void
HyperbolicGapMaterial::positiveIncrement(double dStrain)
{
    // store strain and stress at which unloading first takes place
    if (TstrainMin == Cstrain) {
        TonsetOfUnloadingStrain = TstrainMin;
        TonsetOfUnloadingStress = Cstress;
        TonsetOfReloadingStrain = TstrainMin - TonsetOfUnloadingStress/Kur;
    }

    // unloading
    Tstress = Cstress + Kur*dStrain;
    Ttangent = Kur;

    // check whether zero stress has been achieved
    if (Tstress > 0) {
        Tstress = 0.0;
        Ttangent = 0.0;
    }
}

void
HyperbolicGapMaterial::negativeIncrement(double dStrain)
{
    // reloading
    Tstress = Cstress + Kur*dStrain;
    if (Tstrain < TonsetOfReloadingStrain && Cstress == 0.0)
        Tstress = Kur*(Tstrain-TonsetOfReloadingStrain);
    Ttangent = Kur;

    if (Tstrain > TonsetOfReloadingStrain) {
        Tstress = 0.0;
        Ttangent = 0.0;
    }

}

int 
HyperbolicGapMaterial::commitState(void)
{
    ConsetOfUnloadingStrain = TonsetOfUnloadingStrain;
    ConsetOfReloadingStrain = TonsetOfReloadingStrain;
    Cstress = Tstress;
    Cstrain = Tstrain;
    CstrainMin = TstrainMin;
    return 0;
}

int 
HyperbolicGapMaterial::revertToLastCommit(void)
{
    TonsetOfUnloadingStrain = ConsetOfUnloadingStrain;
    TonsetOfReloadingStrain = ConsetOfReloadingStrain;
    Tstrain = Cstrain;
    Tstress = Cstress;
    TstrainMin = CstrainMin;

    return 0;
}

int 
HyperbolicGapMaterial::revertToStart(void)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = Kmax;
    TstrainMin = 0.0;
    CstrainMin = 0.0; 
    ConsetOfUnloadingStrain = 0.0;
    ConsetOfReloadingStrain = 0.0;
    TonsetOfReloadingStress = 0.0;

    return 0;
}

UniaxialMaterial *
HyperbolicGapMaterial::getCopy(void)
{
    HyperbolicGapMaterial *theCopy = new HyperbolicGapMaterial(this->getTag(),Kmax,Kur,Rf,Fult,gap);

    theCopy->ConsetOfUnloadingStrain = ConsetOfUnloadingStrain;
    theCopy->ConsetOfReloadingStrain = ConsetOfReloadingStrain;
    theCopy->TonsetOfReloadingStress = TonsetOfReloadingStress;
    theCopy->Cstress = Cstress;
    theCopy->Cstrain = Cstrain;
    theCopy->Ttangent = Ttangent;
    theCopy->TstrainMin = TstrainMin;
    theCopy->CstrainMin = TstrainMin;

    return theCopy;
}

int 
HyperbolicGapMaterial::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(15);
    data(0) = this->getTag();
    data(1) = Cstrain;
    data(2) = Kmax;
    data(3) = Kur;
    data(4) = Rf;
    data(5) = Fult;
    data(6) = gap;
    data(7) = ConsetOfUnloadingStrain;
    data(8) = ConsetOfReloadingStrain;
    data(9) = TonsetOfReloadingStress;
    data(10) = Cstress; 
    data(11) = Cstrain;
    data(12) = Ttangent;
    data(13) = TstrainMin;
    data(14) = CstrainMin;

    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "HyperbolicGapMaterial::sendSelf() - failed to send data\n";

    return res;
}

int 
HyperbolicGapMaterial::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(15);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "HyperbolicGapMaterial::recvSelf() - failed to recv data\n";
    else {
        this->setTag((int)data(0));
        Cstrain = data(1);
        Tstrain = Cstrain;
        Kmax = data(2);
        Kur = data(3);
        Rf = data(4);
        Fult = data(5);
        gap = data(6);
        ConsetOfUnloadingStrain = data(7);
        ConsetOfReloadingStrain = data(8);
        TonsetOfReloadingStress = data(9);
        Cstress = data(10); 
        Cstrain = data(11);
        Ttangent = data(12);
        TstrainMin = data(13);
        CstrainMin = data(14);

    }

    return res;
}

void 
HyperbolicGapMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "HyperbolicGapMaterial tag: " << this->getTag() << endln;
        s << "  Kmax: " << Kmax << endln;
        s << "  Kur: " << Kur << endln;
        s << "  Rf: " << Rf << endln;
        s << "  Fult: " << Fult << endln;
        s << "  initial gap: " << gap << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"HyperbolicGapMaterial\", ";
        s << "\"Kmax\": " << Kmax << ", ";
        s << "\"Kur\": " << Kur << ", ";
        s << "\"Rf\": " << Rf << ", ";
        s << "\"Fult\": " << Fult << ", ";
        s << "\"gap\": " << gap << "}";
    }
}

double
HyperbolicGapMaterial::negEnvStress(double strain)
{
    if (strain >= gap)
        return 0.0;
    else
        return (strain-gap)/((1/Kmax)+(Rf*(strain-gap)/Fult));
}

double
HyperbolicGapMaterial::negEnvTangent(double strain)
{
    if (strain > gap)
        return 0.0;
    else
        return 1/(Kmax*pow(((1/Kmax)+(Rf*(strain-gap)/Fult)),2));
}
