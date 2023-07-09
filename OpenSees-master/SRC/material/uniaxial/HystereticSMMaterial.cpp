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

// $Revision: 1.19 $
// $Date: 2008-12-18 23:40:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticSMMaterial.cpp,v $

// Written: MHS
// Created: July 2000
//
// Extended to Multi-Point: SilviaMazzoni
// Started: August 2022
// Updated: November 2022
// NOTE: only the first two segments can have slope >=0. all other segments the slope must be <= 0!
//
// Description: This file contains the implementation of 
// HystereticSMMaterial.  HystereticSMMaterial is
// a one-dimensional HystereticSM model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.
// with 3 arguments, it reverts back to hysteretic


#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <vector>

#include <HystereticSMMaterial.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

#include <MaterialResponse.h> // added by Silvia Mazzoni (silviamazzoni@yahoo.com) on 8/27/2022

static int numHystereticSMMaterials = 0;

void*
OPS_HystereticSMMaterial(void)
{

    if (numHystereticSMMaterials == 0) {
        numHystereticSMMaterials++;
        OPS_Error("HystereticSM: multi-point envelope + DCR recorders  - Code by Silvia Mazzoni, 2022 (silviamazzoni@yahoo.com) \n", 1);
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;


    int numArgs = OPS_GetNumRemainingInputArgs();
    int numData = 1;
    std::vector<double> defoLimitStates;
    int nDefoLimitStates = 0;
    std::vector<double> forceLimitStates;
    int nForceLimitStates = 0;
    int numOptionalArgs = 0;
    int numdata = 1;

    int loc = 2;



    // Read the optional arguments first
    while (OPS_GetNumRemainingInputArgs() > 0) {
        std::string theType = OPS_GetString();
        if (theType == "-defoLimitStates") {
            numOptionalArgs++;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    // OPS_ResetCurrentInputArg(loc);
                    break;
                }
                defoLimitStates.push_back(val);
                //loc++;
                nDefoLimitStates++;
                numOptionalArgs++;
            }
        }
    }

    OPS_ResetCurrentInputArg(-numArgs);
    while (OPS_GetNumRemainingInputArgs() > 0) {
        std::string theType = OPS_GetString();
        if (theType == "-forceLimitStates") {
            numOptionalArgs++;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    // OPS_ResetCurrentInputArg(loc);
                    break;
                }
                forceLimitStates.push_back(val);
                //loc++;
                nForceLimitStates++;
                numOptionalArgs++;
            }
        }
    }


    OPS_ResetCurrentInputArg(-numArgs);
    double degEnvFactor = 0;
    OPS_ResetCurrentInputArg(-numArgs);
    while (OPS_GetNumRemainingInputArgs() > 0) {
        std::string theType = OPS_GetString();
        if (theType == "-degEnvFactor") {
            numOptionalArgs++;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &degEnvFactor) < 0)
                    return 0;
                numOptionalArgs++;
            }
        }
    }

    OPS_ResetCurrentInputArg(-numArgs);

    int numargs0 = numArgs;
    numArgs = numArgs - numOptionalArgs;


    if (numArgs != 34 && numArgs != 33 && numArgs != 25 && numArgs != 26 && numArgs != 18 && numArgs != 17 && numArgs != 14 && numArgs != 13) {
        opserr << "numargs0 HystereticSM " << numargs0 << endln;
        opserr << "numOptionalArgs HystereticSM " << numOptionalArgs << endln;
        opserr << "numArgs HystereticSM " << numArgs << endln;
        opserr << "Want: uniaxialMaterial HystereticSM tag? mom1p? rot1p? mom2p? rot2p? <mom3p? rot3p? mom4p? rot4p? mom5p? rot5p? mom6p? rot6p? mom7p? rot7p?> "
            << "\nmom1n? rot1n? mom2n? rot2n? <mom3n? rot3n? mom4n? rot4n? mom5n? rot5n? mom6n? rot6n? mom7n? rot7n?> pinchX? pinchY? damfc1? damfc2? <beta?> "
            << "\n<-degEnvFactor degEnvFactor?> "
            << "\n<-defoLimitStates lsD1? <lsD2?>...> "
            << "\n<-forceLimitStates lsF1? <lsF2?>...> ";
        return 0;
    }

    int iData[1];
    double dData[33];
    for (int i = 0; i < 33; i++)
        dData[i] = 0.0;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid tag for uniaxialMaterial HystereticSM" << endln;
        return 0;
    }

    numData = numArgs - 1;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for uniaxial HystereticSM " << iData[0] << endln;
        return 0;
    }

    // Parsing was successful, allocate the material
    Vector theLSdefo(&defoLimitStates[0], (int)defoLimitStates.size());
    Vector theLSforce(&forceLimitStates[0], (int)forceLimitStates.size());

    if (numData > 25) {
        theMaterial = new HystereticSMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
            dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
            dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20], dData[21],
            dData[22], dData[23], dData[24], dData[25], dData[26], dData[27], dData[28], dData[29], theLSforce, theLSdefo, degEnvFactor,
            dData[30], dData[31], dData[32]);
    }
    else if (numData > 17)
        theMaterial = new HystereticSMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
            dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
            dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20], dData[21],
            theLSforce, theLSdefo, degEnvFactor,
            dData[22], dData[23], dData[24]);
    else if (numData > 13)
        theMaterial = new HystereticSMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
            dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13],
            theLSforce, theLSdefo, degEnvFactor,
            dData[14], dData[15], dData[16]);
    else
        theMaterial = new HystereticSMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
            dData[6], dData[7], dData[8], dData[9],
            theLSforce, theLSdefo, degEnvFactor,
            dData[10], dData[11], dData[12]);


    //opserr << "defoLimitStates HystereticSM OPS function" << defoLimitStates[0] << endln;
    //opserr << "defoLimitStates HystereticSM OPS function" << defoLimitStates[1] << endln;
    //opserr << "theLSdefo HystereticSM OPS function" << theLSdefo[0] << endln;
    //opserr << "theLSdefo HystereticSM OPS function" << theLSdefo[1] << endln;


    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type HystereticSM\n";
        return 0;
    }

    return theMaterial;
}


HystereticSMMaterial::HystereticSMMaterial(int tag,
    double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
    double m4p, double r4p, double m5p, double r5p, double m6p, double r6p, double m7p, double r7p,
    double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
    double m4n, double r4n, double m5n, double r5n, double m6n, double r6n, double m7n, double r7n,
    double px, double py, const Vector& forceLimitStatesIN, const Vector& defoLimitStatesIN, double degEnvFactorIN, double d1, double d2, double b) :
    UniaxialMaterial(tag, MAT_TAG_HystereticSM),
    pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
    mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
    mom4p(m4p), rot4p(r4p), mom5p(m5p), rot5p(r5p), mom6p(m6p), rot6p(r6p), mom7p(m7p), rot7p(r7p),
    mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n),
    mom4n(m4n), rot4n(r4n), mom5n(m5n), rot5n(r5n), mom6n(m6n), rot6n(r6n), mom7n(m7n), rot7n(r7n),
    forceLimitStates(forceLimitStatesIN), defoLimitStates(defoLimitStatesIN),
    degEnvFactor(degEnvFactorIN)
{

    nDefoLimitStates = defoLimitStates.Size();
    nForceLimitStates = forceLimitStates.Size();


    bool error = false;
    // Positive backbone parameters

    if (rot1p <= 0.0)
        error = true;

    if (rot2p <= rot1p)
        error = true;

    if (rot3p <= rot2p)
        error = true;


    if (rot4p <= rot3p)
        error = true;


    if (rot5p <= rot4p)
        error = true;


    if (rot6p <= rot5p)
        error = true;

    if (rot7p <= rot6p)
        error = true;


    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- POSITIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }
    error = false;
    // Negative backbone parameters
    if (rot1n >= 0.0)
        error = true;

    if (rot2n >= rot1n)
        error = true;

    if (rot3n >= rot2n)
        error = true;

    if (rot4n >= rot3n)
        error = true;

    if (rot5n >= rot4n)
        error = true;

    if (rot6n >= rot5n)
        error = true;

    if (rot7n >= rot6n)
        error = true;

    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- NEGATIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }

    //// ensure there is no zero-tangent line. it causes convergence issues when used in a section

    if (mom3p >= mom2p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3p will be adjusted to be less than m2p. Current value: " << mom3p << "\n";
        mom3p = 0.99999 * mom2p;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3p was adjusted to be less than m2p. New value: " << mom3p << "\n";
    }

    if (mom4p >= mom3p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4p will be adjusted to be less than mom3p. Current value: " << mom4p << "\n";
        mom4p = 0.99999 * mom3p;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4p was adjusted to be less than m3p. New value: " << mom4p << "\n";
    }


    if (mom5p >= mom4p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5p will be adjusted to be less than mom4p. Current value: " << mom5p << "\n";
        mom5p = 0.99999 * mom4p;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5p was adjusted to be less than m4p. New value: " << mom5p << "\n";
    }


    if (mom6p >= mom5p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6p will be adjusted to be less than mom5p. Current value: " << mom6p << "\n";
        mom6p = 0.99999 * mom5p;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6p was adjusted to be less than m5p. New value: " << mom6p << "\n";
    }

    if (mom7p >= mom6p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7p will be adjusted to be less than mom6p. Current value: " << mom7p << "\n";
        mom7p = 0.99999 * mom6p;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7p was adjusted to be less than m6p. New value: " << mom7p << "\n";
    }


    // negative side
    if (mom3n <= mom2n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3n will be adjusted to be less than m2n. Current value: " << mom3n << "\n";
        mom3n = 0.99999 * mom2n;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3n was adjusted to be less than m2n. New value: " << mom3n << "\n";
    }

    if (mom4n <= mom3n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4n will be adjusted to be less than mom3n. Current value: " << mom4n << "\n";
        mom4n = 0.99999 * mom3n;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4n was adjusted to be less than m3n. New value: " << mom4n << "\n";
    }


    if (mom5n <= mom4n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5n will be adjusted to be less than mom4n. Current value: " << mom5n << "\n";
        mom5n = 0.99999 * mom4n;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5n was adjusted to be less than m4n. New value: " << mom5n << "\n";
    }


    if (mom6n <= mom5n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6n will be adjusted to be less than mom5n. Current value: " << mom6n << "\n";
        mom6n = 0.99999 * mom5n;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6n was adjusted to be less than m5n. New value: " << mom6n << "\n";
    }

    if (mom7n <= mom6n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7n will be adjusted to be less than m6n. Current value: " << mom7n << "\n";
        mom7n = 0.99999 * mom6n;
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7n was adjusted to be less than m6n. New value: " << mom7n << "\n";
    }


    // ensure there is negative slope after m2

    if (mom3p > mom2p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3p needs to be less than m2p\n";
        exit(-1);
    }

    if (mom4p > mom3p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4p needs to be less than m3p\n";
        exit(-1);
    }


    if (mom5p > mom4p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5p needs to be less than m4p\n";
        exit(-1);
    }


    if (mom6p > mom5p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6p needs to be less than m5p\n";
        exit(-1);
    }

    if (mom7p > mom6p) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7p needs to be less than m6p\n";
        exit(-1);
    }

    // negative side

    if (mom3n < mom2n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom3n needs to be greater than m2n\n";
        exit(-1);
    }

    if (mom4n < mom3n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom4n needs to be greater than m3n\n";
        exit(-1);
    }


    if (mom5n < mom4n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom5n needs to be greater than m4n\n";
        exit(-1);
    }


    if (mom6n < mom5n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom6n needs to be greater than m5n\n";
        exit(-1);
    }

    if (mom7n < mom6n) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- mom7n needs to be greater than m6n\n";
        exit(-1);
    }


    E1p = mom1p / rot1p;
    E2p = (mom2p - mom1p) / (rot2p - rot1p);
    E3p = (mom3p - mom2p) / (rot3p - rot2p);
    E4p = (mom4p - mom3p) / (rot4p - rot3p);
    E5p = (mom5p - mom4p) / (rot5p - rot4p);
    E6p = (mom6p - mom5p) / (rot6p - rot5p);
    E7p = (mom7p - mom6p) / (rot7p - rot6p);

    E1n = mom1n / rot1n;
    E2n = (mom2n - mom1n) / (rot2n - rot1n);
    E3n = (mom3n - mom2n) / (rot3n - rot2n);
    E4n = (mom4n - mom3n) / (rot4n - rot3n);
    E5n = (mom5n - mom4n) / (rot5n - rot4n);
    E6n = (mom6n - mom5n) / (rot6n - rot5n);
    E7n = (mom7n - mom6n) / (rot7n - rot6n);

    if (E1p <= 0 || E2p <= 0 || E1n <= 0 || E2n <= 0)
        error = true;
    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- The first two segments in positive and negative direction must have positive slope\n";
        exit(-1);
    }
    //  if (E3p <= 0 || E4p <= 0 || E5p <= 0 || E6p <= 0 || E7p <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }
    //  if (E3n <= 0 || E4n <= 0 || E5n <= 0 || E6n <= 0 || E7n <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }




    energyA = 0.5 * (rot1p * mom1p + (rot2p - rot1p) * (mom2p + mom1p) + (rot3p - rot2p) * (mom3p + mom2p) +
        (rot4p - rot3p) * (mom4p + mom3p) + (rot5p - rot4p) * (mom5p + mom4p) + (rot6p - rot5p) * (mom6p + mom5p) + (rot7p - rot6p) * (mom7p + mom6p) +
        rot1n * mom1n + (rot2n - rot1n) * (mom2n + mom1n) + (rot3n - rot2n) * (mom3n + mom2n)) +
        +(rot4n - rot3n) * (mom4n + mom3n) + (rot5n - rot4n) * (mom5n + mom4n) + (rot6n - rot5n) * (mom6n + mom5n) + (rot7n - rot6n) * (mom7n + mom6n);

    //opserr << "energyA  HystereticSM " << energyA << endln;

    // Set envelope slopes
    this->setEnvelope();

    // Initialize history variables
    this->revertToStart();
    this->revertToLastCommit();

}




HystereticSMMaterial::HystereticSMMaterial(int tag,
    double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
    double m4p, double r4p, double m5p, double r5p,
    double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
    double m4n, double r4n, double m5n, double r5n,
    double px, double py, const Vector& forceLimitStatesIN, const Vector& defoLimitStatesIN, double degEnvFactorIN, double d1, double d2, double b) :
    UniaxialMaterial(tag, MAT_TAG_HystereticSM),
    pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
    mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
    mom5p(m4p), rot5p(r4p), mom7p(m5p), rot7p(r5p),
    mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n),
    mom5n(m4n), rot5n(r4n), mom7n(m5n), rot7n(r5n),
    forceLimitStates(forceLimitStatesIN), defoLimitStates(defoLimitStatesIN),
    degEnvFactor(degEnvFactorIN)

{
    bool error = false;
    // Positive backbone parameters

    if (rot1p <= 0.0)
        error = true;

    if (rot2p <= rot1p)
        error = true;

    if (rot3p <= rot2p)
        error = true;


    if (rot5p <= rot3p)
        error = true;


    if (rot7p <= rot5p)
        error = true;




    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- POSITIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }
    error = false;
    // Negative backbone parameters
    if (rot1n >= 0.0)
        error = true;

    if (rot2n >= rot1n)
        error = true;

    if (rot3n >= rot2n)
        error = true;

    if (rot5n >= rot3n)
        error = true;

    if (rot7n >= rot5n)
        error = true;


    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- NEGATIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }


    mom4p = 0.5 * (mom3p + mom5p);
    mom4n = 0.5 * (mom3n + mom5n);
    rot4p = 0.5 * (rot3p + rot5p);
    rot4n = 0.5 * (rot3n + rot5n);


    mom6p = 0.5 * (mom5p + mom7p);
    mom6n = 0.5 * (mom5n + mom7n);
    rot6p = 0.5 * (rot5p + rot7p);
    rot6n = 0.5 * (rot5n + rot7n);



    E1p = mom1p / rot1p;
    E2p = (mom2p - mom1p) / (rot2p - rot1p);
    E3p = (mom3p - mom2p) / (rot3p - rot2p);
    E4p = (mom4p - mom3p) / (rot4p - rot3p);
    E5p = (mom5p - mom4p) / (rot5p - rot4p);


    E1n = mom1n / rot1n;
    E2n = (mom2n - mom1n) / (rot2n - rot1n);
    E3n = (mom3n - mom2n) / (rot3n - rot2n);
    E4n = (mom4n - mom3n) / (rot4n - rot3n);
    E5n = (mom5n - mom4n) / (rot5n - rot4n);


    if (E1p <= 0 || E2p <= 0 || E1n <= 0 || E2n <= 0)
        error = true;
    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- The first two segments in positive and negative direction must have positive slope\n";
        exit(-1);
    }
    //  if (E3p <= 0 || E4p <= 0 || E5p <= 0 || E6p <= 0 || E7p <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }
    //  if (E3n <= 0 || E4n <= 0 || E5n <= 0 || E6n <= 0 || E7n <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }




    energyA = 0.5 * (rot1p * mom1p + (rot2p - rot1p) * (mom2p + mom1p) + (rot3p - rot2p) * (mom3p + mom2p) +
        (rot4p - rot3p) * (mom4p + mom3p) + (rot5p - rot4p) * (mom5p + mom4p) + (rot6p - rot5p) * (mom6p + mom5p) + (rot7p - rot6p) * (mom7p + mom6p) +
        rot1n * mom1n + (rot2n - rot1n) * (mom2n + mom1n) + (rot3n - rot2n) * (mom3n + mom2n)) +
        +(rot4n - rot3n) * (mom4n + mom3n) + (rot5n - rot4n) * (mom5n + mom4n) + (rot6n - rot5n) * (mom6n + mom5n) + (rot7n - rot6n) * (mom7n + mom6n);

    //opserr << "energyA  HystereticSM " << energyA << endln;



    // Set envelope slopes
    this->setEnvelope();

    // Initialize history variables
    this->revertToStart();
    this->revertToLastCommit();

}





HystereticSMMaterial::HystereticSMMaterial(int tag,
    double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
    double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
    double px, double py, const Vector& forceLimitStatesIN, const Vector& defoLimitStatesIN, double degEnvFactorIN, double d1, double d2, double b) :
    UniaxialMaterial(tag, MAT_TAG_HystereticSM),
    pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
    mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom7p(m3p), rot7p(r3p),
    mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom7n(m3n), rot7n(r3n),
    forceLimitStates(forceLimitStatesIN), defoLimitStates(defoLimitStatesIN),
    degEnvFactor(degEnvFactorIN)

{
    bool error = false;
    // Positive backbone parameters

    if (rot1p <= 0.0)
        error = true;

    if (rot2p <= rot1p)
        error = true;


    if (rot7p <= rot2p)
        error = true;




    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- POSITIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }
    error = false;
    // Negative backbone parameters
    if (rot1n >= 0.0)
        error = true;

    if (rot2n >= rot1n)
        error = true;

    if (rot7n >= rot2n)
        error = true;


    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- NEGATIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }


    mom3p = mom2p + 0.15 * (mom7p - mom2p);
    mom3n = mom2n + 0.15 * (mom7n - mom2n);
    rot3p = rot2p + 0.15 * (rot7p - rot2p);
    rot3n = rot2n + 0.15 * (rot7n - rot2n);

    mom4p = mom2p + 0.35 * (mom7p - mom2p);
    mom4n = mom2n + 0.35 * (mom7n - mom2n);
    rot4p = rot2p + 0.35 * (rot7p - rot2p);
    rot4n = rot2n + 0.35 * (rot7n - rot2n);

    mom5p = mom2p + 0.65 * (mom7p - mom2p);
    mom5n = mom2n + 0.65 * (mom7n - mom2n);
    rot5p = rot2p + 0.65 * (rot7p - rot2p);
    rot5n = rot2n + 0.65 * (rot7n - rot2n);

    mom6p = mom2p + 0.85 * (mom7p - mom2p);
    mom6n = mom2n + 0.85 * (mom7n - mom2n);
    rot6p = rot2p + 0.85 * (rot7p - rot2p);
    rot6n = rot2n + 0.85 * (rot7n - rot2n);



    E1p = mom1p / rot1p;
    E2p = (mom2p - mom1p) / (rot2p - rot1p);
    E3p = (mom3p - mom2p) / (rot3p - rot2p);
    E4p = (mom4p - mom3p) / (rot4p - rot3p);
    E5p = (mom5p - mom4p) / (rot5p - rot4p);


    E1n = mom1n / rot1n;
    E2n = (mom2n - mom1n) / (rot2n - rot1n);
    E3n = (mom3n - mom2n) / (rot3n - rot2n);
    E4n = (mom4n - mom3n) / (rot4n - rot3n);
    E5n = (mom5n - mom4n) / (rot5n - rot4n);


    if (E1p <= 0 || E2p <= 0 || E1n <= 0 || E2n <= 0)
        error = true;
    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- The first two segments in positive and negative direction must have positive slope\n";
        exit(-1);
    }
    //  if (E3p <= 0 || E4p <= 0 || E5p <= 0 || E6p <= 0 || E7p <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }
    //  if (E3n <= 0 || E4n <= 0 || E5n <= 0 || E6n <= 0 || E7n <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }




    energyA = 0.5 * (rot1p * mom1p + (rot2p - rot1p) * (mom2p + mom1p) + (rot3p - rot2p) * (mom3p + mom2p) +
        (rot4p - rot3p) * (mom4p + mom3p) + (rot5p - rot4p) * (mom5p + mom4p) + (rot6p - rot5p) * (mom6p + mom5p) + (rot7p - rot6p) * (mom7p + mom6p) +
        rot1n * mom1n + (rot2n - rot1n) * (mom2n + mom1n) + (rot3n - rot2n) * (mom3n + mom2n)) +
        +(rot4n - rot3n) * (mom4n + mom3n) + (rot5n - rot4n) * (mom5n + mom4n) + (rot6n - rot5n) * (mom6n + mom5n) + (rot7n - rot6n) * (mom7n + mom6n);



    //opserr << "energyA  HystereticSM " << energyA << endln;

    // Set envelope slopes
    this->setEnvelope();

    // Initialize history variables
    this->revertToStart();
    this->revertToLastCommit();

}




HystereticSMMaterial::HystereticSMMaterial(int tag,
    double m1p, double r1p, double m2p, double r2p,
    double m1n, double r1n, double m2n, double r2n,
    double px, double py, const Vector& forceLimitStatesIN, const Vector& defoLimitStatesIN, double degEnvFactorIN, double d1, double d2, double b) :
    UniaxialMaterial(tag, MAT_TAG_HystereticSM),
    pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
    mom1p(m1p), rot1p(r1p), mom7p(m2p), rot7p(r2p),
    mom1n(m1n), rot1n(r1n), mom7n(m2n), rot7n(r2n),
    forceLimitStates(forceLimitStatesIN), defoLimitStates(defoLimitStatesIN),
    degEnvFactor(degEnvFactorIN)

{
    bool error = false;
    // Positive backbone parameters

    if (rot1p <= 0.0)
        error = true;


    if (rot7p <= rot1p)
        error = true;




    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- POSITIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }
    error = false;
    // Negative backbone parameters
    if (rot1n >= 0.0)
        error = true;


    if (rot7n >= rot1n)
        error = true;


    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- NEGATIVE input backbone is not unique (one-to-one)\n";
        exit(-1);
    }


    mom2p = mom1p + 0.15 * (mom7p - mom1p);
    mom2n = mom1n + 0.15 * (mom7n - mom1n);
    rot2p = rot1p + 0.15 * (rot7p - rot1p);
    rot2n = rot1n + 0.15 * (rot7n - rot1n);

    mom3p = mom1p + 0.35 * (mom7p - mom1p);
    mom3n = mom1n + 0.35 * (mom7n - mom1n);
    rot3p = rot1p + 0.35 * (rot7p - rot1p);
    rot3n = rot1n + 0.35 * (rot7n - rot1n);

    mom4p = mom1p + 0.5 * (mom7p - mom1p);
    mom4n = mom1n + 0.5 * (mom7n - mom1n);
    rot4p = rot1p + 0.5 * (rot7p - rot1p);
    rot4n = rot1n + 0.5 * (rot7n - rot1n);

    mom5p = mom1p + 0.65 * (mom7p - mom1p);
    mom5n = mom1n + 0.65 * (mom7n - mom1n);
    rot5p = rot1p + 0.65 * (rot7p - rot1p);
    rot5n = rot1n + 0.65 * (rot7n - rot1n);

    mom6p = mom1p + 0.85 * (mom7p - mom1p);
    mom6n = mom1n + 0.85 * (mom7n - mom1n);
    rot6p = rot1p + 0.85 * (rot7p - rot1p);
    rot6n = rot1n + 0.85 * (rot7n - rot1n);



    E1p = mom1p / rot1p;
    E2p = (mom2p - mom1p) / (rot2p - rot1p);
    E3p = (mom3p - mom2p) / (rot3p - rot2p);
    E4p = (mom4p - mom3p) / (rot4p - rot3p);
    E5p = (mom5p - mom4p) / (rot5p - rot4p);


    E1n = mom1n / rot1n;
    E2n = (mom2n - mom1n) / (rot2n - rot1n);
    E3n = (mom3n - mom2n) / (rot3n - rot2n);
    E4n = (mom4n - mom3n) / (rot4n - rot3n);
    E5n = (mom5n - mom4n) / (rot5n - rot4n);


    if (E1p <= 0 || E2p <= 0 || E1n <= 0 || E2n <= 0)
        error = true;
    if (error) {
        opserr << "HystereticSMMaterial::HystereticSMMaterial -- The first two segments in positive and negative direction must have positive slope\n";
        exit(-1);
    }
    //  if (E3p <= 0 || E4p <= 0 || E5p <= 0 || E6p <= 0 || E7p <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }
    //  if (E3n <= 0 || E4n <= 0 || E5n <= 0 || E6n <= 0 || E7n <= 0)
    //      error = true;
    //  if (error) {
    //      opserr << "HystereticSMMaterial::HystereticSMMaterial -- Segments 2-7 must have negative slope\n";
    //      exit(-1);
    //  }




    energyA = 0.5 * (rot1p * mom1p + (rot2p - rot1p) * (mom2p + mom1p) + (rot3p - rot2p) * (mom3p + mom2p) +
        (rot4p - rot3p) * (mom4p + mom3p) + (rot5p - rot4p) * (mom5p + mom4p) + (rot6p - rot5p) * (mom6p + mom5p) + (rot7p - rot6p) * (mom7p + mom6p) +
        rot1n * mom1n + (rot2n - rot1n) * (mom2n + mom1n) + (rot3n - rot2n) * (mom3n + mom2n)) +
        +(rot4n - rot3n) * (mom4n + mom3n) + (rot5n - rot4n) * (mom5n + mom4n) + (rot6n - rot5n) * (mom6n + mom5n) + (rot7n - rot6n) * (mom7n + mom6n);

    //opserr << "energyA  HystereticSM " << energyA << endln;





    // Set envelope slopes
    this->setEnvelope();

    // Initialize history variables
    this->revertToStart();
    this->revertToLastCommit();

}







// WHAT IS THIS??
HystereticSMMaterial::HystereticSMMaterial() :
    UniaxialMaterial(0, MAT_TAG_HystereticSM),
    pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
    mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
    mom4p(0.0), rot4p(0.0), mom5p(0.0), rot5p(0.0), mom6p(0.0), rot6p(0.0), mom7p(0.0), rot7p(0.0),
    mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0),
    mom4n(0.0), rot4n(0.0), mom5n(0.0), rot5n(0.0), mom6n(0.0), rot6n(0.0), mom7n(0.0), rot7n(0.0)

{

}

HystereticSMMaterial::~HystereticSMMaterial()
{
    // Nothing to do
}

int
HystereticSMMaterial::setTrialStrain(double strain, double strainRate)
{
    if (TloadIndicator == 0 && strain == 0.0)
        return 0;

    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TenergyD = CenergyD;
    TrotPu = CrotPu;
    TrotNu = CrotNu;

    Tstrain = strain;
    double dStrain = Tstrain - Cstrain;

    if (fabs(dStrain) < DBL_EPSILON)
        return 0;

    TloadIndicator = CloadIndicator;

    if (TloadIndicator == 0)
        TloadIndicator = (dStrain < 0.0) ? 2 : 1;

    if (Tstrain >= CrotMax) {
        TrotMax = Tstrain;
        Ttangent = posEnvlpTangent(Tstrain);
        Tstress = posEnvlpStress(Tstrain);
        TloadIndicator = 1;
    }
    else if (Tstrain <= CrotMin) {
        TrotMin = Tstrain;
        Ttangent = negEnvlpTangent(Tstrain);
        Tstress = negEnvlpStress(Tstrain);
        TloadIndicator = 2;
    }
    else {
        if (dStrain < 0.0)
            negativeIncrement(dStrain);
        else if (dStrain > 0.0)
            positiveIncrement(dStrain);
    }

    TenergyD = CenergyD + 0.5 * (Cstress + Tstress) * dStrain;


    //  if (this->getTag() == 40)
    //    opserr << "setTrial: " << Tstrain << " " << Ttangent << " " << Tstress << endln;

    return 0;
}


double
HystereticSMMaterial::getStrain(void)
{
    return Tstrain;
}

double
HystereticSMMaterial::getMUy(void)
{
    double MU = 0;

    if (Tstrain > 0)
        MU = Tstrain / rot1p;
    else
        MU = Tstrain / rot1n;
    return MU;
}

double
HystereticSMMaterial::getThetaP(void)
{
    double ThetaP = 0;

    if (Tstrain > rot1p)
        ThetaP = Tstrain - ThetaP;
    else if (Tstrain < rot1n)
        ThetaP = Tstrain - rot1n;
    return ThetaP;
}



double
HystereticSMMaterial::getStress(void)
{
    return Tstress;
}

double
HystereticSMMaterial::getTangent(void)
{
    return Ttangent;
}

void
HystereticSMMaterial::positiveIncrement(double dStrain)
{
    double kn = pow(CrotMin / rot1n, beta);
    kn = (kn < 1.0) ? 1.0 : 1.0 / kn;
    double kp = pow(CrotMax / rot1p, beta);
    kp = (kp < 1.0) ? 1.0 : 1.0 / kp;
    double damfc = 0.0; // 220905sm moved this here
    if (TloadIndicator == 2) {
        TloadIndicator = 1;
        if (Cstress <= 0.0) {
            TrotNu = Cstrain - Cstress / (Eun * kn);
            double energy = CenergyD - 0.5 * Cstress / (Eun * kn) * Cstress;
            // double damfc = 0.0;
            if (CrotMin < rot1n) {
                damfc = damfc2 * energy / energyA;
                damfc += damfc1 * (CrotMin - rot1n) / rot1n;
            }

            TrotMax = CrotMax * (1.0 + damfc);
        }
    }

    TloadIndicator = 1;

    if (TrotMax > POS_INF_STRAIN)
        TrotMax = POS_INF_STRAIN;

    TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;
    double maxmom = posEnvlpStress(TrotMax);

    //double maxmom = posEnvlpStress(TrotMax) * (1.0 - damfc * degEnvFactor); // 220905sm added damfc
    //if (damfc != 0) {
    //    opserr << "damfc positiveIncrement HystereticSM " << damfc << endln;
    //    opserr << "(1.0 - damfc* degEnvFactor) positiveIncrement HystereticSM " << (1.0 - damfc * degEnvFactor) << endln;
    //    opserr << "maxmom positiveIncrement HystereticSM " << maxmom << endln;
    //}

    double rotlim = negEnvlpRotlim(CrotMin);
    double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

    // rotrel = TrotNu;
    // if (negEnvlpStress(CrotMin) >= 0.0)
    //    rotrel = rotlim;

    //	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);

    double rotmp2 = TrotMax - (1.0 - pinchY) * maxmom / (Eup * kp);
    //double rotmp2 = TrotMax-(1-pinchY)*maxmom/Eup;
    //	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
    double rotch = rotrel + (rotmp2 - rotrel) * pinchX;                   // changed on 7/11/2006

    double tmpmo1;
    double tmpmo2;

    if (Tstrain < TrotNu) {
        Ttangent = Eun * kn;
        Tstress = Cstress + Ttangent * dStrain;
        if (Tstress >= 0.0) {
            Tstress = 0.0;
            Ttangent = Eun * 1.0e-9;
        }
    }

    else if (Tstrain >= TrotNu && Tstrain < rotch) {
        if (Tstrain <= rotrel) {
            Tstress = 0.0;
            Ttangent = Eup * 1.0e-9;
        }
        else {
            Ttangent = maxmom * pinchY / (rotch - rotrel);
            tmpmo1 = Cstress + Eup * kp * dStrain;
            tmpmo2 = (Tstrain - rotrel) * Ttangent;
            if (tmpmo1 < tmpmo2) {
                Tstress = tmpmo1;
                Ttangent = Eup * kp;
            }
            else
                Tstress = tmpmo2;
        }
    }

    else {
        Ttangent = (1.0 - pinchY) * maxmom / (TrotMax - rotch);
        tmpmo1 = Cstress + Eup * kp * dStrain;
        tmpmo2 = pinchY * maxmom + (Tstrain - rotch) * Ttangent;
        if (tmpmo1 < tmpmo2) {
            Tstress = tmpmo1;
            Ttangent = Eup * kp;
        }
        else
            Tstress = tmpmo2;
    }

    ////// 220905sm added damfc
    //double TstressRedux = (1.0 - damfc * degEnvFactor);
    //if (TstressRedux < 0) {
    //    TstressRedux = 0;
    //}
    //Tstress = Tstress * TstressRedux;
}

void
HystereticSMMaterial::negativeIncrement(double dStrain)
{
    double kn = pow(CrotMin / rot1n, beta);
    kn = (kn < 1.0) ? 1.0 : 1.0 / kn;
    double kp = pow(CrotMax / rot1p, beta);
    kp = (kp < 1.0) ? 1.0 : 1.0 / kp;

    double damfc = 0.0; // 220905sm moved this here
    if (TloadIndicator == 1) {
        TloadIndicator = 2;
        if (Cstress >= 0.0) {
            TrotPu = Cstrain - Cstress / (Eup * kp);
            double energy = CenergyD - 0.5 * Cstress / (Eup * kp) * Cstress;
            //double damfc = 0.0;
            if (CrotMax > rot1p) {
                damfc = damfc2 * energy / energyA;
                damfc += damfc1 * (CrotMax - rot1p) / rot1p;
            }

            TrotMin = CrotMin * (1.0 + damfc);
        }
    }

    TloadIndicator = 2;

    if (TrotMin < NEG_INF_STRAIN)
        TrotMin = NEG_INF_STRAIN;

    TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;


    double minmom = negEnvlpStress(TrotMin);

    //double minmom = negEnvlpStress(TrotMin) * (1.0 - damfc * degEnvFactor); // 220905sm added damfc
    //if (damfc != 0) {
    //    opserr << "damfc negativeIncrement HystereticSM " << damfc << endln;
    //    opserr << "(1.0 - damfc* degEnvFactor) negativeIncrement HystereticSM " << (1.0 - damfc * degEnvFactor) << endln;
    //    opserr << "minmom negativeIncrement HystereticSM " << minmom << endln;
    //}

    double rotlim = posEnvlpRotlim(CrotMax);
    double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

    //rotrel = TrotPu;
    //if (posEnvlpStress(CrotMax) <= 0.0)
    //  rotrel = rotlim;

    //double rotmp1 = rotrel + pinchY*(TrotMin-rotrel);
    double rotmp2 = TrotMin - (1.0 - pinchY) * minmom / (Eun * kn);
    //double rotmp2 = TrotMin-(1-pinchY)*minmom/Eun;	
    //double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
    double rotch = rotrel + (rotmp2 - rotrel) * pinchX;                   // changed on 7/11/2006

    double tmpmo1;
    double tmpmo2;

    if (Tstrain > TrotPu) {
        Ttangent = Eup * kp;
        Tstress = Cstress + Ttangent * dStrain;
        if (Tstress <= 0.0) {
            Tstress = 0.0;
            Ttangent = Eup * 1.0e-9;
        }
    }

    else if (Tstrain <= TrotPu && Tstrain > rotch) {
        if (Tstrain >= rotrel) {
            Tstress = 0.0;
            Ttangent = Eun * 1.0e-9;
        }
        else {
            Ttangent = minmom * pinchY / (rotch - rotrel);
            tmpmo1 = Cstress + Eun * kn * dStrain;
            tmpmo2 = (Tstrain - rotrel) * Ttangent;
            if (tmpmo1 > tmpmo2) {
                Tstress = tmpmo1;
                Ttangent = Eun * kn;
            }
            else
                Tstress = tmpmo2;
        }
    }

    else {
        Ttangent = (1.0 - pinchY) * minmom / (TrotMin - rotch);
        tmpmo1 = Cstress + Eun * kn * dStrain;
        tmpmo2 = pinchY * minmom + (Tstrain - rotch) * Ttangent;
        if (tmpmo1 > tmpmo2) {
            Tstress = tmpmo1;
            Ttangent = Eun * kn;
        }
        else
            Tstress = tmpmo2;
    }

    //    //// 220905sm added damfc
    //    double TstressRedux = (1.0 - damfc * degEnvFactor);
    //    if (TstressRedux < 0) {
    //        TstressRedux = 0;
    //    }
    //    Tstress = Tstress * TstressRedux;

}

int
HystereticSMMaterial::commitState(void)
{
    CrotMax = TrotMax;
    CrotMin = TrotMin;
    CrotPu = TrotPu;
    CrotNu = TrotNu;
    CenergyD = TenergyD;
    CloadIndicator = TloadIndicator;

    Cstress = Tstress;
    Cstrain = Tstrain;
    return 0;
}

int
HystereticSMMaterial::revertToLastCommit(void)
{
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;

    Tstress = Cstress;
    Tstrain = Cstrain;

    return 0;
}

int
HystereticSMMaterial::revertToStart(void)
{
    CrotMax = 0.0;
    CrotMin = 0.0;
    CrotPu = 0.0;
    CrotNu = 0.0;
    CenergyD = 0.0;
    CloadIndicator = 0;

    Cstress = 0.0;
    Cstrain = 0.0;

    Tstrain = 0;
    Tstress = 0;
    Ttangent = E1p;

    return 0;
}

UniaxialMaterial*
HystereticSMMaterial::getCopy(void)
{
    HystereticSMMaterial* theCopy = new HystereticSMMaterial(this->getTag(),
        mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
        mom4p, rot4p, mom5p, rot5p, mom6p, rot6p, mom7p, rot7p,
        mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
        mom4n, rot4n, mom5n, rot5n, mom6n, rot6n, mom7n, rot7n,

        pinchX, pinchY, forceLimitStates, defoLimitStates, degEnvFactor, damfc1, damfc2, beta);

    theCopy->CrotMax = CrotMax;
    theCopy->CrotMin = CrotMin;
    theCopy->CrotPu = CrotPu;
    theCopy->CrotNu = CrotNu;
    theCopy->CenergyD = CenergyD;
    theCopy->CloadIndicator = CloadIndicator;
    theCopy->Cstress = Cstress;
    theCopy->Cstrain = Cstrain;
    theCopy->Ttangent = Ttangent;

    return theCopy;
}



int
HystereticSMMaterial::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    static Vector data(43);

    data(0) = this->getTag();
    data(1) = mom1p;
    data(2) = rot1p;
    data(3) = mom2p;
    data(4) = rot2p;
    data(5) = mom3p;
    data(6) = rot3p;
    data(7) = mom4p;
    data(8) = rot4p;
    data(9) = mom5p;
    data(10) = rot5p;
    data(11) = mom6p;
    data(12) = rot6p;
    data(13) = mom7p;
    data(14) = rot7p;
    data(15) = mom1n;
    data(16) = rot1n;
    data(17) = mom2n;
    data(18) = rot2n;
    data(19) = mom3n;
    data(20) = rot3n;
    data(21) = mom4n;
    data(22) = rot4n;
    data(23) = mom5n;
    data(24) = rot5n;
    data(25) = mom6n;
    data(26) = rot6n;
    data(27) = mom7n;
    data(28) = rot7n;
    data(29) = pinchX;
    data(30) = pinchY;
    data(31) = damfc1;
    data(32) = damfc2;
    data(33) = beta;
    data(34) = CrotMax;
    data(35) = CrotMin;
    data(36) = CrotPu;
    data(37) = CrotNu;
    data(38) = CenergyD;
    data(39) = CloadIndicator;
    data(40) = Cstress;
    data(41) = Cstrain;
    data(42) = Ttangent;


    res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0)
        opserr << "HystereticSMMaterial::sendSelf() - failed to send data\n";


    return res;
}

int
HystereticSMMaterial::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int res = 0;

    static Vector data(43);
    res = theChannel.recvVector(this->getDbTag(), commitTag, data);

    if (res < 0) {
        opserr << "HystereticSMMaterial::recvSelf() - failed to receive data\n";
        return res;
    }
    else {
        this->setTag((int)data(0));
        mom1p = data(1);
        rot1p = data(2);
        mom2p = data(3);
        rot2p = data(4);
        mom3p = data(5);
        rot3p = data(6);
        mom4p = data(7);
        rot4p = data(8);
        mom5p = data(9);
        rot5p = data(10);
        mom6p = data(11);
        rot6p = data(12);
        mom7p = data(13);
        rot7p = data(14);
        mom1n = data(15);
        rot1n = data(16);
        mom2n = data(17);
        rot2n = data(18);
        mom3n = data(19);
        rot3n = data(20);
        mom4n = data(21);
        rot4n = data(22);
        mom5n = data(23);
        rot5n = data(24);
        mom6n = data(25);
        rot6n = data(26);
        mom7n = data(27);
        rot7n = data(28);
        pinchX = data(29);
        pinchY = data(30);
        damfc1 = data(31);
        damfc2 = data(32);
        beta = data(33);
        CrotMax = data(34);
        CrotMin = data(35);
        CrotPu = data(36);
        CrotNu = data(37);
        CenergyD = data(38);
        CloadIndicator = (int)data(39);
        Cstress = data(40);
        Cstrain = data(41);
        Ttangent = data(42);


        // set the trial values
        TrotMax = CrotMax;
        TrotMin = CrotMin;
        TrotPu = CrotPu;
        TrotNu = CrotNu;
        TenergyD = CenergyD;
        TloadIndicator = CloadIndicator;
        Tstress = Cstress;
        Tstrain = Cstrain;
    }

    // Set envelope slopes
    this->setEnvelope();

    return 0;
}

void
HystereticSMMaterial::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "HHystereticSMMaterial, tag: " << this->getTag() << endln;
        s << "s1p: " << mom1p << endln;
        s << "e1p: " << rot1p << endln;
        s << "E1p: " << E1p << endln;
        s << "s2p: " << mom2p << endln;
        s << "e2p: " << rot2p << endln;
        s << "E2p: " << E2p << endln;
        s << "s3p: " << mom3p << endln;
        s << "e3p: " << rot3p << endln;
        s << "E3p: " << E3p << endln;

        s << "s4p: " << mom4p << endln;
        s << "e4p: " << rot4p << endln;
        s << "E4p: " << E4p << endln;
        s << "s5p: " << mom5p << endln;
        s << "e5p: " << rot5p << endln;
        s << "E5p: " << E5p << endln;
        s << "s6p: " << mom6p << endln;
        s << "e6p: " << rot6p << endln;
        s << "E6p: " << E6p << endln;
        s << "s7p: " << mom7p << endln;
        s << "e7p: " << rot7p << endln;
        s << "E7p: " << E7p << endln;

        s << "s1n: " << mom1n << endln;
        s << "e1n: " << rot1n << endln;
        s << "E1n: " << E1n << endln;
        s << "s2n: " << mom2n << endln;
        s << "e2n: " << rot2n << endln;
        s << "E2n: " << E2n << endln;
        s << "s3n: " << mom3n << endln;
        s << "e3n: " << rot3n << endln;
        s << "E3n: " << E3n << endln;

        s << "s4n: " << mom4n << endln;
        s << "e4n: " << rot4n << endln;
        s << "E4n: " << E4n << endln;
        s << "s5n: " << mom5n << endln;
        s << "e5n: " << rot5n << endln;
        s << "E5n: " << E5n << endln;
        s << "s6n: " << mom6n << endln;
        s << "e6n: " << rot6n << endln;
        s << "E6n: " << E6n << endln;
        s << "s7n: " << mom7n << endln;
        s << "e7n: " << rot7n << endln;
        s << "E7n: " << E7n << endln;

        s << "pinchX: " << pinchX << endln;
        s << "pinchY: " << pinchY << endln;
        s << "damfc1: " << damfc1 << endln;
        s << "damfc2: " << damfc2 << endln;
        s << "energyA: " << energyA << endln;
        s << "beta: " << beta << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"HystereticSMMaterial\", ";
        s << "\"s1p\": " << mom1p << ", ";
        s << "\"e1p\": " << rot1p << ", ";
        s << "\"E1p\": " << E1p << ", ";
        s << "\"s2p\": " << mom2p << ", ";
        s << "\"e2p\": " << rot2p << ", ";
        s << "\"E2p\": " << E2p << ", ";
        s << "\"s3p\": " << mom3p << ", ";
        s << "\"e3p\": " << rot3p << ", ";
        s << "\"E3p\": " << E3p << ", ";

        s << "\"s4p\": " << mom4p << ", ";
        s << "\"e4p\": " << rot4p << ", ";
        s << "\"E4p\": " << E4p << ", ";
        s << "\"s5p\": " << mom5p << ", ";
        s << "\"e5p\": " << rot5p << ", ";
        s << "\"E5p\": " << E5p << ", ";
        s << "\"s6p\": " << mom6p << ", ";
        s << "\"e6p\": " << rot6p << ", ";
        s << "\"E6p\": " << E6p << ", ";
        s << "\"s7p\": " << mom7p << ", ";
        s << "\"e7p\": " << rot7p << ", ";
        s << "\"E7p\": " << E7p << ", ";

        s << "\"s1n\": " << mom1n << ", ";
        s << "\"e1n\": " << rot1n << ", ";
        s << "\"E1n\": " << E1n << ", ";
        s << "\"s2n\": " << mom2n << ", ";
        s << "\"e2n\": " << rot2n << ", ";
        s << "\"E2n\": " << E2n << ", ";
        s << "\"s3n\": " << mom3n << ", ";
        s << "\"e3n\": " << rot3n << ", ";
        s << "\"E3n\": " << E3n << ", ";

        s << "\"s4n\": " << mom4n << ", ";
        s << "\"e4n\": " << rot4n << ", ";
        s << "\"E4n\": " << E4n << ", ";
        s << "\"s5n\": " << mom5n << ", ";
        s << "\"e5n\": " << rot5n << ", ";
        s << "\"E5n\": " << E5n << ", ";
        s << "\"s6n\": " << mom6n << ", ";
        s << "\"e6n\": " << rot6n << ", ";
        s << "\"E6n\": " << E6n << ", ";
        s << "\"s7n\": " << mom7n << ", ";
        s << "\"e7n\": " << rot7n << ", ";
        s << "\"E7n\": " << E7n << ", ";

        s << "\"pinchX\": " << pinchX << ", ";
        s << "\"pinchY\": " << pinchY << ", ";
        s << "\"damfc1\": " << damfc1 << ", ";
        s << "\"damfc2\": " << damfc2 << ", ";
        s << "\"energyA\": " << energyA << ", ";
        s << "\"beta\": " << beta << "}";
    }
}

int
HystereticSMMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    if (strcmp(argv[0], "mom1p") == 0 || strcmp(argv[0], "fy") == 0 || strcmp(argv[0], "Fy") == 0) {
        param.setValue(mom1p);
        return param.addObject(1, this);
    }
    if (strcmp(argv[0], "rot1p") == 0) {
        param.setValue(rot1p);
        return param.addObject(2, this);
    }
    if (strcmp(argv[0], "mom2p") == 0) {
        param.setValue(mom2p);
        return param.addObject(3, this);
    }
    if (strcmp(argv[0], "rot2p") == 0) {
        param.setValue(rot2p);
        return param.addObject(4, this);
    }
    if (strcmp(argv[0], "mom3p") == 0) {
        param.setValue(mom3p);
        return param.addObject(5, this);
    }
    if (strcmp(argv[0], "rot3p") == 0) {
        param.setValue(rot3p);
        return param.addObject(6, this);
    }


    if (strcmp(argv[0], "mom4p") == 0) {
        param.setValue(mom4p);
        return param.addObject(7, this);
    }
    if (strcmp(argv[0], "rot4p") == 0) {
        param.setValue(rot4p);
        return param.addObject(8, this);
    }
    if (strcmp(argv[0], "mom5p") == 0) {
        param.setValue(mom5p);
        return param.addObject(9, this);
    }
    if (strcmp(argv[0], "rot5p") == 0) {
        param.setValue(rot5p);
        return param.addObject(10, this);
    }
    if (strcmp(argv[0], "mom6p") == 0) {
        param.setValue(mom6p);
        return param.addObject(11, this);
    }
    if (strcmp(argv[0], "rot6p") == 0) {
        param.setValue(rot6p);
        return param.addObject(12, this);
    }
    if (strcmp(argv[0], "mom7p") == 0) {
        param.setValue(mom7p);
        return param.addObject(13, this);
    }
    if (strcmp(argv[0], "rot7p") == 0) {
        param.setValue(rot7p);
        return param.addObject(14, this);
    }


    if (strcmp(argv[0], "mom1n") == 0) {
        param.setValue(mom1n);
        return param.addObject(15, this);
    }
    if (strcmp(argv[0], "rot1n") == 0) {
        param.setValue(rot1n);
        return param.addObject(16, this);
    }
    if (strcmp(argv[0], "mom2n") == 0) {
        param.setValue(mom2n);
        return param.addObject(17, this);
    }
    if (strcmp(argv[0], "rot2n") == 0) {
        param.setValue(rot2n);
        return param.addObject(18, this);
    }
    if (strcmp(argv[0], "mom3n") == 0) {
        param.setValue(mom3n);
        return param.addObject(19, this);
    }
    if (strcmp(argv[0], "rot3n") == 0) {
        param.setValue(rot3n);
        return param.addObject(20, this);
    }


    if (strcmp(argv[0], "mom4n") == 0) {
        param.setValue(mom4n);
        return param.addObject(21, this);
    }
    if (strcmp(argv[0], "rot4n") == 0) {
        param.setValue(rot4n);
        return param.addObject(22, this);
    }
    if (strcmp(argv[0], "mom5n") == 0) {
        param.setValue(mom5n);
        return param.addObject(23, this);
    }
    if (strcmp(argv[0], "rot5n") == 0) {
        param.setValue(rot5n);
        return param.addObject(24, this);
    }
    if (strcmp(argv[0], "mom6n") == 0) {
        param.setValue(mom6n);
        return param.addObject(25, this);
    }
    if (strcmp(argv[0], "rot6n") == 0) {
        param.setValue(rot6n);
        return param.addObject(26, this);
    }
    if (strcmp(argv[0], "mom7n") == 0) {
        param.setValue(mom7n);
        return param.addObject(27, this);
    }
    if (strcmp(argv[0], "rot7n") == 0) {
        param.setValue(rot7n);
        return param.addObject(28, this);
    }




    if (strcmp(argv[0], "mom1") == 0) {
        param.setValue(mom1p);
        return param.addObject(29, this);
    }
    if (strcmp(argv[0], "rot1") == 0) {
        param.setValue(rot1p);
        return param.addObject(30, this);
    }
    if (strcmp(argv[0], "mom2") == 0) {
        param.setValue(mom2p);
        return param.addObject(31, this);
    }
    if (strcmp(argv[0], "rot2") == 0) {
        param.setValue(rot2p);
        return param.addObject(32, this);
    }
    if (strcmp(argv[0], "mom3") == 0) {
        param.setValue(mom3p);
        return param.addObject(33, this);
    }
    if (strcmp(argv[0], "rot3") == 0) {
        param.setValue(rot3p);
        return param.addObject(34, this);
    }

    if (strcmp(argv[0], "mom4") == 0) {
        param.setValue(mom4p);
        return param.addObject(35, this);
    }
    if (strcmp(argv[0], "rot4") == 0) {
        param.setValue(rot4p);
        return param.addObject(36, this);
    }
    if (strcmp(argv[0], "mom5") == 0) {
        param.setValue(mom5p);
        return param.addObject(37, this);
    }
    if (strcmp(argv[0], "rot5") == 0) {
        param.setValue(rot5p);
        return param.addObject(38, this);
    }
    if (strcmp(argv[0], "mom6") == 0) {
        param.setValue(mom6p);
        return param.addObject(39, this);
    }
    if (strcmp(argv[0], "rot6") == 0) {
        param.setValue(rot6p);
        return param.addObject(40, this);
    }
    if (strcmp(argv[0], "mom7") == 0) {
        param.setValue(mom7p);
        return param.addObject(41, this);
    }
    if (strcmp(argv[0], "rot7") == 0) {
        param.setValue(rot7p);
        return param.addObject(42, this);
    }


    return -1;
}

int
HystereticSMMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case -1:
        return -1;
    case 1:
        this->mom1p = info.theDouble;
        break;
    case 2:
        this->rot1p = info.theDouble;
        break;
    case 3:
        this->mom2p = info.theDouble;
        break;
    case 4:
        this->rot2p = info.theDouble;
        break;
    case 5:
        this->mom3p = info.theDouble;
        break;
    case 6:
        this->rot3p = info.theDouble;
        break;

    case 7:
        this->mom4p = info.theDouble;
        break;
    case 8:
        this->rot4p = info.theDouble;
        break;
    case 9:
        this->mom5p = info.theDouble;
        break;
    case 10:
        this->rot5p = info.theDouble;
        break;
    case 11:
        this->mom6p = info.theDouble;
        break;
    case 12:
        this->rot6p = info.theDouble;
        break;
    case 13:
        this->mom7p = info.theDouble;
        break;
    case 14:
        this->rot7p = info.theDouble;
        break;

    case 15:
        this->mom1n = info.theDouble;
        break;
    case 16:
        this->rot1n = info.theDouble;
        break;
    case 17:
        this->mom2n = info.theDouble;
        break;
    case 18:
        this->rot2n = info.theDouble;
        break;
    case 19:
        this->mom3n = info.theDouble;
        break;
    case 20:
        this->rot3n = info.theDouble;
        break;

    case 21:
        this->mom4n = info.theDouble;
        break;
    case 22:
        this->rot4n = info.theDouble;
        break;
    case 23:
        this->mom5n = info.theDouble;
        break;
    case 24:
        this->rot5n = info.theDouble;
        break;
    case 25:
        this->mom6n = info.theDouble;
        break;
    case 26:
        this->rot6n = info.theDouble;
        break;
    case 27:
        this->mom7n = info.theDouble;
        break;
    case 28:
        this->rot7n = info.theDouble;
        break;







    case 29:
        this->mom1p = info.theDouble;
        this->mom1n = -mom1p;
        break;
    case 30:
        this->rot1p = info.theDouble;
        this->rot1n = -rot1p;
        break;
    case 31:
        this->mom2p = info.theDouble;
        this->mom2n = -mom2p;
        break;
    case 32:
        this->rot2p = info.theDouble;
        this->rot2n = -rot2p;
        break;
    case 33:
        this->mom3p = info.theDouble;
        this->mom3n = -mom3p;
        break;
    case 34:
        this->rot3p = info.theDouble;
        this->rot3n = -rot3p;
        break;

    case 35:
        this->mom4p = info.theDouble;
        this->mom4n = -mom4p;
        break;
    case 36:
        this->rot4p = info.theDouble;
        this->rot4n = -rot4p;
        break;
    case 37:
        this->mom5p = info.theDouble;
        this->mom5n = -mom5p;
        break;
    case 38:
        this->rot5p = info.theDouble;
        this->rot5n = -rot5p;
        break;
    case 39:
        this->mom6p = info.theDouble;
        this->mom6n = -mom6p;
        break;
    case 40:
        this->rot6p = info.theDouble;
        this->rot6n = -rot6p;
        break;
    case 41:
        this->mom7p = info.theDouble;
        this->mom7n = -mom7p;
        break;
    case 42:
        this->rot7p = info.theDouble;
        this->rot7n = -rot7p;
        break;



    default:
        return -1;
    }

    this->setEnvelope();

    return 0;
}

void
HystereticSMMaterial::setEnvelope(void)
{
    E1p = mom1p / rot1p;
    E2p = (mom2p - mom1p) / (rot2p - rot1p);
    E3p = (mom3p - mom2p) / (rot3p - rot2p);
    E4p = (mom4p - mom3p) / (rot4p - rot3p);
    E5p = (mom5p - mom4p) / (rot5p - rot4p);
    E6p = (mom6p - mom5p) / (rot6p - rot5p);
    E7p = (mom7p - mom6p) / (rot7p - rot6p);

    E1n = mom1n / rot1n;
    E2n = (mom2n - mom1n) / (rot2n - rot1n);
    E3n = (mom3n - mom2n) / (rot3n - rot2n);
    E4n = (mom4n - mom3n) / (rot4n - rot3n);
    E5n = (mom5n - mom4n) / (rot5n - rot4n);
    E6n = (mom6n - mom5n) / (rot6n - rot5n);
    E7n = (mom7n - mom6n) / (rot7n - rot6n);



    Eup = E1p;
    if (E2p > Eup) Eup = E2p;
    if (E3p > Eup) Eup = E3p;
    if (E4p > Eup) Eup = E4p;
    if (E5p > Eup) Eup = E5p;
    if (E6p > Eup) Eup = E6p;
    if (E7p > Eup) Eup = E7p;

    Eun = E1n;
    if (E2n > Eun) Eun = E2n;
    if (E3n > Eun) Eun = E3n;
    if (E4n > Eun) Eun = E4n;
    if (E5n > Eun) Eun = E5n;
    if (E6n > Eun) Eun = E6n;
    if (E7n > Eun) Eun = E7n;
}

double
HystereticSMMaterial::posEnvlpStress(double strain)
{
    if (strain <= 0.0)
        return 0.0;
    else if (strain <= rot1p)
        return E1p * strain;
    else if (strain <= rot2p)
        return mom1p + E2p * (strain - rot1p);
    else if (strain <= rot3p)
        return mom2p + E3p * (strain - rot2p);
    else if (strain <= rot4p)
        return mom3p + E4p * (strain - rot3p);
    else if (strain <= rot5p)
        return mom4p + E5p * (strain - rot4p);
    else if (strain <= rot6p)
        return mom5p + E6p * (strain - rot5p);
    else if (strain <= rot7p || E7p > 0.0)   // removed 
        return mom6p + E7p * (strain - rot6p);
    else
        return mom7p;
}

double
HystereticSMMaterial::negEnvlpStress(double strain)
{
    if (strain >= 0.0)
        return 0.0;
    else if (strain >= rot1n)
        return E1n * strain;
    else if (strain >= rot2n)
        return mom1n + E2n * (strain - rot1n);
    else if (strain >= rot3n)
        return mom2n + E3n * (strain - rot2n);
    else if (strain >= rot4n)
        return mom3n + E4n * (strain - rot3n);
    else if (strain >= rot5n)
        return mom4n + E5n * (strain - rot4n);
    else if (strain >= rot6n)
        return mom5n + E6n * (strain - rot5n);
    else if (strain >= rot7n || E7p > 0.0)  // removed 
        return mom6n + E7n * (strain - rot6n);
    else
        return mom7n;
}

double
HystereticSMMaterial::posEnvlpTangent(double strain)
{
    if (strain < 0.0)
        return E1p * 1.0e-9;
    else if (strain <= rot1p)
        return E1p;
    else if (strain <= rot2p)
        return E2p;
    else if (strain <= rot3p || E3p > 0.0)
        return E3p;
    else if (strain <= rot4p || E4p > 0.0)
        return E4p;
    else if (strain <= rot5p || E5p > 0.0)
        return E5p;
    else if (strain <= rot6p || E6p > 0.0)
        return E6p;
    else if (strain <= rot7p || E7p > 0.0) // removed 
        return E7p;
    else
        return E1p * 1.0e-9;
}

double
HystereticSMMaterial::negEnvlpTangent(double strain)
{
    if (strain > 0.0)
        return E1n * 1.0e-9;
    else if (strain >= rot1n)
        return E1n;
    else if (strain >= rot2n)
        return E2n;
    else if (strain >= rot3n || E3n > 0.0)
        return E3n;
    else if (strain >= rot4n || E4n > 0.0)
        return E4n;
    else if (strain >= rot5n || E5n > 0.0)
        return E5n;
    else if (strain >= rot6n || E6n > 0.0)
        return E6n;
    else if (strain >= rot7n || E7n > 0.0)
        return E7n;
    else
        return E1n * 1.0e-9;
}

double
HystereticSMMaterial::posEnvlpRotlim(double strain)
{
    double strainLimit = POS_INF_STRAIN;

    if (strain <= rot1p)
        return POS_INF_STRAIN;
    if (strain > rot1p && strain <= rot2p)
        strainLimit = rot1p + (mom2p - mom1p) / E2p;
    if (strain > rot2p && strain <= rot3p)
        strainLimit = rot2p + (mom3p - mom2p) / E3p;
    if (strain > rot3p && strain <= rot4p)
        strainLimit = rot3p + (mom4p - mom3p) / E4p;
    if (strain > rot4p && strain <= rot5p)
        strainLimit = rot4p + (mom5p - mom4p) / E5p;
    if (strain > rot5p && strain <= rot6p)
        strainLimit = rot5p + (mom6p - mom5p) / E6p;
    if (strain > rot6p && E7p < 0.0)
        strainLimit = rot6p + (mom7p - mom6p) / E7p;

    //if (strain <= rot1p)
    //    return POS_INF_STRAIN;
    //if (strain > rot1p && strain <= rot2p && E2p < 0.0)
    //    strainLimit = rot1p - mom1p / E2p;
    //if (strain > rot2p && strain <= rot3p && E3p < 0.0)
    //    strainLimit = rot2p - mom2p / E3p;
    //if (strain > rot3p && strain <= rot4p && E4p < 0.0)
    //    strainLimit = rot3p - mom3p / E4p;
    //if (strain > rot4p && strain <= rot5p && E5p < 0.0)
    //    strainLimit = rot4p - mom4p / E5p;
    //if (strain > rot5p && strain <= rot6p && E6p < 0.0)
    //    strainLimit = rot5p - mom5p / E6p;
    //if (strain > rot6p && E7p < 0.0)
    //    strainLimit = rot6p - mom6p / E7p;

    //if (strain <= rot1p)
    //    return POS_INF_STRAIN;
    //if (strain > rot1p && strain <= rot2p)
    //    strainLimit = rot1p - mom1p / E2p;
    //if (strain > rot2p && strain <= rot3p)
    //    strainLimit = rot2p - mom2p / E3p;
    //if (strain > rot3p && strain <= rot4p)
    //    strainLimit = rot3p - mom3p / E4p;
    //if (strain > rot4p && strain <= rot5p)
    //    strainLimit = rot4p - mom4p / E5p;
    //if (strain > rot5p && strain <= rot6p)
    //    strainLimit = rot5p - mom5p / E6p;
    //if (strain > rot6p && E7p < 0.0)
    //    strainLimit = rot6p - mom6p / E7p;

    if (strainLimit == POS_INF_STRAIN)
        return POS_INF_STRAIN;
    else if (posEnvlpStress(strainLimit) > 0)
        return POS_INF_STRAIN;
    else
        return strainLimit;

}

double
HystereticSMMaterial::negEnvlpRotlim(double strain)
{
    double strainLimit = NEG_INF_STRAIN;

    if (strain >= rot1n)
        return NEG_INF_STRAIN;
    if (strain < rot1n && strain >= rot2n)
        strainLimit = rot1n + (mom2n - mom1n) / E2n;
    if (strain < rot2n && strain >= rot3n)
        strainLimit = rot2n + (mom3n - mom2n) / E3n;
    if (strain < rot3n && strain >= rot4n)
        strainLimit = rot3n + (mom4n - mom3n) / E4n;
    if (strain < rot4n && strain >= rot5n)
        strainLimit = rot4n + (mom5n - mom4n) / E5n;
    if (strain < rot5n && strain >= rot6n)
        strainLimit = rot5n + (mom6n - mom5n) / E6n;
    if (strain < rot6n && E7n < 0.0)
        strainLimit = rot6n + (mom7n - mom6n) / E7n;


    if (strainLimit == NEG_INF_STRAIN)
        return NEG_INF_STRAIN;
    else if (negEnvlpStress(strainLimit) < 0)
        return NEG_INF_STRAIN;
    else
        return strainLimit;
}

// All the response options were implemented by Silvia Mazzoni, 2022
Response*
HystereticSMMaterial::setResponse(const char** argv, int argc, OPS_Stream& theOutput)
{

    if (strcmp(argv[0], "MUy") == 0 || strcmp(argv[0], "MU1") == 0) {
        return new MaterialResponse(this, 11, 0.0);
    }

    if (strcmp(argv[0], "defoPlastic") == 0 || strcmp(argv[0], "thetaP") == 0) {
        return new MaterialResponse(this, 21, 0.0);
    }


    // backbone pts DCR
    // deformation
    // current step
    else if (strcmp(argv[0], "defoDCR") == 0) {
        return new MaterialResponse(this, 311, Vector(7));
    }
    // max
    else if (strcmp(argv[0], "defoDCRMax") == 0) {
        return new MaterialResponse(this, 312, Vector(7));
    }

    // user-defined limit states DCR
    // deformation
    else if (strcmp(argv[0], "defoLimitStates") == 0) {
        return new MaterialResponse(this, 96, Vector(nDefoLimitStates));
    }
    else if (strcmp(argv[0], "defoLimitStatesDCR") == 0) {
        return new MaterialResponse(this, 961, Vector(nDefoLimitStates));
    }

    else if (strcmp(argv[0], "defoLimitStatesDCRMax") == 0) {
        return new MaterialResponse(this, 962, Vector(nDefoLimitStates));
    }
    else if (strcmp(argv[0], "defoLimitStatesDCRMaxAbs") == 0) {
        return new MaterialResponse(this, 963, Vector(nDefoLimitStates));
    }

    // user-defined limit states DCR
    // force
    else if (strcmp(argv[0], "forceLimitStates") == 0) {
        return new MaterialResponse(this, 97, Vector(nDefoLimitStates));
    }

    else if (strcmp(argv[0], "forceLimitStatesDCR") == 0) {
        return new MaterialResponse(this, 971, Vector(7));
    }

    else if (strcmp(argv[0], "AllData") == 0) {
        return new MaterialResponse(this, 99, Vector(43));
    }
    //else if (strcmp(argv[0], "AllDataHeader") == 0) {
    //    return new MaterialResponse(this, 991, 0.0);
    //}

    //by default, See if the response is one of the defaults
    Response* res = UniaxialMaterial::setResponse(argv, argc, theOutput);

    if (res != 0)      return res;
    else {
        opserr << "error in HystereticSMMaterial::setResponse" << endln;
        return 0;
    }

}

int
HystereticSMMaterial::getResponse(int responseID, Information& matInfo)
{

    if (responseID == 11) {
        if (Cstrain > 0)
            return matInfo.setDouble(this->Cstrain / rot1p);
        else
            return matInfo.setDouble(this->Cstrain / rot1n);
    }

    // plastic deformation
    else if (responseID == 21) {
        if (Cstrain > 0) {
            if (Cstrain > rot1p) {
                return matInfo.setDouble(this->Cstrain - Cstress / mom1p * rot1p);
            }
        }
        else if (Cstrain < 0) {
            if (Cstrain < rot1n) {
                return matInfo.setDouble(this->Cstrain - Cstress / mom1n * rot1n);
            }
        }
        return matInfo.setDouble(0.0);
    }





    // backbone pts DCR
    // deformation

    // current
    else if (responseID == 311) {
        static Vector data(7);
        if (Cstrain > 0) {
            data(0) = this->Cstrain / rot1p;
            data(1) = this->Cstrain / rot2p;
            data(2) = this->Cstrain / rot3p;
            data(3) = this->Cstrain / rot4p;
            data(4) = this->Cstrain / rot5p;
            data(5) = this->Cstrain / rot6p;
            data(6) = this->Cstrain / rot7p;
            data(0) = this->Cstrain / rot1p;
        }
        else {
            data(0) = this->Cstrain / rot1n;
            data(1) = this->Cstrain / rot2n;
            data(2) = this->Cstrain / rot3n;
            data(3) = this->Cstrain / rot4n;
            data(4) = this->Cstrain / rot5n;
            data(5) = this->Cstrain / rot6n;
            data(6) = this->Cstrain / rot7n;
            data(0) = this->Cstrain / rot1n;
        }

        return matInfo.setVector(data);
    }
    // max
    else if (responseID == 312) {
        static Vector data(14);
        data(0) = this->CrotPu / rot1p;
        data(0 + 7) = this->CrotNu / rot1n;
        data(1) = this->CrotPu / rot2p;
        data(1 + 7) = this->CrotNu / rot2n;
        data(2) = this->CrotPu / rot3p;
        data(2 + 7) = this->CrotNu / rot3n;
        data(3) = this->CrotPu / rot4p;
        data(3 + 7) = this->CrotNu / rot4n;
        data(4) = this->CrotPu / rot5p;
        data(4 + 7) = this->CrotNu / rot5n;
        data(5) = this->CrotPu / rot6p;
        data(5 + 7) = this->CrotNu / rot6n;
        data(6) = this->CrotPu / rot7p;
        data(6 + 7) = this->CrotNu / rot7n;
        return matInfo.setVector(data);
    }




    // user-defined limit states DCR
    // defo
    // input values
    else if (responseID == 96) {
        return matInfo.setVector(defoLimitStates);
    }
    // current step
    else if (responseID == 961) {
        static Vector data(nDefoLimitStates);
        for (int i = 0; i < nDefoLimitStates; i++) {
            data(i) = this->Cstrain / defoLimitStates[i];
        }
        return matInfo.setVector(data);
    }

    // Maximum value
    else if (responseID == 962) {
        static Vector data(nDefoLimitStates);
        for (int i = 0; i < nDefoLimitStates; i++) {
            if (defoLimitStates[i] > 0) {
                data(i) = this->CrotPu / defoLimitStates[i];
            }
            else {
                data(i) = this->CrotNu / defoLimitStates[i];
            }
        }
        return matInfo.setVector(data);
    }

    // Maximum absolute
    else if (responseID == 963) {
        static Vector data(nDefoLimitStates);
        for (int i = 0; i < nDefoLimitStates; i++) {
            data(i) = fabs(CrotNu / defoLimitStates[i]);
            if (fabs(CrotPu / defoLimitStates[i]) > data(i)) {
                data(i) = fabs(CrotPu / defoLimitStates[i]);
            }
        }
        return matInfo.setVector(data);
    }



    // user-defined limit states DCR
    // force
    // input values
    else if (responseID == 97) {
        return matInfo.setVector(forceLimitStates);
    }
    // current step
    else if (responseID == 971) {
        static Vector data(nDefoLimitStates);
        for (int i = 0; i < nDefoLimitStates; i++) {
            data(i) = this->Cstress / forceLimitStates[i];
        }
        return matInfo.setVector(data);
    }




    else if (responseID == 991) {
        // doesn't work, as it is a string, need to work on it...
        return matInfo.setString("matTag,mom1p,rot1p,mom2p,rot2p,mom3p,rot3p,mom4p,rot4p,mom5p,rot5p,mom6p,rot6p,mom7p,rot7p,mom1n,rot1n,mom2n,rot2n,mom3n,rot3n,mom4n,rot4n,mom5n,rot5n,mom6n,rot6n,mom7n,rot7n,pinchX,pinchY,damfc1,damfc2,beta,CrotMax,CrotMin,CrotPu,CrotNu,CenergyD,CloadIndicator,Cstress,Cstrain,Ttangent");

    }

    else if (responseID == 99) {
        // return matInfo.setDouble(this->getThetaP());
        //return matInfo.setVector(this ->getAllData());

        //stressStrainTangent(0) = this->getStress();
        //stressStrainTangent(1) = this->getStrain();
        //stressStrainTangent(2) = this->getTangent();
        //matInfo.setVector(stressStrainTangent);

        static Vector data(43);

        data(0) = this->getTag();
        data(1) = mom1p;
        data(2) = rot1p;
        data(3) = mom2p;
        data(4) = rot2p;
        data(5) = mom3p;
        data(6) = rot3p;
        data(7) = mom4p;
        data(8) = rot4p;
        data(9) = mom5p;
        data(10) = rot5p;
        data(11) = mom6p;
        data(12) = rot6p;
        data(13) = mom7p;
        data(14) = rot7p;
        data(15) = mom1n;
        data(16) = rot1n;
        data(17) = mom2n;
        data(18) = rot2n;
        data(19) = mom3n;
        data(20) = rot3n;
        data(21) = mom4n;
        data(22) = rot4n;
        data(23) = mom5n;
        data(24) = rot5n;
        data(25) = mom6n;
        data(26) = rot6n;
        data(27) = mom7n;
        data(28) = rot7n;
        data(29) = pinchX;
        data(30) = pinchY;
        data(31) = damfc1;
        data(32) = damfc2;
        data(33) = beta;
        data(34) = CrotMax;
        data(35) = CrotMin;
        data(36) = CrotPu;
        data(37) = CrotNu;
        data(38) = CenergyD;
        data(39) = CloadIndicator;
        data(40) = Cstress;
        data(41) = Cstrain;
        data(42) = Ttangent;

        return matInfo.setVector(data);


    }

    else {

        // Just call the base class method ... don't need to define
        // this function, but keeping it here just for clarity
        return UniaxialMaterial::getResponse(responseID, matInfo);

    }

}

