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

// $Revision: 1.03 $
// $Date: 2021/12/13 00:00:00 $
// Written: Hanlin Dong, Xijun Wang, Tongji University, self@hanlindong.com
//
// Description: This file contains the class definition for DowelType.
// DowelType provides the abstraction of a dowel-type timber joint, including nails, screws, and bolts.
// The envelope curve can be chosen from (1) exponential, (2) Bezier, and (3) piecewise linear.
// The hysteretic law considers force intercept variation, strength degradation,
// unloading, pinching, and reloading stiffness degradation, etc.
// 
// Material definition tcl command:
// uniaxialMaterial DowelType $tag
//     $Fi $Kp $Ru $c $beta $gamma $eta $Dy $alpha_p $alpha_u $alpha_r 
//     <-exponential $K0 $R1 $F0 $Dc $Kd <$Du> <$K0N $R1N $F0N $DcN $KdN $DuN> >
//     <-bezier $D1 $F1 $D2 $F2 $Dc $Fc $Kd <$Du> <$D1N $F1N $D2N $F2N $DcN $FcN $KdN $DuN>> >
//     <-piecewise $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>>

#include <elementAPI.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include "DowelType.h"


#define PRECISION 1.0e-12
#define PI 3.14159265358979323846
#define MAX_ITER 2000
#define DEBUG false

static int numDowelType = 0;

void *
OPS_DowelType()
{
    if (numDowelType == 0) {
        opserr << "DowelType v1.03 - Written by Hanlin Dong (self@hanlindong.com) and Xijun Wang ";
        opserr << "from Tongji University, Copyright 2021 - Use at your Own Peril" << endln;
        numDowelType = 1;
    }
    UniaxialMaterial *theMaterial = 0;

    // parse the input line for the material parameters
    // parsing material tag
    int numData;
    int tag;
    numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "ERROR: invalid uniaxialMaterial DowelType tag" << endln;
        return 0;
    }

    // parsing hysteresis data
    double dData[11];
    numData = 11;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "ERROR: expected $Fi $Kp $Ru $c $gamma $eta";
        opserr << "$Dy $alpha_p $alpha_u $alpha_r" << endln;
        return 0;
    }

    // parsing envelope type
    const char* envstr = OPS_GetString();
    if (strcmp(envstr, "-exponential") == 0) {
        numData = OPS_GetNumRemainingInputArgs();
        double eData[12];
        if (OPS_GetDoubleInput(&numData, eData) != 0) {
            opserr << "ERROR: expected exponential envelope parameters: ";
            opserr << "$K0 $R1 $F0 $Dc $Kd <$Du> <$K0N $R1N $F0N $DcN $KdN <$DuN>>" << endln;
            return 0;   
        }
        if (numData == 5) { 
            eData[5] = 0.;
            eData[6] = eData[0];
            eData[7] = eData[1];
            eData[8] = -eData[2];
            eData[9] = -eData[3];
            eData[10] = eData[4];
            eData[11] = 0.;
        } else if (numData == 6) {
            eData[6] = eData[0];
            eData[7] = eData[1];
            eData[8] = -eData[2];
            eData[9] = -eData[3];
            eData[10] = eData[4];
            eData[11] = -eData[5];
        } else if (numData == 10) {
            // leave space for dult
            for (int i = 10; i > 5; i--) {
                eData[i] = eData[i-1];
            }
            eData[5] = 0.;
            eData[11] = 0.;
        } else if (numData != 12) {
            opserr << "ERROR: invailed number of args (should be 5, 6, 10, or 12). Expected:";
            opserr << "$K0 $R1 $F0 $Dc $Kd <$Du> <$K0N $R1N $F0N $DcN $KdN <$DuN>>" << endln;
            return 0;
        }
        theMaterial = new DowelType(tag, 
            dData[0], dData[1], dData[2], dData[3], 
            dData[4], dData[5], dData[6], 
            dData[7], dData[8], dData[9], dData[10],
            eData[0], eData[1], eData[2], eData[3], eData[4], eData[5], 
            eData[6], eData[7], eData[8], eData[9], eData[10], eData[11]);
    } else if (strcmp(envstr, "-bezier") == 0) {
        numData = OPS_GetNumRemainingInputArgs();
        double eData[16];
        if (OPS_GetDoubleInput(&numData, eData) != 0) {
            opserr << "ERROR: expected Bezier envelope parameters: ";
            opserr << "$Db1 $Fb1 $Db2 $Fb2 $Dc $Fc $Kd <$Du> ";
            opserr << "<$Db1N $Fb1N $Db2N $Fb2N $DcN $FcN $KdN <$DuN>>" << endln;
        }
        if (numData == 7) {
            eData[7] = 0.;
            for (int i = 8; i <= 13; i++) {
                eData[i] = -eData[i-8];
            }
            eData[14] = eData[6];
            eData[15] = 0.;
        } else if (numData == 8) {
            for (int i = 8; i <= 13; i++) {
                eData[i] = -eData[i-8];
            }
            eData[14] = eData[6];
            eData[15] = -eData[7];
        } else if (numData == 14) {
            eData[15] = 0.;
            for (int i = 14; i >= 8; i--) {
                eData[i] = eData[i-1];
            }
            eData[7] = 0.;
        } else if (numData != 16) {
            opserr << "ERROR: wrong number of args (should be 7, 8, 14, or 16). Expected: ";
            opserr << "$Db1 $Fb1 $Db2 $Fb2 $Dc $Fc $Kd <$Du> <$Db1N $Fb1N $Db2N $Fb2N $DcN $FcN $KdN <$DuN>>" << endln;
            return 0;
        }
        theMaterial = new DowelType(tag, 
            dData[0], dData[1], dData[2], dData[3], 
            dData[4], dData[5], dData[6], 
            dData[7], dData[8], dData[9], dData[10],
            eData[0], eData[1], eData[2], eData[3], 
            eData[4], eData[5], eData[6], eData[7], 
            eData[8], eData[9], eData[10], eData[11],
            eData[12], eData[13], eData[14], eData[15]);
    } else if (strcmp(envstr, "-piecewise") == 0) {
        numData = OPS_GetNumRemainingInputArgs();
        if (numData < 6) {
            opserr << "ERROR: piecewise envelope: number of definition points must be no less than 3. ";
            opserr << "Expected: $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>" << endln;
            return 0;
        }
        if (numData % 2 != 0) {
            opserr << "ERROR: piecewise envelope: number of definition coordinates must be even. ";
            opserr << "Expected: $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>" << endln;
            return 0;
        }
        // maximum 20 points.
        if (numData > 40) {
            opserr << "ERROR: piecewise envelope: too many definition points (should be no more than 20 points). ";
            opserr << "Expected: $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>" << endln;
            return 0;
        }
        double * eData = new double[numData];
        if (OPS_GetDoubleInput(&numData, eData) != 0) {
            opserr << "ERROR: piecewise envelope: cannot read coordinates.  ";
            opserr << "Expected: $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>" << endln;
            return 0;
        }
        bool hasNeg = false;  // has negative displacement definition
        for (int i = 0; i < numData; i += 2) {
            if (eData[i] < 0) {
                hasNeg = true;
                break;
            }
        }
        int size = hasNeg ? numData / 2 + 1 : numData + 1; // size of env array.
        double * dCoords = new double[size];  // disp array
        double * fCoords = new double[size];  // force array
        dCoords[0] = 0.0;  // add origin
        fCoords[0] = 0.0;  // add origin
        int numSorted = 1; // points sorted (currently the origin only)
        for (int i = 0; i < numData; i += 2) {
            for (int j = numSorted; j >= 0; j--) {
                if (eData[i] > dCoords[j-1]) {
                    dCoords[j] = eData[i];
                    fCoords[j] = eData[i+1];
                    break;
                } else if (j == 0) {
                    dCoords[j] = eData[i];
                    fCoords[j] = eData[i+1];
                } else {
                    dCoords[j] = dCoords[j-1];
                    fCoords[j] = fCoords[j-1];
                }
            }
            numSorted ++;
        }
        // copy symmetric points
        if (!hasNeg) {
            for (int i = size - 1; i >= 0; i--) {
                if (i >= numSorted-1) {
                    dCoords[i] = dCoords[i - numSorted + 1];
                    fCoords[i] = fCoords[i - numSorted + 1];
                } else {
                    dCoords[i] = -dCoords[size - i - 1];
                    fCoords[i] = -fCoords[size - i - 1];
                }
            }
        }
        if (eData) {
            delete [] eData;
        }
        theMaterial = new DowelType(tag, 
                    dData[0], dData[1], dData[2], dData[3], 
                    dData[4], dData[5], dData[6], 
                    dData[7], dData[8], dData[9], dData[10],
                    size, dCoords, fCoords);
    } else {
        opserr << "ERROR: unsupported envelope type. Expected: -exponential, -bezier or -piecewise" << endln;
        return 0;
    }

    // Fail to create the material
    if (theMaterial == 0) {
        opserr << "ERROR: could not create uniaxialMaterial of type DowelType" << endln;
        return 0;
    }

    return theMaterial;
}


DowelType::DowelType(int tag, 
                double Fi, double Kp, double Ru, double c,
                double beta, double gamma, double eta,
                double Dy, double Ap, double Au, double Ar, 
                double K0_p, double R1_p, double F0_p, double Dc_p, double Kd_p, double Du_p,
                double K0_n, double R1_n, double F0_n, double Dc_n, double Kd_n, double Du_n):
UniaxialMaterial(tag, MAT_TAG_DowelType), 
fi(Fi), kp(Kp), ru(Ru), c(c),
beta(beta), gamma(gamma), eta(eta),
dyield(Dy), alpha_p(Ap), alpha_u(Au), alpha_r(Ar),
envType(1), 
k0_p(K0_p), k0_n(K0_n), dcap_p(Dc_p), dcap_n(Dc_n), kdesc_p(Kd_p), kdesc_n(Kd_n), 
k1_p(R1_p*K0_p), k1_n(R1_n*K0_n), f0_p(F0_p), f0_n(F0_n),
denv1_p(0.0), denv1_n(0.0), fenv1_p(0.0), fenv1_n(0.0), 
denv2_p(0.0), denv2_n(0.0), fenv2_p(0.0), fenv2_n(0.0),
envSize(0), envZero(0),
isPHC(true), ePHC_p(0.0), ePHC_n(0.0), eFHC_p(0.0), eFHC_n(0.0),
tStrain(0.0), tStress(0.0), tTangent(K0_p),
cStrain(0.0), cStress(0.0), cTangent(K0_p),
tPath(1), tDmin(0.0), tFdmin(0.0), tDmax(0.0), tFdmax(0.0),
cPath(1), cDmin(0.0), cFdmin(0.0), cDmax(0.0), cFdmax(0.0)
{
    fcap_p = (F0_p + k1_p * Dc_p) * (1 - exp(-1 * K0_p * Dc_p / F0_p));
    fcap_n = (F0_n + k1_n * Dc_n) * (1 - exp(-1 * K0_n * Dc_n / F0_n));
    dult_p = (fabs(Du_p) < DBL_EPSILON) ? dcap_p + fcap_p / kdesc_p : Du_p;
    dult_n = (fabs(Du_n) < DBL_EPSILON) ? dcap_n + fcap_n / kdesc_n : Du_n;
    fyield_p = envelope(Dy);
    fyield_n = envelope(-Dy);
    dinter_p = envIntersection(Kp, Fi);
    dinter_n = envIntersection(Kp, -Fi);
    double x_p = exp(K0_p * Dc_p / F0_p); // integral exponential term
    eMono_p = Dc_p * Dc_p * k1_p / 2 -
        pow(F0_p / K0_p, 2) * (K0_p + k1_p) / x_p * (x_p - 1) +
        Dc_p * (F0_p * k1_p / x_p / K0_p + F0_p) +
        0.5 * (1 + 0.8) * (1 - 0.8) * pow(fcap_p, 2) / kdesc_p;
    double x_n = exp(K0_n * Dc_n / F0_n); // integral exponential term
    eMono_n = Dc_n * Dc_n * k1_n / 2 -
        pow(F0_n / K0_n, 2) * (K0_n + k1_n) / x_n * (x_n - 1) +
        Dc_n * (F0_n * k1_n / x_n / K0_n + F0_n) +
        0.5 * (1 + 0.8) * (1 - 0.8) * pow(fcap_n, 2) / kdesc_n;
    for (int i=0; i < 20; i++) {
        pxs[i] = 0.0;
        pys[i] = 0.0;
    }
    denvs = NULL;
    fenvs = NULL;
}


DowelType::DowelType(int tag, 
                double Fi, double Kp, double Ru, double c,
                double beta, double gamma, double eta,
                double Dy, double Ap, double Au, double Ar, 
                double D1_p, double F1_p, double D2_p, double F2_p,
                double Dc_p, double Fc_p, double Kd_p, double Du_p,
                double D1_n, double F1_n, double D2_n, double F2_n,
                double Dc_n, double Fc_n, double Kd_n, double Du_n):
UniaxialMaterial(tag, MAT_TAG_DowelType), 
fi(Fi), kp(Kp), ru(Ru), c(c),
beta(beta), gamma(gamma), eta(eta),
dyield(Dy), alpha_p(Ap), alpha_u(Au), alpha_r(Ar),
envType(2), 
k0_p(F1_p / D1_p), k0_n(F1_n / D1_n), dcap_p(Dc_p), dcap_n(Dc_n), 
fcap_p(Fc_p), fcap_n(Fc_n), kdesc_p(Kd_p), kdesc_n(Kd_n), 
k1_p(0.0), k1_n(0.0), f0_p(0.0), f0_n(0.0),
denv1_p(D1_p), denv1_n(D1_n), fenv1_p(F1_p), fenv1_n(F1_n), 
denv2_p(D2_p), denv2_n(D2_n), fenv2_p(F2_p), fenv2_n(F2_n),
envSize(0), envZero(0),
isPHC(true), ePHC_p(0.0), ePHC_n(0.0), eFHC_p(0.0), eFHC_n(0.0), 
tStrain(0.0), tStress(0.0), tTangent(k0_p),
cStrain(0.0), cStress(0.0), cTangent(k0_p),
tPath(1), tDmin(0.0), tFdmin(0.0), tDmax(0.0), tFdmax(0.0),
cPath(1), cDmin(0.0), cFdmin(0.0), cDmax(0.0), cFdmax(0.0)
{
    dult_p = fabs(Du_p) < DBL_EPSILON ? Dc_p + fcap_p / kdesc_p : Du_p;
    dult_n = fabs(Du_n) < DBL_EPSILON ? Dc_n + fcap_n / kdesc_n : Du_n;
    fyield_p = envelope(Dy);
    fyield_n = envelope(-Dy);
    dinter_p = envIntersection(Kp, Fi);
    dinter_n = envIntersection(Kp, -Fi);
    // seed 20 segs for monotonic energy
    int seeds = 20;
    eMono_p = 0.0;
    double xprev = 0.0;
    double yprev = 0.0;
    for (int i = 0; i < seeds; i++) {
        double t = i * 0.1 + 0.1;
        double x = 3 * pow(1-t, 2) * t * denv1_p + 
                   3 * (1-t) * pow(t, 2) * denv2_p +
                   pow(t, 3) * dcap_p;
        double y = 3 * pow(1-t, 2) * t * fenv1_p +
                   3 * (1-t) * pow(t, 2) * fenv2_p +
                   pow(t, 3) * fcap_p;
        eMono_p += 0.5 * (y + yprev) * (x - xprev);
        xprev = x;
        yprev = y;
    }
    eMono_n = 0.0;
    xprev = 0.0;
    yprev = 0.0;
    for (int i = 0; i < seeds; i++) {
        double t = i * 0.1 + 0.1;
        double x = 3 * pow(1-t, 2) * t * denv1_n + 
                   3 * (1-t) * pow(t, 2) * denv2_n +
                   pow(t, 3) * dcap_n;
        double y = 3 * pow(1-t, 2) * t * fenv1_n +
                   3 * (1-t) * pow(t, 2) * fenv2_n +
                   pow(t, 3) * fcap_n;
        eMono_n += 0.5 * (y + yprev) * (x - xprev);
        xprev = x;
        yprev = y;
    }

    for (int i = 0; i < 20; i++) {
        pxs[i] = 0.0;
        pys[i] = 0.0;
    }
    denvs = NULL;
    fenvs = NULL;
}

DowelType::DowelType(int tag, 
                double Fi, double Kp, double Ru, double c,
                double beta, double gamma, double eta,
                double Dy, double Ap, double Au, double Ar, 
                int Size, double *Denvs, double *Fenvs):
UniaxialMaterial(tag, MAT_TAG_DowelType), 
fi(Fi), kp(Kp), ru(Ru), c(c),
beta(beta), gamma(gamma), eta(eta),
dyield(Dy), alpha_p(Ap), alpha_u(Au), alpha_r(Ar),
envType(3), 
k1_p(0.0), k1_n(0.0), f0_p(0.0), f0_n(0.0),
denv1_p(0.0), denv1_n(0.0), fenv1_p(0.0), fenv1_n(0.0), 
denv2_p(0.0), denv2_n(0.0), fenv2_p(0.0), fenv2_n(0.0),
envSize(Size),
isPHC(true), ePHC_p(0.0), ePHC_n(0.0), eFHC_p(0.0), eFHC_n(0.0),
tStrain(0.0), tStress(0.0), 
cStrain(0.0), cStress(0.0), 
tPath(1), tDmin(0.0), tFdmin(0.0), tDmax(0.0), tFdmax(0.0),
cPath(1), cDmin(0.0), cFdmin(0.0), cDmax(0.0), cFdmax(0.0)
{
    denvs = Denvs;
    fenvs = Fenvs;
    dcap_p = denvs[0];
    fcap_p = fenvs[0];
    dcap_n = denvs[0];
    fcap_n = denvs[0];
    for (int i = 1; i < Size; i++) {
        if (fenvs[i] > fcap_p) {
            dcap_p = denvs[i];
            fcap_p = fenvs[i];
        }
        if (fenvs[i] < fcap_n) {
            dcap_n = denvs[i];
            fcap_n = fenvs[i];
        }
    }
    k0_p = 0.;
    k0_n = 0.;
    envZero = 0;
    for (int i = 0; i < Size; i++) {
        if (fabs(denvs[i]) < PRECISION) {
            envZero = i;
            break;
        }
    }
    k0_p = fenvs[envZero+1] / denvs[envZero+1];
    k0_n = fenvs[envZero-1] / denvs[envZero-1];
    tTangent = k0_p;
    cTangent = k0_p;

    dult_p = denvs[Size-1];
    dult_n = denvs[0];

    dinter_p = envIntersection(Kp, Fi);
    dinter_n = envIntersection(Kp, -Fi);
    fyield_p = envelope(Dy);
    fyield_n = envelope(-Dy);

    eMono_p = 0.0;
    for (int i = envZero + 1; i < Size; i++) {
        eMono_p += 0.5 * (fenvs[i] + fenvs[i-1]) * (denvs[i] - denvs[i-1]);
    }
    eMono_n = 0.0;
    for (int i = envZero; i > 0; i--) {
        eMono_n += 0.5 * (fenvs[i] + fenvs[i-1]) * (denvs[i] - denvs[i-1]);
    }

    for (int i = 0; i < 20; i++) {
        pxs[i] = 0.0;
        pys[i] = 0.0;
    }
}


DowelType::DowelType():
UniaxialMaterial(0, MAT_TAG_DowelType), 
fi(0.0), kp(0.0), ru(0.0), c(0.0),
beta(0.0), gamma(0.0), eta(0.0),
dyield(0.0), alpha_p(0.0), alpha_u(0.0), alpha_r(0.0),
envType(0), 
k0_p(0.0), k0_n(0.0), dcap_p(0.0), dcap_n(0.0), fcap_p(0.0), fcap_n(0.0), 
fyield_p(0.0), fyield_n(0.0), dult_p(0.0), dult_n(0.0), kdesc_p(0.0), kdesc_n(0.0), 
dinter_p(0.0), dinter_n(0.0), eMono_p(0.0), eMono_n(0.0),
k1_p(0.0), k1_n(0.0), f0_p(0.0), f0_n(0.0),
denv1_p(0.0), denv1_n(0.0), fenv1_p(0.0), fenv1_n(0.0), 
denv2_p(0.0), denv2_n(0.0), fenv2_p(0.0), fenv2_n(0.0),
envSize(0), envZero(0),
isPHC(true), ePHC_p(0.0), ePHC_n(0.0), eFHC_p(0.0), eFHC_n(0.0),
tStrain(0.0), tStress(0.0), tTangent(0.0),
cStrain(0.0), cStress(0.0), cTangent(0.0),
tPath(1), tDmin(0.0), tFdmin(0.0), tDmax(0.0), tFdmax(0.0),
cPath(1), cDmin(0.0), cFdmin(0.0), cDmax(0.0), cFdmax(0.0)
{
    for (int i = 0; i < 20; i++) {
        pxs[i] = 0.0;
        pys[i] = 0.0;
    }
    denvs = NULL;
    fenvs = NULL;
}


DowelType::~DowelType()
{
    if (denvs) {
        delete [] denvs;
        denvs = NULL;
    }
    if (fenvs) {
        delete [] fenvs;
        fenvs = NULL;
    }
}

double DowelType::envelope(double disp)
{
    double force;
    if (envType == 1) {
        force = 
            disp < dult_n ? DBL_EPSILON : 
            disp < dcap_n ? fcap_n - kdesc_n * (disp - dcap_n) :
            disp < 0.0    ? (f0_n + k1_n * disp) * (1 - exp(-k0_n * disp / f0_n)) :
            disp < dcap_p ? (f0_p + k1_p * disp) * (1 - exp(-k0_p * disp / f0_p)) :
            disp < dult_p ? fcap_p - kdesc_p * (disp - dcap_p) :
                            DBL_EPSILON;
    } else if (envType == 2) {
        force = 
            disp < dult_n     ? DBL_EPSILON :
            disp < dcap_n     ? fcap_n - kdesc_n * (disp - dcap_n) :
            disp < 0.0        ? getBezierYK(0.0, denv1_n, denv2_n, dcap_n,
                                            0.0, fenv1_n, fenv2_n, fcap_n,
                                            disp, 0) :
            disp < dcap_p     ? getBezierYK(0.0, denv1_p, denv2_p, dcap_p,
                                            0.0, fenv1_p, fenv2_p, fcap_p,
                                            disp, 0) :
            disp <= dult_p     ? fcap_p - kdesc_p * (disp - dcap_p) :
                                DBL_EPSILON;
    } else if (envType == 3) {
        if (disp < denvs[0] || disp > denvs[envSize-1]) {
            force = DBL_EPSILON;
        } else {
            for (int i = 1; i < envSize; i++) {
                if (disp <= denvs[i]) {
                    force = (fenvs[i] - fenvs[i-1]) / (denvs[i] - denvs[i-1]) * 
                            (disp - denvs[i-1]) + fenvs[i-1];
                    break;
                }
            }
        }
    }
    return force;
}

double DowelType::denvelope(double disp)
{
    double tangent;
    if (envType == 1) {
        tangent = 
            disp < dult_n ? -kdesc_n :
            disp < dcap_n ? -kdesc_n :
            disp < 0.0    ? k1_n + ((k0_n - k1_n) + (disp * k0_n * k1_n) / f0_n) * 
                            exp(-k0_n * disp / f0_n) :
            disp < dcap_p ? k1_p + ((k0_p - k1_p) + (disp * k0_p * k1_p) / f0_p) * 
                            exp(-k0_p * disp / f0_p) :
            disp <= dult_p ? -kdesc_p :
                            -kdesc_p;
    } else if (envType == 2) {
        tangent = 
            disp < dult_n ? -kdesc_n :
            disp < dcap_n ? -kdesc_n :
            disp < 0.0    ? getBezierYK(0.0, denv1_n, denv2_n, dcap_n,
                                        0.0, fenv1_n, fenv2_n, fcap_n,
                                        disp, 0, false) :
            disp < dcap_p ? getBezierYK(0.0, denv1_p, denv2_p, dcap_p,
                                        0.0, fenv1_p, fenv2_p, fcap_p,
                                        disp, 0, false) :
            disp <= dult_p ? -kdesc_p :
                            -kdesc_p;
    } else if (envType == 3) {
        if (disp < denvs[0] || disp > denvs[envSize-1]) {
            return 0.0;
        } else {
            for (int i = 0; i < envSize; i++) {
                if (disp <= denvs[i]) {
                    tangent = (fenvs[i] - fenvs[i-1]) / (denvs[i] - denvs[i-1]);
                    break;
                }
            }
        }
    }
    return tangent;
}

// solve the intersection of the envelope and a line. Provide tangent and intercept.
double DowelType::envIntersection(double k, double b)
{
    double xm = 0.0;
    if (envType == 1 || envType == 2) {
        // use chord section method for curves.
        double x0 = 0.0;
        double offlen = dcap_p + dcap_n > 0 ? dcap_p / 50.0 : -dcap_n / 50.0;  // a small searching region 
        double offset = b > 0 ? offlen : -offlen;
        double x1 = offset;
        double y0 = envelope(x0) - (k * x0 + b);
        double y1 = envelope(x1) - (k * x1 + b);
        while (y0 * y1 > 0 && x1 <= dcap_p && x1 >= dcap_n) {
            if (DEBUG) opserr << "envIntersection: translate x0=" << x0 << " x1=" << x1 << " y0=" << y0 << " y1=" << y1 <<endln;
            y0 = y1;
            x0 = x1;
            x1 += offset;
            y1 = envelope(x1) - (k * x1 + b);
        }
        if (x1 > dcap_p || x1 < dcap_n) {
            x1 = b > 0 ? dcap_p : dcap_n;
        }
        if (fabs(y0) <= PRECISION) {
            return x0;
        } else if (fabs(y1) <= PRECISION) {
            return x1;
        } else if (y0 * y1 > 0) {
            opserr << "ERROR: Pinching path has no intersection with envelope. Check definition." << endln;
            tPath = 4;
            return 0.0;
        }
        double ym;
        int count = 0;
        while (count < MAX_ITER) {
            if (DEBUG) opserr << "envIntersection: iteration x0=" << x0 << " x1=" << x1 << " xm=" << xm << endln;
            xm = x0 - y0 * (x0 - x1) / (y0 - y1);
            ym = envelope(xm) - (k * xm + b);
            if (fabs(ym) < PRECISION || fabs(x1 - x0) < PRECISION) {
                break;
            } else if (ym * y0 < 0 && ym * y1 > 0) {
                x1 = xm;
                y1 = ym;
            } else {
                x0 = xm;
                y0 = ym;
            }
            count ++;
        }
        if (count >= MAX_ITER) {
            opserr << "ERROR: too many iterations when solving envelope and pinching line intersection. Check definitions." << endln;
            tPath = 4;
            return 0.0;
        }
    } else if (envType == 3) {
        double k1;
        if (b > 0) {
            for (int j = envZero + 1; j < envSize; j++) {
                if (fenvs[j] > k * denvs[j] + b) {
                    k1 = (fenvs[j] - fenvs[j-1]) / (denvs[j] - denvs[j-1]);
                    xm = (b - fenvs[j-1] + k1 * denvs[j-1]) / (k1 - k) ;
                    break;
                }
            }
        } else if (b < 0) {
            for (int j = envZero - 1; j >= 0; j--) {
                if (fenvs[j] < k * denvs[j] + b) {
                    k1 = (fenvs[j] - fenvs[j+1]) / (denvs[j] - denvs[j+1]);
                    xm = (b - fenvs[j+1] + k1 * denvs[j+1]) / (k1 - k) ;
                    break;
                }
            }
        }
    }
    return xm;
}

// get a point on the envelope with a specific tangent value.
double DowelType::envWithSlope(double k, bool flag, double xinit)
{
    if (denvelope(xinit) < k) {
        return xinit;
    }
    double xm = xinit;
    if (envType == 1 || envType == 2) {
        double offlen = dcap_p + dcap_n > 0 ? dcap_p / 50.0 : -dcap_n / 50.0;
        double offset = flag ?  -offlen : offlen ;
        double x0 = xinit;
        while (denvelope(x0) > k && x0 > dcap_n && x0 < dcap_p) {
            x0 += offset;
            if (DEBUG) opserr << "Finding slope x0=" << x0 << endln;
        }
        double x1 = x0 - offset;        
        double y0 = denvelope(x0) - k;
        double y1 = denvelope(x1) - k;
        if (fabs(y0) <= PRECISION) {
            return x0;
        } else if (fabs(y1) <= PRECISION) {
            return x1;
        } else if (x0 > dcap_n || x0 < dcap_p) {
            x0 = flag ? dcap_n : dcap_p; 
            if (DEBUG) opserr << "envWithSlope: reloading stiffness is too small." << endln;
            return x0;
        }
        double ym;
        int count = 0;
        while (count < MAX_ITER) {
            if (DEBUG) opserr << "envWithSlope: iteration x0=" << x0 << " x1=" << x1 << " xm=" << xm << endln;
            xm = x0 - y0 * (x0 - x1) / (y0 - y1);
            ym = envelope(xm) - k;
            if (fabs(ym) < PRECISION || fabs(x1 - x0) < PRECISION) {
                break;
            } else if (ym * y0 < 0 && ym * y1 > 0) {
                x1 = xm;
                y1 = ym;
            } else {
                x0 = xm;
                y0 = ym;
            }
            count ++;
        }
        if (count == MAX_ITER) {
            opserr << "WARNING: too many iterations when solving envelope point with a specific slope. Check the definition." << endln;
        }
    } else if (envType == 3) {
        int target = envZero;
        // double kenv;
        if (flag) {
            while (target > 0 && fenvs[target] > envelope(xinit)) {
                target --;
            }
        } else {
            while (target < envSize && fenvs[target] < envelope(xinit)) {
                target ++;
            }
        }
        xm = denvs[target];
    }
    if (DEBUG) opserr << "EnvWithSlope returned " << xm << endln;
    return xm;
}

// Given an x, get the y and k on the reverse pinched path.
void DowelType::getReverseYK(double x, bool flag, double *y, double *k)
{
    int dt = flag ? 0 : 10;  // index offset for Bezier points
    if ((x > pxs[0+dt] && x < pxs[1+dt]) || (x < pxs[0+dt] && x > pxs[1+dt])) 
    {
        *k = (pys[0+dt] - pys[1+dt]) / (pxs[0+dt] - pxs[1+dt]);
        *y = *k * (x - pxs[0+dt]) + pys[0+dt];
    } 
    else if ((x >= pxs[1+dt] && x <= pxs[4+dt]) || (x <= pxs[1+dt] && x >= pxs[4+dt])) 
    {
        *y = getBezierYK(pxs[1+dt], pxs[2+dt], pxs[3+dt], pxs[4+dt],
                         pys[1+dt], pys[2+dt], pys[3+dt], pys[4+dt],
                         x, k);
    } 
    else if ((x > pxs[4+dt] && x < pxs[5+dt]) || (x < pxs[4+dt] && x > pxs[5+dt])) 
    {
        *k = (pys[4+dt] - pys[5+dt]) / (pxs[4+dt] - pxs[5+dt]);
        *y = *k * (x - pxs[4+dt]) + pys[4+dt];
    } 
    else if ((x >= pxs[5+dt] && x <= pxs[8+dt]) || (x <= pxs[5+dt] && x >= pxs[8+dt])) 
    {
        *y = getBezierYK(pxs[5+dt], pxs[6+dt], pxs[7+dt], pxs[8+dt],
                         pys[5+dt], pys[6+dt], pys[7+dt], pys[8+dt],
                         x, k);
    } 
    else if ((x > pxs[8+dt] && x < pxs[9+dt]) || (x < pxs[8+dt] && x > pxs[9+dt]))
    {
        *k = (pys[8+dt] - pys[9+dt]) / (pxs[8+dt] - pxs[9+dt]);
        *y = *k * (x - pxs[8+dt]) + pys[8+dt];
    }
    else
    {
        opserr << "ERROR: x is not on the pinched curve. x=" << x << " Controlling points"<< endln;
        opserr << pxs[0+dt] << " " << pxs[1+dt] << " " << pxs[2+dt] << " " << pxs[3+dt] << " " << pxs[4+dt] ;
        opserr << pxs[5+dt] << " " << pxs[6+dt] << " " << pxs[7+dt] << " " << pxs[8+dt] << " " << pxs[9+dt] <<  endln;
        *y = 0.0;
        *k = 0.0;
        tPath = cPath = 4;
    }
}

// reset the reverse path control points when loading history is changed.
void DowelType::resetReversePoints(double disp, double force, bool flag)
{
    if (DEBUG) opserr << "------> resetReversePoints" << endln;
    // flag: go right to left, true. go left to right, false.
    double dpeak = flag ? cDmin : cDmax; // peak disp on the opposite side
    double dmaxs = flag ? cDmax : cDmin; // maximum disp of the same side
    double Fdmaxs = flag ? cFdmax : cFdmin; // force corresponding to the maximum disp of the same side
    double dmax0 = flag ? cDmax : -cDmin; // absolute maximum disp of the same side
    double dmax1 = flag ? -cDmin : cDmax; // absolute maximum disp of the other side 
    double dmax = cDmax > -cDmin ? cDmax : -cDmin; // absolute maximum disp
    double fb = flag ? -fi : fi;  // intercept with y axis
    double k0 = flag ? k0_p : k0_n; // initial stiffness, same side.
    double k01 = flag ? k0_n : k0_p; // initial stiffness, other side.
    double dmg = flag ? (ePHC_n + eFHC_n) / (eMono_n + eFHC_n) : 
                        (ePHC_p + eFHC_p) / (eMono_p + eFHC_p); // damage index

    double ax, ay, ak;
    double bx, by;
    double cx, cy;
    double dx, dy, dk;
    double mx, my, mk;
    double nx, ny;
    // a, b, c, d: four controlling points of the unloading-pinching-reloading curve
    // m: the force intercept of the pinching path
    // n: the separation point (mid-point) of the pinching curve.
    // ak: unloading stiffness
    // dk: reloading stiffness
    // mk: pinching stiffness

    // A: unloading occurring point
    ax = disp;
    ay = force;
    ak = alpha_u >= 0 ? dmax0 <= dyield ? ru * k0 : ru * k0 * pow(dyield / dmax0, alpha_u) :
                        dmax0 < DBL_EPSILON ? ru * k0 : ru * k0 * pow(Fdmaxs / dmaxs / k0, -alpha_u);
    if (DEBUG) opserr << "ax=" << ax << " ay=" << ay << " ak=" << ak << endln;
    // M: pinch line controlling point: force intercept.
    mx = 0.0;
    my = flag ? dmax0 > dyield ? fabs(Fdmaxs) > fabs(fyield_p) ? fb - eta * (Fdmaxs - fyield_p) : fb :
                                 (dmax0 / dyield) * fb :                   
                dmax0 > dyield ? fabs(Fdmaxs) > fabs(fyield_n) ? fb - eta * (Fdmaxs - fyield_n) : fb :
                                 (dmax0 / dyield) * fb;
    mk = alpha_p >= 0 ? dmax <= dyield ? kp :  kp * pow(dyield / dmax, alpha_p):
                       dmax0 < DBL_EPSILON ? kp :  kp * pow(Fdmaxs / dmaxs / k0, -alpha_p);
    if (DEBUG) opserr << "mx=" << mx << " my=" << my << " mk=" << mk << endln;
    // D: reloading line and envelope intersection
    double energyAmp = pow(gamma, dmg);
    if (DEBUG) opserr << "Energy:" << "ePHC_n=" << ePHC_n << " eFHC_n=" << eFHC_n << " ePHC_p=" << ePHC_p << " eFHC_p=" << eFHC_p << " eMono=" << eMono_n << endln;
    if (DEBUG) opserr << "dmg=" << dmg <<", energyAmp=" << energyAmp << endln;
    dx = dpeak * beta * energyAmp;
    dy = envelope(dx);
    dk = alpha_r >= 0 ? dmax1 <= dyield ? k01 : k01 * pow(dyield / dmax1, alpha_r):
                        dmax0 < DBL_EPSILON ? k01 : k01 * pow(Fdmaxs / dmaxs / k0, -alpha_r);
    if (DEBUG) opserr << "dx=" << dx << " dy="<<dy<< " dk=" << dk << endln;
    // solve B from A and M (intersection)
    bx = ((my - mk * mx) - (ay - ak * ax)) / (ak - mk);
    by = mk * (bx - mx) + my;
    // solve C from D and M (intersection)
    cx = ((my - mk * mx) - (dy - dk * dx)) / (dk - mk);
    cy = mk * (cx - mx) + my;
    if (DEBUG) opserr << "bx=" << bx << " by="<<by<< " cx=" << cx << " cy="<< cy << endln;
    // solve N from B and C (mid point)
    nx = (bx + cx) / 2;
    ny = (by + cy) / 2;
    // special scenario 1
    bool pinching2env = false;
    // No reloading part: goes from pinching to envelope.
    if ((flag && cy <= dy && dx >= envIntersection(mk, my)) || (!flag && cy >= dy && dx <= envIntersection(mk, my))) {
        if (DEBUG) opserr << "pinching scenario 1: no reloading part, goes from pinching to envelope." << endln;
        double delta_cd = cx - dx;
        dx = envIntersection(mk, my);
        dy = envelope(dx);
        dk = mk;
        cx = 0.;
        cy = my;
        if ((flag && disp > dyield) || (!flag && disp < -dyield)) {
            double x = envIntersection(mk, my) + delta_cd;
            dx = x;
            dy = envelope(dx);
            dk = denvelope(dx);
            cx = 0.;
            cy = my;
            pinching2env = true;
            if (DEBUG) opserr << " disp=" << disp << " dx=" << x << " my=" << my << endln;
        }
        nx = (bx + cx) / 2;
        ny = (by + cy) / 2;
    }
    // special scenario 2
    // goes from pinching to the point on the descending part.
    if ((flag && cy <= dy && dx < envIntersection(mk, my)) || (!flag && cy >= dy && dx > envIntersection(mk, my))) {
        if (DEBUG) opserr << "pinching scenario 2: goes from pinching to the point on the descending part" << endln;
        // move pinching line go from unloading end point to target d.
        cx = dx;
        cy = dy;
        // eliminate negative stiffness.
        if ((flag && by < cy) || (!flag && by > cy)) {
            bx = (cy - ay) / ak + ax;
            by = cy;
        }
        nx = (cx + bx) / 2;
        ny = (cy + by) / 2;
        mx = nx;
        my = ny;
        mk = (cy - by) / (cx - bx);
    }
    // special scenario 3
    // BC goes in wrong direction, ignore pinching part. 
    // From unloading to reloading directly.
    if ((flag && bx < cx && bx <= ax) || (!flag && bx > cx && bx >= ax)) {
        if (DEBUG) opserr << "pinching scenario 3: from unloading to reloading directly" << endln;
        // move m, n, b to a
        nx = ax;
        ny = ay;
        mx = ax;
        my = ay;
        mk = ak;
        cx = ((ay - ak * ax) - (dy - dk * dx)) / (dk - ak);
        cy = dk * (cx - dx) + dy;
        bx = ax;
        by = ay;
    }
    // special scenario 4
    // AB goes in wrong direction
    // go directly towards the reloading point C.
    if ((flag && bx > ax) || (!flag && bx < ax)) {
        if (DEBUG) opserr << "pinching scenario 4: go directly towards the reloading point c" << endln;
        bx = ax;
        by = ay;
        nx = ax;
        ny = ay;
        // AC goes in wrong direction
        // go directly towards the target point D.
        if ((flag && (cx > ax)) || (!flag && (cx < ax))) {
            cx = ax;
            cy = ay;
        } 
        // ignore unloading part.
        if ((flag && (cx < ax) && (cy > ay)) || (!flag && (cx > ax) && (cy < ay))) {
            cx = ((ay - mk * ax) - (dy - dk * dx)) / (dk - mk);
            cy = dk * (cx - dx) + dy;
            if ((flag && (cx < dx)) || (!flag && (cx > dx))) {
                cx = ax;
                cy = ay;
            }
        }
    }
    // assign the control points
    int dt = flag ? 0 : 10;
    if (c <= 1.0) {
        pxs[0+dt] = ax;
        pxs[1+dt] = c * (ax - bx) + bx;
        pxs[2+dt] = bx;
        pxs[3+dt] = bx;
        pxs[4+dt] = c * (nx - bx) + bx;
        pxs[5+dt] = c * (nx - cx) + cx;
        pxs[6+dt] = cx;
        pxs[7+dt] = cx;
        pxs[8+dt] = c * (dx - cx) + cx;
        pxs[9+dt] = dx;
        pys[0+dt] = ay;
        pys[1+dt] = c * (ay - by) + by;
        pys[2+dt] = by;
        pys[3+dt] = by;
        pys[4+dt] = c * (ny - by) + by;
        pys[5+dt] = c * (ny - cy) + cy;
        pys[6+dt] = cy;
        pys[7+dt] = cy;
        pys[8+dt] = c * (dy - cy) + cy;
        pys[9+dt] = dy;
    } else if (c < 2.0) {
        pxs[0+dt] = ax;
        pxs[1+dt] = ax;
        pxs[2+dt] = (c - 1) * (ax - bx) + bx;
        pxs[3+dt] = (c - 1) * (nx - bx) + bx;
        pxs[4+dt] = nx;
        pxs[5+dt] = nx;
        pxs[6+dt] = (c - 1) * (nx - cx) + cx;
        pxs[7+dt] = (c - 1) * (dx - cx) + cx;
        pxs[8+dt] = dx;
        pxs[9+dt] = dx;
        pys[0+dt] = ay;
        pys[1+dt] = ay;
        pys[2+dt] = (c - 1) * (ay - by) + by;
        pys[3+dt] = (c - 1) * (ny - by) + by;
        pys[4+dt] = ny;
        pys[5+dt] = ny;
        pys[6+dt] = (c - 1) * (ny - cy) + cy;
        pys[7+dt] = (c - 1) * (dy - cy) + cy;
        pys[8+dt] = dy;
        pys[9+dt] = dy;
    }

    if (pinching2env) {
        if (DEBUG) opserr << "Reloading special scenario: pinching to envelope." << endln;
        pxs[5+dt] = nx;
        pxs[6+dt] = cx;
        pxs[7+dt] = 0;
        pxs[8+dt] = dx;
        pys[5+dt] = ny;
        pys[6+dt] = cy;
        double intercept = dk * -dx + dy;
        if ((flag && intercept > cy) || (!flag && intercept < cy)) {
            intercept = cy;
        }
        pys[7+dt] = intercept;
        pys[8+dt] = dy;
    }
    if (DEBUG) opserr << "resetReversePoints:" << endln;
    if (DEBUG) for (int i = 0; i < 20; i++) {
        opserr << "No." << i <<" px=" << pxs[i] << " py=" << pys[i] << endln;
    }
}

int DowelType::setTrialStrain(double strain, double strainRate)
{
    if (DEBUG) opserr << "------> Set trail strain = " << strain << ", in ";
    if ((fabs(tStrain - strain) < DBL_EPSILON) && (fabs(tStrain) > DBL_EPSILON))
    {
        if (DEBUG)  opserr << "None"<< endln;
        tStrain = strain;
        return 0;
    }
    tStrain = strain;
    tDmin = (cDmin > tStrain) ? tStrain : cDmin;
    tDmax = (cDmax < tStrain) ? tStrain : cDmax;
    if (cPath == 4 || tStrain > dult_p || tStrain < dult_n) 
    {
        // material fail.
        if (DEBUG) opserr << "scenario 4" << endln;
        tStress = DBL_EPSILON;
        tTangent = DBL_EPSILON;
        tPath = 4;
    } 
    else if (cPath == 1)
    {
        // Current committed point on envelope.
        if ((tStrain >= cStrain && tStrain >= pxs[19]) || (tStrain <= cStrain && tStrain <= pxs[9]))
        {
            // scenario 11: loading on backbone curve
            if (DEBUG) opserr << "scenario 11" << endln;
            tStress = envelope(tStrain);
            tTangent = denvelope(tStrain);
            tPath = 1;
        }
        else if (tStrain < cStrain)
        {
            // scenario 12: unloading to pinched curve right->left
            if (DEBUG) opserr << "scenario 12" << endln;
            getReverseYK(tStrain, true, &tStress, &tTangent);
            tPath = 2;
        }
        else if (tStrain > cStrain)
        {
            // path scenario 13: unloading to pinched curve, left->right
            if (DEBUG) opserr << "scenario 13" << endln;
            getReverseYK(tStrain, false, &tStress, &tTangent);
            tPath = 3;
        }
    } 
    else if (cPath == 2) 
    {
        // current point on pinched curve right->left
        if (tStrain <= cStrain && tStrain > pxs[9])
        {
            // path scenario 22: go on with pinched curve right->left
            if (DEBUG)  opserr << "scenario 22 "<< endln;
            getReverseYK(tStrain, true, &tStress, &tTangent);
            tPath = 2;
        }
        else if ((tStrain <= cStrain && tStrain <= pxs[9]) || (tStrain > cStrain && tStrain > pxs[19]))
        {
            // path scenario 21: from pinched curve to envelope
            if (DEBUG) opserr << "scenario 21 "<< endln;
            tStress = envelope(tStrain);
            tTangent = denvelope(tStrain);
            tPath = 1;
        }
        else if (tStrain > cStrain)
        {
            // path scenario 23: pinched from right->left to left->right
            if (DEBUG) opserr << "scenario 23 "<< endln;
            getReverseYK(tStrain, false, &tStress, &tTangent);
            tPath = 3;
        }
    } 
    else if (cPath == 3) 
    {
        // current point on pinched curve left->right
        if (tStrain >= cStrain && tStrain < pxs[19])
        {
            // path scenario 33: go on with pinched curve left->right
            if (DEBUG) opserr << "scenario 33 "<< endln;
            getReverseYK(tStrain, false, &tStress, &tTangent);
            tPath = 3;
        }
        else if ((tStrain >= cStrain && tStrain >= pxs[19]) || (tStrain < cStrain && tStrain < pxs[9]))
        {
            // path scenario 31: from pinched curve to envelope
            if (DEBUG) opserr << "scenario 31 "<< endln;
            tStress = envelope(tStrain);
            tTangent = denvelope(tStrain);
            tPath = 1;
        }
        else if (tStrain < cStrain)
        {
            // path scenario 32: pinched curve left->right to right->left
            if (DEBUG) opserr << "scenario 32 "<< endln;
            getReverseYK(tStrain, true, &tStress, &tTangent);
            tPath = 2;
        }
    } 
    return 0;
}

double DowelType::getStrain(void)
{
    return tStrain;
}

double DowelType::getStress(void)
{
    return tStress;
}

double DowelType::getTangent(void)
{
    return tTangent;
}

double DowelType::getInitialTangent(void)
{
    return k0_p;
}


int DowelType::commitState(void)
{
    double energy = 0.5 * (tStrain - cStrain) * (tStress + cStress);
    cStrain  = tStrain;
    cTangent = tTangent;
    cStress  = tStress;
    cPath = tPath;
    if (DEBUG) opserr << "------> Commit State: " << "cPath=" << cPath << ", cStrain=" << cStrain << ", cStress=" << cStress;
    if (DEBUG) opserr << ", cTangent=" << cTangent << ", cDmin=" << cDmin << ", cDmax=" << cDmax << endln;
    if (cStrain >= cDmax) {
        isPHC = true;
        cDmax = cStrain;
        cFdmax = cStress;
    } 
    if (cStrain <= cDmin) {
        isPHC = true;
        cDmin = cStrain;
        cFdmin = cStress;
    } 
    if (isPHC) {
        cStrain > 0 ? ePHC_p += energy : ePHC_n += energy;
    } else {
        cStrain > 0 ? eFHC_p += energy : eFHC_n += energy;
    }
    if (isPHC && ((cStrain > 0.0 && cStress < 0.0) || (cStrain < 0.0 && cStress > 0.0))) {
        isPHC = false;
    }
    if ((cPath == 1 && cStrain >= 0.0) || cPath == 3) {
        resetReversePoints(cStrain, cStress, true);
    } else if ((cPath == 1 && cStrain < 0.0) || cPath == 2) {
        resetReversePoints(cStrain, cStress, false);
    } 
    return 0;
}


int DowelType::revertToLastCommit(void)
{
    tStrain  = cStrain;
    tTangent = cTangent;
    tStress  = cStress;
    tDmin    = cDmin;
    tFdmin   = cFdmin;
    tDmax    = cDmax;
    tFdmax   = cFdmax;
    tPath    = cPath;
    return 0;
}


int DowelType::revertToStart(void)
{
    tStrain  = cStrain = 0.0;
    tTangent = cTangent = k0_p;
    tStress  = cStress = 0.0;
    tDmin = cDmin = 0.0;
    tFdmin = cFdmin = 0.0;
    tDmax = cDmax = 0.0;
    tFdmax = cFdmax = 0.0;
    tPath = cPath = 1;

    isPHC = true;
    ePHC_p = ePHC_n = 0.0;
    eFHC_p = eFHC_n = 0.0;

    for (int i=0; i < 20; i++) {
      pxs[i] = 0.0;
      pys[i] = 0.0;
    }
    
    return 0;
}


UniaxialMaterial *DowelType::getCopy(void)
{
    DowelType *theCopy = new DowelType();
    theCopy->setTag(this->getTag());
    theCopy->fi = fi;
    theCopy->kp = kp;
    theCopy->ru = ru;
    theCopy->c = c;
    theCopy->beta = beta;
    theCopy->gamma = gamma;
    theCopy->eta = eta;
    theCopy->dyield = dyield; 
    theCopy->alpha_p = alpha_p;
    theCopy->alpha_u = alpha_u;
    theCopy->alpha_r = alpha_r;
    theCopy->envType = envType;
    theCopy->k0_p = k0_p;
    theCopy->k0_n = k0_n;
    theCopy->dcap_p = dcap_p;
    theCopy->dcap_n = dcap_n;
    theCopy->fcap_p = fcap_p;
    theCopy->fcap_n = fcap_n;
    theCopy->fyield_p = fyield_p;
    theCopy->fyield_n = fyield_n;
    theCopy->dult_p = dult_p;
    theCopy->dult_n = dult_n;
    theCopy->kdesc_p = kdesc_p;
    theCopy->kdesc_n = kdesc_n;
    theCopy->dinter_p = dinter_p;
    theCopy->dinter_n = dinter_n;
    theCopy->eMono_p = eMono_p;
    theCopy->eMono_n = eMono_n;
    theCopy->k1_p = k1_p;
    theCopy->k1_n = k1_n;
    theCopy->f0_p = f0_p;
    theCopy->f0_n = f0_n;
    theCopy->denv1_p = denv1_p;
    theCopy->denv1_n = denv1_n;
    theCopy->fenv1_p = fenv1_p;
    theCopy->fenv1_n = fenv1_n;
    theCopy->denv2_p = denv2_p;
    theCopy->denv2_n = denv2_n;
    theCopy->fenv2_p = fenv2_p;
    theCopy->fenv2_n = fenv2_n;
    theCopy->envSize = envSize;
    theCopy->envZero = envZero;
    if (!denvs) {
        theCopy->denvs = NULL;
    } else {
        double * denvsCopy = new double[envSize];
        for (int i = 0; i < envSize; i++) {
            denvsCopy[i] = denvs[i];
        }
        theCopy->denvs = denvsCopy;
    }
    if (!fenvs) {
        theCopy->fenvs = NULL;
    } else {
        double * fenvsCopy = new double[envSize];
        for (int i = 0; i < envSize; i++) {
            fenvsCopy[i] = fenvs[i];
        }
        theCopy->fenvs = fenvsCopy;
    }
    theCopy->isPHC = isPHC;
    theCopy->ePHC_p = ePHC_p;
    theCopy->ePHC_n = ePHC_n;
    theCopy->eFHC_p = eFHC_p;
    theCopy->eFHC_n = eFHC_n;
    for (int i=0; i < 20; i++) {
      theCopy->pxs[i] = pxs[i];
      theCopy->pys[i] = pys[i];
    }
    return theCopy;
}


int DowelType::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    int dbTag = this->getDbTag();
    
    ID idata(3);
    idata(0) = this->getTag();
    idata(1) = envSize;
    idata(2) = envZero;    
    res = theChannel.sendID(dbTag, cTag, idata);
    if (res < 0) {
      opserr << "DowelType::sendSelf() - failed to send ID data" << endln;
      return res;
    }
    
    Vector data(98 + 2 * envSize);
    //data(0)  = this->getTag();
    data(1)  = fi;
    data(2)  = kp;
    data(3)  = ru;
    data(4)  = c;
    data(5)  = beta;
    data(6)  = gamma;
    data(7)  = eta;
    data(8)  = dyield;
    data(9)  = alpha_p;
    data(10) = alpha_u;
    data(11) = alpha_r;
    
    data(12) = envType;
    data(13) = k0_p;
    data(14) = k0_n;
    data(15) = dcap_p;
    data(16) = dcap_n;
    data(17) = fcap_p;
    data(18) = fcap_n;
    data(19) = fyield_p;
    data(20) = fyield_n;
    data(21) = dult_p;
    data(22) = dult_n;
    data(23) = kdesc_p;
    data(24) = kdesc_n;

    data(25) = dinter_p;
    data(26) = dinter_n;
    data(27) = eMono_p;
    data(28) = eMono_n;

    data(29) = k1_p;
    data(30) = k1_n;
    data(31) = f0_p;
    data(32) = f0_n;
    data(33) = denv1_p;
    data(34) = denv1_n;
    data(35) = fenv1_p;
    data(36) = fenv1_n;
    data(37) = denv2_p;
    data(38) = denv2_n;
    data(39) = fenv2_p;
    data(40) = fenv2_n;

    //data(41) = envSize;
    //data(42) = envZero;

    data(43) = isPHC ? 1.0 : -1.0;
    data(44) = ePHC_p;
    data(45) = ePHC_n;
    data(46) = eFHC_p;
    data(47) = eFHC_n;

    for (int i = 0; i < 20; i++) {
      data(48+i) = pxs[i];
      data(68+i) = pys[i];
    }    

    data(90) = cStrain;
    data(91) = cStress;
    data(92) = cTangent;
    data(93) = cPath;   
    data(94) = cDmin;
    data(95) = cDmax;
    data(96) = cFdmin;
    data(97) = cFdmax;
    for (int i = 0; i < envSize; i++) {
        data(98 + i * 2) = denvs[i];
        data(98 + i * 2 + 1) = fenvs[i];
    }
    res = theChannel.sendVector(dbTag, cTag, data);
    if (res < 0) {
      opserr << "DowelType::sendSelf() - failed to send data" << endln;
      return res;
    }
    return 0;
}

int DowelType::recvSelf(int cTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dbTag = this->getDbTag();
    
    ID idata(3);
    if (theChannel.recvID(dbTag, cTag, idata) < 0) {
      opserr << "DowelType::recvSelf() - failed to recv ID data" << endln;
      return -1;
    }

    this->setTag(idata(0));
    envSize = idata(1);
    envZero = idata(2);
    
    Vector data(98 + 2 * envSize);
    res = theChannel.recvVector(dbTag, cTag, data);
    if (res < 0) {
      opserr << "DowelType::recvSelf() - failed to recv data" << endln;
      return res;
    }
    else {
      //this->setTag(data(0));
        fi = data(1);
        kp = data(2);
        ru = data(3);
        c = data(4);
        beta = data(5);
        gamma = data(6);
        eta = data(7);
        dyield = data(8);
        alpha_p = data(9);
        alpha_u = data(10);
        alpha_r = data(11);
	
        envType = data(12);
        k0_p = data(13);
        k0_n = data(14);
        dcap_p = data(15);
        dcap_n = data(16);
        fcap_p = data(17);
        fcap_n = data(18);
        fyield_p = data(19);
        fyield_n = data(20);
        dult_p = data(21);
        dult_n = data(22);
        kdesc_p = data(23);
        kdesc_n = data(24);

	dinter_p = data(25);
        dinter_n = data(26);
        eMono_p = data(27);
        eMono_n = data(28);

	k1_p = data(29);
        k1_n = data(30);
        f0_p = data(31);
        f0_n = data(32);

	denv1_p = data(33);
        denv1_n = data(34);
        fenv1_p = data(35);
        fenv1_n = data(36);
        denv2_p = data(37);
        denv2_n = data(38);
        fenv2_p = data(39);
        fenv2_n = data(40);

	//envSize = data(41);
        //envZero = data(42);

	isPHC = (data(43) > 0.0) ? true : false;
        ePHC_p = data(44);
        ePHC_n = data(45);
        eFHC_p = data(46);
        eFHC_n = data(47);

	for (int i = 0; i < 20; i++) {
	  pxs[i] = data(48+i);
	  pys[i] = data(68+i);
	}

	cStrain = data(90);
        cStress = data(91);
        cTangent = data(92);
        cPath = data(93);   
        cDmin = data(94);
        cDmax = data(95);
	cFdmin = data(96);
	cFdmax = data(97);	
        double * denvsCopy = new double[envSize];
        double * fenvsCopy = new double[envSize];
        for (int i = 0; i < envSize; i++) {
            denvsCopy[i] = data(98 + i * 2);
            fenvsCopy[i] = data(98 + i * 2 + 1);
        }
	if (denvs != 0) delete [] denvs;
	if (fenvs != 0) delete [] fenvs;	
        denvs = denvsCopy;
        fenvs = fenvsCopy;
        this->revertToLastCommit();
    }
    return res;
}

void DowelType::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "DowelType tag: " << this->getTag() << endln;
        s << "Fi=" << fi << ", Kp=" << kp << ", Ru=" << ru << ", c=" << c << endln;
        s << "beta=" << beta << ", gamma=" << gamma << ", eta=" << eta << endln;
        s << "Dy=" << dyield << ", alpha_p=" << alpha_p;
        s << ", alpha_u=" << alpha_u << ", alpha_r" << alpha_r << endln;
        if (envType == 1) {
            s << "Envelope type : exponential. " << endln;
            s << "K0=" << k0_p << ", R1=" << k1_p/k0_p << ", F0=" << f0_p;
            s << ", Dc=" << dcap_p << ", Kd=" << kdesc_p << ", Du=" << dult_p << endln;
            s << "K0N=" << k0_n << ", R1N=" << k1_n/k0_n << ", F0N=" << f0_n;
            s << ", DcN=" << dcap_n << ", KdN=" << kdesc_n << ", DuN=" << dult_n << endln;
        } else if (envType == 2) {
            s << "Envelope type : Bezier. " << endln;
            s << "D1=" << denv1_p << ", F1=" << fenv1_p << ", D2=" << denv2_p << ", F2=" << fenv2_p ;
            s << ", Dc=" << dcap_p << ", Fc=" << fcap_p << ", Kd=" << kdesc_p << ", Du=" << dult_p << endln;
            s << "D1N=" << denv1_n << ", F1N=" << fenv1_n << ", D2N=" << denv2_n << ", F2N=" << fenv2_n ;
            s << ", DcN=" << dcap_n << ", FcN=" << fcap_n << ", KdN=" << kdesc_n << ", DuN=" << dult_n << endln;
        } else if (envType == 3) {
            s << "Envelope type : Piecewise. " << endln;
            for (int i = 0; i < envSize; i++) {
                s << "D" << i << "=" << denvs[i] << ", F" << i << "=" << fenvs[i] << endln;
            }
        }
    }
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{" << endln;
        s << "\t\t\t\t\"tag\": " << this->getTag() << "," << endln;
        s << "\t\t\t\t\"Fi\": " << fi << "," << endln;
        s << "\t\t\t\t\"Kp\": " << kp << "," << endln;
        s << "\t\t\t\t\"Ru\": " << ru << "," << endln;
        s << "\t\t\t\t\"c\": " << c << "," << endln;
        s << "\t\t\t\t\"beta\": " << beta << "," << endln;
        s << "\t\t\t\t\"gamma\": " << gamma << "," << endln;
        s << "\t\t\t\t\"eta\": " << eta << "," << endln;
        s << "\t\t\t\t\"Dy\": " << dyield << "," << endln;
        s << "\t\t\t\t\"alpha_p\": " << alpha_p << "," << endln;
        s << "\t\t\t\t\"alpha_u\": " << alpha_u << "," << endln;
        s << "\t\t\t\t\"alpha_r\": " << alpha_r << "," << endln;
        if (envType == 1) {
            s << "\t\t\t\t\"envelope\": \"exponential\"," << endln;
            s << "\t\t\t\t\"envelope_params\": [" << k0_p << ", " << k1_p/k0_p << ", " << f0_p;
            s << ", " << dcap_p << ", " << kdesc_p << ", " << dult_p << ", " << k0_n;
            s << ", " << k1_n/k0_n << ", " << f0_n << ", " << dcap_n << ", " << kdesc_n;
            s << ", " << dult_n << "]," << endln;
        } else if (envType == 2) {
            s << "\t\t\t\t\"envelope\": \"bezier\"," << endln;
            s << "\t\t\t\t\"envelope_params\": [" << denv1_p << ", " << fenv1_p << ", " << denv2_p;
            s << ", " << fenv2_p << ", " << dcap_p << ", " << fcap_p << ", " << kdesc_p << ", " << dult_p;
            s << ", " << denv1_n << ", " << fenv1_n << ", " << denv2_n << ", " << fenv2_n ;
            s << ", " << dcap_n << ", " << fcap_n << ", " << kdesc_n << ", " << dult_n << "]," endln;
        } else if (envType == 3) {
            s << "\"envelope\": \"-piecewise\"," << endln;
            s << "\"envelope_params\": [";
            for (int i = 0; i < envSize; i++) {
                s << denvs[i] << ", " << fenvs[i];
                if (i != envSize-1) {
                    s << ", ";
                } else {
                    s << "]," << endln;
                }
            }
        }
        s << "\t\t\t}" << endln;
    }
}


double DowelType::getBezierYK(double x1, double x2, double x3, double x4,
                              double y1, double y2, double y3, double y4,
                              double x, double *other, bool returnY)
{
    double PA = -1 * x1 + 3 * x2 - 3 * x3 + x4;
    double PB =  3 * x1 - 6 * x2 + 3 * x3;
    double PC = -3 * x1 + 3 * x2;
    double PD =      x1 - x;
    double t = -1.0;

    if (fabs(x - x1) < PRECISION) {
        t = 0.;
    } else if (fabs(x - x4) < PRECISION) {
        t = 1.;
    } else if (fabs(PA) > DBL_EPSILON) {
        double A = PB / PA;
        double B = PC / PA;
        double C = PD / PA;
        double Q = (3 * B - pow(A, 2)) / 9;
        double R = (9 * A * B - 27 * C - 2 * pow(A, 3)) / 54;
        double D = pow(Q, 3) + pow(R, 2);
        double t1, t2, t3;
        if (D >= 0) {
            double S, T;
            S = (R + sqrt(D) > 0.0) ?
                pow(fabs(R + sqrt(D)), (1.0 / 3)) :
                -1 * pow(fabs(R + sqrt(D)), (1.0 / 3));
            T = (R - sqrt(D) > 0.0) ?
                pow(fabs(R - sqrt(D)), (1.0 / 3)) :
                -1 * pow(fabs(R - sqrt(D)), (1.0 / 3));
            t1 = -1 * A / 3 + S + T;
            if (t1 >= 0.0 && t1 <= 1.0)
                t = t1;
            else if (S == T)
            {
                t2 = -1 * A / 3 - (S + T) / 2;
                if (t2 >= 0.0 && t2 <= 1.0)
                    t = t2;
            }
        } else {
            double th = acos(R / sqrt(-pow(Q, 3)));
            t1 = 2 * sqrt(-Q) * cos(th/3) - A/3;
            if (t1 >= 0.0 && t1 <= 1.0) {
                t = t1;
            } else {
                t2 = 2 * sqrt(-Q) * cos((th + 2*PI)/3) - A/3;
                if (t2 >= 0.0 && t2 <= 1.0) {
                    t = t2;
                } else {
                    t3 = 2 * sqrt(-Q) * cos((th + 4*PI)/3) - A/3;
                    if (t3 >= 0.0 && t3 <= 1.0) {
                        t = t3;
                    }
                }
            }
        }
    } else if (fabs(PB) > DBL_EPSILON) {
        double delta = pow(PC, 2) - 4 * PB * PD;
        if (delta >= 0) {
            double t1 = (-PC + pow(delta, 0.5)) / (2*PB);
            double t2 = (-PC - pow(delta, 0.5)) / (2*PB);
            t = (t1 > 0.0 && t1 < 1.0) ? t1 : t2;
        }
    } else {
        t = -PD / PC;
    }
    if (t < 0 || t > 1) {
        opserr << "ERROR: t is not in [0, 1]" << endln;
        opserr << "xs=" << x1 << " " << x2 << " " << x3 << " " << x4 << endln;
        opserr << "ys=" << y1 << " " << y2 << " " << y3 << " " << y4 << endln;
        opserr << "x=" << x << " t=" << t << endln;
    }
    double dx = (-3 * x1 + 9 * x2 - 9 * x3 + 3 * x4) * t * t + (6 * x1 - 12 * x2 + 6 * x3) * t + (-3 * x1 + 3* x2);
    double dy = (-3 * y1 + 9 * y2 - 9 * y3 + 3 * y4) * t * t + (6 * y1 - 12 * y2 + 6 * y3) * t + (-3 * y1 + 3* y2);
    double y = pow(1-t, 3) * y1 + 
               3 * pow(1-t, 2) * t * y2 +
               3 * (1-t) * pow(t, 2) * y3 +
               pow(t, 3) * y4;
    double k = dx == 0 ? dy : dy / dx;
    if (DEBUG) opserr << "Bezier solver: t=" << t << " x=" << x << " y=" << y << " k=" << k << endln;
    if (returnY) {
        if (other) {
            *other = k;
        }
        // *other = k;
        return y;
    } else {
        if (other) {
            *other = y;
        }
        // *other = y; 
        return k;
    }
}
