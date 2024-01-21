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
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//Modified Ibarra-Medina-Krawinkler with Peak-Oriented Hysteretic Response

//**********************************************************************
// Code Developed by: Ahmed Elkady and Hammad ElJisr
// Last Updated: September 2023
//**********************************************************************

#ifndef IMKPinching_h
#define IMKPinching_h

#include <UniaxialMaterial.h>

class IMKPinching : public UniaxialMaterial
{
public:
    IMKPinching(int tag, double Ke,
        double posUy_0, double posUcap_0, double posUu_0, double posFy_0, double posFcapFy_0, double posFresFy_0,
        double negUy_0, double negUcap_0, double negUu_0, double negFy_0, double negFcapFy_0, double negFresFy_0,
        double LAMBDA_S, double LAMBDA_C, double LAMBDA_A, double LAMBDA_K, double c_S, double c_C, double c_A, double c_K, double D_pos, double D_neg, double kappaF, double kappaD);
    IMKPinching();
    ~IMKPinching();
    const char *getClassType(void) const { return "IMKPinching"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double  getStrain(void);
    double  getStress(void);
    double  getTangent(void);
    double  getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);


protected:

private:
// 25 Fixed input material parameters
    double  Ke;
    double  posUp_0;
    double  posUpc_0;
    double  posUu_0;
    double  posFy_0;
    double  posFcapFy_0;
    double  posFresFy_0;
    double  negUp_0;
    double  negUpc_0;
    double  negUu_0;
    double  negFy_0;
    double  negFcapFy_0;
    double  negFresFy_0;
    double  LAMBDA_S;
    double  LAMBDA_C;
    double  LAMBDA_A;
    double  LAMBDA_K;
    double  c_S;
    double  c_C;
    double  c_A;
    double  c_K;
    double  D_pos;
    double  D_neg;
    double  kappaF;
    double  kappaD;
// 14 Initial Variables
    double  posUy_0,        negUy_0;
    double  posUcap_0,      negUcap_0;
    double  posFcap_0,      negFcap_0;
    double  posKp_0,        negKp_0;
    double  posKpc_0,       negKpc_0;
    double  engRefS;
    double  engRefC;
    double  engRefA;
    double  engRefK;
// History Variables
// 12 Positive U and F
    double  posUy,          cPosUy;
    double  posFy,          cPosFy;
    double  posUcap,        cPosUcap;
    double  posFcap,        cPosFcap;
    double  posUlocal,      cPosUlocal;
    double  posFlocal,      cPosFlocal;
    double  posUglobal,     cPosUglobal;
    double  posFglobal,     cPosFglobal;
    double  posUres,        cPosUres;
    double  posFres,        cPosFres;
    double  posKp,          cPosKp;
    double  posKpc,         cPosKpc;
// 12 Negative U and F
    double  negUy,          cNegUy;
    double  negFy,          cNegFy;
    double  negUcap,        cNegUcap;
    double  negFcap,        cNegFcap;
    double  negUlocal,      cNegUlocal;
    double  negFlocal,      cNegFlocal;
    double  negUglobal,     cNegUglobal;
    double  negFglobal,     cNegFglobal;
    double  negUres,        cNegUres;
    double  negFres,        cNegFres;
    double  negKp,          cNegKp;
    double  negKpc,         cNegKpc;
// 2 Pinching
    double  Fpinch,         cFpinch;
    double  Upinch,         cUpinch;
// 3 State Variables
    double  Ui,             cUi;
    double  Fi,             cFi;
// 3 Stiffness
    double  Kreload,        cKreload, KgetTangent;
    double  Kunload,        cKunload;
// 2 Energy
    double  engAcml,        cEngAcml;
    double  engDspt,        cEngDspt;
// 2 Flag
    bool    Failure_Flag,   cFailure_Flag;
    int     Branch,         cBranch;
};

#endif