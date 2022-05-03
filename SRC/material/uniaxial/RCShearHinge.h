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
// $Date: 2022/05/02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/RCShearHinge.cpp,v $

// Written: Amir Reza Tabkhi Wayghan, MASc, Structural Engineering Graduate, Carleton University
//          & Vahid Sadeghian, PhD, Assistant Professor, Carleton University
// Created: May, 2022
// Revision: A
//
//
// Description: This file contains the implementation for the RCShearHinge model developed based on the following journal papers & thesis:
// 				Tabkhi, A.R. and Sadeghian, V. (accepted) “A Shear Hinge Model for Analysis of Reinforced Concrete Columns,” ACI Structural Journal.
//				Tabkhi, A.R. and Sadeghian, V. (2021) “A Shear Hinge Model for Analysis of Reinforced Concrete Beams,” ACI Structural Journal, Vol. 118, No. 6, pp. 279-291.
//				Tabkhi, A.R. (2021) "Development of Shear Plastic Hinge Models for Analysis of Reinforced Concrete Members," MASc Thesis, Carleton University, 2021.

#ifndef RCShearHinge_h
#define RCShearHinge_h

#include <UniaxialMaterial.h>

class RCShearHinge : public UniaxialMaterial
{
public:
    RCShearHinge(int tag, double b, double h, double d, double pRat, double fc, double Ec, double eps0, double fyt, double fyl, double Es,
        double Ast, double Asl, double Asc, double s, double a, double alpha, double cover, double forcecf, double dispcf, double isDeep = 0);

    RCShearHinge();
    ~RCShearHinge();

    const char* getClassType(void) const { return "RCShearHinge"; };

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStrainRate(void);
    double getStress(void);

    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void) { return 0; };


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
    int tag;
    double tStrain, tStress, tTangent;
    double cStrain, cStress, cTangent;
    double vu, deltau, vy, deltay, vcr, deltacr, vf, deltaf;
};


#endif
#pragma once
