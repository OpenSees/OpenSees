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

// Written: csasj 
// $Revision: 1.0 $
// $Date: 26/07/2023 $
//
// Description: This file contains the class implementation for TzSandCPT material. 
//              Calculates the load-transfer curves that can be used to represent 
//				the axial load-displacement relationship between the pile and the 
//				soil following the unified CPT method for sands (Lehane et al. 2020).

#include <elementAPI.h>
#include "TzSandCPT.h" 

#define _USE_MATH_DEFINES
#include <math.h>
#include <Matrix.h>
#include <Channel.h>
#include <vector>
using namespace std;

// Interface function between the interpreter and the TzSandCPT class. 
void*
OPS_TzSandCPT()
{
    // Checking the number of data 
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 6)
    {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: uniaxialMaterial TzSandCPT tag? qc? sigma? D? t? h? dz?" << endln;
        return 0;
    }

    // parse the input line for the material parameters

    int    iData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) < 0)
    {
        opserr << "WARNING invalid int inputs" << endln;
        return 0;
    }

    double dData[6] = {0,0,0,0,0,0};
    numData = OPS_GetNumRemainingInputArgs();
    if (OPS_GetDoubleInput(&numData, dData) < 0)
    {
        opserr << "WARNING invalid double inputs" << endln;
        return 0;
    }

    // Creating the new material

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;
    theMaterial = new TzSandCPT(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);


    if (theMaterial == 0)
    {
        opserr << "WARNING could not create uniaxialMaterial of type TzSandCPT" << endln;
        return 0;
    }

    // return the material
    return theMaterial;
}

// Full constructor with data 
TzSandCPT::TzSandCPT(int tag, double qc, double Sv, double D, double t, double h, double dz)
    :UniaxialMaterial(tag, 0),
    q_c(qc),
    sigma_vo_eff(Sv),
    diameter(D),
    wall_thickness(t),
    h_dist(h),
    delta_h(dz)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data_c(0, 2);
}

//	Default constructor
TzSandCPT::TzSandCPT()
    :UniaxialMaterial(0, 0),
    q_c(0.0),
    sigma_vo_eff(0.0),
    diameter(0.0),
    wall_thickness(0.0),
    h_dist(0.0),
    delta_h(0.0)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data_c(0, 2);
}

//	Default destructor
TzSandCPT::~TzSandCPT()
{
    // does nothing
}

/**
 * Calculates the ultimate shaft friction and the corresponding "peak" settlement
 * Note that no distinction regarding the loading direction (compression or tension)
 * is not made at this state.
 * @returns 
 */

void TzSandCPT::ultimate_capacity(double qc, double Sv, double D, double t, double h)
{
    // To avoid numerical problems
    if (Sv < 1e-6)
    {
        Sv = 1e-6;
    }
    if (qc < 1e-6)
    {
        qc = 1e-6;
    }
    double diameter_in = D - 2 * t;  // Inner pile diameter (m)
    double d_cpt = 35.7e-3;          // CPT diameter probe (m)
    // Plug length ratio
    double plr = tanh(0.3 * sqrt(diameter_in / d_cpt));
    double a_re = 1 - plr * pow(diameter_in / D, 2);
    // Stationary radial effective stress (kPa)
    double dist;
    if (h/D > 1) {
        dist = h/D;
    }
    else {
        dist = 1;
    }
    double sigma_rc = qc / 44 * pow(a_re, 0.3) * pow(dist, -0.4);
    // Increases in radial effective stress during pile loading (kPa)
    double sigma_rd = qc / 10 * pow(qc / Sv, -0.33) * (d_cpt / D);
    // Ultimate interface friction(degrees)
    double delta_f = 29;  
    // Ultimate shaft capacity (kPa)
    tau_f = (sigma_rc + sigma_rd) * tan(delta_f * M_PI /180);
    // atmospheric pressure (kPa)
    double pa = 100; 
    // Peak settlement
    w_f = sqrt(qc) * pow(Sv, 0.25) * D / pow(pa, 0.75);

    return;
}

// Function to compute the strees and tangent from a trial strain sent by the element
int TzSandCPT::setTrialStrain(double strain, double strainRate)
{
    if (fabs(tStrain - strain) < DBL_EPSILON)
        return 0;
    tStrain = strain;
    tSlope = 0;
    Matrix data;              // stored load transfer function
    if (tStrain <= 0) {
        // Tension
        data = data_t;
    }
    else {
        // Compression
        data = data_c;
    }
    // initial elastic part
    if (fabs(tStrain) <= fabs(data(0, 0))) { 
        tSlope = 0;
        tStress = tStrain * data(0, 2);
        tTangent = data(0, 2);
    }
    // yield
    else if (fabs(tStrain) >= fabs(data(numSlope-1, 0))) {
        tStress = data(numSlope-1, 1);
        tTangent = 0.;
    }
    // searching 
    else { 
        tSlope = 1;
        while (tSlope < numSlope && fabs(tStrain) > fabs(data(tSlope, 0)))
            tSlope++;
        if (tSlope == numSlope)
            tSlope = numSlope - 1;
        tStress = data(tSlope, 1) + (tStrain - data(tSlope, 0)) * data(tSlope, 2);
        tTangent = data(tSlope, 2);
    }

    return 0;
}

double
TzSandCPT::getStrain(void)
{
    return tStrain;
}

double
TzSandCPT::getStress(void)
{
    return tStress;
}

double
TzSandCPT::getTangent(void)
{
    return tTangent;
}

double
TzSandCPT::getInitialTangent(void)
{
    return this->data_c(0, 2);
}

int
TzSandCPT::commitState(void)
{
    cStrain = tStrain;
    cStress = tStress;
    cTangent = tTangent;
    return 0;
}

int
TzSandCPT::revertToLastCommit(void)
{
    tStrain = cStrain;
    tStress = cStress;
    tTangent = cTangent;

    return 0;
}

int
TzSandCPT::sendSelf(int cTag, Channel& theChannel)
{
    return -1;
}

int
TzSandCPT::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    return -1;
}

UniaxialMaterial*
TzSandCPT::getCopy(void)
{
    TzSandCPT* theCopy;                 // pointer to a TzSandCPT class
    theCopy = new TzSandCPT();          // new instance of this class
    *theCopy = *this;                   // theCopy (dereferenced) = this (dereferenced pointer) 
    return theCopy;
}

void
TzSandCPT::Print(OPS_Stream& s, int flag)
{
    s << "TzSandCPT, tag: " << this->getTag() << endln;
    s << "  qc (kPa) : " << q_c << endln;
    s << "  Vertical effective soil stress (kPa): " << sigma_vo_eff << endln;
    s << "  Pile diameter (m): " << diameter << endln;
    s << "  Wall thickness (m): " << wall_thickness << endln;
    s << "  Distance to pile toe (m): " << h_dist << endln;
    s << "  Shaft capacity compression (kPa): " << tau_f << endln;
    s << "  Peak displacement compression (m): " << w_f / 1250 << endln;
    s << "  Shaft capacity tension (kPa): " << 0.75 * tau_f << endln;
    s << "  Peak displacement tension (m): " << w_f / 625 << endln;
}

// Load transfer function initialize here
int
TzSandCPT::revertToStart(void)
{
    // -------- Ultimate shaft friction calculation --------
    ultimate_capacity(q_c, sigma_vo_eff, diameter, wall_thickness, h_dist);
    //-------- Load transfer function --------
    // Points to discretize the curve
    double v_start = 0.125;  // same initial stifness than API
    double delta = 1e-1;
    double v_end = 1 + delta;
    int const v_num = (v_end - v_start) / delta;
    numSlope = v_num + 1;
    vector <double> epsilon_bar_ratio(numSlope, 0.0);
    vector <double> sigma_bar_ratio(numSlope, 0.0);
    // Normalized t-z curve
    for (int i = 0; i < numSlope; i++)
    {
        epsilon_bar_ratio[i] = v_start + i * delta;
        sigma_bar_ratio[i] = 2 * epsilon_bar_ratio[i] * (1 - 0.5 *epsilon_bar_ratio[i]);
    }
    // ------------- t-z curve for compression loading -------------
    double A_C = 1250;
    data_c.resize(numSlope, 4);
    data_c(0, 0) = epsilon_bar_ratio[0] * w_f / A_C;                            // yield strain 
    data_c(0, 1) = sigma_bar_ratio[0] * tau_f * M_PI * diameter * delta_h;      // yield force
    data_c(0, 2) = data_c(0, 1) / data_c(0, 0);                                 // slope
    data_c(0, 3) = data_c(0, 0);                                                // dist 
    for (int i = 1; i < numSlope; i++)
    {
        data_c(i, 0) = epsilon_bar_ratio[i] * w_f / A_C;
        data_c(i, 1) = sigma_bar_ratio[i] * tau_f * M_PI * diameter * delta_h;
        data_c(i, 2) = (data_c(i, 1) - data_c(i - 1, 1)) / (data_c(i, 0) - data_c(i - 1, 0));
        data_c(i, 3) = data_c(i, 0) - data_c(i - 1, 0);
    }
    tStrain = 0.0;
    tStress = 0.0;
    tTangent = data_c(0, 2);

    cStrain = 0.0;
    cStress = 0.0;
    cTangent = tTangent;

    tSlope = 0;

    // ------------- t-z curve for tension loading -------------
    double A_T = 625;
    double f_T = 0.75;  // 75% shaft capacity in tension
    data_t.resize(numSlope, 4);
    data_t(0, 0) = -epsilon_bar_ratio[0] * w_f / A_T;                              // pos yield strain (tension)
    data_t(0, 1) = -f_T * sigma_bar_ratio[0] * tau_f * M_PI * diameter * delta_h;  // pos yield force
    data_t(0, 2) = data_t(0, 1) / data_t(0, 0);                                    // slope
    data_t(0, 3) = data_t(0, 0);                                                   // dist - [0-1)/2
    for (int i = 1; i < numSlope; i++)
    {
        data_t(i, 0) = -epsilon_bar_ratio[i] * w_f / A_T;
        data_t(i, 1) = -f_T * sigma_bar_ratio[i] * tau_f * M_PI * diameter * delta_h;
        data_t(i, 2) = (data_t(i, 1) - data_t(i - 1, 1)) / (data_t(i, 0) - data_t(i - 1, 0));
        data_t(i, 3) = data_t(i, 0) - data_t(i - 1, 0);
    }

    this->commitState();

    return 0;
}

