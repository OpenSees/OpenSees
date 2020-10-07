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
// $Revision: 1.12 $
// $Date: 07/10/2020 $
//
// Description: This file contains the class implementation for HbSandPISA material. 
//				Provide the base shear load curve for sandy soils according to PISA project. 

#include <elementAPI.h>
#include "HbSandPISA.h" 

#define _USE_MATH_DEFINES
#include <math.h>
#include <Matrix.h>
#include <Channel.h>
#include <vector>
using namespace std;

// Interface function between the interpreter and the HbSandPISA class. 
void*
OPS_HbSandPISA()
{
    // Checking the number of data 
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5)
    {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: uniaxialMaterial HbSandPISA tag? diameter? pile embedded lenght? sigma_vo_eff? G0? A0? Au?" << endln;
        return 0;
    }

    // parse the input line for the material parameters

    int    iData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0)
    {
        opserr << "WARNING invalid int inputs" << endln;
        return 0;
    }

    double dData[6] = {0, 0, 0, 0, 1, 1};
    numData = OPS_GetNumRemainingInputArgs();;
    if (OPS_GetDoubleInput(&numData, dData) != 0)
    {
        opserr << "WARNING invalid double inputs" << endln;
        return 0;
    }

    // Creating the new material

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;
    theMaterial = new HbSandPISA(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);


    if (theMaterial == 0)
    {
        opserr << "WARNING could not create uniaxialMaterial of type HbSandPISA" << endln;
        return 0;
    }

    // return the material
    return theMaterial;
}

// Full constructor with data 
HbSandPISA::HbSandPISA(int tag, double D, double L, double sigma, double Go, double A0, double Au)
    :UniaxialMaterial(tag, 0),
    diameter(D),
    embedded_length(L),
    sigma_vo_eff(sigma),
    gmax(Go),
    A0(A0),
    Au(Au)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 4);
}

//	Default constructor
HbSandPISA::HbSandPISA()
    :UniaxialMaterial(0, 0),
    diameter(0.0),
    embedded_length(0.0),
    sigma_vo_eff(0.0),
    gmax(0.0),
    A0(1.0),
    Au(1.0)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 4);
}

//	Default destructor
HbSandPISA::~HbSandPISA()
{
    // does nothing
}

/**
 * Calculate the 4 parameters needed to define the conic expression of PISA project for the p-y curve
 *
 * @returns (Ultimate strain, Initial stiffness,  Curvature, Ultimate reaction)
 */
void HbSandPISA::HB_conic_function(double diameter, double embedded_length, double A0, double Au)
{
    //-----------------------------------------------------
    //HB curve
    //-----------------------------------------------------
    // Ultimate strain 
    double vunorm_coeff1 = -0.29;
    double vunorm_coeff2 = 2.31;
    double vBu = vunorm_coeff1 * (embedded_length / diameter) + vunorm_coeff2;
    epsilon_u = Au / A0 * vBu;
    // Initial stiffness 
    double k0Hb_coeff1 = -0.38;
    double k0Hb_coeff2 = 3.02;
    double k0Hb = k0Hb_coeff1 * (embedded_length / diameter) + k0Hb_coeff2;
    initial_stiffness = A0 * k0Hb;
    // Curvature  
    double nHb_coeff1 = -0.05;
    double nHb_coeff2 = 0.94;
    double nHb = nHb_coeff1 * (embedded_length / diameter) + nHb_coeff2;
    curvature = nHb;
    if (curvature == 0.5) curvature += 1e-4; // Otherwise numerical instability for n = 0.5
    // Ultimate reaction 
    double hbu_coeff1 = -0.07;
    double hbu_coeff2 = 0.62;
    double HBu_norm = hbu_coeff1 * (embedded_length / diameter) + hbu_coeff2;
    sigma_u = Au * HBu_norm;

    return;
}

/**
 * Evaluate conic function from PISA project
 *
 * @returns ratio of generalised normalised stress to ultimate generalised normalised stress.
 */
double HbSandPISA::conic_function(double epsilon_ratio, double k0, double n, double Epsilon_u, double Sigma_u)
{
    // Normalised roots 
    double a = 1 - 2 * n;
    double b = (2 * n * epsilon_ratio) - (1 - n) * (1 + k0 * epsilon_ratio * (Epsilon_u / Sigma_u));
    double c = k0 * (1 - n) * epsilon_ratio * (Epsilon_u / Sigma_u) - n * (pow(epsilon_ratio, 2));
    double sigma_ratio = epsilon_ratio <= 1. ? (2 * c / (-b + sqrt(pow(b, 2) - 4 * a * c))) : 1.0;
    return sigma_ratio;
}

// Function to compute the strees and tangent from a trial strain sent by the element
int HbSandPISA::setTrialStrain(double strain, double strainRate)
{
    if (fabs(tStrain - strain) < DBL_EPSILON)
        return 0;

    tStrain = strain;
    tSlope = 0;

    // elastic
    if (tStrain >= data(0, 0) && tStrain <= data(0, 1)) { // elastic
        tSlope = 0;
        tStress = data(0, 2) + (tStrain - data(0, 0)) * data(0, 4);
        tTangent = data(0, 4);
    }
    else if (tStrain < data(0, 0)) { // search neg of data
        tSlope = 1;
        while (tSlope < numSlope && tStrain < data(tSlope, 0))
            tSlope++;
        if (tSlope == numSlope)
            tSlope = numSlope - 1;
        tStress = data(tSlope, 2) + (tStrain - data(tSlope, 0)) * data(tSlope, 4);
        tTangent = data(tSlope, 4);
    }
    else { // search pos side of data
        tSlope = 1;
        while (tSlope < numSlope && tStrain > data(tSlope, 1))
            tSlope++;
        if (tSlope == numSlope)
            tSlope = numSlope - 1;
        tStress = data(tSlope, 3) + (tStrain - data(tSlope, 1)) * data(tSlope, 4);
        tTangent = data(tSlope, 4);
    }

    return 0;
}

double
HbSandPISA::getStrain(void)
{
    return tStrain;
}

double
HbSandPISA::getStress(void)
{
    return tStress;
}

double
HbSandPISA::getTangent(void)
{
    return tTangent;
}

double
HbSandPISA::getInitialTangent(void)
{
    return this->data(0, 4);
}

int
HbSandPISA::commitState(void)
{
    cStrain = tStrain;
    cStress = tStress;
    cTangent = tTangent;

    return 0;
}

int
HbSandPISA::revertToLastCommit(void)
{
    tStrain = cStrain;
    tStress = cStress;
    tTangent = cTangent;

    return 0;
}

int
HbSandPISA::sendSelf(int cTag, Channel& theChannel)
{
    return -1;
}

int
HbSandPISA::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    return -1;
}

UniaxialMaterial*
HbSandPISA::getCopy(void)
{
    HbSandPISA* theCopy;                 // pointer to a HbSandPISA class
    theCopy = new HbSandPISA();          // new instance of this class
    *theCopy = *this;                    // theCopy (dereferenced) = this (dereferenced pointer) 
    return theCopy;
}

void
HbSandPISA::Print(OPS_Stream& s, int flag)
{
    s << "HbSandPISA, tag: " << this->getTag() << endln;
    s << "  Pile diameter (m) : " << diameter << endln;
    s << "  Embedded pile length (m) : " << embedded_length << endln;
    s << "  Vertical effective soil stress (kPa) : " << sigma_vo_eff << endln;
    s << "  Small-strain shear modulus (kPa) : " << gmax << endln;
    s << "  Scaling factor for initial stiffness (-) : " << A0 << endln;
    s << "  Scaling factor for ultimate stress (-) : " << Au << endln;
}

// Function to create the parameters following the APIrp2GEO method
int
HbSandPISA::revertToStart(void)
{
    // -------- Normalised parameters ------------------------
    HB_conic_function(diameter, embedded_length, A0, Au);

    //------------- Evaluating conic function ---------------
    // Points to discretize the curve
    double v_start = 1e-3;
    double delta = 1e-3;
    double v_end = 1 + delta;
    int const v_num = (v_end - v_start) / delta;
    numSlope = v_num + 1;
    vector <double> epsilon_bar_ratio(numSlope, 0.0);
    vector <double> sigma_bar_ratio(numSlope, 0.0);
    for (int i = 0; i < numSlope; i++)
    {
        epsilon_bar_ratio[i] = v_start + i * delta;
        double sigma_ratio = conic_function(epsilon_bar_ratio[i], initial_stiffness, curvature, epsilon_u, sigma_u);
        sigma_bar_ratio[i] = sigma_ratio;
    }

    //---------- Calculating nor normalized soil curve -------
    vector <double> p_F(numSlope, 0.0);
    vector <double> y(numSlope, 0.0);
    for (int i = 0; i < numSlope; i++)
    {
        // List of P-values for selected points in the P curve
        double hb_norm = sigma_bar_ratio[i] * sigma_u;
        double p_D = hb_norm * sigma_vo_eff * pow(diameter,2);
        p_F[i] = p_D;
        // List of y-values for selected points in the py curve (m)
        double y_norm = epsilon_bar_ratio[i] * epsilon_u;
        y[i] = y_norm * diameter * sigma_vo_eff / gmax;
    }

    //------------- Creating data matrix--------------------
    data.resize(numSlope, 6);

    data(0, 0) = -y[0];             // neg yield strain
    data(0, 1) = y[0];              // pos yield strain
    data(0, 2) = -p_F[0];           // neg yield stress
    data(0, 3) = p_F[0];            // pos yield stress
    data(0, 4) = p_F[0] / y[0];     // slope
    data(0, 5) = y[0];              // dist - [0-1)/2

    for (int i = 1; i < numSlope; i++)
    {
        data(i, 0) = -y[i];
        data(i, 1) = y[i];
        data(i, 2) = -p_F[i];
        data(i, 3) = p_F[i];
        data(i, 4) = (p_F[i] - p_F[i - 1]) / (y[i] - y[i - 1]);
        data(i, 5) = y[i] - y[i - 1];
    }
    tStrain = 0.0;
    tStress = 0.0;
    tTangent = data(0, 4);


    cStrain = 0.0;
    cStress = 0.0;
    cTangent = tTangent;

    tSlope = 0;

    this->commitState();

    return 0;
}

