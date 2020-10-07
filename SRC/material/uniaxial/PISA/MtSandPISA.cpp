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
// Description: This file contains the class implementation for MtSandPISA material. 
//				Provide moment distributed spring for sands according to PISA project. 

#include <elementAPI.h>
#include "MtSandPISA.h" 

#define _USE_MATH_DEFINES
#include <math.h>
#include <Matrix.h>
#include <Channel.h>
#include <vector>
using namespace std;

// Interface function between the interpreter and the MtSandPISA class. 
void*
OPS_MtSandPISA()
{
    // Checking the number of data 
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 8)
    {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: uniaxialMaterial MtSandPISA tag? depth? diameter? pile embedded lenght? sigma_vo_eff? G0? Local distributed lateral load? grid spacing? A0? Au? ?" << endln;
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

    double dData[9] = {0, 0, 0, 0, 0, 0, 0, 1, 1};
    numData = OPS_GetNumRemainingInputArgs();
    if (OPS_GetDoubleInput(&numData, dData) != 0)
    {
        opserr << "WARNING invalid double inputs" << endln;
        return 0;
    }

    // Creating the new material

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;
    theMaterial = new MtSandPISA(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8]);


    if (theMaterial == 0)
    {
        opserr << "WARNING could not create uniaxialMaterial of type MtSandPISA" << endln;
        return 0;
    }

    // return the material
    return theMaterial;
}

// Full constructor with data 
MtSandPISA::MtSandPISA(int tag, double z, double D, double L, double sigma, double Go, double P, double gs, double A0, double Au)
    :UniaxialMaterial(tag, 0),
    depth(z),
    diameter(D),
    embedded_length(L),
    sigma_vo_eff(sigma),
    gmax(Go),
    Pvalue(P),
    grid(gs),
    A0(A0),
    Au(Au)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 4);
}

//	Default constructor
MtSandPISA::MtSandPISA()
    :UniaxialMaterial(0, 0),
    depth(0.0),
    diameter(0.0),
    embedded_length(0.0),
    sigma_vo_eff(0.0),
    gmax(0.0),
    Pvalue(0.0),
    grid(0.0),
    A0(1.0),
    Au(1.0)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 4);
}

//	Default destructor
MtSandPISA::~MtSandPISA()
{
    // does nothing
}

/**
 * Calculate the 4 parameters needed to define the conic expression of PISA project for the p-y curve
 *
 * @returns (Ultimate strain, Initial stiffness,  Curvature, Ultimate reaction)
 */
void MtSandPISA::Mt_conic_function(double depth, double diameter, double embedded_length, double A0, double Au)
{
    //-----------------------------------------------------
    // M-teta curve
    //-----------------------------------------------------
    // Ultimate strain 
    double tetau_norm = 20.0;
    epsilon_u = Au/A0*tetau_norm;
    // Initial stiffness 
    double k0_coeff1 = 20.0;
    double k0m = k0_coeff1;
    initial_stiffness = A0 * k0m;
    // Curvature  
    curvature = 0.0;
    if (curvature == 0.5) curvature += 1e-4; // Otherwise numerical instability for n = 0.5
    // Ultimate reaction 
    double mu_norm_coeff1 = -0.05;
    double mu_norm_coeff2 = 0.21;
    double mu_norm = mu_norm_coeff1 * (depth / embedded_length) + mu_norm_coeff2;
    sigma_u = Au * mu_norm;

    return;
}

/**
 * Evaluate conic function from PISA project
 *
 * @returns ratio of generalised normalised stress to ultimate generalised normalised stress.
 */
double MtSandPISA::conic_function(double epsilon_ratio, double k0, double n, double Epsilon_u, double Sigma_u)
{
    // Normalised roots 
    double a = 1 - 2 * n;
    double b = (2 * n * epsilon_ratio) - (1 - n) * (1 + k0 * epsilon_ratio * (Epsilon_u / Sigma_u));
    double c = k0 * (1 - n) * epsilon_ratio * (Epsilon_u / Sigma_u) - n * (pow(epsilon_ratio, 2));
    double sigma_ratio = epsilon_ratio <= 1. ? (2 * c / (-b + sqrt(pow(b, 2) - 4 * a * c))) : 1.0;
    return sigma_ratio;
}

// Function to compute the strees and tangent from a trial strain sent by the element
int MtSandPISA::setTrialStrain(double strain, double strainRate)
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
MtSandPISA::getStrain(void)
{
    return tStrain;
}

double
MtSandPISA::getStress(void)
{
    return tStress;
}

double
MtSandPISA::getTangent(void)
{
    return tTangent;
}

double
MtSandPISA::getInitialTangent(void)
{
    return this->data(0, 4);
}

int
MtSandPISA::commitState(void)
{
    cStrain = tStrain;
    cStress = tStress;
    cTangent = tTangent;
    return 0;
}

int
MtSandPISA::revertToLastCommit(void)
{
    tStrain = cStrain;
    tStress = cStress;
    tTangent = cTangent;

    return 0;
}

int
MtSandPISA::sendSelf(int cTag, Channel& theChannel)
{
    return -1;
}

int
MtSandPISA::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    return -1;
}

UniaxialMaterial*
MtSandPISA::getCopy(void)
{
    MtSandPISA* theCopy;                 // pointer to a MtSandPISA class
    theCopy = new MtSandPISA();          // new instance of this class
    *theCopy = *this;                    // theCopy (dereferenced) = this (dereferenced pointer) 
    return theCopy;
}

void
MtSandPISA::Print(OPS_Stream& s, int flag)
{
    s << "  MtSandPISA, tag: " << this->getTag() << endln;
    s << "  Depth (m) : " << depth << endln;
    s << "  Pile diameter (m) : " << diameter << endln;
    s << "  Embedded pile length (m) : " << embedded_length << endln;
    s << "  Vertical effective soil stress (kPa) : " << sigma_vo_eff << endln;
    s << "  Small-strain shear modulus (kPa) : " << gmax << endln;
    s << "  Local lateral distributed load (kN/m) at the current depth : " << Pvalue << endln;
    s << "  Scaling factor for initial stiffness (-) : " << A0 << endln;
    s << "  Scaling factor for ultimate stress (-) : " << Au << endln;
}

// Function to create the parameters following the APIrp2GEO method
int
MtSandPISA::revertToStart(void)
{
    // -------- Normalised parameters ------------------------
    if (depth == 0.0)
    {
        depth = 1e-6;
        sigma_vo_eff = 1e-6;
    }

    Mt_conic_function(depth, diameter, embedded_length, A0, Au);

    //------------- Evaluating conic function ---------------
    // Points to discretize the curve
    double v_start = 0.01;
    double v_end = 1.11;
    double delta = 0.01;
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
        // List of distributed moment values for selected points in the mt curve 
        double p_norm = sigma_bar_ratio[i] * sigma_u;
        double p_D = p_norm * abs(Pvalue) * diameter;
        p_F[i] = grid * p_D;
        // List of y-values for selected points in the py curve (m)
        double y_norm = epsilon_bar_ratio[i] * epsilon_u;
        y[i] = y_norm * sigma_vo_eff / gmax;
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

