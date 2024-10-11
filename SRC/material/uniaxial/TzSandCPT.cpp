/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** This file contains the class implementation for TzSandCPT          **
** material. This uniaxial material material represents the           **
** shaft-load transfer curve (t-z spring) according to the new Unifi- **
** ed CPT-based method for driven piles in sands. The formulation     **
** incorporates the maximum skin friction and end-bearing CPT values  **
** based on the new ISO-19901-4.                                      **
**                                                                    **
**                                                                    **
** Code written by:                                                   **
**   Carlos Sastre Jurado (carlos.sastre.jurado@vub.be)               **
**                                                                    **
**                                                                    **
** References:                                                        **
**   - LEHANE, Barry M., et al. A new'unified'CPT-based axial pile    **
**   capacity design method for driven piles in sand. In: 4th Inter-  **
**   national Symposium on Frontiers in Offshore Geotechnics          **
**   (postponed). 2020. p. 462-477.                                   **
**   - LEHANE, B. M.; LI, Lin; BITTAR, E. J. Cone penetration         **
**   test-based load-transfer formulations for driven piles in sand.  **
**   Geotechnique Letters, 2020, 10.4: 568-574.                       **
** ****************************************************************** */

// $Revision: 1.0    $
// $Date: 19/01/2024 $

#include "TzSandCPT.h" 
#include <elementAPI.h>
#include <OPS_Globals.h>


#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>

#include <string.h>

using namespace std;

// Interface function between the interpreter and the TzSandCPT class. 
void*
OPS_TzSandCPT()
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int    iData[1];
    double dData[9];
    int numData = 1;

    // Checking tag material
    if (OPS_GetIntInput(&numData, iData) < 0)
    {
        opserr << "WARNING invalid uniaxialMaterial TzSandCPT tag" << endln;
        return 0;
    }

    // Checking number of arguments
    numData = OPS_GetNumRemainingInputArgs();
    if (numData < 8)
    {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: uniaxialMaterial TzSandCPT "<< iData[0] << "qc? sigma? D? t? h? dz? dcpt? pa?" << endln;
        return 0;
    }
            else if (numData > 9)
    {
        opserr << "WARNING number of arguments exceeded" << endln;
        opserr << "Want: uniaxialMaterial TzSandCPT "<< iData[0] <<"qc? sigma? D? t? h? dz? dcpt? pa? delta_f?"<< endln;
        return 0;
    }

    if (OPS_GetDoubleInput(&numData, dData) != 0) 
    {
        opserr << "Invalid #args, want: uniaxialMaterial TzSandCPT " << iData[0] << " qc? sigma? D? t? h? dz? dcpt? pa? delta_f?" << endln;
        return 0;
    }

    // Default variables
    if (numData == 8)
    {
        dData[8] = IFA_DEFAULT;
    }

    // Parsing was successful, allocate the material
    theMaterial = new TzSandCPT(iData[0], 
        dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], 
        dData[6], dData[7], dData[8]);

    if (theMaterial == 0)
    {
        opserr << "WARNING could not create uniaxialMaterial of type TzSandCPT" << endln;
        return 0;
    }

    return theMaterial;
}

// Full constructor with data 
TzSandCPT::TzSandCPT(int tag, double qc, double Sv, double D, double t, double h, double dz, double dcpt, double pa, double d_f)
    :UniaxialMaterial(tag, 0),
    q_c(qc),
    sigma_vo_eff(Sv),
    diameter(D),
    wall_thickness(t),
    h_dist(h),
    delta_h(dz),
    d_cpt(dcpt),
    p_a(pa), 
    delta_f(d_f)
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
    delta_h(0.0),
    d_cpt(0.0),
    p_a(0.0),
    delta_f(0.0)
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

void TzSandCPT::ultimate_capacity(double qc, double Sv, double D, double t, double h, 
    double dcpt, double pa, double d_f)
{
    // To avoid numerical problems
    if (Sv < DBL_EPSILON)
    {
        Sv = 1e-9;
    }
    if (qc < DBL_EPSILON)
    {
        qc = 1e-9;
    }
    double diameter_in = D - 2 * t;  // Inner pile diameter (m)
    // Plug length ratio
    double plr = tanh(0.3 * sqrt(diameter_in / dcpt));
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
    // Increases in radial effective stress during pile loading 
    double sigma_rd = qc / 10 * pow(qc / Sv, -0.33) * (dcpt / D); 
    // Ultimate shaft capacity
    tau_f = (sigma_rc + sigma_rd) * tan(d_f * M_PI /180);
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
TzSandCPT::sendSelf(int cTag, Channel &theChannel)
{
    static Vector data(15);
    data(0) = this->getTag();

    /***     Input parameters     ***/
    data(1) = q_c;
    data(2) = sigma_vo_eff;
    data(3) = diameter;
    data(4) = wall_thickness;
    data(5) = h_dist;
    data(6) = delta_h;
    data(7) = d_cpt;
    data(8) = p_a;

    /***     Axial calculations    ***/
    data(9) = tau_f;
    data(10) = w_f;

    /*** CONVERGED State Variables ***/
    data(11) = cStrain;
    data(12) = cStress;
    data(13) = cTangent;

    /*** Curve parameters          ***/
    data(14) = numSlope;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)
        opserr << "TzSandCPT::sendSelf() - failed to send data" << endln;

    return 0;
}

int
TzSandCPT::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();
    static Vector data(15);

    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)
    {
        opserr << "TzSandCPT::recvSelf() - failed to receive data" << endln;
        this->setTag(0);
        return -1;
    }

    this->setTag(int(data(0)));

    /***     Input parameters     ***/
    q_c = data(1);
    sigma_vo_eff = data(2);
    diameter = data(3);
    wall_thickness = data(4);
    h_dist = data(5);
    delta_h = data(6);
    d_cpt = data(7);
    p_a = data(8);

    /***     Axial calculations    ***/
    tau_f = data(9);
    w_f = data(10);

    /*** CONVERGED State Variables ***/
    cStrain = data(11);
    cStress = data(12);
    cTangent = data(13);

    /*** Curve parameters          ***/
    numSlope = data(14);

    return 0;
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
    s << "  qc : " << q_c << endln;
    s << "  vertical effective soil stress : " << sigma_vo_eff << endln;
    s << "  pile diameter : " << diameter << endln;
    s << "  wall thickness : " << wall_thickness << endln;
    s << "  distance to pile toe : " << h_dist << endln;
    s << "  diameter CPT probe : " << d_cpt << endln;
    s << "  atmospheric pressure : " << p_a << endln;
    s << "  interface friction angle : " << delta_f << endln;
    s << "  shaft capacity compression : " << tau_f << endln;
    s << "  peak settlement : " << w_f / 1250 << endln;
    s << "  shaft capacity tension : " << 0.75 * tau_f << endln;
    s << "  peak displacement tension : " << w_f / 625 << endln;
}

// Load transfer function initialize here
int
TzSandCPT::revertToStart(void)
{
    // -------- Ultimate shaft friction calculation --------
    ultimate_capacity(q_c, sigma_vo_eff, diameter, wall_thickness, h_dist, d_cpt, p_a, delta_f);
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

