/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** This file contains the class implementation for QbSandCPT          **
** material. This uniaxial material material represents the           **
** shaft-load transfer curve (q-z spring) according to the new Unifi- **
** ed CPT-based method for driven piles in sands. The formulation     **
** incorporates the maximum skin friction and end-bearing CPT values  **
** based on the new ISO-19901-4.                                      **
**                                                                    **
**                                                                    **
** Written:                                                           **
**   csasj (carlos.sastre.jurado@vub.be)                              **
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

// $Revision: 1.1    $
// $Date: 01/02/2024 $

#include <elementAPI.h>
#include "QbSandCPT.h" 
#include <OPS_Globals.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <vector>
#include <Matrix.h>
#include <Channel.h>

#include <string.h>

using namespace std;

// Interface function between the interpreter and the QbSandCPT class. 
void*
OPS_QbSandCPT()
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int    iData[1];
    double dData[4];
    int numData = 1;

    // Checking tag material
    if (OPS_GetIntInput(&numData, iData) < 0)
    {
        opserr << "WARNING invalid uniaxialMaterial QbSandCPT tag" << endln;
        return 0;
    }

    // Checking number of arguments
    numData = OPS_GetNumRemainingInputArgs();
    if (numData < 4)
    {
        opserr << "WARNING insufficient arguments" << endln;
        opserr << "Want: uniaxialMaterial QbSandCPT" << iData[0] << "qp? D? t? dcpt?" << endln;
        return 0;
    }
    else if (numData > 4)
    {
        opserr << "WARNING number of arguments exceeded" << endln;
        opserr << "Want: uniaxialMaterial QbSandCPT " << iData[0] << "qp? D? t? dcpt?" << endln;
        return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();
    if (OPS_GetDoubleInput(&numData, dData) != 0)
    {
        opserr << "Invalid #args, want: uniaxialMaterial QbSandCPT " << iData[0] << "qp? D? t? dcpt?" << endln;
        return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new QbSandCPT(iData[0], dData[0], dData[1], dData[2], dData[3]);

    if (theMaterial == 0)
    {
        opserr << "WARNING could not create uniaxialMaterial of type QbSandCPT" << endln;
        return 0;
    }

    // return the material
    return theMaterial;
}

// Full constructor with data 
QbSandCPT::QbSandCPT(int tag, double qp, double D, double t, double dcpt)
    :UniaxialMaterial(tag, 0),
    q_p(qp),
    diameter(D),
    wall_thickness(t),
    d_cpt(dcpt)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 2);
}

//	Default constructor
QbSandCPT::QbSandCPT()
    :UniaxialMaterial(0, 0),
    q_p(0.0),
    diameter(0.0),
    wall_thickness(0.0),
    d_cpt(0.0)
{
    // Initialize all variables are needed for the material algorithm
    this->revertToStart();
    double initialTangent = data(0, 2);
}

//	Default destructor
QbSandCPT::~QbSandCPT()
{
    // does nothing
}

/**
 * Calculates the ultimate shaft friction and the corresponding "peak" settlement
 * Note that no distinction regarding the loading direction (compression or tension)
 * is not made at this state.
 * @returns 
 */

void QbSandCPT::ultimate_capacity(double qp, double D, double t, double dcpt)
{
    double A_re;
    // Effective area ratio calculation
    if (t < DBL_EPSILON){
        A_re = 1.;
    }
    else {
        double diameter_in = D - 2 * t;  // Inner pile diameter (m)
        // Plug length ratio
        double plr = tanh(0.3 * sqrt(diameter_in / dcpt));
        A_re = 1 - plr * pow(diameter_in / D, 2);
    }
    // qb01
    q_b01 = (0.12 + 0.38 * A_re) * q_p;

    return;
}

// Function to compute the strees and tangent from a trial strain sent by the element
int QbSandCPT::setTrialStrain(double strain, double strainRate)
{
    if (fabs(tStrain - strain) < DBL_EPSILON)
        return 0;
    tStrain = strain;
    tSlope = 0;
    if (tStrain <= 0) {
        // Tension, no end bearing resistance
        tStress = 0.;
        tTangent = 0.;
    }
    else {
        // Compression
        // initial elastic part
        if (fabs(tStrain) <= fabs(data(0, 0))) {
            tSlope = 0;
            tStress = tStrain * data(0, 2);
            tTangent = data(0, 2);
        }
        // yield
        else if (fabs(tStrain) >= fabs(data(numSlope - 1, 0))) {
            tStress = data(numSlope - 1, 1);
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
    }
    return 0;
}

double
QbSandCPT::getStrain(void)
{
    return tStrain;
}

double
QbSandCPT::getStress(void)
{
    return tStress;
}

double
QbSandCPT::getTangent(void)
{
    return tTangent;
}

double
QbSandCPT::getInitialTangent(void)
{
    return this->data(0, 2);
}

int
QbSandCPT::commitState(void)
{
    cStrain = tStrain;
    cStress = tStress;
    cTangent = tTangent;
    return 0;
}

int
QbSandCPT::revertToLastCommit(void)
{
    tStrain = cStrain;
    tStress = cStress;
    tTangent = cTangent;

    return 0;
}

int
QbSandCPT::sendSelf(int cTag, Channel& theChannel)
{
    static Vector data(10);
    data(0) = this->getTag();

    /***     Input parameters     ***/
    data(1) = q_p;
    data(2) = diameter;
    data(3) = wall_thickness;
    data(4) = d_cpt;

    /***     Axial calculations    ***/
    data(5) = q_b01;

    /*** CONVERGED State Variables ***/
    data(6) = cStrain;
    data(7) = cStress;
    data(8) = cTangent;

    /*** Curve parameters          ***/
    data(9) = numSlope;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)
        opserr << "QbSandCPT::sendSelf() - failed to send data" << endln;

    return 0;
}

int
QbSandCPT::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();
    static Vector data(10);

    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)
    {
        opserr << "QbSandCPT::recvSelf() - failed to receive data" << endln;
        this->setTag(0);
        return -1;
    }

    this->setTag(int(data(0)));

    /***     Input parameters     ***/
    q_p = data(1);
    diameter = data(2);
    wall_thickness = data(3);
    d_cpt = data(4);

    /***     Axial calculations    ***/
    q_b01 = data(5);

    /*** CONVERGED State Variables ***/
    cStrain = data(6);
    cStress = data(7);
    cTangent = data(8);

    /*** Curve parameters          ***/
    numSlope = data(9);

    return 0;
}

UniaxialMaterial*
QbSandCPT::getCopy(void)
{
    QbSandCPT* theCopy;                 // pointer to a QbSandCPT class
    theCopy = new QbSandCPT();          // new instance of this class
    *theCopy = *this;                   // theCopy (dereferenced) = this (dereferenced pointer) 
    return theCopy;
}

void
QbSandCPT::Print(OPS_Stream& s, int flag)
{
    s << "QbSandCPT, tag: " << this->getTag() << endln;
    s << "  qp : " << q_p << endln;
    s << "  pile diameter : " << diameter << endln;
    s << "  wall thickness : " << wall_thickness << endln;
    s << "  diameter CPT probe : " << d_cpt << endln;
    s << "  qb01 : " << q_b01 << endln;
}

// Load transfer function initialize here
int
QbSandCPT::revertToStart(void)
{
    // -------- qb01 calculation --------
    ultimate_capacity(q_p, diameter, wall_thickness, d_cpt);
    //-------- Load transfer function --------
    //-------- Load transfer function --------
    // Points to discretize the curve
    double v_start = 0.1;    // similar initial stifness than API
    double delta = 1e-1;
    double v_end = 1 + delta;
    int const v_num = (v_end - v_start) / delta;
    numSlope = v_num;
    vector <double> qb_qb01(numSlope, 0.0);
    vector <double> w_D(numSlope, 0.0);                        // w/D
    // Normalized Qb-wb curve
    for (int i = 0; i < numSlope; i++)
    {
        qb_qb01[i] = v_start + i * delta;
        w_D[i] = 0.01 * qb_qb01[i] / (1 - 0.9 * qb_qb01[i]);
    }
    // ------------- Non normalised curve -------------
    double A_tip = M_PI * pow(diameter, 2) / 4;
    data.resize(numSlope, 4);
    data(0, 0) = w_D[0] * diameter;                             // neg yield displacement (compression)
    data(0, 1) = qb_qb01[0] * q_b01 * A_tip;                    // neg yield force (tip compression)
    data(0, 2) = data(0, 1) / data(0, 0);                       // slope
    data(0, 3) = data(0, 0);                                    // dist 
    for (int i = 1; i < numSlope; i++)
    {
        data(i, 0) = w_D[i] * diameter;
        data(i, 1) = qb_qb01[i] * q_b01 * A_tip;
        data(i, 2) = (data(i, 1) - data(i - 1, 1)) / (data(i, 0) - data(i - 1, 0));
        data(i, 3) = data(i, 0) - data(i - 1, 0);
    }
    tStrain = 0.0;
    tStress = 0.0;
    tTangent = data(0, 2);

    cStrain = 0.0;
    cStress = 0.0;
    cTangent = tTangent;

    tSlope = 0;

    this->commitState();

    return 0;
}

