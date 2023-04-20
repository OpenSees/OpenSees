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

// $Revision: 1.15 $
// $Date: 2008-08-26 16:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FlagShapeMaterial.cpp,v $

// Written: CLLee
// Created: Feb 2021
// Reference: HardeningMaterial
//
// Description: This file contains the class implementation for 
// FlagShapeMaterial. 

#include "FlagShapeMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_FlagShapeMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: uniaxialMaterial FlagShape tag? E? fy? Eh? <beta?>" << endln;
        return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING: failed to read tag\n";
        return 0;
    }

    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data)) {
        opserr << "WARING: failed to read data\n";
        return 0;
    }

    double beta = 0.0;
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
        numdata = 1;
        if (OPS_GetDouble(&numdata, &beta) < 0) {
            opserr << "WARNING: failed to read beta\n";
            return 0;
        }
    }

    UniaxialMaterial* mat = new FlagShapeMaterial(tag, data[0], data[1], data[2], beta);
    if (mat == 0) {
        opserr << "WARNING: failed to create FlagShapeMaterial material\n";
        return 0;
    }

    return mat;
}

FlagShapeMaterial::FlagShapeMaterial(int tag, double e, double s,
    double eh, double b)
    :UniaxialMaterial(tag, MAT_TAG_FlagShapeMaterial),
    E(e), fy(s), Eh(eh), beta(b)
{

    // Initialize variables
    this->revertToStart();
}

FlagShapeMaterial::FlagShapeMaterial()
    :UniaxialMaterial(0, MAT_TAG_FlagShapeMaterial),
    E(0.0), fy(0.0), Eh(0.0), beta(0.0)
{

    // Initialize variables
    this->revertToStart();
}


int
FlagShapeMaterial::setTrialStrain(double strain, double strainRate)
{

    if (fabs(Tstrain - strain) < DBL_EPSILON)
        return 0;

    //Define and compute Kinematic modulus
    double H = Eh / (1 - Eh / E);
    double fyflag = fabs(CbackStress) > 0 ? beta*fy:fy;

    // Set total strain
    Tstrain = strain;

    // Elastic trial stress
    Tstress = E * (Tstrain - CplasticStrain);

    // Compute trial stress relative to committed back stress
    double xsi = Tstress - CbackStress;

    // Compute yield criterion
    double f = fabs(xsi) - fyflag;

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
        //if (f <= 1.0e-8) {
        // Set trial tangent
        Ttangent = E;
    }

    // Plastic step ... perform return mapping algorithm
    else {

        // Compute consistency parameter
        double dGamma = f / (E + H);

        // Find sign of xsi
        int sign = (xsi < 0) ? -1 : 1;

        // Bring trial stress back to yield surface
        Tstress -= dGamma * E * sign;

        // Update plastic strain
        TplasticStrain = CplasticStrain + dGamma * sign;

    	// Overshoot?
    	if (TplasticStrain*CplasticStrain<0)
    	{
            Tstress = E * Tstrain;
            f = fabs(Tstress) - fy;
    		if (f<=-DBL_EPSILON*E)
    		{
                TplasticStrain = 0.0;
                Ttangent = E;
    		}
    		else
            {
                dGamma = f / (E + H);
                sign = (Tstress > 0) - (Tstress < 0);
                TplasticStrain = dGamma * sign;
                Tstress -= dGamma * E * sign;
                Ttangent = Eh;
            }
    	} else
    	{
            Ttangent = Eh;
    	}

        // Update back stress
        if (fabs(TplasticStrain / (fy / E)) > 1.0e-12) {
            int epsp_sign = (TplasticStrain > 0) - (TplasticStrain < 0);
            TbackStress = epsp_sign * (1 - beta) * fy + H * TplasticStrain;
        }
        else {
            TbackStress = 0.0;
        }

        // Set trial tangent
        // Ttangent = Eh;
    }

    if (beta == 0.0) {
        TplasticStrain = 0.0;
        TbackStress = 0.0;
    }

    return 0;
}

double
FlagShapeMaterial::getStress(void)
{
    return Tstress;
}

double
FlagShapeMaterial::getTangent(void)
{
    return Ttangent;
}

double
FlagShapeMaterial::getStrain(void)
{
    return Tstrain;
}

int
FlagShapeMaterial::commitState(void)
{
    // Commit trial history variables
    CplasticStrain = TplasticStrain;
    CbackStress = TbackStress;

    return 0;
}

int
FlagShapeMaterial::revertToLastCommit(void)
{
    return 0;
}

int
FlagShapeMaterial::revertToStart(void)
{
    // Reset committed history variables
    CplasticStrain = 0.0;
    CbackStress = 0.0;

    // Reset trial history variables
    TplasticStrain = 0.0;
    TbackStress = 0.0;

    // Initialize state variables
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = E;


    return 0;
}

UniaxialMaterial*
FlagShapeMaterial::getCopy(void)
{
    return new FlagShapeMaterial(*this);
}

int
FlagShapeMaterial::sendSelf(int cTag, Channel& theChannel)
{
    int res = 0;

    static Vector data(10);

    data(0) = this->getTag();
    data(1) = E;
    data(2) = fy;
    data(3) = Eh;
    data(4) = beta;
    data(5) = CplasticStrain;
    data(6) = CbackStress;
    data(7) = Tstrain;
    data(8) = Tstress;
    data(9) = Ttangent;

    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "FlagShapeMaterial::sendSelf() - failed to send data\n";

    return res;
}

int
FlagShapeMaterial::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int res = 0;

    static Vector data(10);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "FlagShapeMaterial::recvSelf() - failed to receive data\n";
        E = 0;
        this->setTag(0);
    }
    else {
        this->setTag((int)data(0));
        E = data(1);
        fy = data(2);
        Eh = data(3);
        beta = data(4);
        CplasticStrain = data(5);
        CbackStress = data(6);
        Tstrain = data(7);
        Tstress = data(8);
        Ttangent = data(9);
    }

    return res;
}

void
FlagShapeMaterial::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "FlagShapeMaterial, tag: " << this->getTag() << endln;
        s << "  E: " << E << endln;
        s << "  fy: " << fy << endln;
        s << "  Eh: " << Eh << endln;
        s << "  beta: " << beta << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"FlagShapeMaterial\", ";
        s << "\"E\": " << E << ", ";
        s << "\"fy\": " << fy << ", ";
        s << "\"Eh\": " << Eh << ", ";
        s << "\"beta\": " << beta << "}";
    }
}
