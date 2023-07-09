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

// $Revision: 1.2 $
// $Date: 2009-03-27 19:19:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ImpactMaterial.cpp,v $

// File: ~/material/ImpactMaterial.C
//
// Written: md
// Created: 06/2008
//
// Description: This file contains the class implementation for
// ImpactMaterial. This material is based on an approximation to the Hertz contact model proposed by Muthukumar.
//
// References:
// Muthukumar, S., and DesRoches, R. (2006). "A Hertz Contact Model with Non-linear Damping for Pounding Simulation."
//   Earthquake Engineering and Structural Dynamics, 35, 811-828.
// Muthukumar, S. (2003). "A Contact Element Approach with Hysteresis Damping for the Analysis and Design of Pounding
//   in Bridges." PhD Thesis, Georgia Institute of Technology. http://smartech.gatech.edu/
// Nielson, B. (2005). "Analytical Fragility Curves for Highway Bridges in Moderate Seismic Zones." PhD Thesis, Georgia
//  Institute of Technology. http://smartech.gatech.edu/

#include <stdlib.h>

#include <ImpactMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

void *
OPS_ImpactMaterial(void)
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int argc = OPS_GetNumRemainingInputArgs();

    if (argc < 5) {
        opserr << "WARNING incorrect num args want: uniaxialMaterial ImpactMaterial ?tag $K1 $K2 $Delta_y $gap" << endln;
        return 0;
    }

    int    iData[1];
    double dData[4];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial ImpactMaterial tag" << endln;
        return 0;
    }

    numData = 4;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid double data: for ImpactMaterial tag: " << iData[0] << "\n";
        return 0;
    }

    theMaterial = new ImpactMaterial(iData[0], dData[0], dData[1], dData[2], dData[3]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type ImpactMaterial\n";
        return 0;
    }

    return theMaterial;
}


ImpactMaterial::ImpactMaterial(int tag, double k1, double k2, double delta_y, double gap0)
    :UniaxialMaterial(tag,MAT_TAG_ImpactMaterial),
    K1(k1), K2(k2), Delta_y(delta_y), gap(gap0)
{
    if (gap>=0) {
        opserr << "ImpactMaterial::ImpactMaterial -- Initial gap size must be negative for compression-only material\n";
        exit(-1);
    }
    if (Delta_y>=0) {
        opserr << "ImpactMaterial::ImpactMaterial -- Yield displacement must be negative for compression-only material\n";
        exit(-1);
    }

    // Initialize history variables
    this->revertToStart();
    this->revertToLastCommit();
}

ImpactMaterial::ImpactMaterial()
    :UniaxialMaterial(0,MAT_TAG_ImpactMaterial),
    K1(0.0), K2(0.0), Delta_y(0.0), gap(0.0)
{
    // does nothing
}

ImpactMaterial::~ImpactMaterial()
{
    // does nothing
}

int 
ImpactMaterial::setTrialStrain(double strain, double strainRate)
{
    // set the trial strain
    Tstrain = strain;

    // set the strain increment
    dStrain = Tstrain - Cstrain;

    // no gap closure
    if (Tstrain >= gap) {
        Tstress = 0.0;
        Ttangent = 0.0;
    } else {
        // loading or reloading
        if (dStrain < 0.0) {
            // loading with initial stiffness
            Tstress = Cstress+(K1*dStrain);
            Ttangent = K1;
            // loading with secondary stiffness
            if ( (Cstress+(K1*dStrain)) < ((K1*Delta_y)+K2*(Tstrain-gap-Delta_y)) ) {
                Tstress = (K1*Delta_y)+K2*(Tstrain-gap-Delta_y);
                Ttangent = K2;
            }
        }
        // unloading
        else if (dStrain > 0.0) {
            // unloading with initial stiffness
            Tstress = Cstress+(K1*dStrain);
            Ttangent = K1;
            // unloading with secondary stiffness
            if ( (Cstress+(K1*dStrain)) > (K2*(Tstrain-gap)) ) {
                Tstress = K2*(Tstrain-gap);	
                Ttangent = K2;	
            }
        }
    }
    return 0;
}

double 
ImpactMaterial::getStrain(void)
{
    return Tstrain;
}

double 
ImpactMaterial::getStress(void)
{
    return Tstress;
}

double 
ImpactMaterial::getTangent(void)
{
    return Ttangent;
}

double 
ImpactMaterial::getInitialTangent(void)
{
    return K1; 
}

int 
ImpactMaterial::commitState(void)
{
    Cstress = Tstress;
    Cstrain = Tstrain;
    return 0;
}

int 
ImpactMaterial::revertToLastCommit(void)
{
    Tstrain = Cstrain;
    Tstress = Cstress;
    return 0;
}

int 
ImpactMaterial::revertToStart(void)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = K1;
    return 0;
}

UniaxialMaterial *
ImpactMaterial::getCopy(void)
{
    ImpactMaterial *theCopy = new ImpactMaterial(this->getTag(),K1,K2,Delta_y,gap);
    theCopy->Cstress = Cstress;
    theCopy->Cstrain = Cstrain;
    theCopy->Ttangent = Ttangent;
    return theCopy;
}

int 
ImpactMaterial::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(8);
    data(0) = this->getTag();
    data(1) = K1;
    data(2) = K2;
    data(3) = Delta_y;
    data(4) = gap;
    data(5) = Cstress; 
    data(6) = Cstrain;
    data(7) = Ttangent;
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "ImpactMaterial::sendSelf() - failed to send data\n";
    return res;
}

int 
ImpactMaterial::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(8);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "ImpactMaterial::recvSelf() - failed to recv data\n";
    else {
        this->setTag((int)data(0));
        K1 = data(1);
        K2 = data(2);
        Delta_y = data(3);
        gap = data(4);
        Cstress = data(5); 
        Cstrain = data(6);
        Ttangent = data(7);
        Tstress = Cstress;
        Tstrain = Cstrain;
    }
    return res;
}

void 
ImpactMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "ImpactMaterial tag: " << this->getTag() << endln;
        s << "  K1: " << K1 << endln;
        s << "  K2: " << K2 << endln;
        s << "  Delta_y: " << Delta_y << endln;
        s << "  initial gap: " << gap << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ImpactMaterial\", ";
        s << "\"K1\": " << K1 << ", ";
        s << "\"K2\": " << K2 << ", ";
        s << "\"deltaY\": " << Delta_y << ", ";
        s << "\"gap\": " << gap << "}";
    }
}
