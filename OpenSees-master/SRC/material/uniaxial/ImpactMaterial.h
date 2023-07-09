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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ImpactMaterial.h,v $

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


#ifndef ImpactMaterial_h
#define ImpactMaterial_h

#include <UniaxialMaterial.h>

class ImpactMaterial : public UniaxialMaterial
{
public:
    ImpactMaterial(int tag, double K1, double K2, double Delta_y, double gap);
    ImpactMaterial();
    ~ImpactMaterial();

    const char *getClassType(void) const {return "ImpactMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);    

    void Print(OPS_Stream &s, int flag =0);

protected:

private:
    double K1;
    double K2;
    double Delta_y;
    double gap;
    double dStrain;
    double Tstress;
    double Tstrain;
    double Ttangent;
    double Cstrain;
    double Cstress;
};

#endif
