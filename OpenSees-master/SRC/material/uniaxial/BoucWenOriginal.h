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

// $Revision: $
// $Date: $
// $URL: $

#ifndef BoucWenOriginal_h
#define BoucWenOriginal_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 07/17
// Revision: A
//
// Description: This file contains the class definition for BoucWenOriginal.
// uniaxial material class. BoucWenOriginal is based on the original Bouc-Wen
// hysteretic material model without strength degradation and without pinching
// effects.

#include <UniaxialMaterial.h>
#include <Matrix.h>

class BoucWenOriginal : public UniaxialMaterial
{
public:
    // constructors
    BoucWenOriginal(int tag,
        double Ei,
        double fy,
        double alphaL,
        double alphaNL = 0.0,
        double mu = 2.0,
        double eta = 1.0,
        double beta = 0.5,
        double gamma = 0.5,
        double tol = 1E-12,
        int maxIter = 25);
    BoucWenOriginal();

    // destructor
    ~BoucWenOriginal();

    // method to get class type
    const char *getClassType(void) const {return "BoucWenOriginal";};

    // public methods to set the state of the material
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain();
    double getStress();
    double getTangent();
    double getInitialTangent();
    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    
    UniaxialMaterial *getCopy();
    
    // public methods for material output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    // private methods
    double sgn(double x);

    // material parameters
    double Ei;          // initial stiffness of material
    double fy;          // yield stress
    double alphaL;      // stiffness ratio of linear elastic component
    double alphaNL;     // stiffness ratio of nonlinear elastic component
    double mu;          // exponent of nonlinear elastic component
    double eta;         // yielding exponent (sharpness of hysteresis loop corners)
    double beta;        // hysteretic shape parameter
    double gamma;       // hysteretic shape parameter
    double tol;         // tolerance for convergence criterion
    int maxIter;        // maximum number of iterations

    // state variables
    double eps;         // trial strain
    double z;           // hysteretic evolution parameter
    double sig;         // trial stress
    double Et;          // tangent stiffness

    // committed history variables
    double epsC;        // strain
    double zC;          // hysteretic evolution parameter
};

#endif
