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

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 07/17
// Revision: A
//
// Description: This file contains the implementation of the
// BoucWenOriginal uniaxial material class.

#include <BoucWenOriginal.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

void* OPS_BoucWenOriginal()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
        opserr << "WARNING: Insufficient arguments\n";
        opserr << "Want: uniaxialMaterial BoucWenOriginal tag E fy alphaL" << endln;
        return 0;
    }
    
    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid tag\n";
        return 0;
    }
    
    double data[9] = { 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.5, 0.5, 1.0E-8 };
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 9) {
        numdata = 9;
    }
    if (OPS_GetDoubleInput(&numdata, data)) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
    }
    
    int maxIter = 25;
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }
    }
    
    UniaxialMaterial* mat = new BoucWenOriginal(tag, data[0], data[1], data[2],
        data[3], data[4], data[5], data[6], data[7], data[8], maxIter);
    if (mat == 0) {
        opserr << "WARNING: failed to create BoucWenOriginal material\n";
        return 0;
    }
    
    return mat;
}


BoucWenOriginal::BoucWenOriginal(int tag,
    double ei, double _fy, double alphal, double alphanl,
    double _mu, double _eta, double _beta, double _gamma,
    double _tol, int maxiter)
    : UniaxialMaterial(tag, MAT_TAG_BoucWenOriginal),
    Ei(ei), fy(_fy), alphaL(alphal), alphaNL(alphanl),
    mu(_mu), eta(_eta), beta(_beta), gamma(_gamma),
    tol(_tol), maxIter(maxiter),
    eps(0.0), z(0.0), sig(0.0), Et(ei),
    epsC(0.0), zC(0.0)
{
    
}


BoucWenOriginal::BoucWenOriginal()
    : UniaxialMaterial(0, MAT_TAG_BoucWenOriginal),
    Ei(0.0), fy(0.0), alphaL(0.0), alphaNL(0.0),
    mu(2.0), eta(1.0), beta(0.5), gamma(0.5),
    tol(1E-12), maxIter(25),
    eps(0.0), z(0.0), sig(0.0), Et(0.0),
    epsC(0.0), zC(0.0)
{
    
}


BoucWenOriginal::~BoucWenOriginal()
{

}


int BoucWenOriginal::setTrialStrain(double strain, double strainRate)
{
    eps = strain;
    double delta_eps = eps - epsC;
    if (fabs(delta_eps) > 0.0) {
        
        // initialize stiffnesses
        double k2 = alphaL*Ei;
        double k3 = alphaNL*Ei;
        double k0 = Ei - k2;

        // get yield strain
        double epsy = fy / Ei;

        // get yield force of hysteretic component
        double qd = fy - k2*epsy - k3*pow(epsy, mu);
        
        // calculate hysteretic evolution parameter z using Newton-Raphson
        int iter = 0;
        double zAbs, tmp1, f, Df, delta_z;
        do {
            zAbs = fabs(z);
            if (zAbs == 0.0)    // check because of negative exponents
                zAbs = DBL_EPSILON;
            tmp1 = gamma + beta*sgn(z*delta_eps);
            
            // function and derivative
            f = z - zC - delta_eps / epsy*(1.0 - pow(zAbs, eta)*tmp1);
            Df = 1.0 + delta_eps / epsy*eta*pow(zAbs, eta - 1.0)*sgn(z)*tmp1;
            
            // issue warning if derivative Df is zero
            if (fabs(Df) <= DBL_EPSILON) {
                opserr << "WARNING: BoucWenOriginal::setTrialStrain() - "
                    << "zero derivative in Newton-Raphson scheme for "
                    << "hysteretic evolution parameter z.\n";
                return -1;
            }
            
            // advance one step
            delta_z = f / Df;
            z -= delta_z;
            iter++;
        } while ((fabs(delta_z) >= tol) && (iter < maxIter));
        
        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter) {
            opserr << "WARNING: BoucWenOriginal::setTrialStrain() - "
                << "did not find the hysteretic evolution parameter z after "
                << iter << " iterations and norm: " << fabs(delta_z) << endln;
            return -2;
        }
        
        // get derivative of hysteretic evolution parameter * epsy
        double dzdeps = 1.0 - pow(fabs(z), eta)*(gamma + beta*sgn(z*delta_eps));
        // set stress
        sig = qd*z + k2*eps + k3*sgn(eps)*pow(fabs(eps), mu);
        // set tangent stiffness
        Et = k0*dzdeps + k2 + k3*mu*pow(fabs(eps), mu - 1.0);
    }
    
    return 0;
}


double BoucWenOriginal::getStress()
{
    return sig;
}


double BoucWenOriginal::getInitialTangent()
{
    return Ei;
}


double BoucWenOriginal::getTangent()
{
    return Et;
}


double BoucWenOriginal::getStrain()
{
    return eps;
}


int BoucWenOriginal::commitState()
{
    // commit trial history variables
    epsC = eps;
    zC = z;
    
    return 0;
}


int BoucWenOriginal::revertToLastCommit()
{
    // revert trial history variables
    eps = epsC;
    z = zC;
    
    return 0;
}


int BoucWenOriginal::revertToStart()
{
    // reset trial history variables
    eps = 0.0;
    z = 0.0;
    sig = 0.0;
    
    // reset tangent stiffness
    Et = Ei;
    
    // reset committed history variables
    epsC = 0.0;
    zC = 0.0;
    
    return 0;
}


UniaxialMaterial *BoucWenOriginal::getCopy()
{
    BoucWenOriginal *theCopy =
        new BoucWenOriginal(this->getTag(), Ei, fy, alphaL, alphaNL,
            mu, eta, beta, gamma, tol, maxIter);
    
    theCopy->eps = eps;
    theCopy->z = z;
    theCopy->sig = sig;
    theCopy->Et = Et;
    theCopy->epsC = epsC;
    theCopy->zC = zC;
    
    return theCopy;
}


int BoucWenOriginal::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    
    static Vector data(13);
    data(0) = this->getTag();
    data(1) = Ei;
    data(2) = fy;
    data(3) = alphaL;
    data(4) = alphaNL;
    data(5) = mu;
    data(6) = eta;
    data(7) = beta;
    data(8) = gamma;
    data(9) = tol;
    data(10) = maxIter;
    data(11) = epsC;
    data(12) = zC;
    
    // Data is only sent after convergence, so no trial variables
    // need to be sent through data vector
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "BoucWenOriginal::sendSelf() - failed to send data\n";
    
    return res;
}


int BoucWenOriginal::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(11);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    
    if (res < 0) {
        opserr << "BoucWenOriginal::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        this->setTag(int(data(0)));
        Ei = data(1);
        fy = data(2);
        alphaL = data(3);
        alphaNL = data(4);
        mu = data(5);
        eta = data(6);
        beta = data(7);
        gamma = data(8);
        tol = data(9);
        maxIter = int(data(10));
        epsC = data(11);
        zC = data(12);

        this->revertToLastCommit();
    }
    
    return res;
}


void BoucWenOriginal::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "BoucWenOriginal, tag: " << this->getTag() << endln;
        s << "  E: " << Ei << endln;
        s << "  fy: " << fy << endln;
        s << "  alphaL: " << alphaL << endln;
        s << "  alphaNL: " << alphaNL << endln;
        s << "  mu: " << mu << endln;
        s << "  eta: " << eta << endln;
        s << "  beta: " << beta << endln;
        s << "  gamma: " << gamma << endln;
        s << "  tol: " << tol << endln;
        s << "  maxIter: " << maxIter << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"BoucWenOriginal\", ";
        s << "\"E\": " << Ei << ", ";
        s << "\"fy\": " << fy << ", ";
        s << "\"alphaL\": " << alphaL << ", ";
        s << "\"alphaNL\": " << alphaNL << ", ";
        s << "\"mu\": " << mu << ", ";
        s << "\"eta\": " << eta << ", ";
        s << "\"beta\": " << beta << ", ";
        s << "\"gamma\": " << gamma << ", ";
        s << "\"tol\": " << tol << ", ";
        s << "\"maxIter\": " << maxIter << "}";
    }
}


double BoucWenOriginal::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
