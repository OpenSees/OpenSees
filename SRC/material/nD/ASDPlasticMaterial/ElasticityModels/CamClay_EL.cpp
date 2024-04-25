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
                                                                        
// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial
//
// Fully general templated material class for plasticity modeling

#include "CamClay_EL.h"
#include "Vector.h"
#include "../ASDPlasticMaterialGlobals.h"

VoigtMatrix CamClay_EL::Ee(3, 3, 3, 3, 0.0);


CamClay_EL::CamClay_EL(double e0_, double kappa_, double nu_) : ElasticityBase<CamClay_EL>::ElasticityBase(),  // Note the full-qualification of ElasticityBase through the scope resolution operator (::)
    e0(e0_),
    kappa(kappa_),
    nu(nu_)
{

}


VoigtMatrix& CamClay_EL::operator()(const VoigtVector& stress) //See note on base class
{
    using namespace ASDPlasticMaterialGlobals;

    Ee *= 0; //Zero it. It may have values from another instance with different parameters;
    double p = -stress(i, i) / 3;
    double e = e0 - kappa * log(p);
    double K = (1 + e) * p / kappa;
    double lambda = 3 * K * nu / (1 + nu);
    double mu = 3 * K * (1 - 2 * nu) / (2 * (1 + nu));

    Ee( 0, 0, 0, 0 ) = lambda + 2 * mu;
    Ee( 0, 0, 1, 1 ) = lambda;
    Ee( 0, 0, 2, 2 ) = lambda;
    Ee( 0, 1, 0, 1 ) = mu;
    Ee( 0, 1, 1, 0 ) = mu;
    Ee( 0, 2, 0, 2 ) = mu;
    Ee( 0, 2, 2, 0 ) = mu;
    Ee( 1, 0, 0, 1 ) = mu;
    Ee( 1, 0, 1, 0 ) = mu;
    Ee( 1, 1, 0, 0 ) = lambda;
    Ee( 1, 1, 1, 1 ) = lambda + 2 * mu;
    Ee( 1, 1, 2, 2 ) = lambda;
    Ee( 1, 2, 1, 2 ) = mu;
    Ee( 1, 2, 2, 1 ) = mu;
    Ee( 2, 0, 0, 2 ) = mu;
    Ee( 2, 0, 2, 0 ) = mu;
    Ee( 2, 1, 1, 2 ) = mu;
    Ee( 2, 1, 2, 1 ) = mu;
    Ee( 2, 2, 0, 0 ) = lambda;
    Ee( 2, 2, 1, 1 ) = lambda;
    Ee( 2, 2, 2, 2 ) = lambda + 2 * mu;

    return Ee;
}

int CamClay_EL::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(3);
    data(0) = e0;
    data(1) = kappa;
    data(2) = nu;

    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "CamClay_EL::sendSelf() - Failed to send data. " << endl;
        return -1;
    }

    return 0;
}

int CamClay_EL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(3);

    if (theChannel.receiveVector(0, commitTag, data) != 0)
    {
        cerr << "CamClay_EL::recvSelf() - Failed to receive data. " << endl;
        return -1;
    }

    e0 = data(0);
    kappa = data(1);
    nu = data(2);

    return 0;
}
