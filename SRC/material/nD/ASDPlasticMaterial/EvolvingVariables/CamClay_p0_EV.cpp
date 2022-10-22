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

#include "CamClay_p0_EV.h"
#include "Vector.h"

double CamClay_p0_EV::derivative = 0.0;


CamClay_p0_EV::CamClay_p0_EV( double M_, double lambda_, double kappa_, double e0_, double p0_)
    :
    EvolvingVariable(p0_),
    M(M_),
    lambda(lambda_),
    kappa(kappa_),
    e0(e0_)
{}

const double& CamClay_p0_EV::getDerivative(const VoigtVector &depsilon,
        const VoigtVector &m,
        const VoigtVector& stress) const
{
    using namespace ASDPlasticMaterialGlobals;
    // Clear the static variables
    derivative = 0;

    // double e = e0 + (1.0 + e0) * stress(i, i);
    double p = -stress(i, i) / 3;
    double e = e0 - kappa * log(p);

    const double &p0 = this->getVariableConstReference();
    derivative = (1.0 + e) * p0 / (lambda - kappa) * m(i, i);
    return derivative ;
}

// double const& CamClay_p0_EV::getHardeningType() const
// {
//     return H;
// }


int CamClay_p0_EV::sendSelf(int commitTag, Channel &theChannel)
{
    //Shove all data into single vector for sending
    static Vector data(7);
    const double &p0 = this->getVariableConstReference();
    const double &p0_committed = this->getVariableConstReference();

    data(0) = M;
    data(1) = lambda;
    data(2) = kappa;
    data(3) = e0;
    data(5) = p0;
    data(6) = p0_committed;

    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "CamClay_p0_EV::sendSelf() - Failed sending data" << endl;
        return -1;
    }

    return 0;
}

int CamClay_p0_EV::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(7);
    if (theChannel.receiveVector(0, commitTag, data) != 0)
    {
        cerr << "CamClay_p0_EV::recvSelf() - Failed recieving data" << endl;
        return -1;
    }

    //Extract data from vector
    double tmp_p0;
    double tmp_p0_committed;
    M = data(0) ;
    lambda = data(1) ;
    kappa = data(2) ;
    e0 = data(3) ;
    tmp_p0 = data(5);
    tmp_p0_committed = data(6);

    this->setVar(tmp_p0);
    this->setCommittedVar(tmp_p0_committed);

    return 0;
}