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

#include "LinearHardeningScalar_EV.h"
#include "Vector.h"

double LinearHardeningScalar_EV::derivative = 0.0;

LinearHardeningScalar_EV::LinearHardeningScalar_EV( double H_) : EvolvingVariable(0.0), H(H_)
{}

LinearHardeningScalar_EV::LinearHardeningScalar_EV( double H_, double k0) : EvolvingVariable(k0), H(H_)
{}

const double& LinearHardeningScalar_EV::getDerivative(const VoightTensor6 &depsilon,
        const VoightTensor6 &m,
        const VoightTensor6& stress) const
{
    using namespace ASDPlasticMaterialGlobals;
    // Clear the static variables
    derivative = 0;

    //Now compute the equivalent m
    double m_eq = sqrt((2 * m.dot(m)) / 3);

    //Compute the derivative.
    derivative = H * m_eq;
    return derivative ;
}


int LinearHardeningScalar_EV::sendSelf(int commitTag, Channel &theChannel)
{
    //Shove all data into single vector for sending
    static Vector data(3);
    const double &a = this->getVariableConstReference();
    const double &a_committed = this->getVariableConstReference();

    data(0) = H;
    data(1) = a;
    data(2) = a_committed;

    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "LinearHardeningScalar_EV::sendSelf() - Failed sending data" << endl;
        return -1;
    }

    return 0;
}

int LinearHardeningScalar_EV::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(3);
    if (theChannel.recvVector(0, commitTag, data) != 0)
    {
        cerr << "LinearHardeningScalar_EV::recvSelf() - Failed recieving data" << endl;
        return -1;
    }

    //Extract data from vector
    double tmp_a;
    double tmp_a_committed;
    H = data(0);
    tmp_a = data(1);
    tmp_a_committed = data(2);

    this->setVar(tmp_a);
    this->setCommittedVar(tmp_a_committed);

    return 0;
}