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


#include "LinearHardeningTensor_EV.h"
#include <Vector.h>

VoigtVector LinearHardeningTensor_EV::derivative(3, 3, 0.0);

LinearHardeningTensor_EV::LinearHardeningTensor_EV( double H_) : EvolvingVariable(VoigtVector(3, 3, 0.0)), H(H_)
{}

LinearHardeningTensor_EV::LinearHardeningTensor_EV( double H_, VoigtVector& alpha0) : EvolvingVariable(alpha0), H(H_)
{}

const VoigtVector& LinearHardeningTensor_EV::getDerivative(const VoigtVector &depsilon,
        const VoigtVector &m,
        const VoigtVector& stress) const
{
    using namespace ASDPlasticMaterialGlobals;
    //Zero de static variable
    derivative *= 0;

    //Compute the derivative (hardening function)
    double mkk = - (m(k, k) / 3);
    derivative(i, j) =  H * (m(i, j) + mkk * kronecker_delta(i, j));

    return derivative;
}


int LinearHardeningTensor_EV::sendSelf(int commitTag, Channel &theChannel)
{
    //Shove all data into single vector for sending
    static Vector data(9 + 9 + 1);
    const VoigtVector &a = this->getVariableConstReference();
    const VoigtVector &a_committed = this->getVariableConstReference();
    int pos = 0;

    data(pos++) = H;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            data(pos++) = a(i, j);
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            data(pos++) = a_committed(i, j);
        }
    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "LinearHardeningTensor_EV::sendSelf() - Failed sending data" << endl;
        return -1;
    }

    return 0;
}

int LinearHardeningTensor_EV::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(9 + 9 + 1);
    if (theChannel.receiveVector(0, commitTag, data) != 0)
    {
        cerr << "LinearHardeningTensor_EV::recvSelf() - Failed recieving data" << endl;
        return -1;
    }

    //Extract data from vector
    int pos = 0;
    static VoigtVector tmp_a(3, 3, 0);
    static VoigtVector tmp_a_committed(3, 3, 0);
    H = data(pos++);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp_a(i, j) = data(pos++);
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp_a_committed(i, j) = data(pos++) ;
        }

    this->setVar(tmp_a);
    this->setCommittedVar(tmp_a_committed);

    return 0;
}