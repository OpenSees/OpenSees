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


#ifndef Linear_HardeningScalarV_H
#define Linear_HardeningScalarV_H


#include "../EvolvingVariable.h"
#include "../ASDPlasticMaterialGlobals.h" // Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.


class LinearHardeningScalar_EV : public EvolvingVariable<VoigtScalar, LinearHardeningScalar_EV> //CRTP on LinearHardeningScalar_EV
{
public:

    LinearHardeningScalar_EV() : EvolvingVariable(0.0), H(0.0);

    LinearHardeningScalar_EV( VoigtScalar H_) : EvolvingVariable(0.0), H(H_) {};

    LinearHardeningScalar_EV( VoigtScalar H_, VoigtScalar k0) : EvolvingVariable(k0), H(H_) {};

    const VoigtScalar& getDerivative(const VoigtVector &depsilon,
                                const VoigtVector &m,
                                const VoigtVector& stress) const
    {
        derivative = H * sqrt((2 * m.dot(m)) / 3);
        return derivative;
    }


private:
    VoigtScalar H;
    static VoigtScalar derivative; //Must return a reference.
};

 VoigtScalar LinearHardeningScalar_EV::derivative;

#endif //Linear_HardeningScalarV_H