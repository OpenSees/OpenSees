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
#include "../ASDPlasticMaterialGlobals.h"
#include "../AllASDModelParameterTypes.h"


class LinearHardeningScalar_EV : public EvolvingVariable<VoigtScalar, LinearHardeningScalar_EV> //CRTP on LinearHardeningScalar_EV
{
public:

    LinearHardeningScalar_EV() : EvolvingVariable(0.0) {}
    LinearHardeningScalar_EV(VoigtScalar k0) : EvolvingVariable(k0) {};

    template<class ParameterStorageType>
    const VoigtScalar& getDerivative(const VoigtVector &depsilon,
                                const VoigtVector &m,
                                const VoigtVector& stress,
        const ParameterStorageType& parameters_storage) const
    {
        double H = parameters_storage.template get<LinearHardeningForScalar> ().value;
        derivative = H * sqrt((2 * m.dot(m)) / 3);
        return derivative;
    }

    using parameters_t = std::tuple<LinearHardeningForScalar>;

private:
    static VoigtScalar derivative; //Must return a reference.
};

 VoigtScalar LinearHardeningScalar_EV::derivative;

#endif //Linear_HardeningScalarV_H