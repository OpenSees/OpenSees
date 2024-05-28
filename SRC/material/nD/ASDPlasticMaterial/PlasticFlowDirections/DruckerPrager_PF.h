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

#ifndef DruckerPrager_PF_H
#define DruckerPrager_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterialGlobals.h"

#include <cmath>
#include <typeinfo>


template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_PF : public PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_PF";


    DruckerPrager_PF( ):
        PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, KHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        double p = -sigma.meanStress();
        auto s = sigma.deviator();
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        auto r = s / p;

        double den = (SQRT_2_over_3 * k).value();
        auto n = (r - alpha) / den;
        double nr = n.dot(r);
        result = n - nr * kronecker_delta() / 3;

        return result;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;
    using parameters_t = std::tuple<>;

private:

    static VoigtVector result; //For returning VoigtVectors

};


template<class AlphaHardeningType, class KHardeningType>
VoigtVector DruckerPrager_PF<AlphaHardeningType, KHardeningType  >::result;

#endif
