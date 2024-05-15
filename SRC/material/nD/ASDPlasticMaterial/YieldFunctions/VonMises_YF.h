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

#ifndef VonMises_YF_H
#define VonMises_YF_H

// #include "../EvolvingVariable.h"
#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

#include "../ASDPlasticMaterialGlobals.h"
using namespace ASDPlasticMaterialGlobals;




template<class AlphaHardeningType, class KHardeningType>
class VonMises_YF : public YieldFunctionBase<VonMises_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "VonMises_YF";

    VonMises_YF( ):
        YieldFunctionBase<VonMises_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase() // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                {}

    YIELD_FUNCTION
    {
        auto s = sigma.deviator();

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        double tmp = (s - alpha).dot(s - alpha);
        tmp = tmp > 0 ? tmp : 0;
        return std::sqrt( tmp ) - SQRT_2_over_3 * k.value() ;  // This one assumes p positive in tension
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        result = sigma.deviator() - alpha;

        double den = sqrt(result.dot(result));
        if (abs(den) > 100*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
            result = result / den;

        return result;

    }
    
    YIELD_FUNCTION_XI_STAR_H_STAR
    {
        double dbl_result = 0.0;

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        //Zero the stress deviator
        auto s = sigma.deviator();

        // This is for the hardening of k
        double df_dk = -SQRT_2_over_3;
        dbl_result +=  (df_dk * GET_INTERNAL_VARIABLE_HARDENING(KHardeningType)).value();

        //This is for the hardening of alpha
        double den = sqrt((s - alpha).dot(s - alpha));

        if (abs(den) < sqrt(s.dot(s))*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
        {
            return dbl_result;
        }

        auto df_dalpha = -(s - alpha) / den;
        dbl_result +=  df_dalpha.dot(GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType));

        return dbl_result;
    }

    bool hasCorner() const{
        return false;
    }

    bool in_Apex(VoigtVector const& TrialStress)
    {
        std::cout<<"von Mises yield surface does not have a corner. This function should never be callled!"<<std::endl;
        return false;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;

    using parameters_t = std::tuple<>;


private:

    static VoigtVector result; //For returning VoigtVector's

};

template <class AlphaHardeningType,  class KHardeningType>
VoigtVector VonMises_YF<AlphaHardeningType, KHardeningType>::result;


#endif