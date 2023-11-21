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

#ifndef DruckerPrager_YF_H
#define DruckerPrager_YF_H

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>


template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_YF";


    DruckerPrager_YF( ):
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase() // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
        {}

    YIELD_FUNCTION 
    {
        double p = -sigma.meanStress();
        auto s = sigma.deviator();
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        double tmp = (s - p*alpha).dot(s - p*alpha);

        double yf = sqrt( tmp ) - (SQRT_2_over_3 * k * p).value(); // This one assumes p positive in tension

        return yf;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
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

    YIELD_FUNCTION_XI_STAR_H_STAR
    {
        double dbl_result = 0.0;

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);
        
        double p = -sigma.meanStress();
        auto s = sigma.deviator();

        double den = (SQRT_2_over_3 * k).value();

        // This is for the hardening of k
        double df_dk = -SQRT_2_over_3 * p;
        dbl_result +=  (df_dk * GET_INTERNAL_VARIABLE_HARDENING(KHardeningType)).value();

        //This is for the hardening of alpha
        auto df_dalpha = ( p * alpha - s) / den;
        dbl_result +=  df_dalpha.dot(GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType));

        return dbl_result;
    }

  
    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;

    using parameters_t = std::tuple<>;

private:

    static VoigtVector result; //For returning VoigtVector's
};

template <class AlphaHardeningType,  class KHardeningType>
VoigtVector DruckerPrager_YF<AlphaHardeningType, KHardeningType>::result;

//Declares this YF as featuring an apex
template<class AlphaHardeningType, class KHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> : std::true_type {};

#endif