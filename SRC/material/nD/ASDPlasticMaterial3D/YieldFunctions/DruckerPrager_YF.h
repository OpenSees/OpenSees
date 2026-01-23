#ifndef DruckerPrager_YF_H
#define DruckerPrager_YF_H

#include "../YieldFunctionBase.h"
#include <cmath>
#include <iostream>
#include "../AllASDModelParameterTypes.h"

template<class AlphaHardeningType, class CohesionHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_YF";

    DruckerPrager_YF():
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>>::YieldFunctionBase()
    {}

    YIELD_FUNCTION
    {
        auto s = sigma.deviator();
        double p = sigma.meanStress();  // mean stress (positive in compression in geomechanics convention)
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        // auto eta = GET_TRIAL_INTERNAL_VARIABLE(CohesionHardeningType);
        
        // Get friction parameter (eta) and cohesion (c) from parameters
        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);      // cohesion
        double eta = GET_PARAMETER_VALUE(DP_eta);      // friction slope
        
        double tmp = tensor_dot_stress_like(s - alpha, s - alpha);
        tmp = tmp > 0 ? tmp : 0;
        double sqrt_J2 = std::sqrt(0.5*tmp);
        
        // Drucker-Prager yield function: sqrt(J2) + eta * p - (c + k)
        // where:
        // - 0.5*sqrt(J2) is the deviatoric stress magnitude
        // - eta is the friction parameter 
        // - p is the mean stress (positive in compression)
        // - xi_c is the adjusted cohesion
        return sqrt_J2 + eta * p - (xi_c);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        auto s = sigma.deviator();
        double p = sigma.meanStress();
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);

        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);      // cohesion
        double eta = GET_PARAMETER_VALUE(DP_eta);      // friction slope

        // Derivative with respect to deviatoric stress
        VoigtVector dev_part = s - alpha;
        double den = sqrt(0.5*tensor_dot_stress_like(dev_part, dev_part));
        
        if (abs(den) > 100*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            dev_part = dev_part / den;
        else
            dev_part *= 0.0;  // Zero out if denominator too small
            
        // Add pressure-dependent part: eta * dp/dsigma = eta/3 * I
        VoigtVector pressure_part;
        pressure_part *= 0.0;  // Initialize to zero
        pressure_part(0) = eta / 3.0;  // sigma_xx component
        pressure_part(1) = eta / 3.0;  // sigma_yy component  
        pressure_part(2) = eta / 3.0;  // sigma_zz component
        // shear components remain zero
        
        vv_out = dev_part + pressure_part;
        
        return vv_out;
    }
    
    YIELD_FUNCTION_HARDENING
    {
        double dbl_result = 0.0;
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto eta = GET_TRIAL_INTERNAL_VARIABLE(CohesionHardeningType);
        
        auto s = sigma.deviator();
        
        // Hardening contribution from eta (isotropic hardening)
        double df_deta = -1.0;  // derivative of f with respect to eta
        dbl_result += (df_deta * GET_INTERNAL_VARIABLE_HARDENING(CohesionHardeningType)).value();
        
        // Hardening contribution from alpha (kinematic hardening)
        double den = sqrt(0.5*tensor_dot_stress_like(s - alpha, s - alpha));
        
        if (abs(den) < sqrt(0.5*tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
        {
            return dbl_result;
        }
        
        auto df_dalpha = -(s - alpha) / den;
        VoigtVector hh = GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType);
        dbl_result += tensor_dot_stress_like(df_dalpha, hh);
        
        return dbl_result;
    }


    CHECK_APEX_REGION
    {

        // Implement!!! 

        return false;
    }

    APEX_STRESS
    {

        // Implement!!! 
        vv_out = VoigtVector(0,0,0,0,0,0);

        return vv_out;

    }


    using internal_variables_t = std::tuple<AlphaHardeningType, CohesionHardeningType>;
    using parameters_t         = std::tuple<DP_xi_c, DP_eta>;

private:
    static VoigtVector vv_out; // For returning VoigtVectors
};

// Static member definition
template <class AlphaHardeningType, class CohesionHardeningType>
VoigtVector DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>::vv_out;

// Declares this YF as featuring an apex
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};

#endif
