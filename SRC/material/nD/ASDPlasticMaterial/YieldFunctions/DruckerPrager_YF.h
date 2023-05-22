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

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ASDPlasticMaterialGlobals.h"
using namespace ASDPlasticMaterialGlobals;




template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_YF";


    DruckerPrager_YF( ):
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase() // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
        {}

    template <typename IVStorageType, typename ParameterStorageType>
    double operator()( const VoigtVector& sigma, 
        const IVStorageType& internal_variables_storage,
        const ParameterStorageType& parameters_storage) const 
    {
        double p = -sigma.meanStress();
        auto s = sigma.deviator();
        
        // auto alpha = alpha_.getVariableConstReference();
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();
        auto alpha = AHT.trial_value;

        // auto k = k_.getVariableConstReference();
        const KHardeningType& KHT = 
        internal_variables_storage.template get<KHardeningType> ();
        auto k = KHT.trial_value;

        double tmp = (s - p*alpha).dot(s - p*alpha);

        double yf = sqrt( tmp ) - (SQRT_2_over_3 * k * p).value(); // This one assumes p positive in tension

        return yf;
    }

    template <typename IVStorageType, typename ParameterStorageType>
    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma,
        const IVStorageType& internal_variables_storage,
        const ParameterStorageType& parameters_storage)
    {  


        double p = -sigma.meanStress();
        auto s = sigma.deviator();
        
        // auto alpha = alpha_.getVariableConstReference();
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();
        auto alpha = AHT.trial_value;

        // auto k = k_.getVariableConstReference();
        const KHardeningType& KHT = 
        internal_variables_storage.template get<KHardeningType> ();
        auto k = KHT.trial_value;

        auto r = s / p;

        double den = (SQRT_2_over_3 * k).value();
        auto n = (r - alpha) / den;
        double nr = n.dot(r);
        result = n - nr * kronecker_delta() / 3;

        return result;
    }

    template <typename IVStorageType, typename ParameterStorageType>
    double xi_star_h_star(const VoigtVector& depsilon, 
        const VoigtVector& m, 
        const VoigtVector& sigma,
        const IVStorageType& internal_variables_storage,
        const ParameterStorageType& parameters_storage)
    {
        double dbl_result = 0.0;

        // const VoigtVector &alpha = alpha_.getVariableConstReference();
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();
        auto alpha = AHT.trial_value;

        // const VoigtScalar &k = k_.getVariableConstReference();
        const KHardeningType& KHT = 
        internal_variables_storage.template get<KHardeningType> ();
        auto k = KHT.trial_value;
        
        double p = -sigma.meanStress();
        auto s = sigma.deviator();

        double den = (SQRT_2_over_3 * k).value();

        // This is for the hardening of k
        double df_dk = -SQRT_2_over_3 * p;
        dbl_result +=  (df_dk * KHT.hardening_function(depsilon, m, sigma, parameters_storage)).value();

        //This is for the hardening of alpha
        auto df_dalpha = ( p * alpha - s) / den;
        dbl_result +=  df_dalpha.dot(AHT.hardening_function(depsilon, m, sigma, parameters_storage));

        return dbl_result;
    }

  
    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;

    using parameters_t = std::tuple<>;

private:

    static VoigtVector result; //For returning VoigtVector's
};

template <class AlphaHardeningType,  class KHardeningType>
VoigtVector DruckerPrager_YF<AlphaHardeningType, KHardeningType>::result;






  // bool hasCorner() const
    // {
    //     return true;
    // }

    // bool in_Apex(VoigtVector const& TrialStress)
    // {
    //     using namespace ASDPlasticMaterialGlobals;
    //     double I1 = TrialStress(i, i);
    //     double sigma_m  = I1 / 3.0;
    //     double DP_p = -sigma_m;
    //     static VoigtVector DP_s(3, 3, 0.0);
    //     DP_s(i, j) = TrialStress(i, j) - kronecker_delta(i, j) * sigma_m;
    //     double DP_k = k_.getVariableConstReference();
    //     static VoigtVector DP_alpha(3, 3, 0.0);
    //     DP_alpha    =  alpha_.getVariableConstReference();
    //     static VoigtVector relative_s(3, 3, 0.0);
    //     // Remove the influence of kinematic hardening:
    //     relative_s(i, j) = DP_s(i, j) - DP_p * DP_alpha(i, j);
    //     double J2   = 0.5 * relative_s(i, j) * relative_s(i, j) ;
    //     double DP_q = sqrt(3 * J2);

    //     if (DP_k < MACHINE_EPSILON)
    //     {
    //         cout << "DP_k (denominator) = 0! Error in ASDPlasticMaterial-->requires_DruckerPrager_Apex_Tension_check_\n ";
    //     }

    //     // In this range (condition 1 && 2), DP stress state cannot be returned correctly.
    //     // So the the trick (return to small_stress) will be applied.
    //     bool condition1 = (  (1.0 / DP_k * DP_p + DP_q) < 0  );
    //     bool condition2 = (  (1.0 / DP_k * DP_p - DP_q) < 0  );

    //     if (condition1 && condition2)
    //     {
    //         // cout<<"PredictorStress in DruckerPrager apex region. Return to a small stress.\n";
    //         return true;
    //     }
    //     return false;
    // }

    // double get_k() const
    // {
    //     return k_.getVariableConstReference();
    // }
    // VoigtVector const& get_alpha() const
    // {
    //     return alpha_.getVariableConstReference();
    // }
#endif