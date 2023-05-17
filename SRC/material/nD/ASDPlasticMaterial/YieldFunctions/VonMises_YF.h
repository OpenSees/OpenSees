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

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
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

    template <typename IVStorageType, typename ParameterStorageType>
    double operator()( const VoigtVector& sigma, 
    	const IVStorageType& internal_variables_storage,
        const ParameterStorageType& parameters_storage) const 
    {
        auto s = sigma.deviator();

        // auto alpha = alpha_.getVariableConstReference();
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();
        auto alpha = AHT.trial_value;

        // auto k = k_.getVariableConstReference();
        const KHardeningType& KHT = 
        internal_variables_storage.template get<KHardeningType> ();
        auto k = KHT.trial_value;

        double tmp = (s - alpha).dot(s - alpha);

        // printTensor("s", s);
        // printTensor("alpha", alpha);
        // // printTensor("k", k);
        // cout << "k = " << k(0) << endl;

        return std::sqrt( tmp ) - SQRT_2_over_3 * k.value() ;  // This one assumes p positive in tension
    }

    template <typename IVStorageType, typename ParameterStorageType>
    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma,
    	const IVStorageType& internal_variables_storage,
        const ParameterStorageType& parameters_storage)
    {
        // VoigtVector alpha = alpha_.getVariable();
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();
        auto alpha = AHT.trial_value;

        result = sigma.deviator() - alpha;

        double den = sqrt(result.dot(result));
        if (abs(den) > 100*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
            result = result / den;

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

        //Zero the stress deviator
        auto s = sigma.deviator();

        //Denominator of the expression

        // This is for the hardening of k
        constexpr double df_dk = -SQRT_2_over_3;
        dbl_result +=  (df_dk * KHT.hardening_function(depsilon, m, sigma, parameters_storage)).value();

        //This is for the hardening of alpha
        double den = sqrt((s - alpha).dot(s - alpha));

        if (abs(den) < 100*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
        {
            return dbl_result;
        }

        auto df_dalpha = -(s - alpha) / den;
        dbl_result +=  df_dalpha.dot(AHT.hardening_function(depsilon, m, sigma, parameters_storage));

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

    // using iv_parameters_t = std_tuple_concat_Type <
    //                         typename AlphaHardeningType::parameters_t,
    //                         typename KHardeningType::parameters_t >;

    // using more_parameters_t = std::tuple<>;

    // using parameters_t = std_tuple_concat_Type<iv_parameters_t, more_parameters_t>;

    using parameters_t = std::tuple<>;


private:

    static VoigtVector s; //Stress deviator
    static VoigtVector result; //For returning VoigtVector's

};

template <class AlphaHardeningType,  class KHardeningType>
VoigtVector VonMises_YF<AlphaHardeningType, KHardeningType>::s;
template <class AlphaHardeningType,  class KHardeningType>
VoigtVector VonMises_YF<AlphaHardeningType, KHardeningType>::result;


#endif