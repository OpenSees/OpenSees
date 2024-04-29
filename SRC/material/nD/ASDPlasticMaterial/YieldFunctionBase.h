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

#ifndef YieldFunctionBase_H
#define YieldFunctionBase_H

#include "EigenAPI.h"
#include <Channel.h>
#include <type_traits>


// Trait for enabling special algorithms in YF that have apexes

// Base template assuming no apex member exists
template <typename T>
struct yf_has_apex : std::false_type {};



// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct yf_has_parameters_t : std::false_type {};

template <typename T>
struct yf_has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

// Helper template to check if a class has a internal_variables_t type alias
template <typename T, typename = void>
struct yf_has_internal_variables_t : std::false_type {};

template <typename T>
struct yf_has_internal_variables_t<T, typename std::enable_if<!std::is_same<typename T::internal_variables_t, void>::value>::type> : std::true_type {};


#define YIELD_FUNCTION template <typename IVStorageType, typename ParameterStorageType> \
    double operator()( const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage) const 

#define YIELD_FUNCTION_STRESS_DERIVATIVE template <typename IVStorageType, typename ParameterStorageType> \
    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage)

#define YIELD_FUNCTION_XI_STAR_H_STAR template <typename IVStorageType, typename ParameterStorageType> \
    double xi_star_h_star(const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma,\
        const IVStorageType& internal_variables_storage,\
        const ParameterStorageType& parameters_storage)

#define GET_INTERNAL_VARIABLE_HARDENING(type) \
    internal_variables_storage.template get<type> ().hardening_function(depsilon, m, sigma, parameters_storage)

#define GET_TRIAL_INTERNAL_VARIABLE(type) \
    ((internal_variables_storage).template get<type>().trial_value)

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

template <class T>
class YieldFunctionBase
{
public:
    YieldFunctionBase() { 
        static_assert(yf_has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
        static_assert(yf_has_internal_variables_t<T>::value, "Derived class must have a 'internal_variables_t' type alias.");
    }

    YIELD_FUNCTION
    {
        return static_cast<T*>(this)->operator()(sigma, internal_variables_storage, parameters_storage);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        return static_cast<T*>(this)->df_dsigma_ij(sigma, internal_variables_storage, parameters_storage);
    }

    YIELD_FUNCTION_XI_STAR_H_STAR
    {
        return static_cast<T*>(this)->df_dxi_star_h_star(depsilon, m , sigma, internal_variables_storage, parameters_storage);
    }

    inline const char* getName() const { return static_cast<T*>(this)->NAME; }

};

#endif
