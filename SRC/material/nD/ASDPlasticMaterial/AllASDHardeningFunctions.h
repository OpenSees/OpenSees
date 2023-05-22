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


#ifndef _AllASDHardeningFunctions
#define _AllASDHardeningFunctions

#include "ASDPlasticMaterialGlobals.h" 
#include "AllASDModelParameterTypes.h" 




// Hardening policies
struct LinearHardeningForTensorPolicy {
    static constexpr const char* NAME = "TensorLinearHardeningFunction";
	template <class ParameterStorageType>
    static VoigtVector f(
    	const VoigtVector& depsilon,
        const VoigtVector& m,
        const VoigtVector& sigma,
        const ParameterStorageType& parameters) {
        double H = parameters.template get<TensorLinearHardeningParameter>().value;
        auto h = H * m.deviator();
        return h;
    }
    using parameters_t = tuple<TensorLinearHardeningParameter>;
};

struct LinearHardeningForScalarPolicy {
    static constexpr const char* NAME = "ScalarLinearHardeningFunction";
	template <class ParameterStorageType>
    static double f(
    	const VoigtVector& depsilon,
        const VoigtVector& m,
        const VoigtVector& sigma,
        const ParameterStorageType& parameters) {
        double H = parameters.template get<ScalarLinearHardeningParameter>().value;
        double h = H * sqrt((2 * m.dot(m)) / 3);
        return h;
    }
    using parameters_t = tuple<ScalarLinearHardeningParameter>;
};

// Function wrapper base class
template <typename T, class HardeningPolicy>
struct HardeningFunction {
	template< class ParameterStorageType>
    static T f(
        const VoigtVector& depsilon,
        const VoigtVector& m,
        const VoigtVector& sigma,
        const ParameterStorageType& parameters) 
    {
        return HardeningPolicy::f(depsilon, m, sigma, parameters);
    }
    static constexpr const char* NAME = HardeningPolicy::NAME;
    using parameters_t = typename HardeningPolicy::parameters_t;
};


// Aliases for HardeningFunction with specific hardening policies
using TensorLinearHardeningFunction = HardeningFunction<VoigtVector, LinearHardeningForTensorPolicy>;
using ScalarLinearHardeningFunction = HardeningFunction<VoigtScalar, LinearHardeningForScalarPolicy>;


template <typename T, class HardeningPolicy>
std::ostream& operator<<(std::ostream& os, const HardeningFunction<T, HardeningPolicy>& obj) {
    os << "HardeningFunction<" << typeid(T).name() << ", " << typeid(HardeningPolicy).name() << ">";
    return os;
}


#endif //not defined _AllASDHardeningFunctions