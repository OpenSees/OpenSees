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
template <class ParameterStorageType>
struct LinearHardeningForTensorPolicy {
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

template <class ParameterStorageType>
struct LinearHardeningForScalarPolicy {
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
template <typename T, template <class> class HardeningPolicy>
struct HardeningFunction {
	template< class ParameterStorageType>
    static T f(
        const VoigtVector& depsilon,
        const VoigtVector& m,
        const VoigtVector& sigma,
        const ParameterStorageType& parameters) 
    {
        return HardeningPolicy<ParameterStorageType>::f(depsilon, m, sigma, parameters);
    }

	template< class ParameterStorageType>
    using parameters_t = typename HardeningPolicy<ParameterStorageType>::parameters_t;
};


// Aliases for HardeningFunction with specific hardening policies
using TensorLinearHardeningFunction = HardeningFunction<VoigtVector, LinearHardeningForTensorPolicy>;
using ScalarLinearHardeningFunction = HardeningFunction<VoigtScalar, LinearHardeningForScalarPolicy>;


template <typename T, template <class> class HardeningPolicy>
std::ostream& operator<<(std::ostream& os, const HardeningFunction<T, HardeningPolicy>& obj) {
    os << "HardeningFunction<" << typeid(T).name() << ", " << typeid(HardeningPolicy<void>).name() << ">";
    return os;
}



// // Function wrapper base class
// template <typename T, class ParameterStorageType>
// struct HardeningFunctionWrapper {
//     static T hardening(
//         const T& variable,
//         const VoigtVector &depsilon,
//         const VoigtVector &m,
//         const VoigtVector& sigma,
//         const ParameterStorageType& parameters);
// };

// // Template specialization for VoigtVector and ParameterStorageType
// template <class ParameterStorageType>
// struct HardeningFunctionWrapper<VoigtVector, ParameterStorageType> {
//     static VoigtVector hardening(
//         const VoigtVector& variable,
//         const VoigtVector& depsilon,
//         const VoigtVector& m,
//         const VoigtVector& sigma,
//         const ParameterStorageType& parameters) 
//     {
//         double H = parameters.template get<LinearHardeningForTensor>().value;
//         auto derivative = H * m.deviator();
//         return derivative;
//     }

//     using parameters_t = std::tuple<LinearHardeningForTensor>;
// };

// // Alias for HardeningFunctionWrapper with VoigtVector and specific ParameterStorageType
// template <class ParameterStorageType>
// using VoigtVectorHardening = HardeningFunctionWrapper<VoigtVector, ParameterStorageType>;


// // Template specialization for VoigtScalar and ParameterStorageType
// template <class ParameterStorageType>
// struct HardeningFunctionWrapper<VoigtScalar, ParameterStorageType> {
//     static VoigtScalar hardening(
//         const VoigtScalar& variable,
//         const VoigtVector& depsilon,
//         const VoigtVector& m,
//         const VoigtVector& sigma,
//         const ParameterStorageType& parameters) 
//     {
//         double H = parameters.template get<LinearHardeningForScalar>().value;
//         auto derivative = H * sqrt((2 * m.dot(m)) / 3);
//         return derivative;
//     }

//     using parameters_t = std::tuple<LinearHardeningForScalar>;
// };

// // Alias for HardeningFunctionWrapper with VoigtScalar and specific ParameterStorageType
// template <class ParameterStorageType>
// using ScalarLinearHardening = HardeningFunctionWrapper<VoigtScalar, ParameterStorageType>;


// // Function wrapper base class
// template <typename T, class ParameterStorageType>
// struct HardeningFunctionWrapper {
//     static T hardening(
//                 const T& variable,
//                 const VoigtVector &depsilon,
//                 const VoigtVector &m,
//                 const VoigtVector& sigma,
//                 const ParameterStorageType& parameters);
// };

// template <class ParameterStorageType>
// struct HardeningFunctionWrapper<VoigtVector, ParameterStorageType> {
//     static VoigtVector function(const VoigtVector& variable,
//                 const VoigtVector &depsilon,
//                 const VoigtVector &m,
//                 const VoigtVector& sigma,
//                 const ParameterStorageType& parameters) 
//     {
//         double H = parameters.template get<LinearHardeningForTensor> ().value;
//         auto derivative = H * m.deviator();
//         return derivative;
//     }

//     using parameters_t = std::tuple<LinearHardeningForTensor>;
// } ;

// template <class ParameterStorageType>
// struct HardeningFunctionWrapper<VoigtScalar, ParameterStorageType> {
//     static VoigtScalar function(const VoigtScalar& variable,
//                 const VoigtVector &depsilon,
//                 const VoigtVector &m,
//                 const VoigtVector& sigma,
//                 const ParameterStorageType& parameters) 
//     {
//         double H = parameters.template get<LinearHardeningForScalar> ().value;
//         auto derivative = H * sqrt((2 * m.dot(m)) / 3);
//         return derivative;
//     }

//     using parameters_t = std::tuple<LinearHardeningForScalar>;
// } ;


#endif //not defined _AllASDHardeningFunctions