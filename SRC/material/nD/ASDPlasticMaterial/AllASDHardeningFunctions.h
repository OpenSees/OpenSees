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
#include "HardeningFunction.h" 


// Hardening policies
struct LinearHardeningForTensorPolicy {
    static constexpr const char* NAME = "TensorLinearHardeningFunction";
    HARDENING_FUNCTION_DEFINITION
    {
        double H = GET_PARAMETER_VALUE(TensorLinearHardeningParameter);
        auto h = H * m.deviator();
        return h;
    }
    using parameters_t = tuple<TensorLinearHardeningParameter>;
};

struct LinearHardeningForScalarPolicy {
    static constexpr const char* NAME = "ScalarLinearHardeningFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);
        double h = H * sqrt((2 * m.dot(m)) / 3);
        return h;
    }
    using parameters_t = tuple<ScalarLinearHardeningParameter>;
};

// struct ScalarExponentialLinear {
//     static constexpr const char* NAME = "ScalarExponentialLinear";
//     HARDENING_FUNCTION_DEFINITION 
//     {
//         double Sigma0 = GET_PARAMETER_VALUE(ScalarExponentialLinear_Sigma0);
//         double SigmaInf = GET_PARAMETER_VALUE(ScalarExponentialLinear_SigmaInf);
//         double delta = GET_PARAMETER_VALUE(ScalarExponentialLinear_delta);
//         double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);

//         // This is the same as the hardening of the "q(xi)" in terms of the internal variable "xi" in J2Plasticity,
//         // xi is root23 accumulated plastic which we here compute
//         double dxi = QRT_2_over_3 * dLambda;
        
//         // Because of how ASDPlasticMaterial works, instead of using the explicit form of 
//         // q(xi) used in J2 plasticity, we here use the incremental form. Thus, this
//         // function returns dq / dxi 
//         // double dq_over_dxi = -delta * (SigmaInf - Sigma0) * exp(-delta*)
        
//     }
//     using parameters_t = tuple<ScalarExponentialLinear_Sigma0,ScalarExponentialLinear_SigmaInf, ScalarExponentialLinear_delta, ScalarLinearHardeningParameter>;
// };

struct ArmstrongFrederickPolicy {
    static constexpr const char* NAME = "ArmstrongFrederickHardeningFunction";
    HARDENING_FUNCTION_DEFINITION
    {
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);

        auto alpha = current_value;


        auto mdev = m.deviator();
        //Compute the derivative (hardening function)
        // const VoigtVector &alpha = this->getVariableConstReference();
        // static VoigtVector mdev(3, 3, 0);
        // mdev *= 0;
        // mdev(i, j) = m(i, j) - m(k, k) / 3 * kronecker_delta(i, j);
        // derivative(i, j) =  (2. / 3.) * ha * m(i, j) - cr * sqrt((2. / 3.) * m(k, l) * m(k, l)) * alpha(i, j);
        auto derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha;

        return derivative;

    }
    using parameters_t = tuple<AF_ha, AF_cr>;
};


// Aliases for HardeningFunction with specific hardening policies
using TensorLinearHardeningFunction = HardeningFunction<VoigtVector, LinearHardeningForTensorPolicy>;
using ScalarLinearHardeningFunction = HardeningFunction<VoigtScalar, LinearHardeningForScalarPolicy>;
using ArmstrongFrederickHardeningFunction = HardeningFunction<VoigtVector, ArmstrongFrederickPolicy>;


#endif //not defined _AllASDHardeningFunctions
