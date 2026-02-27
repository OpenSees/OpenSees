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
// ASDPlasticMaterial3D
//
// Fully general templated material class for plasticity modeling


#ifndef _AllASDHardeningFunctions
#define _AllASDHardeningFunctions

#include "ASDPlasticMaterial3DGlobals.h" 
#include "AllASDModelParameterTypes.h" 
#include "HardeningFunction.h" 


// Hardening policies
struct LinearHardeningForTensorPolicy {
    static constexpr const char* NAME = "TensorLinearHardeningFunction";
    HARDENING_FUNCTION_DEFINITION
    {
        double H = GET_PARAMETER_VALUE(TensorLinearHardeningParameter);
        VoigtVector h = H * m.deviator();  // best not to use 'auto' here
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

struct NullHardeningScalarPolicy {
    static constexpr const char* NAME = "NullHardeningScalarFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        double zero=0;
        return zero;
    }
    using parameters_t = tuple<>;
};

struct NullHardeningTensorPolicy {
    static constexpr const char* NAME = "NullHardeningTensorFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        VoigtVector zero;
        return zero;
    }
    using parameters_t = tuple<>;
};

// struct ExponentialHardeningForScalarPolicy {
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
        
//         // Because of how ASDPlasticMaterial3D works, instead of using the explicit form of 
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
        auto alpha_dev = alpha.deviator();


        auto eq_norm = [](VoigtVector v){  return sqrt((2./3.)*v.squaredNorm());};

        auto mdev = m.deviator();
        double mdev_eq = eq_norm(mdev);
        auto dEPS_dev = depsilon.deviator();
        double dEPS_dev_eq = eq_norm(dEPS_dev);
        double alpha_norm = eq_norm(alpha_dev);
        // double alpha_norm = alpha_dev.norm();
        // double alpha_limit = sqrt(2. / 3.) * ha / cr;
        double alpha_limit =  ha / cr;

        cout << "depsilon    = " << depsilon.transpose() << endl;
        cout << "m           = " << m.transpose() << endl;
        cout << "mdev        = " << mdev.transpose() << endl;
        cout << "alpha       = " << alpha.transpose() << endl;
        cout << "alpha_norm  = " << alpha_norm << "  <= alpha_limit = " << alpha_limit <<  endl;
        VoigtVector derivative;

        //Compute the derivative (hardening function)
        if (alpha_norm >= alpha_limit)
        {
            cout << "Saturation!" << endl;
            derivative *= 0;  // Take care of the saturation limit in case of overshooting
        }
        else
        {
            derivative =   ha * mdev - cr * mdev_eq * alpha_dev;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha_dev;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha;
            // derivative =  ha * dEPS_dev - cr * dEPS_dev_eq * alpha_dev;
        }
        cout << "----> derivative = " << derivative.transpose() << endl;

        return derivative;

    }
    using parameters_t = tuple<AF_ha, AF_cr>;
};


// Aliases for HardeningFunction with specific hardening policies
using TensorLinearHardeningFunction = HardeningFunction<VoigtVector, LinearHardeningForTensorPolicy>;
using ScalarLinearHardeningFunction = HardeningFunction<VoigtScalar, LinearHardeningForScalarPolicy>;
// using ScalarExponentialHardeningFunction = HardeningFunction<VoigtScalar, ExponentialHardeningForScalarPolicy>;
using ArmstrongFrederickHardeningFunction = HardeningFunction<VoigtVector, ArmstrongFrederickPolicy>;
using NullHardeningScalarFunction = HardeningFunction<VoigtScalar, NullHardeningScalarPolicy>;
using NullHardeningTensorFunction = HardeningFunction<VoigtVector, NullHardeningTensorPolicy>;


#endif //not defined _AllASDHardeningFunctions
