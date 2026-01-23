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

#ifndef MohrCoulomb_YF_H
#define MohrCoulomb_YF_H

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"


template<class NO_HARDENING>
class MohrCoulomb_YF : public YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_YF";


    MohrCoulomb_YF( ):
        YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>>::YieldFunctionBase() 
        {}

    YIELD_FUNCTION 
    {
        using namespace std;

        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);

        double rThresh = c * cos(phi);
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        double yf = rEquivalentStress - rThresh;

        return yf;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);

        using namespace std;

        double sigma_norm = sigma.norm();

        // Perturbation to smooth the YF
        ds = std::max(ds, ds*sigma_norm);

        // Helper lambda for numerical differentiation
        auto computeNumericalDerivative = [this, &internal_variables_storage, &parameters_storage](const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double yf1 = YF(SIG1);
                double yf2 = YF(SIG2);

                result(i) = (yf1 - yf2) / (2*perturbation);
            }
            return result;
        };

        // If perturbation is set, use numerical differentiation
        if (ds > 0) {
            vv_out = computeNumericalDerivative(sigma, ds);
        } 
        else {
            // Try analytical solution with fallback to numerical
            bool useNumerical = false;
            
            try {
                double J2 = sigma.getJ2();
                
                // Check for numerical issues - use simplified approach for hydrostatic states
                if (J2 < 1e-15) {
                    vv_out = std::sin(phi) / 3.0 * calculate_first_vector();
                } else {
                    VoigtVector first_vector = calculate_first_vector();
                    VoigtVector second_vector = calculate_second_vector(sigma);
                    VoigtVector third_vector = calculate_third_vector(sigma);

                    double lode_angle = sigma.lodeAngle();
                    double c1, c2, c3;
                    double checker = std::abs(lode_angle * 180.0 / M_PI);

                    if (std::abs(checker) < 29.0) { // Regular case
                        c1 = std::sin(phi) / 3.0;
                        
                        double denominator = 2.0 * J2 * std::cos(3.0 * lode_angle);
                        if (std::abs(denominator) < 1e-15) {
                            useNumerical = true; // Division by zero risk
                        } else {
                            c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(phi) * std::cos(lode_angle)) / denominator;
                            c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                                std::sin(phi) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
                        }
                    } else { // Edge smoothing with Drucker-Prager
                        c1 = 3.0 * (2.0 * std::sin(phi) / (std::sqrt(3.0) * (3.0 - std::sin(phi))));
                        c2 = 1.0;
                        c3 = 0.0;
                    }

                    if (!useNumerical) {
                        vv_out = c1 * first_vector + c2 * second_vector + c3 * third_vector;

                        // Validate result - check for NaN/Inf
                        for (int i = 0; i < 6; ++i) {
                            if (!std::isfinite(vv_out(i))) {
                                useNumerical = true;
                                break;
                            }
                        }
                    }
                }
            }
            catch (...) {
                useNumerical = true;
            }

            // Fallback to numerical if analytical failed
            if (useNumerical) {
                vv_out = computeNumericalDerivative(sigma, 1e-6);
            }
        }

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // This model does not support hardening 
        return 0.0;
    }


    CHECK_APEX_REGION
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double I1 = sigma.getI1();
        double p = -I1 / 3;

        if (p < - c / tan(phi))
        {
            return true;
        }

        return false;
    }

    APEX_STRESS
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double p_apex = c / tan(phi);

        vv_out = VoigtVector(p_apex,p_apex,p_apex,0,0,0);

        return vv_out;
    }

  
    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds>;

private:


    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class NO_HARDENING>
VoigtVector MohrCoulomb_YF<NO_HARDENING>::vv_out;

//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<MohrCoulomb_YF<NO_HARDENING>> : std::true_type {};

#endif