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

#ifndef TensionCutoff_YF_H
#define TensionCutoff_YF_H


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"


template<class NO_HARDENING>
class TensionCutoff_YF : public YieldFunctionBase<TensionCutoff_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "TensionCutoff_YF";


    TensionCutoff_YF( ):
        YieldFunctionBase<TensionCutoff_YF<NO_HARDENING>>::YieldFunctionBase() 
        {}

    YIELD_FUNCTION 
    {
        using namespace std;
        double min_stress = GET_PARAMETER_VALUE(TC_min_stress);

        double s1,s2,s3;
        std::tie(s1,s2,s3) =  sigma.principalStresses();  //Ordered

        double yf = s3-min_stress;  // If the smallest principal stress is less than 0, is elastic


        return yf;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        double ds = GET_PARAMETER_VALUE(MC_ds);

        using namespace std;

        double sigma_norm = sigma.norm();

        // const double DL = 1e-8;
        ds = ds == 0 ? 1e-8 : ds;


        // If the perturbation is set greater than zero, use numerical differentiation to get normal to YF

        for (int i = 0; i < 6; ++i) {
            VoigtVector SIG1 = sigma;
            VoigtVector SIG2 = sigma;
            
            // Increment SIG at index i by a small amount for numerical differentiation
            SIG1(i) += ds;
            SIG2(i) -= ds;

            // Compute the yield function at the perturbed state
            double yf1 = YF(SIG1);
            double yf2 = YF(SIG2);

            // Calculate the derivative
            vv_out(i) = (yf1 - yf2) / (2*ds);
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
        // double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        // double c = GET_PARAMETER_VALUE(MC_c);
        // double I1 = sigma.getI1();
        // double p = -I1 / 3;

        // if (p < - c / tan(phi))
        // {
        //     return true;
        // }

        // Implement!!! 

        return false;
    }

    APEX_STRESS
    {

        // Implement!!! 
        vv_out = VoigtVector(0,0,0,0,0,0);

        return vv_out;

    }

  
    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_ds, TC_min_stress>;

private:


    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class NO_HARDENING>
VoigtVector TensionCutoff_YF<NO_HARDENING>::vv_out;

//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<TensionCutoff_YF<NO_HARDENING>> : std::true_type {};

#endif



