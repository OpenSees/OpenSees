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

#include "../../../ltensor/LTensor.h"
#include "../EvolvingVariable.h"
#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ClassicElastoplasticityGlobals.h"
using namespace ClassicElastoplasticityGlobals;




template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<DTensor2, AlphaHardeningType> AlphaType;
    typedef EvolvingVariable<double, KHardeningType> KType;


    DruckerPrager_YF( AlphaType &alpha_in, KType& k_in):
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                alpha_(alpha_in), k_(k_in)
    {
    }

    double operator()(const DTensor2& sigma) const
    {
        double p;
        static DTensor2 s(3, 3, 0.0);
        const DTensor2 &alpha = alpha_.getVariableConstReference();
        const double &k = k_.getVariableConstReference();
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p = -p;

        double yf = sqrt( (s(i, j) - p * alpha(i, j)) * (s(i, j) - p * alpha(i, j)) ) - SQRT_2_over_3 * k * p; // This one assumes p positive in tension

        return yf;
    }

    const DTensor2& df_dsigma_ij(const DTensor2& sigma)
    {

        static DTensor2 r(3, 3, 0.0);
        static DTensor2 n(3, 3, 0.0);
        const DTensor2 &alpha = alpha_.getVariableConstReference();
        const double &k = k_.getVariableConstReference();

        //Zero these tensors
        s *= 0;
        r *= 0;
        n *= 0;
        result *= 0;

        double p;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p = -p;

        r(i, j) = s(i, j) / p;

        double den = SQRT_2_over_3 * k;
        n(i, j) = (r(i, j) - alpha(i, j)) / den;
        double nr = n(i, j) * r(i, j);
        result(i, j) = n(i, j) - nr * kronecker_delta(i, j) / 3;

        // static DTensor2 n(3, 3, 0.0);
        // static DTensor2 sbar(3, 3, 0.0);
        // const DTensor2 &alpha = alpha_.getVariableConstReference();
        // const double &k = k_.getVariableConstReference();

        // n *= 0;
        // sbar *= 0;

        // double p;
        // result *= 0;
        // sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        // p = -p;

        // sbar(i, j) = s(i, j) - p * alpha(i, j);
        // double s_norm = sqrt(sbar(i, j) * sbar(i, j));

        // n(i, j) = sbar(i, j) / s_norm;

        // double a_n = alpha(i, j) * n(i, j);

        // result(i, j) = n(i, j) + kronecker_delta(i, j) * (k * sqrt(2.0 / 27.) + a_n / 3);

        return result;
    }

    double xi_star_h_star(const DTensor2& depsilon, const DTensor2& m, const DTensor2& sigma)
    {
        double dbl_result = 0.0;

        const DTensor2 &alpha = alpha_.getVariableConstReference();
        const double &k = k_.getVariableConstReference();

        //Zero the stress deviator
        s *= 0;

        //Compute stress deviator (s) and mean pressure (p)
        double p;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension, so flip the sign
        p = -p;

        double den = SQRT_2_over_3 * k;

        // This is for the hardening of k
        double df_dk = -SQRT_2_over_3 * p;
        dbl_result +=  df_dk * k_.getDerivative(depsilon, m, sigma);

        //This is for the hardening of alpha
        dbl_result +=  (( p * alpha(i, j) - s(i, j)) / den) * alpha_.getDerivative(depsilon, m, sigma)(i, j);
        // static DTensor2 sbar(3, 3, 0);
        // sbar *= 0;
        // sbar(i, j) = s(i, j) - p * alpha(i, j);
        // double s_norm = sqrt(sbar(i, j) * sbar(i, j));
        // dbl_result +=  sbar(i, j) / s_norm * (-p) * alpha_.getDerivative(depsilon, m, sigma)(i, j);

        return dbl_result;
    }

    bool hasCorner() const
    {
        return true;
    }
    bool in_Apex(DTensor2 const& TrialStress)
    {
        using namespace ClassicElastoplasticityGlobals;
        double I1 = TrialStress(i, i);
        double sigma_m  = I1 / 3.0;
        double DP_p = -sigma_m;
        static DTensor2 DP_s(3, 3, 0.0);
        DP_s(i, j) = TrialStress(i, j) - kronecker_delta(i, j) * sigma_m;
        double DP_k = k_.getVariableConstReference();
        static DTensor2 DP_alpha(3, 3, 0.0);
        DP_alpha    =  alpha_.getVariableConstReference();
        static DTensor2 relative_s(3, 3, 0.0);
        // Remove the influence of kinematic hardening:
        relative_s(i, j) = DP_s(i, j) - DP_p * DP_alpha(i, j);
        double J2   = 0.5 * relative_s(i, j) * relative_s(i, j) ;
        double DP_q = sqrt(3 * J2);

        if (DP_k < MACHINE_EPSILON)
        {
            cout << "DP_k (denominator) = 0! Error in ClassicElastoplasticMaterial-->requires_DruckerPrager_Apex_Tension_check_\n ";
        }

        // In this range (condition 1 && 2), DP stress state cannot be returned correctly.
        // So the the trick (return to small_stress) will be applied.
        bool condition1 = (  (1.0 / DP_k * DP_p + DP_q) < 0  );
        bool condition2 = (  (1.0 / DP_k * DP_p - DP_q) < 0  );

        if (condition1 && condition2)
        {
            // cout<<"PredictorStress in DruckerPrager apex region. Return to a small stress.\n";
            return true;
        }
        return false;
    }

    double get_k() const
    {
        return k_.getVariableConstReference();
    }
    DTensor2 const& get_alpha() const
    {
        return alpha_.getVariableConstReference();
    }

private:

    AlphaType &alpha_;
    KType &k_;
    static DTensor2 s; //Stress deviator
    static DTensor2 result; //For returning Dtensor2's
};

template <class AlphaHardeningType,  class KHardeningType>
DTensor2 DruckerPrager_YF<AlphaHardeningType, KHardeningType>::s(3, 3, 0.0);
template <class AlphaHardeningType,  class KHardeningType>
DTensor2 DruckerPrager_YF<AlphaHardeningType, KHardeningType>::result(3, 3, 0.0);


#endif