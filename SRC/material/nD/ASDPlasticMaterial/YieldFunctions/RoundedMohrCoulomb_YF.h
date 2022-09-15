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

#ifndef RoundedMohrCoulomb_YF_H
#define RoundedMohrCoulomb_YF_H

#include "../../../ltensor/LTensor.h"
#include "../EvolvingVariable.h"
#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ClassicElastoplasticityGlobals.h"
using namespace ClassicElastoplasticityGlobals;



template<class ETAHardeningType>
class RoundedMohrCoulomb_YF : public YieldFunctionBase<RoundedMohrCoulomb_YF< ETAHardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<double, ETAHardeningType> ETAType;

    RoundedMohrCoulomb_YF( double m_in, double qa_in, double pc_in, double e_in, ETAType& eta_in):
        YieldFunctionBase<RoundedMohrCoulomb_YF<ETAHardeningType>>::YieldFunctionBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                m(m_in), qa(qa_in), pc(pc_in), e(e_in),
                eta_(eta_in)
    {
        // std::cout << "k_in = " << &k_in << std::endl;
    }

    double operator()(const DTensor2& sigma) const
    {
        double p, q, theta;
        std::tie(p, q, theta) = getpqtheta(sigma);

        theta *= M_PI / 180;  // Convert theta back to radians! :/

        const double &eta = eta_.getVariableConstReference();

        double yf = q * pow(1 + q / qa, m)  - eta * (p - pc) / g(theta);

        return yf;
    }

    const DTensor2& df_dsigma_ij(const DTensor2& sigma)
    {
        const double &eta = eta_.getVariableConstReference();

        double p, q, theta;

        std::tie(p, q, theta) = getpqtheta(sigma);
        theta *= M_PI / 180;  // Convert theta back to radians! :/
        double gg = g(theta);
        double dgg = dg_dlode(theta);

        static DTensor2 dq_dsij(3, 3, 0);
        static DTensor2 dtheta_dsij(3, 3, 0);
        dq_dsij *= 0;
        dtheta_dsij *= 0;
        dq_dsigma_ij(sigma, dq_dsij);
        dtheta_dsigma_ij(sigma, dtheta_dsij);

        result(i, j) = dq_dsij(i, j) * pow(1 + 1 / qa, m - 1) * (1 + q / qa + q * m) -
                       eta * (-kronecker_delta(i, j) / 3 / gg - (p - pc) / (gg * gg) * dgg * dtheta_dsij(i, j) );

        // cout << "sigma  = " << sigma  << endl;
        // cout << "p = " << p << endl;
        // cout << "q = " << q << endl;
        // cout << "theta = " << theta << endl;
        // cout << "gg = " << gg << endl;
        // cout << "dq_dsij = " << dq_dsij << endl;
        // cout << "dtheta_dsij = " << dtheta_dsij << endl;
        // cout << "df_dsigma_ij = " << result << endl;
        // cout << "a = " << a << endl;

        return result;
    }

    double xi_star_h_star(const DTensor2& depsilon, const DTensor2& m, const DTensor2& sigma)
    {
        // const double &k = k_.getVariableConstReference();

        double dbl_result = 0;
        double p, q, theta;
        std::tie(p, q, theta) = getpqtheta(sigma);
        theta *= M_PI / 180;
        dbl_result += -(p - pc) / g(theta) * eta_.getDerivative(depsilon, m, sigma) ;

        // cout << "RMC -> theta = " << theta << endl;
        // cout << "RMC -> g(theta) = " << g(theta) << endl;
        // cout << "RMC -> (p - pc) = " << (p - pc) << endl;
        // cout << "RMC -> dbl_result = " << dbl_result << endl;

        return dbl_result;
    }
    bool hasCorner() const
    {
        return false;
    }

    double get_k() const
    {
        return 0;
    }

    DTensor2 & get_alpha() const
    {
        return s;
    }

    bool in_Apex(DTensor2 const& TrialStress)
    {

        return false;
    }


private:

    //================================================================================
    double g(double lode) const // Willam-Warnke function
    {
        double coslode = cos(lode);
        double coslode2 = coslode * coslode;
        double num = 4 * (1 - e * e) * coslode2  + (2 * e - 1) * (2 * e - 1);
        double den = 2 * (1 - e * e) * coslode + (2 * e - 1) * sqrt(4 * (1. - e * e) * coslode2 + 5 * e * e - 4 * e);
        return num / den;
    }

    double dg_dlode(double theta) const //Derivative of above.
    {
        // === Python code to generate this function
        // import sympy
        //e = sympy.S("e")
        //theta = sympy.S("theta")
        //costheta = sympy.cos(theta)
        //num = 4 * (1 - e * e) * costheta**2  + (2 * e - 1) * (2 * e - 1)
        //den = 2 * (1 - e * e) * costheta + (2 * e - 1) * sympy.sqrt(4 * (1 - e * e) * costheta**2 + 5 * e * e - 4 * e)
        //g = num/den
        //dg = sympy.simplify(sympy.diff(g, theta))
        //print sympy.ccode(dg)
        double dg = (2 * pow(e, 2) - 2) * (4 * ((2 * e - 1) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) + (-2 * pow(e, 2) + 2) * cos(theta)) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) * cos(theta) - ((4 * e - 2) * cos(theta) + sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2))) * (pow(2 * e - 1, 2) + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2))) * sin(theta) / (pow((2 * e - 1) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) + (-2 * pow(e, 2) + 2) * cos(theta), 2) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)));
        return dg;
    }

    double m, qa, pc, e;
    ETAType &eta_;
    static DTensor2 s; //Stress deviator
    static DTensor2 result; //For returning Dtensor2's
};

template < class ETAHardeningType>
DTensor2 RoundedMohrCoulomb_YF< ETAHardeningType>::s(3, 3, 0.0);
template <  class ETAHardeningType>
DTensor2 RoundedMohrCoulomb_YF< ETAHardeningType>::result(3, 3, 0.0);




//  Jose (Sep 27 2016)
//  Unfinished implementation trying to use kinematic hardening... this was folly.
// The reason is that RMC should be used to set the ultimate strength of the material
//  and not initial yielding. The purpose of KH is to set the ultimate through the
// hardening law.
//
// It makes more sense to use RMC as a bounding surface than as a YF.
//
//
// template<class AlphaHardeningType, class ETAHardeningType>
// class RoundedMohrCoulomb_YF : public YieldFunctionBase<RoundedMohrCoulomb_YF<AlphaHardeningType, ETAHardeningType>> // CRTP
// {
// public:

//     typedef EvolvingVariable<DTensor2, AlphaHardeningType> AlphaType;
//     typedef EvolvingVariable<double, ETAHardeningType> ETAType;


//     RoundedMohrCoulomb_YF( AlphaType &alpha_in, ETAType& k_in):
//         YieldFunctionBase<RoundedMohrCoulomb_YF<AlphaHardeningType, ETAHardeningType>>::YieldFunctionBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
//                 alpha_(alpha_in), k_(k_in)
//     {
//         // std::cout << "k_in = " << &k_in << std::endl;
//     }

//     double operator()(const DTensor2& sigma) const
//     {
//         double p, q, theta;
//         std::tie(p, q, theta) = getpqtheta(sigma);

//         theta *= M_PI / 180;  // Convert theta back to radians! :/

//         const DTensor2 &alpha = alpha_.getVariableConstReference();
//         const double &k = k_.getVariableConstReference();

//         double yf = q * pow(1 + q / qa, m)  - k * (p - pc) / g(theta);

//         return yf;
//     }

//     const DTensor2& df_dsigma_ij(const DTensor2& sigma)
//     {
//         const DTensor2 &alpha = alpha_.getVariableConstReference();
//         const double &k = k_.getVariableConstReference();

//         double p, q, theta;
//         static DTensor2 sigma_minus_alpha(3, 3, 0.0);
//         sigma_minus_alpha *= 0;
//         sigma_minus_alpha(i, j) = sigma(i, j) - sigma(i, i) / 3 * alpha(i, j);

//         std::tie(p, q, theta) = getpqtheta(sigma_minus_alpha);
//         theta *= M_PI / 180;  // Convert theta back to radians! :/
//         double gg = g(theta);
//         double dgg = dg_dlode(theta);

//         static DTensor2 dq_dsij(3, 3, 0);
//         static DTensor2 dtheta_dsij(3, 3, 0);
//         dq_dsij * = 0;
//         dtheta_dsij * = 0;
//         dq_dsigma_ij(sigma_minus_alpha, dq_dsij);
//         dtheta_dsigma_ij(sigma_minus_alpha, dtheta_dsij);

//         result(i, j) = dq_dsij(i, j) * pow(1 + 1 / qa, m - 1) * (1 + q / qa + q * m) -
//                        k * (-kronecker_delta(i, j) / 3 / gg - (p - pc) / (gg * gg) * dgg * dtheta_dsij(i, j) );

//         return result;
//     }

//     double xi_star_h_star(const DTensor2& depsilon, const DTensor2& m, const DTensor2& sigma)
//     {

//         return dbl_result;
//     }

//     bool hasCorner() const
//     {
//         return true;
//     }
//     double get_k() const
//     {
//         return k_.getVariableConstReference();
//     }
//     DTensor2 const& get_alpha() const
//     {
//         return alpha_.getVariableConstReference();
//     }
// private:

//     //================================================================================
//     double g(lode)
//     {
//         double coslode = cos(lode);
//         double coslode2 = coslode * coslode;
//         double num = 4 * (1 - e * e) * coslode2  + (2 * e - 1) * (2 * e - 1);
//         double den = 2 * (1 - e * e) * coslode + (2 * e - 1) * sqrt(4 * (1 - e * e) * coslode2 + 5 * e * e - 4 * e);
//         return num / den;
//     }

//     double dg_dlode(lode)
//     {
//         // === Python code to generate this function
//         // import sympy
//         //e = sympy.S("e")
//         //theta = sympy.S("theta")
//         //costheta = sympy.cos(theta)
//         //num = 4 * (1 - e * e) * costheta**2  + (2 * e - 1) * (2 * e - 1)
//         //den = 2 * (1 - e * e) * costheta + (2 * e - 1) * sympy.sqrt(4 * (1 - e * e) * costheta**2 + 5 * e * e - 4 * e)
//         //g = num/den
//         //dg = sympy.simplify(sympy.diff(g, theta))
//         //print sympy.ccode(dg)
//         double dg = (2 * pow(e, 2) - 2) * (4 * ((2 * e - 1) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) + (-2 * pow(e, 2) + 2) * cos(theta)) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) * cos(theta) - ((4 * e - 2) * cos(theta) + sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2))) * (pow(2 * e - 1, 2) + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2))) * sin(theta) / (pow((2 * e - 1) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)) + (-2 * pow(e, 2) + 2) * cos(theta), 2) * sqrt(5 * pow(e, 2) - 4 * e + (-4 * pow(e, 2) + 4) * pow(cos(theta), 2)));
//         return dg;
//     }

//     double m, qa, pc, e;
//     AlphaType &alpha_;
//     ETAType &k_;
//     static DTensor2 s; //Stress deviator
//     static DTensor2 result; //For returning Dtensor2's
// };

// template <class AlphaHardeningType,  class ETAHardeningType>
// DTensor2 RoundedMohrCoulomb_YF<AlphaHardeningType, ETAHardeningType>::s(3, 3, 0.0);
// template <class AlphaHardeningType,  class ETAHardeningType>
// DTensor2 RoundedMohrCoulomb_YF<AlphaHardeningType, ETAHardeningType>::result(3, 3, 0.0);


#endif