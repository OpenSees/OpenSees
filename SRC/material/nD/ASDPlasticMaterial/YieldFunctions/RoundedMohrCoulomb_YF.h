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

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ASDPlasticMaterialGlobals.h"
using namespace ASDPlasticMaterialGlobals;



template<class ETAHardeningType>
class RoundedMohrCoulomb_YF : public YieldFunctionBase<RoundedMohrCoulomb_YF< ETAHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "RoundedMohrCoulomb_YF";

    RoundedMohrCoulomb_YF( double m_in, double qa_in, double pc_in, double e_in, ETAType& eta_in):
        YieldFunctionBase<RoundedMohrCoulomb_YF<ETAHardeningType>>::YieldFunctionBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                m(m_in), qa(qa_in), pc(pc_in), e(e_in),
                eta_(eta_in)
    {
        // std::cout << "k_in = " << &k_in << std::endl;
    }

    double operator()(const VoigtVector& sigma) const
    {
        double p, q, theta;
        std::tie(p, q, theta) = getpqtheta(sigma);

        theta *= M_PI / 180;  // Convert theta back to radians! :/

        const double &eta = eta_.getVariableConstReference();

        double yf = q * pow(1 + q / qa, m)  - eta * (p - pc) / g(theta);

        return yf;
    }

    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma)
    {
        const double &eta = eta_.getVariableConstReference();

        double p, q, theta;

        std::tie(p, q, theta) = getpqtheta(sigma);
        theta *= M_PI / 180;  // Convert theta back to radians! :/
        double gg = g(theta);
        double dgg = dg_dlode(theta);

        static VoigtVector dq_dsij(3, 3, 0);
        static VoigtVector dtheta_dsij(3, 3, 0);
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

    double xi_star_h_star(const VoigtVector& depsilon, const VoigtVector& m, const VoigtVector& sigma)
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

    VoigtVector & get_alpha() const
    {
        return s;
    }

    bool in_Apex(VoigtVector const& TrialStress)
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
    using internal_variables_t = std::tuple<ETAType>;
    using parameters_t = std::tuple<RMC_m,RMC_qa,RMC_pc,RMC_e>;
    
    // double m, qa, pc, e;
    // ETAType &eta_;
    // static VoigtVector s; //Stress deviator
    static VoigtVector result; //For returning VoigtVector's
};


template <  class ETAHardeningType>
VoigtVector RoundedMohrCoulomb_YF< ETAHardeningType>::result;



#endif