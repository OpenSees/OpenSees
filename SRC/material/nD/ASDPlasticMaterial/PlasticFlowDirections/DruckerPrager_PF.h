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

#ifndef DruckerPrager_PF_H
#define DruckerPrager_PF_H

#include "../../../ltensor/LTensor.h"
#include "../PlasticFlowBase.h"

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ClassicElastoplasticityGlobals.h"
#include "../EvolvingVariable.h"
#include <cmath>



template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_PF : public PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<DTensor2, AlphaHardeningType> AlphaType;
    typedef EvolvingVariable<double, KHardeningType> KType;


    // PlasticFlowBase<DruckerPrager_PF<HardeningType>>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
    DruckerPrager_PF( AlphaType &alpha_in, KType &k_in):
        PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType , KHardeningType >>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                alpha_(alpha_in), k_(k_in)
    {
        // std::cout << "k_in = " << &k_in << std::endl;
    }


    const DTensor2& operator()(const DTensor2 &depsilon, const DTensor2& sigma)
    {
        using namespace ClassicElastoplasticityGlobals;
        //Identical to derivative of DruckerPrager_YF wrt sigma (a.k.a nij)
        const DTensor2 &alpha = alpha_.getVariableConstReference();
        const double &k = k_.getVariableConstReference();

        // cout << "     --> PF alpha = " << alpha << endl;
        // cout << "     --> PF k     = " << k << endl;


        //Zero these tensors
        s *= 0;
        result *= 0;

        double p;

        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p = -p;

        // cout << "     --> PF p     = " << p << endl;
        // cout << "     --> PF s     = " << s << endl;


        double den = sqrt((s(i, j) - p * alpha(i, j)) * (s(i, j) - p * alpha(i, j)));
        // cout << "     --> PF den     = " << den << endl;

        result(i, j) =
            (
                (s(i, j) - p * alpha(i, j)) + alpha(m, n) * kronecker_delta(i, j) * (s(m, n) - p * alpha(m, n)) / 3.0
            )
            / den;
        result(i, j) += SQRT_2_over_27 * k * kronecker_delta(i, j);

        // DTensor2 bshit(3, 3, 0.0);
        // bshit(i, j) = SQRT_2_over_27 * k * kronecker_delta(i, j);

        // cout << "     --> PF result     = " << result << endl;
        // cout << "     --> PF bshit     = " << bshit << endl;
        return result;
    }

    DTensor4 const& dm_over_dsigma(DTensor2 const& sigma)
    {
        static DTensor2 s(3, 3, 0.0);
        const DTensor2 &alpha = alpha_.getVariableConstReference();
        // const double &k = k_.getVariableConstReference();
        double p = 0.0;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p = -p;
        static DTensor2 s_minus_p_alpha(3, 3, 0.0);
        s_minus_p_alpha(i, j) = s(i, j) - p * alpha(i, j);
        double s_minus_p_alpha_square = s_minus_p_alpha(i, j) * s_minus_p_alpha(i, j) ;

        // // ======================================================================
        // //  Backup . LTensor does not accept this. So change to the naive for-loop.
        // // ======================================================================
        // dm__dsigma(i,j,m,n) =
        //     (
        //         (kronecker_delta(m,i)*kronecker_delta(n,j) - 1./3.0 * kronecker_delta(m,n) * kronecker_delta(i,j)
        //             + 1./3.0 * kronecker_delta(m,n)*alpha(i,j) ) + 1./3.0 * alpha(p,q)*kronecker_delta(i,j) *
        //         (kronecker_delta(m,p)*kronecker_delta(n,q) - 1.0/3.0*kronecker_delta(m,n)*kronecker_delta(p,q)
        //             + 1./3.0 * kronecker_delta(m,n)*alpha(p,q) )
        //     ) * pow(s_minus_p_alpha_square, -0.5)  -
        //     (
        //          (s(i,j)-p*alpha(i,j) + 1./3.0 *alpha(p,q) * kronecker_delta(i,j) * (s(p,q) - p*alpha(p,q))) *
        //          (kronecker_delta(m,r)*kronecker_delta(n,s) - 1./3.0*kronecker_delta(m,n)*kronecker_delta(r,s)
        //             +1./3.0 * kronecker_delta(m,n) * alpha(r,s)) *
        //          (s(r,s)-p*alpha(r,s))
        //     ) * pow(s_minus_p_alpha_square, -1.5);
        // // ======================================================================
        // =========================================
        // The minimal failed example
        // Possible reasons: i,j,m,n are free indices but ?...
        // =========================================
        // test(i,j,m,n)=kronecker_delta(i,m) * kronecker_delta(j,n) - 1.0/3.0 * kronecker_delta(i,j) * kronecker_delta(m,n) ;
        // =========================================
        // This also failed:
        // test(i,j,m,n)=IdentityTensor4(i,m,j,n) - 1.0/3.0 * IdentityTensor4(i,j,m,n) ;
        // =========================================


        // =========================================
        // the naive for-loop:
        // =========================================
        // pre-computation for the dummy indices.
        double alpha_square = alpha(i, j) * alpha(i, j);
        double alpha_volume = alpha(i, i);
        double alpha_times_s = alpha(i, j) * s(i, j);
        // =========================================

        static DTensor4 dm__dsigma(3, 3, 3, 3, 0.0);
        dm__dsigma *= 0;
        // Four free indices
        for (int ig = 0; ig < 3; ++ig)
            for (int jg = 0; jg < 3; ++jg)
                for (int mg = 0; mg < 3; ++mg)
                    for (int ng = 0; ng < 3; ++ng)
                    {
                        dm__dsigma(ig, jg, mg, ng) =
                            (
                                (
                                    kronecker_delta(mg, ig) * kronecker_delta(ng, jg)
                                    - 1. / 3.0 * kronecker_delta(mg, ng) * kronecker_delta(ig, jg)
                                    + 1. / 3.0 * kronecker_delta(mg, ng) * alpha(ig, jg)
                                )
                                +
                                1. / 3.0 * kronecker_delta(ig, jg) *
                                (
                                    alpha(mg, ng)
                                    - 1. / 3.0 * kronecker_delta(mg, ng) * alpha_volume
                                    + 1. / 3.0 * kronecker_delta(mg, ng) * alpha_square
                                )
                            ) * pow(s_minus_p_alpha_square, -0.5)
                            -
                            (
                                (
                                    s_minus_p_alpha(ig, jg)
                                    + 1. / 3.0  * kronecker_delta(ig, jg) *
                                    (
                                        alpha_times_s - p * alpha_square
                                    )
                                )
                                *
                                (
                                    s_minus_p_alpha(mg, ng)
                                    - 1. / 3.0 * kronecker_delta(mg, ng) * (-p * alpha_volume)
                                    + 1. / 3.0 * kronecker_delta(mg, ng) * (alpha_times_s - p * alpha_square)
                                )
                            ) * pow(s_minus_p_alpha_square, -1.5);
                    }


        return dm__dsigma;
    }



    DTensor2 const& dm_over_dq_start_h_star(DTensor2 const& depsilon, DTensor2 const& pf_m, const DTensor2& stress){
        static DTensor2 s(3, 3, 0.0);
        const DTensor2 &alpha = alpha_.getVariableConstReference();
        // const double &_k_ = k_.getVariableConstReference();
        double p=0;
        stress.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p=-p;

        static DTensor4 IdentityTensor4(3,3,3,3, 0); //optimize this to global later.
        IdentityTensor4(i,j,k,l)=kronecker_delta(i, j)*kronecker_delta(k,l);
        // (1) isotropic hardening part. 
        static DTensor2 dm_dk(3,3,0.0);
        dm_dk(i,j) = SQRT_2_over_27 * kronecker_delta(i, j) ; 

        // (2) kinematic hardening part: dm_dalpha
        static DTensor4 dm_dalpha(3,3,3,3,0.0);
        static DTensor2 s_minus_p_alpha(3,3,0.0);
        s_minus_p_alpha(i,j) = s(i,j) - p * alpha(i,j);
        double s_minus_p_alpha_square = s_minus_p_alpha(i,j) * s_minus_p_alpha(i,j) ; 
        double alpha_times_s_minus_p_alpha = s_minus_p_alpha(i,j) * alpha(i,j);
        for (int ig = 0; ig < 3; ++ig)
            for (int jg = 0; jg < 3; ++jg)
                for (int kg = 0; kg < 3; ++kg)
                    for (int lg = 0; lg < 3; ++lg){
                        dm_dalpha(ig,jg,kg,lg) = 
                            (
                                - p * IdentityTensor4(kg,ig,lg,jg) 
                                + 1./3. * kronecker_delta(ig, jg) *  s_minus_p_alpha(kg,lg)  
                                - 1./3. * p * kronecker_delta(ig, jg) *  alpha(kg,lg) 
                            ) * pow(s_minus_p_alpha_square,-0.5) 
                            -
                            (
                                (s_minus_p_alpha(ig,jg) + 1./3.*kronecker_delta(ig,jg) * alpha_times_s_minus_p_alpha)
                                *(- p *  s_minus_p_alpha(kg,lg) )
                                * pow(s_minus_p_alpha_square,-1.5)
                            );
                        }

        static DTensor2 ret(3,3,0.0);
        ret(i,j) = dm_dalpha(i,j,m,n) * alpha_.getDerivative(depsilon, pf_m, stress)(m,n);
        ret(i,j) += dm_dk(i,j) * k_.getDerivative(depsilon, pf_m, stress) ;

        return ret;
    }


private:

    AlphaType &alpha_;
    KType &k_;

    static DTensor2 s; //sigma deviator
    static DTensor2 result; //For returning Dtensor2s
    // static DTensor4 dm__dsigma; //For returning dm_over_dsigma
    // static DTensor4 dm__dalpha; //For returning dm_over_dsigma

};

template<class AlphaHardeningType, class KHardeningType>
DTensor2 DruckerPrager_PF<AlphaHardeningType , KHardeningType >::s(3, 3, 0.0);
template<class AlphaHardeningType, class KHardeningType>
DTensor2 DruckerPrager_PF<AlphaHardeningType , KHardeningType >::result(3, 3, 0.0);
// template<typename AlphaHardeningType, typename KHardeningType>
// DTensor4 DruckerPrager_PF<AlphaHardeningType , KHardeningType >::dm__dsigma(3, 3, 3, 3, 0.0);
// template<typename AlphaHardeningType, typename KHardeningType>
// DTensor4 DruckerPrager_PF<AlphaHardeningType , KHardeningType >::dm__dalpha(3, 3, 3, 3, 0.0);

#endif


// ============================
// legacy
// should be removed later
// ============================

// DTensor4 const& dm_over_dalpha(DTensor2 const& sigma){
//     static DTensor2 s(3, 3, 0.0);
//     const DTensor2 &alpha = alpha_.getVariableConstReference();
//     // const double &k = k_.getVariableConstReference();
//     double p=0.0;
//     sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
//     p=-p;
//     static DTensor2 s_minus_palpha(3,3,0.0);
//     s_minus_palpha(i,j) = s(i,j) - p*alpha(i,j);
//     double s_minus_p_alpha_square = s_minus_palpha(i,j) * s_minus_palpha(i,j) ;
//     static DTensor4 dm__dalpha(3,3,3,3,0.0);
//     dm__dalpha*=0;
//     for (int ig = 0; ig < 3; ++ig)
//         for (int mg = 0; mg < 3; ++mg)
//             for (int jg = 0; jg < 3; ++jg)
//                 for (int ng = 0; ng < 3; ++ng)
//                     for (int pg = 0; pg < 3; ++pg)
//                         for (int qg = 0; qg < 3; ++qg)
//                             for (int rg = 0; rg < 3; ++rg)
//                                 for (int sg = 0; sg < 3; ++sg){
//                                     dm__dalpha(ig,jg,mg,ng) +=
//                                         (
//                                             -p*kronecker_delta(mg,ig)*kronecker_delta(ng,jg) + 1./3.0 * kronecker_delta(mg,pg) * kronecker_delta(ig,jg)

//                                         ) * pow(s_minus_p_alpha_square, -0.5)  ;
//                                         // Not finished yet!
//                                         // -
//                                         // (
//                                         //      (s(ig,jg)-p*alpha(ig,jg) + 1./3.0 *alpha(pg,qg) * kronecker_delta(ig,jg) * (s(pg,qg) - p*alpha(pg,qg))) *
//                                         //      (kronecker_delta(mg,rg)*kronecker_delta(ng,sg) - 1./3.0*kronecker_delta(mg,ng)*kronecker_delta(rg,sg)
//                                         //         +1./3.0 * kronecker_delta(mg,ng) * alpha(rg,sg)) *
//                                         //      (s(rg,sg)-p*alpha(rg,sg))
//                                         // ) * pow(s_minus_p_alpha_square, -1.5);
//                                 }


//     return dm__dalpha;
// }