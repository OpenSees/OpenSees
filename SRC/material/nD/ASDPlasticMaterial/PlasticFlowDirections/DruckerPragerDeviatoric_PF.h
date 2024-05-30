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

#ifndef DruckerPragerDeviatoric_PF_H
#define DruckerPragerDeviatoric_PF_H

#include "../../../ltensor/LTensor.h"
#include "../PlasticFlowBase.h"
#include "../EvolvingVariable.h"
// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ASDPlasticMaterialGlobals.h"





template<class AlphaHardeningType, class KHardeningType>
class DruckerPragerDeviatoric_PF : public PlasticFlowBase<DruckerPragerDeviatoric_PF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<VoigtVector, AlphaHardeningType> AlphaType;
    typedef EvolvingVariable<double, KHardeningType> KType;


    // PlasticFlowBase<DruckerPragerDeviatoric_PF<HardeningType>>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
    DruckerPragerDeviatoric_PF( AlphaType &alpha_in, KType &k_in):
        PlasticFlowBase<DruckerPragerDeviatoric_PF<AlphaHardeningType , KHardeningType >>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                alpha_(alpha_in), k_(k_in)
    {

    }


    const VoigtVector& operator()(const VoigtVector &depsilon, const VoigtVector& sigma)
    {
        //Identical to derivative of VonMises_YF wrt sigma (a.k.a nij)
        const VoigtVector &alpha = alpha_.getVariableConstReference();

        //Zero these tensors
        s *= 0;
        result *= 0;

        double p = -sigma(i, i) / 3;

        s(i, j) = sigma(i, j) + p * kronecker_delta(i, j);
        result(i, j) = s(i, j) - p * alpha(i, j);
        double den = sqrt(result(i, j) * result(i, j));

        if (den == 0)
        {
            return result;
        }
        else
        {
            result(i, j) = result(i, j) / den;
        }
        // cout << "m = [";
        // for (int ii = 0; ii < 3; ii++)
        //     for (int jj = 0; jj < 3; jj++)
        //     {
        //         cout << result(ii, jj) << " ";
        //     }
        // cout << "]\n";
        // cout << "alpha = [";
        // for (int ii = 0; ii < 3; ii++)
        //     for (int jj = 0; jj < 3; jj++)
        //     {
        //         cout << alpha(ii, jj) << " ";
        //     }
        // cout << "]\n";

        return result;
    }
    VoigtMatrix const& dm_over_dsigma(VoigtVector const& sigma){
        static VoigtVector s(3, 3, 0.0);
        const VoigtVector &alpha = alpha_.getVariableConstReference();
        // const double &k = k_.getVariableConstReference();
        double p=0.0;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p=-p;
        // if(p<MACHINE_EPSILON){
        //     cout<<"DruckerPragerNonAssociate_PF:: pressuse < 0 ! The Drucker-Prager has tensile force.\n";
        // }
        static VoigtVector s_minus_palpha(3,3,0.0);
        s_minus_palpha(i,j) = s(i,j) - p*alpha(i,j);
        double s_minus_p_alpha_square = s_minus_palpha(i,j) * s_minus_palpha(i,j) ; 
        double s_square = s(i,j) * s(i,j) ; 
        // =========================================
        // the naive for-loop:
        // =========================================
        // pre-computation for the dummy indices. 
        double alpha_square = alpha(i,j) * alpha(i,j);
        double alpha_volume = alpha(i,i);
        double alpha_times_s= alpha(i,j) * s(i,j);
        // =========================================
        static VoigtMatrix dm__dsigma(3,3,3,3,0.0);
        dm__dsigma*=0;
        // Four free indices
        for (int ig = 0; ig < 3; ++ig)
            for (int jg = 0; jg < 3; ++jg)
                for (int mg = 0; mg < 3; ++mg)
                    for (int ng = 0; ng < 3; ++ng){
                        dm__dsigma(ig,jg,mg,ng) = 
                            ( 
                                kronecker_delta(mg,ig)*kronecker_delta(ng,jg) 
                                - 1./3.0 * kronecker_delta(mg,ng) * kronecker_delta(ig,jg) 
                                + 1./3.0 * kronecker_delta(mg,ng)*alpha(ig,jg)   
                            ) * pow(s_minus_p_alpha_square, -0.5)  
                            -
                            (
                                s_minus_palpha(ig,jg) 
                                *
                                (
                                    s_minus_palpha(mg,ng)
                                    - 1./3.0 * kronecker_delta(mg,ng) * (-p*alpha_volume)
                                    + 1./3.0 * kronecker_delta(mg,ng) * (alpha_times_s-p*alpha_square)
                                ) 
                                * 
                                pow(s_minus_p_alpha_square, -1.5)
                            )   
                            + 
                            (
                                1./3. * kronecker_delta(ig,jg) * s(mg,ng) * pow(p,-1) 
                            ) * pow(s_square,-0.5) 
                            +
                            (
                                1./9. * kronecker_delta(ig,jg) * kronecker_delta(mg,ng) * pow(p,-2)
                            ) * pow(s_square,0.5); 
                    }


        return dm__dsigma;
    }
    

    VoigtVector const& dm_over_dq_start_h_star(VoigtVector const& depsilon, VoigtVector const& pf_m, const VoigtVector& stress){
        static VoigtVector s(3, 3, 0.0);
        const VoigtVector &alpha = alpha_.getVariableConstReference();
        const double &k = k_.getVariableConstReference();
        double p=0;
        stress.compute_deviatoric_tensor(s, p); // here p is positive if in tension
        p=-p;

        static VoigtMatrix IdentityTensor4(3,3,3,3, 0); 
        IdentityTensor4(i,j,k,l)=kronecker_delta(i, j)*kronecker_delta(k,l);
        // (1) isotropic hardening part. 
        static VoigtVector dm_dk(3,3,0.0);
        dm_dk(i,j) = SQRT_2_over_27 * kronecker_delta(i, j) ; 

        // (2) kinematic hardening part
        static VoigtMatrix dm_dalpha(3,3,3,3,0.0);
        static VoigtVector s_minus_p_alpha(3,3,0.0);
        s_minus_p_alpha(i,j) = s(i,j) - p * alpha(i,j);
        double s_minus_p_alpha_square = s_minus_p_alpha(i,j) * s_minus_p_alpha(i,j) ; 

        for (int ig = 0; ig < 3; ++ig)
            for (int jg = 0; jg < 3; ++jg)
                for (int kg = 0; kg < 3; ++kg)
                    for (int lg = 0; lg < 3; ++lg){
                        dm_dalpha(ig,jg,kg,lg) = 
                            (
                                - p * IdentityTensor4(kg,ig,lg,jg) 
                            ) * pow(s_minus_p_alpha_square,-0.5) 
                            -
                            (
                                s_minus_p_alpha(ig,jg) 
                                *(- p *  s_minus_p_alpha(kg,lg) )
                                * pow(s_minus_p_alpha_square,-1.5)
                            );
                        }

        static VoigtVector ret(3,3,0.0);
        ret(i,j) = dm_dalpha(i,j,m,n) * alpha_.getDerivative(depsilon, pf_m, stress)(m,n);
        ret(i,j) += dm_dk(i,j) * k_.getDerivative(depsilon, pf_m, stress) ;
                
        return ret;
    }
private:

    AlphaType &alpha_;
    KType &k_;

    static VoigtVector s; //sigma deviator
    static VoigtVector result; //For returning VoigtVectors

};


template<class AlphaHardeningType, class KHardeningType>
VoigtVector DruckerPragerDeviatoric_PF<AlphaHardeningType , KHardeningType >::s(3, 3, 0.0);
template<class AlphaHardeningType, class KHardeningType>
VoigtVector DruckerPragerDeviatoric_PF<AlphaHardeningType , KHardeningType >::result(3, 3, 0.0);

#endif



    // VoigtMatrix const& dm_over_dalpha(VoigtVector const& sigma, VoigtVector const& m){
    //     static VoigtMatrix placeholder(3,3,3,3,0.0);
    //     return placeholder;
    // }
    // 