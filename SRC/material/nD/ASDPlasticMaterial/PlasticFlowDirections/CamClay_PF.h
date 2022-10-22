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

#ifndef CamClay_PF_H
#define CamClay_PF_H

#include "../../../ltensor/LTensor.h"
#include "../PlasticFlowBase.h"

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ASDPlasticMaterialGlobals.h"

#include <cmath>



template<class p0HardeningType>
class CamClay_PF : public PlasticFlowBase<CamClay_PF<p0HardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<double, p0HardeningType> p0Type;

    // PlasticFlowBase<CamClay_PF<HardeningType>>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
    CamClay_PF( double M_in, p0Type& p0_in):
        PlasticFlowBase<CamClay_PF<p0HardeningType >>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                M(M_in), p0_(p0_in)
    {    }


    const VoigtVector& operator()(const VoigtVector &depsilon, const VoigtVector& sigma)
    {
        //Zero these tensors
        s *= 0;
        result *= 0;

        double p;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension, so flip the sign
        p = -p;
        const double &p0 = p0_.getVariableConstReference();


        double scalar1 = M * M * (p0 - 2 * p) / 3;

        result(i, j) = s(i, j) * 3 + kronecker_delta(i, j) * scalar1;
        return result;
    }

    VoigtMatrix const& dm_over_dsigma(VoigtVector const& sigma)
    {
        static VoigtMatrix dm__dsigma(3, 3, 3, 3, 0.0);
        static VoigtMatrix delta_imjn(3, 3, 3, 3, 0.0)        ;
        static VoigtMatrix delta_ijmn(3, 3, 3, 3, 0.0)        ;
        delta_imjn(i, j, m, n) = kronecker_delta(i, m) * kronecker_delta(j, n);
        delta_ijmn(i, j, m, n) = kronecker_delta(i, j) * kronecker_delta(m, n);

        s *= 0;
        dm__dsigma *= 0;
        double scalar = (2 * M * M) / 9 - 1 ;

        dm__dsigma(i, j, m, n) = 3 * delta_imjn(i, j, m, n) + scalar * delta_ijmn(i, j, m, n);

        return dm__dsigma;
    }


    VoigtVector const& dm_over_dq_start_h_star(const VoigtVector&dumb1, const VoigtVector&dumb2, const VoigtVector& stress)
    {
        result *= 0;
        // const double &p0 = p0_.getVariableConstReference();

        result(i, j) = (M * M) / 3 * kronecker_delta(i, j) * p0_.getDerivative(stress, stress, stress);  // Will only use last entry, which is sigma

        return result;
    }


private:

    double M;
    p0Type &p0_;

    static VoigtVector s; //sigma deviator
    static VoigtVector result; //For returning VoigtVectors
    // static VoigtMatrix dm__dsigma; //For returning dm_over_dsigma
    // static VoigtMatrix dm__dalpha; //For returning dm_over_dsigma

};

template<class p0HardeningType>
VoigtVector CamClay_PF<p0HardeningType >::s(3, 3, 0.0);
template<class p0HardeningType>
VoigtVector CamClay_PF<p0HardeningType >::result(3, 3, 0.0);

#endif
