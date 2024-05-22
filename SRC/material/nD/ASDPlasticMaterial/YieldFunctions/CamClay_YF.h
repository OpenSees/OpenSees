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

#ifndef CamClay_YF_H
#define CamClay_YF_H

#include "../../../ltensor/LTensor.h"
#include "../EvolvingVariable.h"
#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

// Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.
#include "../ASDPlasticMaterialGlobals.h"
using namespace ASDPlasticMaterialGlobals;




template<class p0HardeningType>
class CamClay_YF : public YieldFunctionBase<CamClay_YF<p0HardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<double, p0HardeningType> p0Type;


    CamClay_YF( double M_in, p0Type& p0_in):
        YieldFunctionBase<CamClay_YF<p0HardeningType>>::YieldFunctionBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                M(M_in), p0_(p0_in)
    {    }

    double operator()(const VoigtVector& sigma) const
    {
        double p, q, theta;
        std::tie(p, q, theta) = getpqtheta(sigma);
        const double &p0 = p0_.getVariableConstReference();

        double yf = q * q - M * M * p * (p0 - p);

        return yf;
    }

    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma)
    {
        //Zero these tensors
        s *= 0;
        result *= 0;

        double p;
        sigma.compute_deviatoric_tensor(s, p); // here p is positive if in tension, so flip the sign
        p = -p;
        const double &p0 = p0_.getVariableConstReference();


        double scalar1 = M * M * (p0 - 2 * p) / 3;

        result(i, j) = 3 * s(i, j)  + kronecker_delta(i, j) * scalar1;
        return result;
    }

    double xi_star_h_star(const VoigtVector& depsilon, const VoigtVector& m, const VoigtVector& sigma)
    {
        double dbl_result = 0.0;

        // const double &p0 = p0_.getVariableConstReference();

        //Compute stress deviator (s) and mean pressure (p)
        double p = -sigma(i, i) / 3;

        // This is for the hardening of k
        double df_dp0 = M * M * p;
        dbl_result +=  df_dp0 * p0_.getDerivative(depsilon, m, sigma);

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

    double M;
    p0Type &p0_;
    static VoigtVector s; //Stress deviator
    static VoigtVector result; //For returning VoigtVector's
};

template <class p0HardeningType>
VoigtVector CamClay_YF<p0HardeningType>::s(3, 3, 0.0);
template <class p0HardeningType>
VoigtVector CamClay_YF<p0HardeningType>::result(3, 3, 0.0);


#endif