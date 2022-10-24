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

#ifndef YieldFunctionBase_H
#define YieldFunctionBase_H

#include "EigenAPI.h"
#include "EvolvingVariable.h"
#include <Channel.h>


template <class T>
class YieldFunctionBase
{
public:
    YieldFunctionBase()
    {
    }

    double operator()( const VoigtVector& sigma) const
    {
        return static_cast<T*>(this)->operator()(sigma);
    }

    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma)
    {
        return static_cast<T*>(this)->df_dsigma_ij(sigma);
    }

    double xi_star_h_star(const VoigtVector& depsilon, const VoigtVector& depsilon_pl, const VoigtVector& sigma)
    {
        return static_cast<T*>(this)->df_dxi_star_h_star(depsilon, depsilon_pl , sigma);
    }

private:

};




#endif