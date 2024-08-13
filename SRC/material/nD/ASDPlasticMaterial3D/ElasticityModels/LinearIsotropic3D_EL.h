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

#ifndef LinearIsotropic3D_EL_H
#define LinearIsotropic3D_EL_H

#include "../ElasticityBase.h"
#include "../AllASDModelParameterTypes.h"

#include <iostream>


class LinearIsotropic3D_EL : public ElasticityBase<LinearIsotropic3D_EL> // CRTP on ElasticityBase
{
public:

    static constexpr const char* NAME = "LinearIsotropic3D_EL";

    LinearIsotropic3D_EL(): ElasticityBase<LinearIsotropic3D_EL>::ElasticityBase()  // Note the full-qualification of ElasticityBase through the scope resolution operator (::)
    {

    }

    ELASTICITY_MATRIX
    {

        double E = GET_PARAMETER_VALUE(YoungsModulus);
        double nu = GET_PARAMETER_VALUE(PoissonsRatio);

        ELAST.setZero(); //Zero it. It may have values from another instance with different parameters;
        double lambda = ( nu * E ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double mu = E / ( 2.0 * ( 1.0 + nu ) );

        ELAST(0, 0) = ELAST(1, 1) = ELAST(2, 2) = 2*mu + lambda;
        ELAST(0, 1) = ELAST(1, 0) = ELAST(0, 2) = ELAST(2, 0) = ELAST(1, 2) = ELAST(2, 1) = lambda;
        ELAST(3, 3) = mu;
        ELAST(4, 4) = mu;
        ELAST(5, 5) = mu;

        return ELAST;
    }

    using parameters_t = std::tuple<YoungsModulus, PoissonsRatio>;

private:

    static VoigtMatrix ELAST;  

};

VoigtMatrix LinearIsotropic3D_EL::ELAST;

#endif
