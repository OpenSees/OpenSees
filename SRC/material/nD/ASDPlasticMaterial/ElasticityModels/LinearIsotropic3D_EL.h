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

        Ee.setZero(); //Zero it. It may have values from another instance with different parameters;
        double lambda = ( nu * E ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double mu = E / ( 2.0 * ( 1.0 + nu ) );
        // const double mu2 = mu * mu;

        Ee(0, 0) = Ee(1, 1) = Ee(2, 2) = 2*mu + lambda;
        Ee(0, 1) = Ee(1, 0) = Ee(0, 2) = Ee(2, 0) = Ee(1, 2) = Ee(2, 1) = lambda;
        Ee(3, 3) = mu;
        Ee(4, 4) = mu;
        Ee(5, 5) = mu;

        return Ee;
    }

    using parameters_t = std::tuple<YoungsModulus, PoissonsRatio>;

private:

    static VoigtMatrix Ee;  //Provides class-wide storage, which avoids mallocs and allows const returning a const & to this object.

};

VoigtMatrix LinearIsotropic3D_EL::Ee;

#endif
