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

#ifndef DuncanChang_EL_H
#define DuncanChang_EL_H

#include "../ElasticityBase.h"
#include "../AllASDModelParameterTypes.h"

#include <iostream>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

class DuncanChang_EL : public ElasticityBase<DuncanChang_EL> // CRTP on ElasticityBase
{
public:

    static constexpr const char* NAME = "DuncanChang_EL";

    DuncanChang_EL(): ElasticityBase<DuncanChang_EL>::ElasticityBase()  // Note the full-qualification of ElasticityBase through the scope resolution operator (::)
    {

    }

    ELASTICITY_MATRIX
    {

        double Eref = GET_PARAMETER_VALUE(ReferenceYoungsModulus);
        double nu = GET_PARAMETER_VALUE(PoissonsRatio);
        double pa = GET_PARAMETER_VALUE(ReferencePressure);
        double sigma3_max = GET_PARAMETER_VALUE(DuncanChang_MaxSigma3);
        double n = GET_PARAMETER_VALUE(DuncanChang_n);


        double p = -stress.meanStress();
        double q = stress.stressDeviatorQ();
        double theta = stress.lodeAngle();

        double sigma3 = -p + (2 * q) / 3 * cos(theta * M_PI / 180 + (2 * M_PI) / 3);

        if (sigma3 > sigma3_max)
        {
            sigma3 = sigma3_max;
        }

        double E = Eref * pa * pow(abs(sigma3) / pa, n);
        double lambda = ( nu * E ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double mu = E / ( 2.0 * ( 1.0 + nu ) );

        EE_MATRIX *= 0; //Use the base-class static elasticity matrix for return value

        EE_MATRIX(0, 0) = EE_MATRIX(1, 1) = EE_MATRIX(2, 2) = 2*mu + lambda;
        EE_MATRIX(0, 1) = EE_MATRIX(1, 0) = EE_MATRIX(0, 2) = EE_MATRIX(2, 0) = EE_MATRIX(1, 2) = EE_MATRIX(2, 1) = lambda;
        EE_MATRIX(3, 3) = mu;
        EE_MATRIX(4, 4) = mu;
        EE_MATRIX(5, 5) = mu;

        return EE_MATRIX;
    }

    using parameters_t = std::tuple<ReferenceYoungsModulus,PoissonsRatio,ReferencePressure,DuncanChang_MaxSigma3,DuncanChang_n>;


};


#endif


