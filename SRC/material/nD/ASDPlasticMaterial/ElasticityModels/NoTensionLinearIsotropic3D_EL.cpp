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

#include "NoTensionLinearIsotropic3D_EL.h"
#include "Vector.h"

VoigtMatrix NoTensionLinearIsotropic3D_EL::Ee(3, 3, 3, 3, 0.0);


NoTensionLinearIsotropic3D_EL::NoTensionLinearIsotropic3D_EL(double E, double nu) : ElasticityBase<NoTensionLinearIsotropic3D_EL>::ElasticityBase()  // Note the full-qualification of ElasticityBase through the scope resolution operator (::)
{
    lambda = ( nu * E ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    mu = E / ( 2.0 * ( 1.0 + nu ) );
    // std::cout << "E  = " << E << std::endl;
    // std::cout << "nu = " << nu << std::endl;
}


VoigtMatrix& NoTensionLinearIsotropic3D_EL::operator()(const VoigtVector& stress) //See note on base class
{
    Ee *= 0; //Zero it. It may have values from another instance with different parameters;

    Ee( 0, 0, 0, 0 ) = lambda + 2 * mu;
    Ee( 0, 0, 1, 1 ) = lambda;
    Ee( 0, 0, 2, 2 ) = lambda;
    Ee( 0, 1, 0, 1 ) = mu;
    Ee( 0, 1, 1, 0 ) = mu;
    Ee( 0, 2, 0, 2 ) = mu;
    Ee( 0, 2, 2, 0 ) = mu;
    Ee( 1, 0, 0, 1 ) = mu;
    Ee( 1, 0, 1, 0 ) = mu;
    Ee( 1, 1, 0, 0 ) = lambda;
    Ee( 1, 1, 1, 1 ) = lambda + 2 * mu;
    Ee( 1, 1, 2, 2 ) = lambda;
    Ee( 1, 2, 1, 2 ) = mu;
    Ee( 1, 2, 2, 1 ) = mu;
    Ee( 2, 0, 0, 2 ) = mu;
    Ee( 2, 0, 2, 0 ) = mu;
    Ee( 2, 1, 1, 2 ) = mu;
    Ee( 2, 1, 2, 1 ) = mu;
    Ee( 2, 2, 0, 0 ) = lambda;
    Ee( 2, 2, 1, 1 ) = lambda;
    Ee( 2, 2, 2, 2 ) = lambda + 2 * mu;

    double p = -(stress(0, 0) + stress(1, 1) + stress(2, 2)) / 3;

    if (p < 0)
    {
        Ee *= 1e-5;
    }

    return Ee;
}

int NoTensionLinearIsotropic3D_EL::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(2);
    data(0) = lambda;
    data(1) = mu;

    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "NoTensionLinearIsotropic3D_EL::sendSelf() - Failed to send data. " << endl;
        return -1;
    }

    return 0;
}

int NoTensionLinearIsotropic3D_EL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(2);

    if (theChannel.receiveVector(0, commitTag, data) != 0)
    {
        cerr << "NoTensionLinearIsotropic3D_EL::recvSelf() - Failed to receive data. " << endl;
        return -1;
    }

    lambda = data(0);
    mu = data(1);

    return 0;
}
