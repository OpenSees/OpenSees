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

#ifndef NoTensionLinearIsotropic3D_EL_H
#define NoTensionLinearIsotropic3D_EL_H

#include "../../../ltensor/LTensor.h"
#include "../EvolvingVariable.h"
#include "../ElasticityBase.h"

#include <iostream>


class NoTensionLinearIsotropic3D_EL : public ElasticityBase<NoTensionLinearIsotropic3D_EL> // CRTP on ElasticityBase
{
public:
    NoTensionLinearIsotropic3D_EL(double E, double nu);

    VoigtMatrix& operator()(const VoigtVector& stress); //See note on base class

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:

    double lambda;
    double mu;
    static VoigtMatrix Ee;  //Provides class-wide storage, which avoids mallocs and allows const returning a const & to this object.

};



#endif
