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


#ifndef Linear_HardeningScalarV_H
#define Linear_HardeningScalarV_H


#include "../EvolvingVariable.h"
#include "../ASDPlasticMaterialGlobals.h" // Defines indices i,j,k,l,m,n,p,q,r,s and the kronecker_delta.


class LinearHardeningScalar_EV : public EvolvingVariable<double, LinearHardeningScalar_EV> //CRTP on LinearHardeningScalar_EV
{
public:

    LinearHardeningScalar_EV( double H_);

    LinearHardeningScalar_EV( double H_, double k0);

    const double& getDerivative(const VoightTensor6 &depsilon,
                                const VoightTensor6 &m,
                                const VoightTensor6& stress) const;
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);


private:
    double H;
    static double derivative; //Must return a reference.
};

 
#endif //Linear_HardeningScalarV_H