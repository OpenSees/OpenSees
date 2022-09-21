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


#ifndef ASDPlasticMaterialGlobals_H
#define ASDPlasticMaterialGlobals_H

#include "EigenAPI.h"
#include <limits>
#include <tuple>
#include <string>
namespace ASDPlasticMaterialGlobals
{


// constexpr double SQRT_2_over_3  = 0.816496580928;//sqrt(2.0 / 3.0); //in case unsupported by compiler
constexpr double SQRT_2_over_3  = sqrt(2.0 / 3.0); 
// constexpr double SQRT_2_over_27 = 0.272165526976;//sqrt(2.0 / 27.0); //in case unsupported by compiler
constexpr double SQRT_2_over_27 = sqrt(2.0 / 27.0);

constexpr double MACHINE_EPSILON = std::numeric_limits<double>::epsilon();


// void printTensor(std::string const& name, DTensor2 const& v);
// void printTensor4(std::string const& name, DTensor4 const& v);
// std::tuple<double, double, double> getpqtheta(const DTensor2 &mystress);
// std::tuple<double, double, double> getI1J2J3(const DTensor2 &mystress);
// bool inverse4thTensor(DTensor4 const& rhs, DTensor4& ret);

// void dJ2_dsigma_ij(const DTensor2& sigma, DTensor2 &result);   // Stress derivative of second deviatoric stress invariant
// void dJ3_dsigma_ij(const DTensor2& sigma, DTensor2 &result);   // Stress derivative of third deviatoric stress invariant
// void dq_dsigma_ij(const DTensor2& sigma, DTensor2 &result);   // Stress derivative of deviatoric stress q
// void dtheta_dsigma_ij(const DTensor2& sigma, DTensor2 &result);   // Stress derivative of Lode angle


// Macaulay Bracket < >  operator. (Integral of Heaviside function)
inline double macaulay_bracket(double x)
{
    return x > 0 ? x : 0;
}
// static int Nsteps = 0 ;  // Removes warning
}

#endif