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



#ifndef ClassicElastoplasticMaterial_H
#define ClassicElastoplasticMaterial_H


#include "NDMaterialLT.h"
#include <G3Globals.h>
#include <iostream>
#include "../../ltensor/LTensor.h"
#include <Channel.h>
#include <limits>

#include "CEPTraits.h"
#include <type_traits>

#include "ClassicElastoplasticityGlobals.h"
// for print p, q, theta
#include <tuple>
// for debugging printing
#include <fstream>

#define ClassicElastoplasticMaterial_MAXITER_BRENT 20
#define TOLERANCE1 1e-6
/*

YieldFunctionType, ElasticityType, PlasticFlowType, HardeningLawType are the user-defined classes
that provide the different components of this elastoplastic model. Called Base Elastoplastic Template
Classes (BET Classes)

* All BET Classes  must provide:

C++ "Rule of 5"
 + An empty constructor (default constructor)
 + A copy constructor
 + A copy assignment operator
 + A move constructor
 + A move assign operator
...and
 + A name member function that returns a std::string. For example.
    SomeYieldFunctionType yf(...parameters...);
    std::string name = yf.name();
 + A getObjectSize() function. Returnig size in bytes of the object in memory.


* YieldFunctionType must additionally provide



* ElasticityType must additionally provide



* PlasticFlowType must additionally provide





*/



// This anonymous namespace avoids multiple-definition name clashes.
// namespace
// {

// using namespace std;

// void printTensor(std::string const& name, DTensor2 const& v)
// {

//     // All 9 elements in one line.
//     // std::cout << name << " = [ " ;
//     // std::cout << v(0, 0) << " "
//     //           << v(1, 1) << " "
//     //           << v(2, 2) << " "
//     //           << v(0, 1) << " "
//     //           << v(0, 2) << " "
//     //           << v(1, 2) << " "
//     //           << v(1, 0) << " "
//     //           << v(2, 0) << " "
//     //           << v(2, 1) << " ]" << std::endl;

//     // This is in good format but take 3 lines.
//     // stderr will print immediately, not like cout (may be reordered by CPU).
//     fprintf(stderr, "%s = \n", name.c_str());
//     fprintf(stderr, "[%16.8f \t %16.8f \t %16.8f \n",   v(0, 0), v(0, 1), v(0, 2));
//     fprintf(stderr, " %16.8f \t %16.8f \t %16.8f \n",   v(1, 0), v(1, 1), v(1, 2));
//     fprintf(stderr, " %16.8f \t %16.8f \t %16.8f] \n",  v(2, 0), v(2, 1), v(2, 2));
// }

// void printTensor4(std::string const& name, DTensor4 const& v)
// {


//     std::cout << name << " = [ " ;
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//             for (int kk = 0; kk < 3; kk++)
//                 for (int ll = 0; ll < 3; ll++)
//                 {
//                     std::cout << v(ii, jj, kk, ll ) << " ";
//                 }
//     std::cout << "]" << std::endl;
// }
// // ------------------------------------------------------------
// // refer to https://en.wikipedia.org/wiki/Yield_surface
// // p = 1/3 * I1
// // q = sqrt(3* J2)
// // cos(3*theta) = 3/2 * sqrt(3) * J3 / J2^(3/2)
// // ------------------------------------------------------------
// std::tuple<double, double, double> getpqtheta(const DTensor2 &mystress)
// {
//     // ------------------------------------------------------------
//     // preliminary
//     const double I1 = mystress(0, 0) + mystress(1, 1) + mystress(2, 2);
//     const double sigma_m = I1 / 3.0;
//     DTensor2 s = mystress;
//     s(0, 0) -= sigma_m;
//     s(1, 1) -= sigma_m;
//     s(2, 2) -= sigma_m;
//     // J2=0.5*s(i,j)*s(i,j)
//     const double J2 = 0.5 * (
//                           s(0, 0) * s(0, 0) +  s(0, 1) * s(0, 1)  + s(0, 2) * s(0, 2)
//                           +   s(1, 0) * s(1, 0) +  s(1, 1) * s(1, 1)  + s(1, 2) * s(1, 2)
//                           +   s(2, 0) * s(2, 0) +  s(2, 1) * s(2, 1)  + s(2, 2) * s(2, 2)   );
//     // 3by3 Determinant: Refer to http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Tensors/Tensors.htm
//     const double J3 = s(0, 0) * (s(1, 1) * s(2, 2) - s(1, 2) * s(2, 1))
//                       +  s(0, 1) * (s(1, 2) * s(2, 0) - s(1, 0) * s(2, 2))
//                       +  s(0, 2) * (s(1, 0) * s(2, 1) - s(2, 0) * s(1, 1));
//     // ------------------------------------------------------------

//     // (1) calculate p
//     const double p = -I1 / 3;
//     // (2) calculate q
//     const double q = sqrt(3 * J2);
//     // (3) calculate theta
//     const double cos3theta = 1.5 * sqrt(3) * J3 / pow(J2, 1.5);
//     // const static double PI = 3.14159265358979323846; //20 digits from Wiki.
//     const double theta = acos(cos3theta) / 3.0 * 180 / M_PI; //theta is in degree.

//     return std::make_tuple(p, q, theta);
//     // return result;
// }

// bool inverse4thTensor(DTensor4 const& rhs, DTensor4& ret)
// {
//     using namespace ClassicElastoplasticityGlobals;
//     static DTensor2 intermediate_matrix(9, 9, 0.0);
//     intermediate_matrix *= 0;
//     // static DTensor4 ret(3,3,3,3,0.0);
//     int m41 = 0,  m42 = 0;
//     // (1). convert 4th order Tensor to matrix
//     for ( int c44 = 1 ; c44 <= 3 ; c44++ )
//         for ( int c43 = 1 ; c43 <= 3 ; c43++ )
//             for ( int c42 = 1 ; c42 <= 3 ; c42++ )
//                 for ( int c41 = 1 ; c41 <= 3 ; c41++ )
//                 {
//                     m41 = 3 * (c41 - 1) + c42;
//                     m42 = 3 * (c43 - 1) + c44;

//                     intermediate_matrix( m41 - 1, m42 - 1 ) = rhs( c41 - 1, c42 - 1, c43 - 1, c44 - 1 );

//                 }
//     // (2). Inverse the matrix .
//     double det = intermediate_matrix.compute_Determinant();

//     // fprintf(stderr, "----> det: %f\n", det);

//     if (det < MACHINE_EPSILON)
//     {
//         cout << "ClassicElastoplasticMaterial matrix T is not invertible to get the consistent stiffness tensor. Use the inconsistent stiffness instead! " << endl;
//         return false;
//     }
//     else
//     {
//         static DTensor2 inv_matrix(3, 3, 0.0);
//         inv_matrix = intermediate_matrix.Inv();
//         // (3). convert Matrix to 4th order tensor
//         for ( int c44 = 1 ; c44 <= 3 ; c44++ )
//             for ( int c43 = 1 ; c43 <= 3 ; c43++ )
//                 for ( int c42 = 1 ; c42 <= 3 ; c42++ )
//                     for ( int c41 = 1 ; c41 <= 3 ; c41++ )
//                     {
//                         m41 = 3 * (c41 - 1) + c42;
//                         m42 = 3 * (c43 - 1) + c44;

//                         ret(c41 - 1, c42 - 1, c43 - 1, c44 - 1) = inv_matrix(m41 - 1, m42 - 1);
//                     }
//         return true;
//     }

// }

// static int Nsteps = 0 ;
// }




template < class ElasticityType, class YieldFunctionType, class PlasticFlowType, class MaterialInternalVariablesType, int thisClassTag, class T >
class ClassicElastoplasticMaterial : public NDMaterialLT
{

public:



//==================================================================================================
//  Void constructor
//==================================================================================================
//Empty constructor is needed for parallel, sets all the pointers to zero. Resources will be
//    allocated later
// ClassicElastoplasticMaterial( )
//     : NDMaterialLT(0, thisClassTag),
//       NDMaterialLT(tag, thisClassTag),
//       TrialStress(3, 3, 0.0),
//       TrialStrain(3, 3, 0.0),
//       TrialPlastic_Strain(3, 3, 0.0),
//       CommitStress(3, 3, 0.0),
//       CommitStrain(3, 3, 0.0),
//       CommitPlastic_Strain(3, 3, 0.0),
//       Stiffness(3, 3, 3, 3, 0.0),
//       yf(),
//       et(),
//       pf(),
//       internal_variables()
// {

// }



//==================================================================================================
//  Main constructor
//==================================================================================================
// Constructor... invokes copy constructor for base elastoplastic template parameters
    ClassicElastoplasticMaterial(int tag,
                                 double rho_,
                                 double p0,
                                 const YieldFunctionType& yf_in,
                                 const ElasticityType& et_in,
                                 const PlasticFlowType& pf_in,
                                 const MaterialInternalVariablesType& internal_variables_in
                                )
        : NDMaterialLT(tag, thisClassTag),
          rho(rho_),
          TrialStrain(3, 3, 0.0),
          TrialStress(3, 3, 0.0),
          TrialPlastic_Strain(3, 3, 0.0),
          CommitStress(3, 3, 0.0),
          CommitStrain(3, 3, 0.0),
          CommitPlastic_Strain(3, 3, 0.0),
          yf(yf_in),
          et(et_in),
          pf(pf_in),
          internal_variables(internal_variables_in)
    {
        // cout << "Copy ctor in CEP" << endl;
        //Set initial stress to some value. Needed for models
        // like Drucker-Prager which blow up at sigma_ii = 0 :)
        TrialStress(0, 0) =  - p0;
        TrialStress(1, 1) =  - p0;
        TrialStress(2, 2) = - p0;

        CommitStress(0, 0) = - p0;
        CommitStress(1, 1) = - p0;
        CommitStress(2, 2) = - p0;

        first_step = true;
    }



//==================================================================================================
//  Destructor
//==================================================================================================
    ~ClassicElastoplasticMaterial(void)
    {

    }




//==================================================================================================
//  Class type function
//==================================================================================================
// Returns a null-terminated character string with the name of the class
    const char *getClassType(void) const
    {
        std::string name("ClassicElastoplasticMaterial");

        return name.c_str();
    };

    double getRho(void)
    {
        return rho;
    }

    double getPressure(void)
    {
        using namespace ClassicElastoplasticityGlobals;
        return -CommitStress(i, i) / 3;
    }

//==================================================================================================
//  Setters for components (invoke assign copy on BET Classes)
//==================================================================================================

    // void setYieldFunction(YieldFunctionType& yf_)
    // {
    //     yf = std::move(yf_);
    // }

    // void setElasticity(ElasticityType& et_)
    // {
    //     et = std::move(et_);
    // }

    // void setPlasticFlow(PlasticFlowType& pf_)
    // {
    //     pf = std::move(pf_);
    // }


//==================================================================================================
//  Set Trial strain and trial strain increment
//==================================================================================================
    // For total strain-based elements.
    // Receives the current total strain at a GP.
    // This function then computes the incremental strain (subtracting from the committed one)
    // and sets the increment.
    // Returns a success flag from the call to setTrialStrainIncr
    int setTrialStrain( const DTensor2 &v )
    {
        using namespace ClassicElastoplasticityGlobals;
        static DTensor2 result( 3, 3, 0.0 );
        result *= 0;

        TrialStrain(i, j) = v(i, j);
        result( i, j ) = v( i, j ) - CommitStrain( i, j );

        return setTrialStrainIncr( result );
    }

    // Directly sets the trial strain increment and does an explicit or implicit step.
    // Returns a flag depending on the result of the step.
    int setTrialStrainIncr( const DTensor2 &strain_increment )
    {
        using namespace ClassicElastoplasticityGlobals;

        int exitflag = 0;

        // ==========================================
        // Check for quick return:
        // If strain_increment are all zeroes. Return
        // ==========================================
        // The quick return is useful for displacement_control.
        // In displacement_control, we have to do double update.
        // ==========================================
        double max_component = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (max_component < abs(strain_increment(i, j)))
                {
                    max_component = abs(strain_increment(i, j));
                }
        // ==========================================
        // We do not use below(eps_norm) because 1E-8 will go to 1E-16. --> False return.
        // double eps_norm = strain_increment(i, j) * strain_increment(i, j);
        // ==========================================
        if (max_component == 0)
        {
            return exitflag;
        }
        // ==========================================

        switch (this->constitutive_integration_method)
        {
        case NDMaterialLT_Constitutive_Integration_Method::Not_Set :
            exitflag = -1;
            cerr << "CEP::setTrialStrainIncr - Integration method not set!\n" ;
            break;
        case NDMaterialLT_Constitutive_Integration_Method::Forward_Euler :
            exitflag = this->Forward_Euler(strain_increment);
            break;
        case NDMaterialLT_Constitutive_Integration_Method::Forward_Euler_Subincrement :
            exitflag = this->Forward_Euler_Subincrement(strain_increment);
            break;
        case NDMaterialLT_Constitutive_Integration_Method::Backward_Euler :
            exitflag = this->Backward_Euler(strain_increment);;
            break;
        case NDMaterialLT_Constitutive_Integration_Method::Backward_Euler_ddlambda :
            exitflag = this->Backward_Euler_ddlambda(strain_increment);;
            break;
        case NDMaterialLT_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement :
            exitflag = this->Backward_Euler_ddlambda_Subincrement(strain_increment);;
            break;
        // case NDMaterialLT_Constitutive_Integration_Method::Forward_Euler_Crisfield :
        //     exitflag = this->Forward_Euler(strain_increment, true);
        //     break;
        // case NDMaterialLT_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield :
        //     exitflag = this->Multistep_Forward_Euler(strain_increment, true);
        //     break;
        // case NDMaterialLT_Constitutive_Integration_Method::Modified_Euler_Error_Control :
        //     exitflag = this->Modified_Euler_Error_Control(strain_increment);
        //     break;
        // case NDMaterialLT_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control :
        //     // exitflag = this->Runge_Kutta_45_Error_Control(strain_increment);;
        //     break;
        default:
            cerr << "ClassicElastoplasticMaterial::setTrialStrainIncr - Integration method not available!\n" ;
            exitflag = -1;
        }

        return exitflag;
    }

//==================================================================================================
//  Getters
//==================================================================================================

    const DTensor4 &getTangentTensor( void )
    {
        using namespace ClassicElastoplasticityGlobals;

        // double yf_val = yf(TrialStress);

        // DTensor4& Eelastic = et(TrialStress);


        // if (yf_val <= 0.0)
        // {
        //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        // }
        // else
        // {
        //     const DTensor2& n = yf.df_dsigma_ij(intersection_stress);
        //     const DTensor2& m = pf(depsilon_elpl, intersection_stress);

        //     double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, depsilon_elpl,  intersection_stress);
        //     double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

        //     //Compute tangent stiffness
        //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / den;
        // }

        if (first_step)
        {
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
            first_step = false;
        }

        return Stiffness;
    }

    const DTensor2 &getStressTensor( void )
    {
        return TrialStress;
    }

    const DTensor2 &getStrainTensor( void )
    {
        return TrialStrain;
    }

    const DTensor2 &getPlasticStrainTensor( void )
    {
        return TrialPlastic_Strain;
    }

    const DTensor2  &getCommittedStressTensor(void)
    {
        return CommitStress;
    }

    const DTensor2 &getCommittedStrainTensor(void)
    {
        return CommitStrain;
    }

    const DTensor2 &getCommittedPlasticStrainTensor(void)
    {
        return CommitPlastic_Strain;
    }



//==================================================================================================
//  State commiting and reversion
//==================================================================================================

    // Forwards the trial variables to commited variables and calls commit state on the
    // BET classes
    int commitState(void)
    {

        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;
        CommitStress(i, j) = TrialStress(i, j);
        CommitStrain(i, j) = TrialStrain(i, j);
        CommitPlastic_Strain(i, j) = TrialPlastic_Strain(i, j);

        internal_variables.commit();

        // if (first_step)
        // {
        //     first_step = false;
        // }

        return errorcode;
    }

    //Reverts the commited variables to the trials and calls revert on BET Classes.
    int revertToLastCommit(void)
    {
        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;
        TrialStress(i, j) = CommitStress(i, j);
        TrialStrain(i, j) = CommitStrain(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        internal_variables.revert();

        return errorcode;
    }

    int revertToStart(void)
    {
        first_step = true;
        return 0;
    }

    NDMaterialLT *getCopy(void)
    {
        ClassicElastoplasticMaterial< ElasticityType,
                                      YieldFunctionType,
                                      PlasticFlowType,
                                      MaterialInternalVariablesType,
                                      thisClassTag,
                                      T > * newmaterial = new T( this->getTag(), rho, this->getPressure(), yf, et, pf, internal_variables);
        newmaterial->internal_variables.setVars(this->internal_variables);
        newmaterial->TrialStrain = (this->TrialStrain);
        newmaterial->TrialStress = (this->TrialStress);
        newmaterial->TrialPlastic_Strain = (this->TrialPlastic_Strain);
        newmaterial->CommitStress = (this->CommitStress);
        newmaterial->CommitStrain = (this->CommitStrain);
        newmaterial->CommitPlastic_Strain = (this->CommitPlastic_Strain);
        // newmaterial->Stiffness = (this->Stiffness);
        return newmaterial;
    }

    // NDMaterialLT *getCopy(const char *code);
    // const char *getType(void) const;

    int sendSelf(int commitTag, Channel &theChannel)
    {
        int pos = 0;
        static Vector data(1 + 9 * 6 + 4); // rho and all the DTensors2 get packed into one vector

        // double rho;
        data(pos++) = rho;
        // DTensor2 TrialStrain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = TrialStrain(i, j);
            }
        // DTensor2 TrialStress;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = TrialStress(i, j);
            }
        // DTensor2 TrialPlastic_Strain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = TrialPlastic_Strain(i, j);
            }
        // DTensor2 CommitStress;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = CommitStress(i, j);
            }
        // DTensor2 CommitStrain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = CommitStrain(i, j);
            }
        // DTensor2 CommitPlastic_Strain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                data(pos++) = CommitPlastic_Strain(i, j);
            }

        data(pos++)  = NDMaterialLT::f_relative_tol        ;
        data(pos++)  = NDMaterialLT::stress_relative_tol   ;
        data(pos++)  = NDMaterialLT::n_max_iterations      ;
        data(pos++)  = (double) NDMaterialLT::constitutive_integration_method ;





        // cout << "Sending data" << endl;
        if (theChannel.sendVector(0, commitTag, data) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::sendSelf() - Failed sending data" << endl;
            return -1;
        }

        // cout << "Sending elasticity" << endl;
        // ElasticityType    et;
        if (et.sendSelf(commitTag, theChannel) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::sendSelf() - Failed sending elasticity data" << endl;
            return -1;
        }

        // cout << "Sending variables" << endl;
        // MaterialInternalVariablesType internal_variables;
        if (internal_variables.sendSelf(commitTag, theChannel) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::sendSelf() - Failed sending internal variables data" << endl;
            return -1;
        }

        return 0;
    }
    int receiveSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
    {
        // cout << "Receiving data" << endl;
        static Vector data(1 + 9 * 6 + 4); // rho and all the DTensors2 get packed into one vector
        if (theChannel.receiveVector(0, commitTag, data) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::receiveSelf() - Failed receiving data" << endl;
            return -1;
        }

        int pos = 0;
        // double rho;
        rho = data(pos++);
        // DTensor2 TrialStrain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TrialStrain(i, j) = data(pos++);
            }
        // DTensor2 TrialStress;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TrialStress(i, j) = data(pos++);
            }
        // DTensor2 TrialPlastic_Strain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                TrialPlastic_Strain(i, j) = data(pos++);
            }
        // DTensor2 CommitStress;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                CommitStress(i, j) = data(pos++);
            }
        // DTensor2 CommitStrain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                CommitStrain(i, j) = data(pos++);
            }
        // DTensor2 CommitPlastic_Strain;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                CommitPlastic_Strain(i, j) = data(pos++);
            }

        NDMaterialLT::f_relative_tol         = data(pos++) ;
        NDMaterialLT::stress_relative_tol    = data(pos++) ;
        NDMaterialLT::n_max_iterations       = data(pos++) ;
        NDMaterialLT::constitutive_integration_method       = (NDMaterialLT_Constitutive_Integration_Method) data(pos++) ;



        // cout << "Receiving elasticity" << endl;
        // ElasticityType    et;
        if (et.receiveSelf(commitTag, theChannel, theBroker) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::receiveSelf() - Failed receiving elasticity data" << endl;
            return -1;
        }
        // cout << "Receiving variables" << endl;
        // MaterialInternalVariablesType internal_variables;
        if (internal_variables.receiveSelf(commitTag, theChannel, theBroker) != 0)
        {
            cerr << "ClassicElastoplasticMaterial::receiveSelf() - Failed receiving internal variables data" << endl;
            return -1;
        }

        return 0;
    }
    void Print(ostream &s, int flag = 0)
    {
        using namespace ClassicElastoplasticityGlobals;
        static DTensor2 zeroTensor(3, 3, 0);
        s << "TrialStress =  " <<  TrialStress << endl;
        s << "CommitStress = " <<  CommitStress << endl;
        s << "yf  = " << yf(CommitStress) << endl;
        const DTensor2& n = yf.df_dsigma_ij(CommitStress);
        const DTensor2& m = pf(zeroTensor, CommitStress);
        s << "n = " <<  n  << endl;
        s << "m = " <<  m  << endl;
        s << "xi_star_h_star  = " << yf.xi_star_h_star(zeroTensor, m, CommitStress) << endl;
    }

    int getObjectSize()
    {
        int size = 0;

        // 6 3x3 DTensor2s and 1 DTensor4 (3x3x3x3)
        size += (3 * 3 * 6 + 3 * 3 * 3 * 3) * sizeof(double);

        //Four pointers
        size += 4 * sizeof(YieldFunctionType*);

        //Whatever the base components size is
        size += sizeof(yf);//yf->getObjectSize();
        size += sizeof(et);//et->getObjectSize();
        size += sizeof(pf);//pf->getObjectSize();
        size += sizeof(internal_variables);//internal_variables->getObjectSize();
        size += sizeof(NDMaterialLT);

        // size += static_cast<T*>(this)->getObjectSize();

        return size;
    }

    void setTrialStress(const DTensor2& stress)
    {
        using namespace ClassicElastoplasticityGlobals;
        TrialStress(i, j) = stress(i, j);
    }

protected:

    void setTrialPlastic_Strain(const DTensor2& strain)
    {
        using namespace ClassicElastoplasticityGlobals;
        TrialPlastic_Strain(i, j) = strain(i, j);
    }

    void setCommitStress(const DTensor2& stress)
    {
        using namespace ClassicElastoplasticityGlobals;
        CommitStress(i, j) = stress(i, j);
    }

    void setCommitStrain(const DTensor2& strain)
    {
        using namespace ClassicElastoplasticityGlobals;
        CommitStrain(i, j) = strain(i, j);
    }

    void setCommitPlastic_Strain(const DTensor2& strain)
    {
        using namespace ClassicElastoplasticityGlobals;
        CommitPlastic_Strain(i, j) = strain(i, j);
    }

    void setStiffness(const DTensor4& stiff)
    {
        using namespace ClassicElastoplasticityGlobals;
        Stiffness = stiff;
    }


    void setStressTensor(DTensor2 &stress)
    {
        using namespace ClassicElastoplasticityGlobals;
        CommitStress(i, j) = stress(i, j);
        TrialStress(i, j) = stress(i, j);
        return;
    }

private:


    // int Forward_Euler(const DTensor2 &strain_incr, bool const& with`_return2yield_surface)
    int Forward_Euler(const DTensor2 &strain_incr)
    {
        using namespace ClassicElastoplasticityGlobals;
        // ----------------------------------------------------------------
        // Print p, q, theta for debug--------------------------------------
        // ----------------------------------------------------------------
        //double enter_yf = yf(CommitStress);
        //double enter_p,enter_q,enter_theta;
        //std::tie(enter_p,enter_q,enter_theta) = getpqtheta(CommitStress);
        //fprintf(stderr, "--------------------------------------------\n ");
        //fprintf(stderr, "--------------------------------------------\n ");
        //fprintf(stderr, "-----------start iteration step %d----------\n",Nsteps );
        //fprintf(stderr, "When Enter Euler Step\n");
        //fprintf(stderr, "yf start (:<0) = %16.8f \n" , enter_yf );
        //fprintf(stderr, "    enter_p    = %16.8f \n"  , enter_p );
        //fprintf(stderr, "    enter_q    = %16.8f \n"  , enter_q );
        //fprintf(stderr, "enter_theta    = %16.8f \n"  , enter_theta );
        // ----------------------------------------------------------------

        int errorcode = 0;
        static DTensor2 depsilon(3, 3, 0);
        depsilon *= 0;
        depsilon(i, j) = strain_incr(i, j);

        const DTensor2& sigma = CommitStress;
        const DTensor2& epsilon = CommitStrain;

        internal_variables.revert();
        internal_variables.commit_tmp();

        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        DTensor4& Eelastic = et(sigma);
        Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

        dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

        TrialStress(i, j) = sigma(i, j) + dsigma(i, j);
        TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        double yf_val_start = yf(sigma);
        double yf_val_end = yf(TrialStress);

        DTensor2& start_stress = CommitStress;
        DTensor2& end_stress = TrialStress;

        intersection_stress(i, j) = start_stress(i, j);

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        }
        else  //Plasticity
        {
            depsilon_elpl(i, j) = depsilon(i, j);
            if (yf_val_start < 0)
            {
                double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );
                intersection_stress(i, j) = start_stress(i, j) * (1 - intersection_factor) + end_stress(i, j) * intersection_factor;
                intersection_strain(i, j) = epsilon(i, j)  + depsilon(i, j) * intersection_factor;
                depsilon_elpl(i, j) = (1 - intersection_factor) * depsilon(i, j);
            }

            TrialStress(i, j) = intersection_stress(i, j);

            Eelastic = et(intersection_stress);
            TrialStress(i, j)  += Eelastic(i, j, k, l) * depsilon_elpl(k, l);

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const DTensor2& n = yf.df_dsigma_ij(intersection_stress);
            const DTensor2& m = pf(depsilon_elpl, intersection_stress);

            double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, m,  intersection_stress);
            double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

            //Compute the plastic multiplier
            // Yuan and Boris 10Sep2016, change from "==0"
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                printTensor("m", m);
                printTensor("n", n);
                cout << "xi_star_h_star" << xi_star_h_star << endl;
                cout << "den" << den << endl;
                printTensor("depsilon_elpl", depsilon_elpl);
                return -1;
            }
            double dLambda =  n(i, j) * Eelastic(i, j, k, l) * depsilon_elpl(k, l);
            dLambda /= den;

            // if (dLambda <= 0)
            // {
            //     cout << "CEP - dLambda = " << dLambda << " <= 0\n";
            //     printTensor("m", m);
            //     printTensor("n", n);
            //     cout << "xi_star_h_star = " << xi_star_h_star << endl;
            //     cout << "den = " << den << endl;
            //     printTensor("depsilon_elpl", depsilon_elpl);
            //     // return -1;
            // }

            // Update the trial plastic strain.
            TrialPlastic_Strain(i, j) += dLambda * m(i, j);
            // Update the internal variables (k and alpha)
            internal_variables.evolve(dLambda, depsilon_elpl, m, intersection_stress);
            internal_variables.commit_tmp();
            // vonMises does NOT enter this part.
            // DruckerPrager requries this part.
            if (yf.hasCorner() && yf.in_Apex(TrialStress))
            {
                static DTensor2 small_stress(3, 3, 0.0);
                small_stress *= 0;
                // The small value 50*Pa refers to the lowest confinement test:
                // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
                double DP_k = yf.get_k();
                double DP_p = 50 ;
                // To make it on the yield surface, the q is equal to k*p.
                double DP_q = DP_k * DP_p ;
                // Assume the triaxial conditions sigma_2 = sigma_3.
                small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
                small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
                small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;
                static DTensor2 dstress(3, 3, 0.0);
                dstress(i, j) = small_stress(i, j) - sigma(i, j);
                static DTensor2 depsilon_Inv(3, 3, 0.0);
                depsilon_Inv = depsilon.Inv();

                // Return results (member variables) :
                Stiffness(i, j, k, l) = dstress(i, j) * depsilon_Inv(k, l);
                TrialStress(i, j) = small_stress(i, j);
                // plastic_strain and internal variables are already updated.
                return 0;
            }

            //Correct the trial stress
            TrialStress(i, j) = TrialStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);

            // Calculate the stiffness for the global iteration:
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / den;


            // // ============================================================================================
            // // Add the additional step: returning to the yield surface.
            // // This algorithm is based on Crisfield(1996). Page 171. Section 6.6.3
            // // After this step, the TrialStress(solution), TrialPlastic_Strain, and Stiffness will be updated to the yield surface.
            // // ============================================================================================
            // if(with_return2yield_surface){
            // if (true)
            // {
            //     // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
            //     // Make surface the internal variables are already updated. And then, return to the yield surface.
            //     double yf_val_after_corrector = yf(TrialStress);
            //     const DTensor2& n_after_corrector = yf.df_dsigma_ij(TrialStress);
            //     const DTensor2& m_after_corrector = pf(depsilon_elpl, TrialStress);
            //     // In the function below, depsilon_elpl is actually not used at all in xi_star_h_star
            //     double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m_after_corrector,  TrialStress);
            //     double dLambda_after_corrector = yf_val_after_corrector / (
            //                                          n_after_corrector(i, j) * Eelastic(i, j, k, l) * m_after_corrector(k, l) - xi_star_h_star_after_corrector
            //                                      );
            //     TrialStress(i, j) = TrialStress(i, j) - dLambda_after_corrector * Eelastic(i, j, k, l) * m_after_corrector(k, l);
            //     TrialPlastic_Strain(i, j) += dLambda_after_corrector * m_after_corrector(i, j);

            //     double den_after_corrector = n_after_corrector(p, q) * Eelastic(p, q, r, s) * m_after_corrector(r, s) - xi_star_h_star_after_corrector;
            //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m_after_corrector(p, q)) * (n_after_corrector(r, s) * Eelastic(r, s, k, l) ) / den_after_corrector;
            // }
            // // ============================================================================================
            // // ============================================================================================

            double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "Numeric error!\n";
                printTensor("TrialStress = " , TrialStress);
                printTensor("CommitStress = " , CommitStress);
                printTensor("depsilon = " , depsilon);
                printTensor("dsigma   = " , dsigma);
                printTensor("intersection_stress = " , intersection_stress);
                printTensor4("Eelastic = " , Eelastic);
                printTensor4("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor("n = " , n );
                printTensor("m = " , m );
                cout << "xi_star_h_star  = " << xi_star_h_star << endl;
                cout << "den = " << den << endl;
                cout << "dLambda = " << dLambda << endl;

                errorcode = -1;
            }

        }

        return errorcode;
    }



    // int Forward_Euler_Subincrement(const DTensor2 &strain_incr, bool const& with_return2yield_surface)
    int Forward_Euler_Subincrement(const DTensor2 &strain_incr)
    {
        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;

        static DTensor2 depsilon(3, 3, 0);
        depsilon *= 0;
        depsilon(i, j) = strain_incr(i, j);

        const DTensor2& sigma = CommitStress;
        const DTensor2& epsilon = CommitStrain;
        internal_variables.revert();
        internal_variables.commit_tmp();

        dsigma *= 0;//Zero-out the stress increment tensor
        intersection_stress *= 0;
        intersection_strain *= 0;

        DTensor4& Eelastic = et(sigma);
        Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

        dsigma(i, j) += Eelastic(i, j, k, l) * depsilon(k, l);

        TrialStress(i, j) = sigma(i, j) + dsigma(i, j);
        TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        double yf_val_start = yf(sigma);
        double yf_val_end = yf(TrialStress);

        DTensor2& start_stress = CommitStress;
        DTensor2& end_stress = TrialStress;

        intersection_stress(i, j) = start_stress(i, j);


        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        }
        else  //Plasticity
        {
            depsilon_elpl(i, j) = depsilon(i, j);
            if (yf_val_start < 0)
            {
                double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );
                // cout << " intersection_factor = " << intersection_factor << endl;
                intersection_stress(i, j) = start_stress(i, j) * (1 - intersection_factor) + end_stress(i, j) * intersection_factor;
                intersection_strain(i, j) = epsilon(i, j)  + depsilon(i, j) * intersection_factor;
                depsilon_elpl(i, j) = (1 - intersection_factor) * depsilon(i, j);
            }

            TrialStress(i, j) = intersection_stress(i, j);


            int Nsubsteps = this-> n_max_iterations;
            static DTensor2 sub_depsilon_elpl(3, 3, 0);

            sub_depsilon_elpl *= 0;
            sub_depsilon_elpl(i, j) = depsilon_elpl(i, j) / Nsubsteps;

            // ====================================
            // Initialize the stiffness to be zero.
            // In multiStep, this entire stiffness should be the average of the substep
            // stiffness to match the real stiffness in one global step.
            Stiffness *= 0;
            // ====================================
            for (int iteration = 0; iteration < Nsubsteps; iteration++)
            {
                internal_variables.revert_tmp();
                Eelastic = et(TrialStress);

                //Compute normal to YF (n) and Plastic Flow direction (m) at the starting point.
                // Now TrialStress is at the starting point.
                const DTensor2& n = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m = pf(sub_depsilon_elpl, TrialStress);

                double xi_star_h_star = yf.xi_star_h_star( sub_depsilon_elpl, m,  TrialStress);
                double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

                //Compute the plastic multiplier
                // if (abs(den) < MACHINE_EPSILON)
                // {
                //     cout << "CEP - den = 0\n";
                //     printTensor("m", m);
                //     printTensor("n", n);
                //     cout << "xi_star_h_star" << xi_star_h_star << endl;
                //     cout << "den" << den << endl;
                //     printTensor("sub_depsilon_elpl", sub_depsilon_elpl);
                //     return -1;
                // }
                double dLambda =  n(i, j) * Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
                dLambda /= den;

                // Add the elastic predictor.
                TrialStress(i, j)  += Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l) ;


                static DTensor4 Stiffness_substep(3, 3, 3, 3, 0.0);
                Stiffness_substep *= 0;
                // vonMises does NOT enter this part.
                // DruckerPrager requries this part.
                if (yf.hasCorner() && yf.in_Apex(TrialStress))
                {
                    static DTensor2 small_stress(3, 3, 0.0);
                    small_stress *= 0;
                    // The small value 50*Pa refers to the lowest confinement test:
                    // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
                    double DP_k = yf.get_k();
                    double DP_p = 50 ;
                    // To make it on the yield surface, the q is equal to k*p.
                    double DP_q = DP_k * DP_p ;
                    // Assume the triaxial conditions sigma_2 = sigma_3.
                    small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
                    small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
                    small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

                    static DTensor2 Predictor_Stress(3, 3, 0.0);
                    Predictor_Stress(i, j) = TrialStress(i, j);
                    // (1) Update the trial stress
                    TrialStress(i, j) = small_stress(i, j);

                    // (2) Update the trial plastic strain
                    const DTensor2& n = yf.df_dsigma_ij(small_stress);
                    const DTensor2& m = pf(depsilon, small_stress);
                    const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
                    double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                    double dLambda = yf(Predictor_Stress) / denominator;
                    TrialPlastic_Strain(i, j) += dLambda * m(i, j);

                    // (3) Update the internal variables
                    internal_variables.evolve(dLambda, depsilon, m, TrialStress);
                    internal_variables.commit_tmp();

                    // (4) Update the stiffness
                    static DTensor2 dstress(3, 3, 0.0);
                    dstress(i, j) = small_stress(i, j) - sigma(i, j);
                    static DTensor2 depsilon_Inv(3, 3, 0.0);
                    depsilon_Inv = depsilon.Inv();
                    Stiffness_substep(i, j, k, l) = dstress(i, j) * depsilon_Inv(k, l);

                    // Go to the next sub-step directly:
                    continue;
                }

                // Update the trial plastic strain.
                TrialPlastic_Strain(i, j) += dLambda * m(i, j);
                // Update the internal variables
                internal_variables.evolve(dLambda, sub_depsilon_elpl, m, TrialStress);
                internal_variables.commit_tmp();

                //Correct the trial stress
                TrialStress(i, j) = TrialStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
                Stiffness_substep(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / den;
                Stiffness(i, j, k, l) += Stiffness_substep(i, j, k, l) / Nsubsteps;

                // // ===========================================================
                // // According to Crisfield(1996) Page 172. Section 6.6.4.
                // // The stiffness for the entire step should be the average of the substeps'
                // // stiffness
                // // ===========================================================
                // if(!with_return2yield_surface){
                //     Stiffness_substep(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / den;
                //     Stiffness(i, j, k, l) += Stiffness_substep(i, j, k, l)/Nsubsteps;
                // }else{
                //     // ==============================================================
                //     // Add the additional step: returning to the yield surface.
                //     // This algorithm is based on Crisfield(1996). Page 171. Section 6.6.3
                //     // After this step, the TrialStress(solution), TrialPlastic_Strain, and Stiffness will be updated again
                //     // to the yield surface.
                //     // In addition, each substep will have this behavior of returning to yield surface.
                //     // ==============================================================
                //     // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                //     // Make surface the internal variables are already updated. And then, return to the yield surface.
                //     double yf_val_after_corrector = yf(TrialStress);
                //     const DTensor2& n_after_corrector = yf.df_dsigma_ij(TrialStress);
                //     const DTensor2& m_after_corrector = pf(depsilon_elpl, TrialStress);
                //     // In the function below, depsilon_elpl is actually not used at all in xi_star_h_star
                //     double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m_after_corrector,  TrialStress);
                //     double dLambda_after_corrector = yf_val_after_corrector / (
                //         n_after_corrector(i,j)*Eelastic(i,j,k,l)*m_after_corrector(k,l) - xi_star_h_star_after_corrector
                //         );
                //     TrialStress(i, j) = TrialStress(i, j) - dLambda_after_corrector * Eelastic(i, j, k, l) * m_after_corrector(k, l);
                //     TrialPlastic_Strain(i, j) += dLambda_after_corrector * m_after_corrector(i, j);

                // double den_after_corrector = n_after_corrector(p, q) * Eelastic(p, q, r, s) * m_after_corrector(r, s) - xi_star_h_star_after_corrector;
                // Stiffness_substep(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m_after_corrector(p, q)) * (n_after_corrector(r, s) * Eelastic(r, s, k, l) ) / den_after_corrector;
                // Stiffness(i, j, k, l) += Stiffness_substep(i, j, k, l)/Nsubsteps;
                // }
                // // ===================================================================

                double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
                if (norm_trial_stress != norm_trial_stress )//denf <= 0 ) //check for nan
                {
                    cout << "Numeric error!\n";
                    printTensor("TrialStress = " , TrialStress);
                    printTensor("CommitStress = " , CommitStress);
                    printTensor("depsilon = " , depsilon);
                    printTensor("intersection_stress = " , intersection_stress);
                    printTensor4("Eelastic = " , Eelastic);
                    printTensor4("Stiffness = " , Stiffness);
                    cout << "yf_val_start = " << yf_val_start << endl;
                    cout << "yf_val_end = " << yf_val_end << endl;
                    printTensor("n = " , n );
                    printTensor("m = " , m );
                    cout << "xi_star_h_star  = " << xi_star_h_star << endl;
                    cout << "den = " << den << endl;
                    cout << "dLambda = " << dLambda << endl;

                    errorcode = -1;
                }
            }
        } // end plasticity

        return errorcode;
    }

    int Backward_Euler(const DTensor2 &strain_incr, bool debugrun = false, int NSteps = 1)
    {
        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;
        // // =====Sub-increments =====
        // // Problem is in the stiffness
        // double max_comp = Max_abs_Component(strain_incr)/ (double)NSteps;
        // double allowed_increment_magnitude = 1E-5;
        // if(max_comp > allowed_increment_magnitude){
        //     int required_Nstep = (int)(max_comp/ allowed_increment_magnitude)+1 ;
        //     Backward_Euler(strain_incr, false, required_Nstep);
        //     return 0;
        // }
        // // =====Sub-increments =====END
        internal_variables.revert();
        internal_variables.commit_tmp();
        static DTensor2 depsilon(3, 3, 0);
        static DTensor2 PredictorStress(3, 3, 0);
        const DTensor2& sigma = CommitStress;
        depsilon(i, j) = strain_incr(i, j) / NSteps;
        TrialStress(i, j) = sigma(i, j);
        TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        // Uncomment next line to enable the subincrement in Backward_Euler.
        // for (int step = 0; step < NSteps; step++)
        {
            DTensor4& Eelastic = et(sigma);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

            dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

            TrialStress(i, j) +=  dsigma(i, j);
            PredictorStress(i, j) = TrialStress(i, j);

            double yf_val_start = yf(sigma);
            double yf_val_end = yf(PredictorStress);

            // printTensor (" CommitStress       " , CommitStress);
            // printTensor (" depsilon           " , depsilon);
            // printTensor (" dsigma             " , dsigma);
            // printTensor (" PredictorStress    " , PredictorStress);
            // cout <<      " yf_val_start        = " << yf_val_start << endl;
            // cout <<      " yf_val_end          = " << yf_val_end << endl;
            // printTensor (" TrialStress        " , TrialStress);

            if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
            {
                DTensor4& Eelastic = et(TrialStress);
                Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
            }
            else  //Plasticity
            {
                static DTensor2 ResidualStress(3, 3, 0);
                static DTensor2 TrialStress_prev(3, 3, 0);
                ResidualStress *= 0;
                TrialStress_prev *= 0;
                double normResidualStress = -1;
                double stress_relative_error = this->stress_relative_tol * 10; //
                double dLambda = 0;
                double denominator = 0;
                int iteration_count = 0;

                double yf_PredictorStress = yf(PredictorStress);
                double yf_TrialStress = yf(TrialStress);

                // vonMises does NOT enter this part.
                // DruckerPrager requries this part.
                if (yf.hasCorner() && yf.in_Apex(PredictorStress))
                {
                    static DTensor2 small_stress(3, 3, 0.0);
                    small_stress *= 0;
                    // The small value 50*Pa refers to the lowest confinement test:
                    // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
                    double DP_k = yf.get_k();
                    double DP_p = 50 ;
                    // To make it on the yield surface, the q is equal to k*p.
                    double DP_q = DP_k * DP_p ;
                    // Assume the triaxial conditions sigma_2 = sigma_3.
                    small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
                    small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
                    small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

                    // (1) Update the trial stress
                    TrialStress(i, j) = small_stress(i, j);
                    // (2) Update the trial plastic strain
                    const DTensor2& n = yf.df_dsigma_ij(small_stress);
                    const DTensor2& m = pf(depsilon, small_stress);
                    const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
                    double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                    double dLambda = yf_PredictorStress / denominator;
                    TrialPlastic_Strain(i, j) += dLambda * m(i, j);
                    // (3) Update the internal variables
                    internal_variables.evolve(dLambda, depsilon, m, TrialStress);
                    internal_variables.commit_tmp();
                    // (4) Update the stiffness.
                    // Backward_Euler use the inconsistent stiffness.
                    // Full_Backward_Euler use the consistent stiffness.
                    Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

                    // // ===================================================================
                    // // Working on consistent stiffness
                    // // ===================================================================
                    // // Construct the Tensor T
                    static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
                    IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                    static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                    dm_dsigma = pf.dm_over_dsigma(TrialStress);
                    static DTensor4 Ts(3, 3, 3, 3, 0.0);
                    Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                    // // ===================================================================
                    // // Construct the Tensor H
                    static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                    dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
                    static DTensor2 H(3, 3, 0.0);
                    H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

                    static DTensor4 invT(3, 3, 3, 3, 0.0);
                    // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
                    if ( !inverse4thTensor(Ts, invT) )
                    {
                        return 0;
                    }
                    static DTensor4 R(3, 3, 3, 3, 0.0);
                    R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

                    // // ===================================================================
                    // // Construct the consistent stiffness
                    double denomin{0.0};
                    denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
                    Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;
                    // // ===================================================================
                    // // consistent stiffness  END
                    // // ===================================================================

                    // Skip the following iteration procedures.
                    return 0;
                }

                double f_relative_error = this->f_relative_tol * 10;

                while ((stress_relative_error > this-> stress_relative_tol ||
                        f_relative_error > this->f_relative_tol) &&
                        iteration_count < this->n_max_iterations
                      )
                {
                    iteration_count++;
                    // In this approach, we always start from
                    // PredictorStress and the initial committed internal variables.
                    internal_variables.revert_tmp();

                    const DTensor2& n = yf.df_dsigma_ij(TrialStress);
                    const DTensor2& m = pf(depsilon, TrialStress);
                    double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);
                    denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

                    //Compute the plastic multiplier
                    // if (abs(denominator) < MACHINE_EPSILON)
                    // {
                    //     cout << "CEP - denominator of plastic multiplier n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star = 0\n";
                    //     printTensor("m", m);
                    //     printTensor("n", n);
                    //     cout << "xi_star_h_star" << xi_star_h_star << endl;
                    //     cout << "denominator" << denominator << endl;
                    //     printTensor("depsilon", depsilon);
                    //     return -1;
                    // }
                    dLambda =  yf_PredictorStress / denominator;

                    TrialStress_prev(i, j) = TrialStress(i, j);

                    TrialStress(i, j) = PredictorStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
                    internal_variables.evolve(dLambda, depsilon, m, TrialStress);
                    // // ============================================================
                    // // Line search in the constitutive level
                    // // ============================================================
                    double subdLambda = dLambda;
                    while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
                    {
                        internal_variables.revert_tmp();
                        subdLambda *= 0.5;
                        TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
                        // Internal variable should be updated with subdLambda.
                        // The condition in while-loop will use the internal variables.
                        internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
                    }
                    // ============================================================
                    dLambda = subdLambda ;
                    yf_TrialStress = yf(TrialStress);

                    ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

                    normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

                    // update the stress and f relative error.
                    stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
                    f_relative_error = abs(yf_TrialStress / yf_PredictorStress);

                    double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
                    if (norm_trial_stress != norm_trial_stress || debugrun) //check for nan
                    {
                        cout << "=============================================================" << endl;
                        cout << "\nIteration # " << iteration_count <<  endl;
                        printTensor (" CommitStress       " , CommitStress);
                        double p, q, theta;
                        std::tie(p, q, theta) = getpqtheta(CommitStress);
                        fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
                        fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
                        fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
                        printTensor (" depsilon           " , depsilon);
                        printTensor (" dsigma             " , dsigma);
                        printTensor (" PredictorStress    " , PredictorStress);
                        std::tie(p, q, theta) = getpqtheta(PredictorStress);
                        fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
                        fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
                        fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
                        cout <<      " yf_val_start        = " << yf_val_start << endl;
                        cout <<      " yf_val_end          = " << yf_val_end << endl;
                        printTensor (" n                  " , n );
                        printTensor (" m                  " , m );
                        cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
                        cout <<      " denominator         = " << denominator << endl;
                        cout <<      " dLambda             = " << dLambda << endl;
                        printTensor (" TrialStress        " , TrialStress);
                        std::tie(p, q, theta) = getpqtheta(TrialStress);
                        fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
                        fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
                        fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

                        printTensor (" TrialStress_prev   " , TrialStress_prev);
                        std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
                        fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
                        fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
                        fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

                        printTensor (" ResidualStress     " , ResidualStress );
                        cout <<      " normResidualStress  = " << normResidualStress << endl;
                        cout <<      " stress_relative_error      = " << stress_relative_error << endl;
                        cout <<      " f_relative_error      = " << f_relative_error << endl;
                        cout <<      " iteration_count     = " << iteration_count << endl;
                        cout <<      " .........................................................." << endl;
                        cout <<      "  Internal variables:" << endl;
                        printTensor (" back_stress (alpha)      " , yf.get_alpha() );
                        cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
                        cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
                        cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
                        cout << "===============================================END===========" << endl;
                        cout << endl << endl;
                    }

                } // while for el-pl this step

                //Update the trial plastic strain.
                const DTensor2& n = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m = pf(depsilon, TrialStress);
                const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);

                denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                TrialPlastic_Strain(i, j) += dLambda * m(i, j);
                internal_variables.commit_tmp();

                Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

                // ===================================================================
                // Working on consistent stiffness
                // ===================================================================
                // Update the m and n:
                static DTensor2 nf(3, 3, 0.0);
                nf = yf.df_dsigma_ij(TrialStress);
                static DTensor2 mf(3, 3, 0.0);
                mf = pf(depsilon, TrialStress);
                // Construct the Tensor T
                static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
                IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                dm_dsigma = pf.dm_over_dsigma(TrialStress);
                static DTensor4 Ts(3, 3, 3, 3, 0.0);
                Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                // ===================================================================
                // Construct the Tensor H
                static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
                static DTensor2 H(3, 3, 0.0);
                H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

                static DTensor4 invT(3, 3, 3, 3, 0.0);
                // If cannot inverse T, return directly.
                // So "Stiffness" is the inconsistent stiffness tensor.
                if ( !inverse4thTensor(Ts, invT) )
                {
                    return 0;
                }
                static DTensor4 R(3, 3, 3, 3, 0.0);
                R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

                // ===================================================================
                // Construct the consistent stiffness
                double denomin{0.0};
                double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
                denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

                Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;

                // ===================================================================
                // Working on consistent stiffness  END
                // ===================================================================

            } //else  //Plasticity
            internal_variables.commit_tmp();
        } //End the loop for (int step = 0; step < NSteps; step++)
        return errorcode;
    }

    int Backward_Euler_ddlambda(const DTensor2 &strain_incr, bool debugrun = false)
    {
        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;

        internal_variables.revert();
        internal_variables.commit_tmp();
        static DTensor2 depsilon(3, 3, 0);
        static DTensor2 PredictorStress(3, 3, 0);
        const DTensor2& sigma = CommitStress;
        depsilon(i, j) = strain_incr(i, j) ;
        TrialStress(i, j) = sigma(i, j);
        TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        DTensor4& Eelastic = et(sigma);
        Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

        dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

        TrialStress(i, j) +=  dsigma(i, j);
        PredictorStress(i, j) = TrialStress(i, j);

        double yf_val_start = yf(sigma);
        double yf_val_end = yf(PredictorStress);

        // printTensor (" CommitStress       " , CommitStress);
        // printTensor (" depsilon           " , depsilon);
        // printTensor (" dsigma             " , dsigma);
        // printTensor (" PredictorStress    " , PredictorStress);
        // cout <<      " yf_val_start        = " << yf_val_start << endl;
        // cout <<      " yf_val_end          = " << yf_val_end << endl;
        // printTensor (" TrialStress        " , TrialStress);

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        }
        else  //Plasticity
        {
            static DTensor2 ResidualStress(3, 3, 0);
            static DTensor2 TrialStress_prev(3, 3, 0);
            ResidualStress *= 0;
            TrialStress_prev *= 0;
            double normResidualStress = -1;
            double stress_relative_error = this->stress_relative_tol * 10; //
            double dLambda = 0;
            double denominator = 0;
            int iteration_count = 0;

            double yf_PredictorStress = yf(PredictorStress);
            double yf_TrialStress = yf(TrialStress);
            // vonMises does NOT enter this part.
            // DruckerPrager requries this part.
            if (yf.hasCorner() && yf.in_Apex(PredictorStress))
            {
                static DTensor2 small_stress(3, 3, 0.0);
                small_stress *= 0;
                // The small value 50*Pa refers to the lowest confinement test:
                // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
                double DP_k = yf.get_k();
                double DP_p = 50 ;
                // To make it on the yield surface, the q is equal to k*p.
                double DP_q = DP_k * DP_p ;
                // Assume the triaxial conditions sigma_2 = sigma_3.
                small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
                small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
                small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

                // (1) Update the trial stress
                TrialStress(i, j) = small_stress(i, j);
                // (2) Update the trial plastic strain
                const DTensor2& n = yf.df_dsigma_ij(small_stress);
                const DTensor2& m = pf(depsilon, small_stress);
                const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
                double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                double dLambda = yf_PredictorStress / denominator;
                TrialPlastic_Strain(i, j) += dLambda * m(i, j);
                // (3) Update the internal variables
                internal_variables.evolve(dLambda, depsilon, m, TrialStress);
                internal_variables.commit_tmp();
                // (4) Update the stiffness.
                Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

                // // =========================================================
                // // Working on consistent stiffness
                // // ========================================================
                // // Construct the Tensor T
                static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
                IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                dm_dsigma = pf.dm_over_dsigma(TrialStress);
                static DTensor4 Ts(3, 3, 3, 3, 0.0);
                Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                // // Construct the Tensor H
                static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
                static DTensor2 H(3, 3, 0.0);
                H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

                static DTensor4 invT(3, 3, 3, 3, 0.0);
                // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
                if ( !inverse4thTensor(Ts, invT) )
                {
                    return 0;
                }
                static DTensor4 R(3, 3, 3, 3, 0.0);
                R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

                // // =======================================================
                // // Construct the consistent stiffness
                double denomin{0.0};
                denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
                Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;

                // // =====================================================
                // // Working on consistent stiffness  END
                // // =====================================================

                // skip the following iteration
                return 0;
            }//END yield surface APEX check

            double f_relative_error = this->f_relative_tol * 10;

            while ((stress_relative_error > this-> stress_relative_tol ||
                    f_relative_error > this->f_relative_tol) &&
                    iteration_count < this->n_max_iterations
                  )
            {
                iteration_count++;
                internal_variables.revert_tmp();

                const DTensor2& n = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m = pf(depsilon, TrialStress);
                double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);
                denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                if (abs(denominator) < MACHINE_EPSILON)
                {
                    cout << "CEP - denominator of plastic multiplier  \
                               n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star = 0\n";
                    printTensor("m", m);
                    printTensor("n", n);
                    cout << "xi_star_h_star" << xi_star_h_star << endl;
                    cout << "denominator" << denominator << endl;
                    printTensor("depsilon", depsilon);
                    return -1;
                }
                // The first dLambda is based on the elastic predictor. The next dLambda's are updated by ddlambda
                TrialStress_prev(i, j) = TrialStress(i, j);

                dLambda =  yf_PredictorStress / denominator;
                TrialStress(i, j) = TrialStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
                internal_variables.evolve(dLambda, depsilon, m, TrialStress);

                // In the first step, TrialStress is the predictor stress.
                // =======Line Search ====
                // // Zeta is the cooefficent for line search.
                double subdLambda = dLambda;
                while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
                {
                    internal_variables.revert_tmp();
                    subdLambda *= 0.5;
                    TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
                    // update the internal variables with subdLambda and commit to update yield surface.
                    internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
                }
                dLambda = subdLambda;
                ResidualStress(i, j) =  TrialStress(i, j) - TrialStress_prev(i, j);
                // ==========================================
                // Works to calculate the ddlambda
                // ==========================================
                // Update the m and n:
                static DTensor2 n_new(3, 3, 0.0);
                n_new = yf.df_dsigma_ij(TrialStress);
                static DTensor2 m_new(3, 3, 0.0);
                m_new = pf(depsilon, TrialStress);
                // Construct the Tensor T
                static DTensor4 IdentityTensor4(3, 3, 3, 3, 0.0);
                IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                dm_dsigma = pf.dm_over_dsigma(TrialStress);
                static DTensor4 T_new(3, 3, 3, 3, 0.0);
                T_new(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                static DTensor4 invTnew(3, 3, 3, 3, 0.0);
                if ( !inverse4thTensor(T_new, invTnew) )
                {
                    // cout<<"Singularity. Cannot inverse tensor T! "<<endl;
                    return 0;
                }

                static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m_new, TrialStress);
                static DTensor2 H_new(3, 3, 0.0);
                H_new(k, l) = m_new(k, l) + dLambda * dm_dq_star_h_star(k, l);

                double const& xi_star_h_star_new = yf.xi_star_h_star( depsilon, m_new,  TrialStress);

                double denomina_ddlambda = 0.0;
                static DTensor2 E_times_H (3, 3, 0.0);
                Eelastic = et(TrialStress);
                E_times_H(i, j) = Eelastic(i, j, k, l) * H_new(k, l);
                denomina_ddlambda = n_new(r, s) * (E_times_H(p, q) * invTnew(p, q, r, s)) - xi_star_h_star_new ;
                double numerator_ddlambda = 0.0;
                numerator_ddlambda = yf(TrialStress) - n_new(r, s) * (ResidualStress(i, j) * invTnew(i, j, r, s)) ;
                double ddlambda =   numerator_ddlambda / denomina_ddlambda ;
                static DTensor2 dSigma_(3, 3, 0.0);
                dSigma_(r, s) = -
                                (
                                    ResidualStress(i, j) + ddlambda * Eelastic(i, j, p, q) * H_new(p, q)
                                )
                                * invTnew(i, j, r, s) ;
                TrialStress_prev(i, j) = TrialStress(i, j) ;
                TrialStress(i, j) =  TrialStress(i, j) + dSigma_(i, j) ;
                double yf_TrialStress_prev = yf(TrialStress_prev) ;
                internal_variables.evolve(ddlambda, depsilon, m_new, TrialStress);

                // =======Line Search ddlambda (double delta ) =======
                double subddlambda = ddlambda;
                double zeta = 1.0;
                while (yf(TrialStress) > yf_TrialStress_prev && subddlambda > MACHINE_EPSILON * 10)
                {
                    internal_variables.revert_tmp();
                    subddlambda *= 0.5 ;
                    zeta *= 0.5 ;
                    TrialStress(i, j) = TrialStress(i, j) + zeta * dSigma_(i, j) ;
                    // Internal variable should be updated with subddlambda.
                    // The condition in while-loop will use the internal variables.
                    internal_variables.evolve(subddlambda, depsilon, m, TrialStress);
                }

                dLambda = dLambda + subddlambda ;

                yf_TrialStress = yf(TrialStress);

                ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

                normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

                stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
                f_relative_error = abs(yf_TrialStress / yf_PredictorStress);
                debugrun = true;
                if (normResidualStress != normResidualStress || debugrun)  // check for NAN
                {
                    cout << "=============================================================" << endl;
                    cout << "\nIteration # " << iteration_count <<  endl;
                    printTensor (" CommitStress       " , CommitStress);
                    double p, q, theta;
                    std::tie(p, q, theta) = getpqtheta(CommitStress);
                    fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
                    printTensor (" depsilon           " , depsilon);
                    printTensor (" dsigma             " , dsigma);
                    printTensor (" PredictorStress    " , PredictorStress);
                    std::tie(p, q, theta) = getpqtheta(PredictorStress);
                    fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
                    cout <<      " yf_val_start        = " << yf_val_start << endl;
                    cout <<      " yf_val_end          = " << yf_val_end << endl;
                    printTensor (" n                  " , n );
                    printTensor (" m                  " , m );
                    cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
                    cout <<      " denominator         = " << denominator << endl;
                    cout <<      " dLambda             = " << dLambda << endl;
                    printTensor (" TrialStress        " , TrialStress);
                    std::tie(p, q, theta) = getpqtheta(TrialStress);
                    fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

                    printTensor (" TrialStress_prev   " , TrialStress_prev);
                    std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
                    fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
                    fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
                    fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

                    printTensor (" ResidualStress     " , ResidualStress );
                    cout <<      " normResidualStress  = " << normResidualStress << endl;
                    cout <<      " stress_relative_error      = " << stress_relative_error << endl;
                    cout <<      " f_relative_error      = " << f_relative_error << endl;
                    cout <<      " iteration_count     = " << iteration_count << endl;
                    cout <<      " .........................................................." << endl;
                    cout <<      "  Internal variables:" << endl;
                    printTensor (" back_stress (alpha)      " , yf.get_alpha() );
                    cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
                    cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
                    cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
                    cout << "===============================================END===========" << endl;
                    cout << endl << endl;
                }

            } // while for el-pl this step


            // ===============================================
            // Start working on the stiffness
            // ===============================================
            // Inconsistent stiffness
            //Update the trial plastic strain. and internal variables
            const DTensor2& n = yf.df_dsigma_ij(TrialStress);
            const DTensor2& m = pf(depsilon, TrialStress);
            const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);

            denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
            TrialPlastic_Strain(i, j) += dLambda * m(i, j);

            Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

            // ===================================================================
            // Working on consistent stiffness
            // ===================================================================
            // Update the m and n:
            static DTensor2 nf(3, 3, 0.0);
            nf = yf.df_dsigma_ij(TrialStress);
            static DTensor2 mf(3, 3, 0.0);
            mf = pf(depsilon, TrialStress);
            // Construct the Tensor T
            static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
            IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

            static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
            dm_dsigma = pf.dm_over_dsigma(TrialStress);
            static DTensor4 Ts(3, 3, 3, 3, 0.0);
            Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

            // Construct the Tensor H
            static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
            dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
            static DTensor2 H(3, 3, 0.0);
            H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

            static DTensor4 invT(3, 3, 3, 3, 0.0);
            // If cannot inverse T, return directly.
            // So "Stiffness" is the inconsistent stiffness tensor.
            if ( !inverse4thTensor(Ts, invT) )
            {
                return 0;
            }
            static DTensor4 R(3, 3, 3, 3, 0.0);
            R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

            // ===================================================================
            // Construct the consistent stiffness
            double denomin{0.0};
            double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
            denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

            Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;
            // ===================================================================
            // Working on consistent stiffness  END
            // ===================================================================

        } //else  //Plasticity
        internal_variables.commit();
        return errorcode;
    }


    int Backward_Euler_ddlambda_Subincrement(const DTensor2 &strain_incr, bool debugrun = false)
    {
        using namespace ClassicElastoplasticityGlobals;
        int errorcode = 0;

        internal_variables.revert();
        internal_variables.commit_tmp();
        static DTensor2 depsilon(3, 3, 0);
        static DTensor2 PredictorStress(3, 3, 0);
        const DTensor2& sigma = CommitStress;
        depsilon(i, j) = strain_incr(i, j) ;
        TrialStress(i, j) = sigma(i, j);
        TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        DTensor4& Eelastic = et(sigma);
        Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

        dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

        TrialStress(i, j) +=  dsigma(i, j);
        PredictorStress(i, j) = TrialStress(i, j);

        double yf_val_start = yf(sigma);
        double yf_val_end = yf(PredictorStress);

        // printTensor (" CommitStress       " , CommitStress);
        // printTensor (" depsilon           " , depsilon);
        // printTensor (" dsigma             " , dsigma);
        // printTensor (" PredictorStress    " , PredictorStress);
        // cout <<      " yf_val_start        = " << yf_val_start << endl;
        // cout <<      " yf_val_end          = " << yf_val_end << endl;
        // printTensor (" TrialStress        " , TrialStress);

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        }
        else  //Plasticity
        {
            static DTensor2 ResidualStress(3, 3, 0);
            static DTensor2 TrialStress_prev(3, 3, 0);
            ResidualStress *= 0;
            TrialStress_prev *= 0;
            double normResidualStress = -1;
            double stress_relative_error = this->stress_relative_tol * 10; //
            double dLambda = 0;
            double denominator = 0;
            int iteration_count = 0;

            double yf_PredictorStress = yf(PredictorStress);
            double yf_TrialStress = yf(TrialStress);
            // vonMises does NOT enter this part.
            // DruckerPrager requries this part.
            if (yf.hasCorner() && yf.in_Apex(PredictorStress))
            {
                static DTensor2 small_stress(3, 3, 0.0);
                small_stress *= 0;
                // The small value 50*Pa refers to the lowest confinement test:
                // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
                double DP_k = yf.get_k();
                double DP_p = 50 ;
                // To make it on the yield surface, the q is equal to k*p.
                double DP_q = DP_k * DP_p ;
                // Assume the triaxial conditions sigma_2 = sigma_3.
                small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
                small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
                small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

                // (1) Update the trial stress
                TrialStress(i, j) = small_stress(i, j);
                // (2) Update the trial plastic strain
                const DTensor2& n = yf.df_dsigma_ij(small_stress);
                const DTensor2& m = pf(depsilon, small_stress);
                const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
                double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
                double dLambda = yf_PredictorStress / denominator;
                TrialPlastic_Strain(i, j) += dLambda * m(i, j);
                // (3) Update the internal variables
                internal_variables.evolve(dLambda, depsilon, m, TrialStress);
                internal_variables.commit_tmp();
                // (4) Update the stiffness.
                Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

                // // =========================================================
                // // Working on consistent stiffness
                // // ========================================================
                // // Construct the Tensor T
                static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
                IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                dm_dsigma = pf.dm_over_dsigma(TrialStress);
                static DTensor4 Ts(3, 3, 3, 3, 0.0);
                Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                // // Construct the Tensor H
                static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
                static DTensor2 H(3, 3, 0.0);
                H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

                static DTensor4 invT(3, 3, 3, 3, 0.0);
                // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
                if ( !inverse4thTensor(Ts, invT) )
                {
                    return 0;
                }
                static DTensor4 R(3, 3, 3, 3, 0.0);
                R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

                // // =======================================================
                // // Construct the consistent stiffness
                double denomin{0.0};
                denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
                Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;

                // // =====================================================
                // // Working on consistent stiffness  END
                // // =====================================================

                // skip the following iterations
                return 0;
            }//END yield surface APEX check

            // ==================================================================
            // From the start point(elastic predictor) to the first plastic corrector.
            // ==================================================================
            yf_PredictorStress = yf(PredictorStress);
            const DTensor2& n = yf.df_dsigma_ij(PredictorStress);
            const DTensor2& m = pf(depsilon, PredictorStress);
            double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  PredictorStress);
            denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

            dLambda =  yf_PredictorStress / denominator;

            TrialStress_prev(i, j) = TrialStress(i, j);

            TrialStress(i, j) = PredictorStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
            internal_variables.evolve(dLambda, depsilon, m, TrialStress);
            // // ============================================================
            // // Line search in dlambda
            // // ============================================================
            double subdLambda = dLambda;
            while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
            {
                internal_variables.revert_tmp();
                subdLambda *= 0.5;
                TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
                // Internal variable should be updated with subdLambda.
                // The condition in while-loop will use the internal variables.
                internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
            }
            // ============================================================
            dLambda = subdLambda;
            ResidualStress(i, j) =  TrialStress(i, j) - TrialStress_prev(i, j);

            double f_relative_error = this->f_relative_tol * 10;

            // Iteration over the ddlambda
            while ((stress_relative_error > this-> stress_relative_tol ||
                    f_relative_error > this->f_relative_tol) &&
                    iteration_count < this->n_max_iterations
                  )
            {
                iteration_count++;

                // ==========================================
                // Works to calculate the ddlambda
                // ==========================================
                // Update the m and n:
                static DTensor2 n_new(3, 3, 0.0);
                n_new = yf.df_dsigma_ij(TrialStress);
                static DTensor2 m_new(3, 3, 0.0);
                m_new = pf(depsilon, TrialStress);
                // Construct the Tensor T
                static DTensor4 IdentityTensor4(3, 3, 3, 3, 0.0);
                IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

                static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
                dm_dsigma = pf.dm_over_dsigma(TrialStress);
                static DTensor4 T_new(3, 3, 3, 3, 0.0);
                T_new(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

                static DTensor4 invTnew(3, 3, 3, 3, 0.0);
                if ( !inverse4thTensor(T_new, invTnew) )
                {
                    // cout<<"Singularity. Cannot inverse tensor T! "<<endl;
                    return 0;
                }

                static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
                dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m_new, TrialStress);
                static DTensor2 H_new(3, 3, 0.0);
                H_new(k, l) = m_new(k, l) + dLambda * dm_dq_star_h_star(k, l);

                double const& xi_star_h_star_new = yf.xi_star_h_star( depsilon, m_new,  TrialStress);

                double denomina_ddlambda = 0.0;
                static DTensor2 E_times_H (3, 3, 0.0);
                Eelastic = et(TrialStress);
                E_times_H(i, j) = Eelastic(i, j, k, l) * H_new(k, l);
                denomina_ddlambda = n_new(r, s) * (E_times_H(p, q) * invTnew(p, q, r, s)) - xi_star_h_star_new ;
                double numerator_ddlambda = 0.0;
                numerator_ddlambda = yf(TrialStress) - n_new(r, s) * (ResidualStress(i, j) * invTnew(i, j, r, s)) ;
                double ddlambda =   numerator_ddlambda / denomina_ddlambda ;

                static DTensor2 dSigma_(3, 3, 0.0);
                dSigma_(r, s) = -
                                (
                                    ResidualStress(i, j) + ddlambda * Eelastic(i, j, p, q) * H_new(p, q)
                                )
                                * invTnew(i, j, r, s) ;

                TrialStress_prev(i, j) = TrialStress(i, j) ;
                TrialStress(i, j) =  TrialStress(i, j) + dSigma_(i, j) ;
                double yf_TrialStress_prev = yf(TrialStress_prev) ;
                internal_variables.evolve(ddlambda, depsilon, m_new, TrialStress);

                // =======Line Search ddlambda (double delta ) =======
                double subddlambda = ddlambda;
                double zeta = 1.0;
                while (yf(TrialStress) > yf_TrialStress_prev && subddlambda > MACHINE_EPSILON * 10)
                {
                    internal_variables.revert_tmp();
                    subddlambda *= 0.5 ;
                    zeta *= 0.5 ;
                    TrialStress(i, j) = TrialStress(i, j) + zeta * dSigma_(i, j) ;
                    // Internal variable should be updated with subddlambda.
                    // The condition in while-loop will use the internal variables.
                    internal_variables.evolve(subddlambda, depsilon, m, TrialStress);
                }

                dLambda = dLambda + subddlambda ;
                yf_TrialStress = yf(TrialStress);

                ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

                normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

                stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
                f_relative_error = abs(yf_TrialStress / yf_PredictorStress);
                // debugrun = true;
                if (normResidualStress != normResidualStress || debugrun)  // check for NAN
                {
                    cout << "=============================================================" << endl;
                    cout << "\nIteration # " << iteration_count <<  endl;
                    printTensor (" CommitStress       " , CommitStress);
                    double p, q, theta;
                    std::tie(p, q, theta) = getpqtheta(CommitStress);
                    fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
                    printTensor (" depsilon           " , depsilon);
                    printTensor (" dsigma             " , dsigma);
                    printTensor (" PredictorStress    " , PredictorStress);
                    std::tie(p, q, theta) = getpqtheta(PredictorStress);
                    fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
                    cout <<      " yf_val_start        = " << yf_val_start << endl;
                    cout <<      " yf_val_end          = " << yf_val_end << endl;
                    printTensor (" n                  " , n );
                    printTensor (" m                  " , m );
                    cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
                    cout <<      " denominator         = " << denominator << endl;
                    cout <<      " dLambda             = " << dLambda << endl;
                    printTensor (" TrialStress        " , TrialStress);
                    std::tie(p, q, theta) = getpqtheta(TrialStress);
                    fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
                    fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
                    fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

                    printTensor (" TrialStress_prev   " , TrialStress_prev);
                    std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
                    fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
                    fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
                    fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

                    printTensor (" ResidualStress     " , ResidualStress );
                    cout <<      " normResidualStress  = " << normResidualStress << endl;
                    cout <<      " stress_relative_error      = " << stress_relative_error << endl;
                    cout <<      " f_relative_error      = " << f_relative_error << endl;
                    cout <<      " iteration_count     = " << iteration_count << endl;
                    cout <<      " .........................................................." << endl;
                    cout <<      "  Internal variables:" << endl;
                    printTensor (" back_stress (alpha)      " , yf.get_alpha() );
                    cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
                    cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
                    cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
                    cout << "===============================================END===========" << endl;
                    cout << endl << endl;
                }

            } // while for el-pl this step


            // ===============================================
            // Start working on the stiffness
            // ===============================================
            // Inconsistent stiffness
            //Update the trial plastic strain. and internal variables
            const DTensor2& n_ = yf.df_dsigma_ij(TrialStress);
            const DTensor2& m_ = pf(depsilon, TrialStress);
            const double xi_star_h_star_ = yf.xi_star_h_star( depsilon, m_,  TrialStress);

            denominator = n_(p, q) * Eelastic(p, q, r, s) * m_(r, s) - xi_star_h_star_;

            Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m_(p, q)) * (n_(r, s) * Eelastic(r, s, k, l) ) / denominator;

            // ===================================================================
            // Working on consistent stiffness
            // ===================================================================
            // Update the m and n:
            static DTensor2 nf(3, 3, 0.0);
            nf = yf.df_dsigma_ij(TrialStress);
            static DTensor2 mf(3, 3, 0.0);
            mf = pf(depsilon, TrialStress);
            // Construct the Tensor T
            static DTensor4 IdentityTensor4(3, 3, 3, 3, 0);
            IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

            static DTensor4 dm_dsigma(3, 3, 3, 3, 0.0);
            dm_dsigma = pf.dm_over_dsigma(TrialStress);
            static DTensor4 Ts(3, 3, 3, 3, 0.0);
            Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

            // Construct the Tensor H
            static DTensor2 dm_dq_star_h_star(3, 3, 0.0);
            dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
            static DTensor2 H(3, 3, 0.0);
            H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

            static DTensor4 invT(3, 3, 3, 3, 0.0);
            // If cannot inverse T, return directly.
            // So "Stiffness" is the inconsistent stiffness tensor.
            if ( !inverse4thTensor(Ts, invT) )
            {
                return 0;
            }
            static DTensor4 R(3, 3, 3, 3, 0.0);
            R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

            // ===================================================================
            // Construct the consistent stiffness
            double denomin{0.0};
            double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
            denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

            Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;
            // ===================================================================
            // Working on consistent stiffness  END
            // ===================================================================

        } //else  //Plasticity
        internal_variables.commit();
        return errorcode;
    }


    // The algorithm below Modified_Euler_Error_Control was implemented by Jose. 18Sep2016
    int Modified_Euler_Error_Control(const DTensor2 &strain_incr)
    {
        using namespace ClassicElastoplasticityGlobals;  // Brings indexes i,j,k,l,m,n,p,q into current scope

        int errorcode = 0;

        static DTensor2 depsilon(3, 3, 0);
        static DTensor2 start_stress(3, 3, 0.0);
        static DTensor2 end_stress(3, 3, 0.0);
        // const DTensor2& sigma = CommitStress;

        //Zero out static variables to ensure nothing is there.
        depsilon *= 0;
        dsigma *= 0;
        start_stress *= 0;
        end_stress *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        depsilon(i, j) = strain_incr(i, j);         //This is needed because some function don't take const references

        //This makes sure we start with last committed values for internal_variables
        internal_variables.revert();
        TrialStress(i, j) = CommitStress(i, j);
        TrialStrain(i, j) = CommitStrain(i, j);
        TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

        //Compute elasticity for current (committed) stress level and assign stiffness to be this.
        DTensor4& Eelastic = et(TrialStress);
        Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

        //Compute the elastic stress increment
        dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

        //Increment the trial stresses and strains
        TrialStress(i, j) += dsigma(i, j);
        TrialStrain(i, j) += depsilon(i, j);

        //Compute the values of yield function before and after elastic increment
        double yf_val_start = yf(CommitStress);
        double yf_val_end = yf(TrialStress);

        //These are the start and end stresses for the brent algorithm which
        //finds the yield surface intersection stress (taken from old NewTemplate3dEP classes)
        start_stress(i, j) = CommitStress(i, j);
        end_stress(i, j) = TrialStress(i, j);

        //Set the intersection stress at the start
        intersection_stress(i, j) = start_stress(i, j);

        //This pre_integration_callback_(...) provides an opportunity to exit before integration
        //in case its needed. For example when some materials go in tension.
        bool returns = false;
        int retval = pre_integration_callback_(depsilon, dsigma, TrialStress, Stiffness, yf_val_start, yf_val_end,  returns);
        if (returns)
        {
            return retval;
        }


        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            //For an elastic step, compute the elasticity at the incremented stress,
            //for use in next step.
            DTensor4& Eelastic = et(TrialStress);
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        }
        else  //Plasticity
        {
            //Some temporary variables used herein
            static DTensor2 sub_depsilon_elpl(3, 3, 0);
            static DTensor2 dsigma1(3, 3, 0);
            static DTensor2 dsigma2(3, 3, 0);
            static DTensor2 sigma_error(3, 3, 0);
            static DTensor2 m(3, 3, 0);

            sub_depsilon_elpl *= 0;
            dsigma1 *= 0;
            dsigma2 *= 0;
            sigma_error *= 0;
            m *= 0;

            //If the material was elastic but now is plastic, compute the intersection point
            //ie. the point at which the trial stress path intersects the yield surface.
            //Also compute the intersection strain increment and compute the plastic strain increment
            //from that.
            depsilon_elpl(i, j) = depsilon(i, j);
            if (yf_val_start < 0)
            {
                double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );
                intersection_stress(i, j) = start_stress(i, j) * (1 - intersection_factor) + end_stress(i, j) * intersection_factor;
                depsilon_elpl(i, j) = (1 - intersection_factor) * depsilon(i, j);
            }

            //Error controlled iteration after Sloan
            int count = 0;      //Number of iterations done
            double Tau = 0;     //Integration pseudo-time
            double dT = 1;      //Current integration pseudo-time step
            double q = 0;       //A factor tu multiply dT with for next step

            //Set both the start stress and the trial stress at the intersection stress
            //if trial fails, will return to start_stress.
            start_stress(i, j) = intersection_stress(i, j);
            TrialStress(i, j) = intersection_stress(i, j);

            internal_variables.commit_tmp(); // Stores temporary internal variables in static storage space

            while (Tau < 1)
            {
                //Compute current plastic strain sub-increment.
                sub_depsilon_elpl(i, j) = dT * depsilon_elpl(i, j);

                //Compute elasticity at current trial stress and determine the first elastic stress
                //increment based on this.
                Eelastic = et(TrialStress);
                dsigma1(i, j)  = Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);

                printTensor("  dsigma_elastic", dsigma1);

                //Compute normal to YF (n) and Plastic Flow direction (m)
                const DTensor2& n1 = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m1 = pf(sub_depsilon_elpl, TrialStress);
                printTensor("  n1", n1);
                printTensor("  m1", m1);

                //Compute hardening parameter
                double xi_star_h_star1 = yf.xi_star_h_star( sub_depsilon_elpl, m1,  TrialStress);

                //Compute denominator of the plastic multiplier
                double den1 = n1(p, q) * Eelastic(p, q, r, s) * m1(r, s) - xi_star_h_star1;

                //Compute the plastic multiplier
                if (den1 == 0)
                {
                    cerr << "CEP - den = 0\n";
                    return -1;
                }
                double dLambda1 =  n1(i, j) * Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
                dLambda1 /= den1;

                //Correct the trial stress increment
                dsigma1(i, j) +=  - dLambda1 * Eelastic(i, j, k, l) * m1(k, l);

                //Now evolve material state
                internal_variables.evolve(dLambda1, sub_depsilon_elpl, m1, TrialStress);
                TrialStress(i, j) += dsigma1(i, j);

                //Now evaluate everything at the new trial stress
                Eelastic = et(TrialStress);
                dsigma2(i, j)  = Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
                const DTensor2& n2 = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m2 = pf(sub_depsilon_elpl, TrialStress);
                double xi_star_h_star2 = yf.xi_star_h_star( sub_depsilon_elpl, m2,  TrialStress);
                double den2 = n2(p, q) * Eelastic(p, q, r, s) * m2(r, s) - xi_star_h_star2;
                if (den2 == 0)
                {
                    cerr << "CEP - den = 0\n";
                    return -1;
                }
                double dLambda2 =  n2(i, j) * Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
                dLambda2 /= den2;
                dsigma2(i, j) += - dLambda2 * Eelastic(i, j, k, l) * m2(k, l);


                //Compute a better estimate of TrialStress by averaging the two increments
                TrialStress(i, j) = start_stress(i, j);
                sigma_error(i, j) = (dsigma2(i, j) - dsigma1(i, j)) / 2;
                TrialStress(i, j) += (dsigma1(i, j) + dsigma2(i, j)) / 2;

                //Compute the relative error based on the stress error and the trial stress
                double Relative_Error = sqrt(sigma_error(i, j) * sigma_error(i, j)) / sqrt(TrialStress(i, j) * TrialStress(i, j));

                double dLambda = 0;

                cout << "  den1 = " << den1 << ", dLambda1 = " << dLambda1 << endl;
                cout << "  den2 = " << den2 << ", dLambda2 = " << dLambda2 << endl;
                printTensor("  depsilon_elpl", depsilon_elpl);
                printTensor("  sub_depsilon_elpl", sub_depsilon_elpl);
                printTensor("  dsigma1", dsigma1);
                printTensor("  dsigma2", dsigma2);
                cout << "  dT = " << dT << " (T = " << Tau << " q = " << q << ")\n";
                cout << "  RE = " << Relative_Error << " (tol=" << this-> stress_relative_tol << ") ->";

                //Revert the variables to the start (use the temporary storage)
                internal_variables.revert_tmp();
                TrialStress(i, j) = start_stress(i, j);
                if (Relative_Error < this-> stress_relative_tol) // Accept this step
                {
                    cout << "Accept!\n";
                    Tau += dT;
                    q = fmin(0.8 * sqrt(this-> stress_relative_tol / Relative_Error), 2.0);
                    dT = q * dT;
                    dT = fmin(dT, 1.0 - Tau);

                    dLambda = dLambda1;//(dLambda1 + dLambda2) / 2;
                    internal_variables.evolve(dLambda, sub_depsilon_elpl, m1, TrialStress);
                    internal_variables.commit_tmp();
                    TrialStress(i, j) += dsigma1(i, j);//(dsigma1(i, j) + dsigma2(i, j)) / 2;


                    start_stress(i, j) = TrialStress(i, j);

                    //Update the trial plastic strain.
                    TrialPlastic_Strain(i, j) += dLambda * m1(i, j);
                }
                else //Reject step, take a smaller one!
                {
                    cout << "Reject!\n";
                    // cout << "Reject! R = " <<  Relative_Error << endl;
                    q = fmax(0.8 * sqrt(this-> stress_relative_tol / Relative_Error), 0.1);
                    dT = q * dT;
                }

                count += 1;
                if (count > this-> n_max_iterations)
                {
                    cerr << "ClassicElastoplasticMaterial -- Modified euler failed to converge after " << count << " iterations with : \n"
                         << "  Relative_Error = " << Relative_Error  << endl
                         << "  Relative Tolerance = " << this-> stress_relative_tol  << endl
                         << "  q = " << q  << endl
                         << "  dT = " << dT  << endl
                         << "  dLambda = " << dLambda  << endl
                         << "  T = " << Tau  << endl;
                    printTensor("TrialStress = " , TrialStress);
                    printTensor("CommitStress = " , CommitStress);
                    printTensor("depsilon = " , depsilon);
                    printTensor("depsilon_elpl = " , depsilon_elpl);
                    printTensor("sub_depsilon_elpl = " , sub_depsilon_elpl);
                    printTensor("dsigma1 = " , dsigma1);
                    printTensor("dsigma2 = " , dsigma2);
                    printTensor("sigma_error = " , sigma_error);
                    printTensor("start_stress = " , start_stress);
                    printTensor("intersection_stress = " , intersection_stress);
                    return -1;
                }
            }

            //Stress correction back to the yield surface
            double yc = yf(TrialStress);
            int ys_correction_count = 0;
            printTensor("  TrialStress before correction", TrialStress);
            while (yc > this-> f_relative_tol && ys_correction_count < this-> n_max_iterations)
            {
                cout << "    +++ yf() = " << yc << endl;
                const DTensor2& n = yf.df_dsigma_ij(TrialStress);
                const DTensor2& m = pf(depsilon_elpl, TrialStress);
                double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s);
                double dLambda_correction = yc / den;
                TrialStress(i, j) = TrialStress(i, j) - dLambda_correction * Eelastic(i, j, k, l) * m(k, l);
                printTensor("    +++ TrialStress correction", TrialStress);
                yc = yf(TrialStress);
                ys_correction_count++;
            }
            printTensor("  TrialStress after correction", TrialStress);


            // Once done compute the stiffness for the next step
            Eelastic = et(TrialStress);
            const DTensor2& nf = yf.df_dsigma_ij(TrialStress);
            const DTensor2& mf = pf(depsilon_elpl, TrialStress);
            double xi_star_h_star_f = yf.xi_star_h_star( depsilon_elpl, mf,  TrialStress);
            double denf = nf(p, q) * Eelastic(p, q, r, s) * mf(r, s) - xi_star_h_star_f;
            Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * mf(p, q)) * (nf(r, s) * Eelastic(r, s, k, l) ) / denf;



            //Detect a NAN result and print diagnostics if one is found
            // double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
            if (true)//norm_trial_stress != norm_trial_stress) //check for nan
            {
                // cout << "Nan Detected in Modified_Euler !\n";
                printTensor("TrialStress = " , TrialStress);
                printTensor("CommitStress = " , CommitStress);
                printTensor("depsilon = " , depsilon);
                printTensor("intersection_stress = " , intersection_stress);
                printTensor4("Eelastic = " , Eelastic);
                printTensor4("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor("n = " , nf );
                printTensor("m = " , mf );
                cout << "xi_star_h_star_f  = " << xi_star_h_star_f << endl;
                cout << "denf = " << denf << endl;
                // cout << "dLambda = " << dLambda << endl;
                // cout << "count = " << count << endl;
                // cout << "ys_correction_count = " << ys_correction_count << endl;

                // errorcode = -1;
            }

        }

        return errorcode;
    }

    template <typename U = T>
    typename std::enable_if < !supports_pre_integration_callback<U>::value, int >::type
// typename std::enable_if < !std::is_base_of<defines_pre_integration_callback, U>::value, int >::type
    pre_integration_callback_(const DTensor2 &depsilon, const DTensor2 &dsigma,  const DTensor2 &TrialStress, const DTensor4 &Stiffness, double yf1, double yf2, bool & returns)
    {
        // cout << "pre_integration_callback_ disabled\n";
        returns = false;
        return 0;
    }

    template <typename U = T>
    typename std::enable_if<supports_pre_integration_callback<U>::value, int>::type
// typename std::enable_if<std::is_base_of<defines_pre_integration_callback, U>::value, int>::type
    pre_integration_callback_(const DTensor2 &depsilon, const DTensor2 &dsigma, const DTensor2 &TrialStress, const DTensor4 &Stiffness, double yf1, double yf2, bool & returns)
    {
        // cout << "pre_integration_callback_ enabled\n";
        return static_cast<U*>(this)->pre_integration_callback(depsilon, dsigma, TrialStress, Stiffness,  yf1,  yf2, returns);
    }

private:
// Routine used by yield_surface_cross to find the stresstensor at cross point
//================================================================================
    double zbrentstress(const DTensor2& start_stress,
                        const DTensor2& end_stress,
                        double x1, double x2, double tol) const
    {
        using namespace ClassicElastoplasticityGlobals;
        double EPS = numeric_limits<double>::epsilon();

        int iter;
        double a = x1;
        double b = x2;
        double c = 0.0;
        double d = 0.0;
        double e = 0.0;
        double min1 = 0.0;
        double min2 = 0.0;
        double fc = 0.0;
        double p = 0.0;
        double q = 0.0;
        double r = 0.0;
        double s = 0.0;
        double tol1 = 0.0;
        double xm = 0.0;

        // double fa = func(start_stress, end_stress, *ptr_material_parameter, a);
        // double fb = func(start_stress, end_stress, *ptr_material_parameter, b);

        static DTensor2 sigma_a(3, 3, 0.0);
        static DTensor2 sigma_b(3, 3, 0.0);

        sigma_a(i, j) = start_stress(i, j) * (1 - a)  + end_stress(i, j) * a;
        sigma_b(i, j) = start_stress(i, j) * (1 - b)  + end_stress(i, j) * b;

        double fa = yf(sigma_a);
        double fb = yf(sigma_b);

        // cout << "   brent fa = " << fa << " fb = " << fb << endl;


        if ( (fb * fa) > 0.0)
        {
            std::cout << "\a\n Root must be bracketed in ZBRENTstress " << std::endl;
            exit(1);
        }

        fc = fb;

        for ( iter = 1; iter <= ClassicElastoplasticMaterial_MAXITER_BRENT; iter++ )
        {
            if ( (fb * fc) > 0.0)
            {
                c = a;
                fc = fa;
                e = d = b - a;
            }

            if ( fabs(fc) < fabs(fb) )
            {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);

            if ( fabs(xm) <= tol1 || fb == 0.0 )
            {
                return b;
            }

            if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) )
            {
                s = fb / fa;

                if (a == c)
                {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else
                {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * ( 2.0 * xm * q * (q - r) - (b - a) * (r - 1.0) );
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0)
                {
                    q = -q;
                }

                p = fabs(p);
                min1 = 3.0 * xm * q - fabs(tol1 * q);
                min2 = fabs(e * q);

                if (2.0 * p < (min1 < min2 ? min1 : min2))
                {
                    e = d;
                    d = p / q;
                }
                else
                {
                    d = xm;
                    e = d;
                }
            }
            else
            {
                d = xm;
                e = d;
            }

            a = b;
            fa = fb;

            if (fabs(d) > tol1)
            {
                b += d;
            }
            else
            {
                b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
            }

            // fb = func(start_stress, end_stress, *ptr_material_parameter, b);
            sigma_b(i, j) = start_stress(i, j) * (1 - b)  + end_stress(i, j) * b;
            fb = yf(sigma_b);
        }

        return 0.0;
    }

    double Max_abs_Component(DTensor2 const& matrix3by3)
    {
        double ret = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (ret < abs(matrix3by3(i, j)))
                {
                    ret = abs(matrix3by3(i, j));
                }
        return ret;
    }

private:

    double rho;

    DTensor2 TrialStrain;
    DTensor2 TrialStress;
    DTensor2 TrialPlastic_Strain;

    DTensor2 CommitStress;
    DTensor2 CommitStrain;
    DTensor2 CommitPlastic_Strain;


    YieldFunctionType yf;
    ElasticityType    et;
    PlasticFlowType   pf;
    MaterialInternalVariablesType internal_variables;

    bool first_step;

    static DTensor2 dsigma;
    static DTensor2 depsilon_elpl;    //Elastoplastic strain increment : For a strain increment that causes first yield, the step is divided into an elastic one (until yield) and an elastoplastic one.
    static DTensor2 intersection_stress;
    static DTensor2 intersection_strain;
    static DTensor4 Stiffness;
// static DTensor2 m;


};

// int ClassicElastoplasticMaterial< class ElasticityType, class YieldFunctionType, , class PlasticFlowType, class HardeningLawType, int thisClassTag >::constitutive_integration_method = 0;

template < class E, class Y, class P, class M, int tag, class T >
DTensor2 ClassicElastoplasticMaterial< E,  Y,  P,  M,  tag,  T >::dsigma(3, 3, 0.0);
template < class E, class Y, class P, class M, int tag, class T >
DTensor2 ClassicElastoplasticMaterial< E,  Y,  P,  M,  tag,  T >::depsilon_elpl(3, 3, 0.0);  //Used to compute the yield surface intersection.
template < class E, class Y, class P, class M, int tag, class T >
DTensor2 ClassicElastoplasticMaterial< E,  Y,  P,  M,  tag,  T >::intersection_stress(3, 3, 0.0);  //Used to compute the yield surface intersection.
template < class E, class Y, class P, class M, int tag, class T >
DTensor2 ClassicElastoplasticMaterial< E,  Y,  P,  M,  tag,  T >::intersection_strain(3, 3, 0.0);  //Used to compute the yield surface intersection.
template < class E, class Y, class P, class M, int tag, class T >
DTensor4 ClassicElastoplasticMaterial< E,  Y,  P,  M,  tag,  T >::Stiffness(3, 3, 3, 3, 0.0);  //Used to compute the yield surface intersection.


#endif
