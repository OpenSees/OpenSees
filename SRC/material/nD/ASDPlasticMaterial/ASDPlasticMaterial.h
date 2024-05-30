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



#ifndef ASDPlasticMaterial_H
#define ASDPlasticMaterial_H


#include "NDMaterial.h"
#include <G3Globals.h>
#include <iostream>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>

#include "ASDPlasticMaterialTraits.h"

#include "ASDPlasticMaterialGlobals.h"


// #include "MaterialInternalVariables.h"
#include "YieldFunctions/AllYieldFunctions.h"
#include "PlasticFlowDirections/AllPlasticFlowDirections.h"
#include "ElasticityModels/AllElasticityModels.h"
#include "AllASDModelParameterTypes.h"
#include "AllASDInternalVariableTypes.h"
#include "AllASDHardeningFunctions.h"

#include "utuple_storage.h"

// for print p, q, theta
#include <tuple>
#include <utility> // For std::pair
#include <map> // For std::pair
#include <limits>
#include <type_traits>

#include "std_tuple_concat.h"

// for debugging printing
#include <fstream>

#define ASDPlasticMaterial_MAXITER_BRENT 50
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



using namespace ASDPlasticMaterialGlobals;


// template <
//     class ElasticityType,
//     class YieldFunctionType,
//     class PlasticFlowType,
//     int thisClassTag,
//     class T >
template <
    class ElasticityType,
    class YieldFunctionType,
    class PlasticFlowType,
    int thisClassTag >
class ASDPlasticMaterial : public NDMaterial
{

public:

    // Concatenate the internal varibles into the storage
    using iv_concat_types = utuple_concat_unique_type <
                            typename YieldFunctionType::internal_variables_t,
                            typename PlasticFlowType::internal_variables_t >;
    using iv_storage_t = utuple_storage<iv_concat_types>;

    // Concatenate the model parameters into the parameters storage
    using extracted_parameters_t = utuple_concat_unique_type <
                                   ExtractNestedParameterTypes_t<typename YieldFunctionType::internal_variables_t>,
                                   ExtractNestedParameterTypes_t<typename PlasticFlowType::internal_variables_t>
                                   >;

    using parameters_concat_types = utuple_concat_type <
                                    typename YieldFunctionType::parameters_t,
                                    typename PlasticFlowType::parameters_t,
                                    typename ElasticityType::parameters_t,
                                    extracted_parameters_t,
                                    std::tuple<MassDensity>,
                                    std::tuple<InitialP0>
                                    >;
    using parameters_storage_t = utuple_storage<parameters_concat_types>;


    //==================================================================================================
    //  Constructors
    //==================================================================================================

    ASDPlasticMaterial( )
        : NDMaterial(0, thisClassTag)
    {

    }


    ASDPlasticMaterial(int tag)
        : NDMaterial(tag, thisClassTag)
    {


        TrialStress *= 0;
        CommitStress *= 0;
        TrialStrain *= 0;
        CommitStrain *= 0;
        TrialPlastic_Strain *= 0;
        CommitPlastic_Strain *= 0;

        first_step = true;
    }


    ~ASDPlasticMaterial(void)
    {

    }

    //==================================================================================================
    // To set internal variables values for the model
    //==================================================================================================
    auto getInternalVariablesNames() const
    {
        return iv_storage.getParameterNames();
    }

    int getInternalVariableSizeByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableSizeByName(iv_name);
    }

    int getInternalVariableIndexByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableIndexByName(iv_name);
    }

    auto setInternalVariableByName(const char * iv_name, int iv_size, double* iv_values)
    {
        cout << "  --->  Setting " << iv_name << " = ";
        for (int i = 0; i < iv_size; ++i)
        {
            cout << iv_values[i] << " ";
        }
        cout << endl;
        return iv_storage.setInternalVariableByName(iv_name, iv_size, iv_values);
    }



    //==================================================================================================
    // To set parameter values for the model
    //==================================================================================================
    auto getParameterNames() const
    {
        return parameters_storage.getParameterNames();
    }

    auto setParameterByName(const char * param_name, double param_value)
    {
        cout << "  --->  Setting " << param_name << " = " << param_value << endl;
        return parameters_storage.setParameterByName(param_name, param_value);
    }


    //==================================================================================================
    //  Class type function
    //==================================================================================================
    const char *getClassType(void) const
    {
        std::string name("ASDPlasticMaterial");

        return name.c_str();
    };

    double getRho(void)
    {
        return parameters_storage.template get<MassDensity>().value;
    }

    double getPressure(void)
    {
        using namespace ASDPlasticMaterialGlobals;
        return -CommitStress.trace() / 3;
    }

    std::string getYFName() const {return YieldFunctionType::NAME;}
    std::string getPFName() const {return PlasticFlowType::NAME;}
    std::string getELName() const {return ElasticityType::NAME;}
    std::string getIVName() const {return iv_storage.getVariableNamesAndHardeningLaws();}


    //==================================================================================================
    //  Set Trial strain and trial strain increment
    //==================================================================================================
    // For total strain-based elements.
    // Receives the current total strain at a GP.
    // This function then computes the incremental strain (subtracting from the committed one)
    // and sets the increment.
    // Returns a success flag from the call to setTrialStrainIncr

    int setTrialStrain(const Vector &v)
    {

        if (first_step)
        {
            double p0 = parameters_storage.template get<InitialP0>().value;
            TrialStress(0) = p0;
            TrialStress(1) = p0;
            TrialStress(2) = p0;
            CommitStress(0) = p0;
            CommitStress(1) = p0;
            CommitStress(2) = p0;
        }


        TrialStrain = VoigtVector::fromStrain(v);
        return setTrialStrainIncr( TrialStrain - CommitStrain );
    }

    // Directly sets the trial strain increment and does an explicit or implicit step.
    // Returns a flag depending on the result of the step.
    int setTrialStrainIncr( const VoigtVector &strain_increment )
    {
        // using namespace ASDPlasticMaterialGlobals;

        int exitflag = -1;

        // ==========================================
        // Check for quick return:
        // If strain_increment are all zeroes. Return
        // ==========================================
        // The quick return is useful for displacement_control.
        // In displacement_control, we have to do double update.
        // ==========================================

        // cout << "strain_increment = " << strain_increment << endl;

        double max_component = 0.0;
        for (int i = 0; i < 6; ++i)
            if (max_component < abs(strain_increment(i)))
            {
                max_component = abs(strain_increment(i));
            }
        // ==========================================
        // We do not use below(eps_norm) because 1E-8 will go to 1E-16. --> False return.
        // double eps_norm = strain_increment(i, j) * strain_increment(i, j);
        // ==========================================
        if (max_component == 0)
        {
            exitflag = 0;
            return exitflag;
        }
        // ==========================================

        switch (INT_OPT_constitutive_integration_method[this->getTag()])
        {
        case ASDPlasticMaterial_Constitutive_Integration_Method::Not_Set :
            exitflag = -1;
            cerr << "CEP::setTrialStrainIncr - Integration method not set!\n" ;
            break;
        case ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler :
            exitflag = this->Forward_Euler(strain_increment);
            break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler_Subincrement :
        //     exitflag = this->Forward_Euler_Subincrement(strain_increment);
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler :
        //     exitflag = this->Backward_Euler(strain_increment);;
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler_ddlambda :
        //     exitflag = this->Backward_Euler_ddlambda(strain_increment);;
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement :
        //     exitflag = this->Backward_Euler_ddlambda_Subincrement(strain_increment);;
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler_Crisfield :
        //     exitflag = this->Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield :
        //     exitflag = this->Multistep_Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial_Constitutive_Integration_Method::Modified_Euler_Error_Control :
        //     exitflag = this->Modified_Euler_Error_Control(strain_increment);
        //     break;
        case ASDPlasticMaterial_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control :
            exitflag = this->Runge_Kutta_45_Error_Control(strain_increment);;
            break;
        default:
            cerr << "ASDPlasticMaterial::setTrialStrainIncr - Integration method not available!\n" ;
            exitflag = -1;
        }

        return exitflag;
    }

//==================================================================================================
//  Getters
//==================================================================================================

    const VoigtMatrix &getTangentTensor( void )
    {
        // using namespace ASDPlasticMaterialGlobals;

        // // double yf_val = yf(TrialStress);

        // // VoigtMatrix& Eelastic = et(TrialStress);


        // // if (yf_val <= 0.0)
        // // {
        // //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        // // }
        // // else
        // // {
        // //     const VoigtVector& n = yf.df_dsigma_ij(intersection_stress);
        // //     const VoigtVector& m = pf(depsilon_elpl, intersection_stress);

        // //     double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, depsilon_elpl,  intersection_stress);
        // //     double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

        // //     //Compute tangent stiffness
        // //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / den;
        // // }

        // if (first_step)
        // {
        //     VoigtMatrix& Eelastic = et(TrialStress);
        //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
        //     first_step = false;
        // }

        return Stiffness;
    }

    const Vector &getStress(void)
    {
        static Vector result(6);
        TrialStress.toStress(result);
        return result;
    }

    const Vector &getStrain(void)
    {
        static Vector result(6);
        TrialStrain.toStrain(result);
        return result;
    }

    const Vector &getPstrain(void)
    {
        static Vector result(6);
        TrialPlastic_Strain.toStress(result);
        return result;
    }

    const Vector &getInternalVariableByPos(int pos)
    {
        static Vector return_vector(6);
        int find_pos = 0;

        // Note by J. Abell on Wed 06 Dec 2023 11:44:24
        //
        // The lambda capture of return_vector which is declared static above triggers
        // the following warning in g++ 11.4.0
        //
        // warning: capture of variable ‘return_vector’ with non-automatic storage duration
        //
        // Explanation:
        //   Because usually the lambda functions are used for threaded applications
        //   this warning is there to help with race conditions on the return_vector.
        //   In this case the warning is benign because opensees is not threaded
        //   at this level.
        //   We want to keep the static allocation of return_vector to avoid
        //   many calls to malloc (new) every time this function is called
        //   for performance reasons, so we have to live with the warning.
        iv_storage.apply([&pos, &find_pos, this](auto & internal_variable)
        {
            if (pos == find_pos)
            {
                auto &iv = internal_variable.trial_value;
                int iv_size = iv.size();
                return_vector.resize(iv_size);
                for (int i = 0; i < iv_size; ++i)
                {
                    return_vector(i) = iv(i);
                }
            }
            find_pos ++;
        });

        return return_vector;
    }

    const VoigtVector &getStressTensor( void )
    {
        return TrialStress;
    }

    const VoigtVector &getStrainTensor( void )
    {
        return TrialStrain;
    }

    const VoigtVector &getPlasticStrainTensor( void )
    {
        return TrialPlastic_Strain;
    }

    const VoigtVector  &getCommittedStressTensor(void)
    {
        return CommitStress;
    }

    const VoigtVector &getCommittedStrainTensor(void)
    {
        return CommitStrain;
    }

    const VoigtVector &getCommittedPlasticStrainTensor(void)
    {
        return CommitPlastic_Strain;
    }

    void ComputeTangentStiffness()
    {
        if (INT_OPT_tangent_operator_type[this->getTag()] == ASDPlasticMaterial_Tangent_Operator_Type::Elastic)
        {
            VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
            Stiffness = Eelastic;
        }
        // else if (INT_OPT_tangent_operator_type[this->getTag()] == ASDPlasticMaterial_Tangent_Operator_Type::Continuum)
        // {

        //     VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
        //     const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
        //     const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

        //     double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

        //     double den_after_corrector = n.transpose() * Eelastic * m - xi_star_h_star;

        //     VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
        //     Stiffness = Econtinuum;
        // }
        // else if (INT_OPT_tangent_operator_type[this->getTag()] == ASDPlasticMaterial_Tangent_Operator_Type::Secant)
        // {

        //     VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
        //     const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
        //     const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

        //     double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

        //     double den_after_corrector = n.transpose() * Eelastic * m - xi_star_h_star;

        //     VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
        //     Stiffness = (Econtinuum + Eelastic)/2;
        // }
    }

    const Matrix& getTangent()
    {
        static Matrix return_matrix(6, 6);

        copyToMatrixReference(Stiffness, return_matrix);

        return return_matrix;
    }


    const Matrix& getInitialTangent()
    {
        static Matrix return_matrix(6, 6);

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
        Stiffness = Eelastic;

        copyToMatrixReference(this->Stiffness, return_matrix);

        return return_matrix;
    }



//==================================================================================================
//  State commiting and reversion
//==================================================================================================

    // Forwards the trial variables to commited variables and calls commit state on the
    // BET classes
    int commitState(void)
    {
        // cout << "Committing!" << endl;

        // cout <<  "CommitStress = " << CommitStress.transpose() << endl;
        // cout <<  "CommitStrain = " << CommitStrain.transpose() << endl;
        // cout <<  "TrialStress = " << TrialStress.transpose() << endl;
        // cout <<  "TrialStrain = " << TrialStrain.transpose() << endl;

        // cout << "Committed!" << endl;

        CommitStress = TrialStress;
        CommitStrain = TrialStrain;
        CommitPlastic_Strain = TrialPlastic_Strain;

        // cout <<  "CommitStress = " << CommitStress.transpose() << endl;
        // cout <<  "CommitStrain = " << CommitStrain.transpose() << endl;
        // cout <<  "TrialStress = " << TrialStress.transpose() << endl;
        // cout <<  "TrialStrain = " << TrialStrain.transpose() << endl;


        iv_storage.commit_all();

        // cout << "Internal model info:";
        // Print(cout);
        // cout << "\n" << endl;

        if (first_step)
        {
            first_step = false;
        }

        return 0;
    }

    //Reverts the commited variables to the trials and calls revert on BET Classes.
    int revertToLastCommit(void)
    {

        TrialStress = CommitStress;
        TrialStrain = CommitStrain;
        TrialPlastic_Strain = CommitPlastic_Strain;

        iv_storage.revert_all();

        return 0;
    }

    int revertToStart(void)
    {
        cerr << "ASDPlasticMaterial::revertToStart - not implemented!!!\n" ;
        return -1;
    }

    NDMaterial *getCopy(void)
    {

        ASDPlasticMaterial <ElasticityType,
                           YieldFunctionType,
                           PlasticFlowType,
                           thisClassTag> *newmaterial = new ASDPlasticMaterial<ElasticityType,
        YieldFunctionType,
        PlasticFlowType,
        thisClassTag>(this->getTag());
        newmaterial->TrialStrain = this->TrialStrain;
        newmaterial->TrialStress = this->TrialStress;
        newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
        newmaterial->CommitStress = this->CommitStress;
        newmaterial->CommitStrain = this->CommitStrain;
        newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
        newmaterial->iv_storage = this->iv_storage;
        newmaterial->parameters_storage = this->parameters_storage;

        return newmaterial;
    }


    NDMaterial *getCopy(const char *type) {
        if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
            ASDPlasticMaterial <ElasticityType,
                               YieldFunctionType,
                               PlasticFlowType,
                               thisClassTag> *newmaterial = new ASDPlasticMaterial<ElasticityType,
            YieldFunctionType,
            PlasticFlowType,
            thisClassTag>(this->getTag());
            newmaterial->TrialStrain = this->TrialStrain;
            newmaterial->TrialStress = this->TrialStress;
            newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
            newmaterial->CommitStress = this->CommitStress;
            newmaterial->CommitStrain = this->CommitStrain;
            newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
            newmaterial->iv_storage = this->iv_storage;
            newmaterial->parameters_storage = this->parameters_storage;

            return newmaterial;
        } else
        {
            cout << "ASDPlasticMaterial::getCopy(const char *type) - Only 3D is currently supported. Use 3D elements!" << endl;
        }
        return 0;
    }


    Response *setResponse (const char **argv, int argc,
                           OPS_Stream & s)
    {

        static Vector return_vector(6);

        if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
            return new MaterialResponse(this, 1, this->getStress());
        else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
            return new MaterialResponse(this, 2, this->getStrain());
        else if (strcmp(argv[0], "pstrain") == 0 || strcmp(argv[0], "pstrains") == 0)
            return new MaterialResponse(this, 3, this->getPstrain());
        else
        {
            const char *iv_name = argv[0];

            int iv_size = this->getInternalVariableSizeByName(iv_name);
            int pos = this->getInternalVariableIndexByName(iv_name);

            return new MaterialResponse(this, 1000 + pos, Vector(iv_size));
        }

        return 0;
    }


    int getResponse (int responseID, Information & matInformation)
    {


        if (matInformation.theVector == 0)
            return 0;

        if (responseID == -1)
        {
            return -1;
        }
        else if (responseID == 1)
        {
            *(matInformation.theVector) = getStress();
        }
        else if (responseID == 2)
            *(matInformation.theVector) = getStrain();
        else if (responseID == 3)
            *(matInformation.theVector) = getPstrain();
        else if (responseID >= 1000)
        {
            int pos = responseID - 1000;
            *(matInformation.theVector) = getInternalVariableByPos(pos);
        }

        return 0;
    }

    // NDMaterial *getCopy(const char *code);
    const char *getType(void) const {return "ThreeDimensional";}

    int sendSelf(int commitTag, Channel & theChannel)
    {
        cerr << "ASDPlasticMaterial::sendSelf - not implemented!!!\n" ;


        return 0;
    }

    int recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
    {
        cerr << "ASDPlasticMaterial::recvSelf - not implemented!!!\n" ;


        return 0;
    }

    void Print(OPS_Stream & s, int flag = 0) {
        s << "ASDPlasticMaterial" << endln;
        s << "  Yield Function          : " << yf.NAME << endln;
        s << "  Plastic flow direction  : " << pf.NAME << endln;
        s << "  Elasticity              : " << et.NAME << endln;
        s << "  # of Internal variables : " << (int) iv_storage.size() <<  endln;
        iv_storage.print_components();
        s << "  # of Parameters         : " << (int) parameters_storage.size() <<  endln;
        parameters_storage.print_components();
    }

    void Print(ostream & s, int flag = 0)
    {
        using namespace ASDPlasticMaterialGlobals;

        s << "ASDPlasticMaterial" << endl;
        s << "  Yield Function          : " << yf.NAME << endl;
        s << "  Plastic flow direction  : " << pf.NAME << endl;
        s << "  Elasticity              : " << et.NAME << endl;
        s << "  # of Internal variables : " << iv_storage.size() <<  endl;
        iv_storage.print_components();
        s << "  # of Parameters         : " << parameters_storage.size() <<  endl;
        parameters_storage.print_components();

    }

    int getObjectSize()
    {
        int size = 0;

        // 6 3x3 VoigtVectors and 1 VoigtMatrix (3x3x3x3)
        size += (3 * 3 * 6 + 3 * 3 * 3 * 3) * sizeof(double);

        //Four pointers
        size += 4 * sizeof(YieldFunctionType*);

        //Whatever the base components size is
        size += sizeof(yf);//yf->getObjectSize();
        size += sizeof(et);//et->getObjectSize();
        size += sizeof(pf);//pf->getObjectSize();
        // size += sizeof(internal_variables);//internal_variables->getObjectSize();
        size += sizeof(NDMaterial);

        // size += static_cast<T*>(this)->getObjectSize();

        return size;
    }

    void setTrialStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterialGlobals;
        TrialStress = stress;
    }

    bool set_constitutive_integration_method(int method, int tangent, double f_relative_tol, double stress_relative_tol, int n_max_iterations, int return_to_yield_surface)
    {
        if ( method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Not_Set
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler_Crisfield
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Multistep_Forward_Euler
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Modified_Euler_Error_Control
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler_Subincrement
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler_ddlambda
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement
                || method == (int) ASDPlasticMaterial_Constitutive_Integration_Method::Full_Backward_Euler)
        {
            int tag = this->getTag();
            INT_OPT_constitutive_integration_method[tag] = (ASDPlasticMaterial_Constitutive_Integration_Method) method ;
            INT_OPT_tangent_operator_type[tag] = (ASDPlasticMaterial_Tangent_Operator_Type) tangent ;
            INT_OPT_f_relative_tol[tag] = f_relative_tol ;
            INT_OPT_stress_relative_tol[tag] = stress_relative_tol ;
            INT_OPT_n_max_iterations[tag] = n_max_iterations ;
            INT_OPT_return_to_yield_surface[tag] = return_to_yield_surface ;

            cout << "set_constitutive_integration_method:" << endl;
            cout << "   method = " << method << endl;
            cout << "   f_relative_tol = " << f_relative_tol << endl;
            cout << "   stress_relative_tol = " << stress_relative_tol << endl;
            cout << "   n_max_iterations = " << n_max_iterations << endl;
            cout << "   return_to_yield_surface = " << return_to_yield_surface << endl;

            return true;
        }
        else
        {
            cerr << "ASDPlasticMaterial::set_constitutive_integration_method - Unknown constitutive_integration_method\n";
            return false;
        }
    }

protected:

    void setTrialPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterialGlobals;
        TrialPlastic_Strain = strain;
    }

    void setCommitStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterialGlobals;
        CommitStress = stress;
    }

    void setCommitStrain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterialGlobals;
        CommitStrain = strain;
    }

    void setCommitPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterialGlobals;
        CommitPlastic_Strain = strain;
    }

    void setStiffness(const VoigtMatrix & stiff)
    {
        using namespace ASDPlasticMaterialGlobals;
        Stiffness = stiff;
    }


    void setStressTensor(VoigtVector & stress)
    {
        using namespace ASDPlasticMaterialGlobals;
        CommitStress = stress;
        TrialStress = stress;
        return;
    }




private:


    int Forward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterialGlobals;



        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon *= 0;
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        dsigma = Eelastic * depsilon;

        TrialStress = sigma + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(sigma, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = epsilon  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;

            Eelastic = et(intersection_stress, parameters_storage);
            TrialStress  += Eelastic * depsilon_elpl;

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

            double xi_star_h_star = yf.xi_star_h_star( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
            double den = n.transpose() * Eelastic * m - xi_star_h_star;

            //Compute the plastic multiplier
            // Yuan and Boris 10Sep2016, change from "==0"
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor("m", m);
                printTensor("n", n);
                cout << "xi_star_h_star = " << xi_star_h_star << endl;
                cout << "den = " << den << endl;
                printTensor("depsilon_elpl", depsilon_elpl);
                return -1;
            }

            double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
            dLambda /= den;

            if (dLambda <= 0)
            {
                cout << "CEP - dLambda = " << dLambda << " <= 0\n";
                printTensor("m", m);
                printTensor("n", n);
                cout << "xi_star_h_star = " << xi_star_h_star << endl;
                cout << "den = " << den << endl;
                printTensor("depsilon_elpl", depsilon_elpl);
            }

            // Update the trial plastic strain.
            TrialPlastic_Strain += dLambda * m;

            // This code iterates internal variables and updates the trial values
            iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                internal_variable.trial_value += dLambda * h;
            });


            // Deal with the APEX if needed
            // if constexpr (yf_has_apex<YieldFunctionType>::value) {
            // {
            //     static VoigtVector small_stress(3, 3, 0.0);
            //     small_stress *= 0;
            //     // The small value 50*Pa refers to the lowest confinement test:
            //     // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
            //     double DP_k = yf.get_k();
            //     double DP_p = 50 ;
            //     // To make it on the yield surface, the q is equal to k*p.
            //     double DP_q = DP_k * DP_p ;
            //     // Assume the triaxial conditions sigma_2 = sigma_3.
            //     small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
            //     small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
            //     small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;
            //     static VoigtVector dstress(3, 3, 0.0);
            //     dstress(i, j) = small_stress(i, j) - sigma(i, j);
            //     static VoigtVector depsilon_Inv(3, 3, 0.0);
            //     depsilon_Inv = depsilon.Inv();

            //     // Return results (member variables) :
            //     Stiffness(i, j, k, l) = dstress(i, j) * depsilon_Inv(k, l);
            //     TrialStress(i, j) = small_stress(i, j);
            //     // plastic_strain and internal variables are already updated.
            //     return 0;
            // }

            //Correct the trial stress
            TrialStress = TrialStress - dLambda * Eelastic * m;


            // ============================================================================================
            // Add the additional step: returning to the yield surface.
            // This algorithm is based on Crisfield(1996). Page 171. Section 6.6.3
            // After this step, the TrialStress(solution), TrialPlastic_Strain, and Stiffness will be updated to the yield surface.
            // ============================================================================================
            if (INT_OPT_return_to_yield_surface[this->getTag()])
            {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                // In the function below, depsilon_elpl is actually not used at all in xi_star_h_star
                // double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m_after_corrector,  TrialStress);
                double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);
                double dLambda_after_corrector = yf_val_after_corrector / (
                                                     n_after_corrector.transpose() * Eelastic * m_after_corrector - xi_star_h_star_after_corrector
                                                 );
                TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;

                //double den_after_corrector = n_after_corrector.transpose() * Eelastic * m_after_corrector - xi_star_h_star_after_corrector;
                // Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m_after_corrector(p, q)) * (n_after_corrector(r, s) * Eelastic(r, s, k, l) ) / den_after_corrector;
            }
            else
            {
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
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
            else
            {
                errorcode = 0;
            }

            ComputeTangentStiffness();

        }

        return errorcode;
    }



    // int Forward_Euler_Subincrement(const VoigtVector &strain_incr, bool const& with_return2yield_surface)
    int Forward_Euler_Subincrement(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterialGlobals;
        int errorcode = -1;



        return errorcode;
    }

    std::pair<double, VoigtVector> CalculateLambdaM(
        const VoigtVector & thisSigma,
        const VoigtVector & depsilon_elpl,
        const parameters_storage_t &this_parameters_storage,
        const iv_storage_t &this_iv_storage)
    {
        VoigtMatrix Eelastic = et(thisSigma, this_parameters_storage);
        const VoigtVector& n = yf.df_dsigma_ij(thisSigma, this_iv_storage, this_parameters_storage);
        VoigtVector m = pf(depsilon_elpl, thisSigma, this_iv_storage, this_parameters_storage);
        double xi_star_h_star = yf.xi_star_h_star(depsilon_elpl, m, thisSigma, this_iv_storage, this_parameters_storage);

        double den = n.transpose() * Eelastic * m - xi_star_h_star;
        double dLambda = (den != 0) ? (n.transpose() * Eelastic * depsilon_elpl).value() / den : 0;

        if (dLambda != dLambda)
        {
            cout << "CalculateLambdaM error" << endl;
            cout << "yf = " << yf(thisSigma, this_iv_storage, this_parameters_storage) << endl;
            cout << "thisSigma = " << thisSigma.transpose() << endl;
            cout << "depsilon_elpl = " << depsilon_elpl.transpose() << endl;
            cout << "n = " << n.transpose() << endl;
            cout << "m = " << m.transpose() << endl;
            cout << "xi_star_h_star = " << xi_star_h_star << endl;
            cout << "den = " << den << endl;
            cout << "Eelastic = " << Eelastic << endl;

            cout << "IVSTORAGE" << endln;
            this_iv_storage.print_components();

            cout << "PARAMSTORAGE" << endln;
            this_parameters_storage.print_components();
        }


        return std::make_pair(dLambda, m);
    }

    int Runge_Kutta_45_Error_Control(const VoigtVector & strain_incr)
    {
        // cout << "Runge_Kutta_45_Error_Control" << endl;

        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon *= 0;
        depsilon = strain_incr;


        iv_storage.revert_all();


        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);

        dsigma = Eelastic * depsilon;


        TrialStress = CommitStress + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(CommitStress, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = 1e-3, TolE = this->INT_OPT_stress_relative_tol[this->getTag()];


            VoigtVector next_Sigma = TrialStress;
            iv_storage_t next_iv_storage = iv_storage;
            VoigtVector next_EpsilonPl = CommitPlastic_Strain;

            VoigtVector prev_Sigma = TrialStress;
            VoigtVector prev_EpsilonPl = CommitPlastic_Strain;
            iv_storage_t prev_iv_storage = iv_storage;

            VoigtVector dSigma, dSigma1, dSigma2, dSigma3, dSigma4, dSigma5, dSigma6;
            VoigtVector dEpsilonPl, dEpsilonPl1, dEpsilonPl2, dEpsilonPl3, dEpsilonPl4, dEpsilonPl5, dEpsilonPl6;
            iv_storage_t iv_storage1 = iv_storage;
            iv_storage_t iv_storage2 = iv_storage;
            iv_storage_t iv_storage3 = iv_storage;
            iv_storage_t iv_storage4 = iv_storage;
            iv_storage_t iv_storage5 = iv_storage;


            int niter = 0;
            double maxStepError = 0;
            while (T < 1.0)
            {
                niter ++;

                VoigtVector dEPS = dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                //Delta 1
                std::tie(dLambda, m) = CalculateLambdaM(prev_Sigma, dEPS, parameters_storage, prev_iv_storage);
                dSigma1  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl1 = dLambda * m;
                iv_storage1.apply([&m, &dLambda, &prev_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, prev_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + 0.5 * dLambda * h;
                });
                next_Sigma =  prev_Sigma + 0.5 * dSigma1;

                //Delta 2
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage1);
                dSigma2  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl2 = dLambda * m;
                iv_storage2.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, this](auto & iv2)
                {
                    using VT1 = std::decay_t<decltype(iv2)>;
                    const VT1 &iv1 = iv_storage1.template get<VT1>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;

                    auto h2 = iv2.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH2 = dLambda * h2;

                    iv2.trial_value = iv2.committed_value + 0.25 * (DH1 + DH2);
                });
                next_Sigma =  prev_Sigma + 0.25 * (dSigma1 + dSigma2);

                //Delta 3
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage2);
                dSigma3  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl3 = dLambda * m;
                iv_storage3.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, this](auto & iv3)
                {
                    using VT = std::decay_t<decltype(iv3)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;

                    auto h3 = iv3.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH3 = dLambda * h3;

                    iv3.trial_value = iv3.committed_value +  -DH2 + 2 * DH3;
                });
                next_Sigma =  prev_Sigma - dSigma2 + 2 * dSigma3;

                //Delta 4
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage3);
                dSigma4  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl4 = dLambda * m;
                iv_storage4.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, this](auto & iv4)
                {
                    using VT = std::decay_t<decltype(iv4)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;

                    auto h4 = iv4.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH4 = dLambda * h4;

                    iv4.trial_value = iv4.committed_value +  (7 * DH1 + 10 * DH2 + DH4) / 27;
                });
                next_Sigma =  prev_Sigma + (7 * dSigma1 + 10 * dSigma2 + dSigma4) / 27;

                //Delta 5
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage4);
                dSigma5  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl5 = dLambda * m;
                iv_storage5.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, this](auto & iv5)
                {
                    using VT = std::decay_t<decltype(iv5)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;

                    auto h5 = iv5.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH5 = dLambda * h5;

                    iv5.trial_value = iv5.committed_value +  (28 * DH1 - 125 * DH2 + 546 * DH3 + 54 * DH4 - 378 * DH5) / 625;
                });
                next_Sigma =  prev_Sigma + (28 * dSigma1 - 125 * dSigma2 + 546 * dSigma3 + 54 * dSigma4 - 378 * dSigma5) / 625;

                //Delta 6 - The final predictor
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage5);
                dSigma6  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl6 = dLambda * m;
                next_iv_storage.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, &iv_storage5, this](auto & niv)
                {
                    using VT = std::decay_t<decltype(niv)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;
                    const VT &iv5 = iv_storage5.template get<VT>();
                    auto DH5 = iv5.trial_value - iv5.committed_value;

                    auto h6 = niv.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH6 = dLambda * h6;

                    niv.trial_value = niv.committed_value +  ( DH1 +  4 * DH3 + DH4) / 6;
                });
                dSigma =  (dSigma1 + 4 * dSigma3 + dSigma4) / 6;
                dEpsilonPl =  (dEpsilonPl1 + 4 * dEpsilonPl3 + dEpsilonPl4) / 6;
                next_Sigma =  prev_Sigma + dSigma;
                next_EpsilonPl = prev_EpsilonPl + dEpsilonPl;


                if (dSigma != dSigma)
                {
                    cout << "ASDPlasticMaterial::Runge_Kutta_45_Error_Control Integration error" << endl;
                    cout << "T = " << T << endl;
                    cout << "dT = " << dT << endl;
                    cout << "dSigma1 = " << dSigma1.transpose() << endl;
                    cout << "dSigma2 = " << dSigma2.transpose() << endl;
                    cout << "dSigma3 = " << dSigma3.transpose() << endl;
                    cout << "dSigma4 = " << dSigma4.transpose() << endl;
                    cout << "dSigma5 = " << dSigma5.transpose() << endl;
                    cout << "dSigma6 = " << dSigma6.transpose() << endl;
                    cout << "m = " << m.transpose() << endl;
                    cout << "dEpsilonPl = " << dEpsilonPl.transpose() << endl;
                    cout << "TrialStress = " << TrialStress.transpose() << endl;
                    cout << "dEPS = " << dEPS.transpose() << endl;
                    exit(-1);
                }




                //Stress norm and stress error
                double stressNorm = next_Sigma.norm();
                double curStepError1 = (-42 * dSigma1 - 224 * dSigma3 - 21 * dSigma4 + 162 * dSigma5 + 125 * dSigma6).norm() / 336;
                if (stressNorm >= 0.5) { curStepError1 /= (2 * stressNorm); }
                //Internal variables norm and internal variables error (TODO)

                double curStepError = curStepError1; //fmax(curStepError1, curStepError2);

                //Check convergence and adjust integration timestep
                if (curStepError > TolE)
                {
                    double q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.1);

                    if (dT == dT_min) {

                        prev_Sigma = next_Sigma;
                        prev_EpsilonPl = next_EpsilonPl;
                        prev_iv_storage = next_iv_storage;

                        T += dT;
                    }
                    dT = fmax(q * dT, dT_min);
                }
                else {

                    prev_Sigma = next_Sigma;
                    prev_EpsilonPl = next_EpsilonPl;
                    prev_iv_storage = next_iv_storage;

                    double q = fmin(0.8 * pow(TolE / curStepError, 0.2), 2.0);
                    T += dT;
                    dT = fmax(q * dT, dT_min);
                    dT = fmin(dT, 1 - T);

                    maxStepError = max(maxStepError, curStepError);
                }

            }

            TrialStress = next_Sigma;
            TrialPlastic_Strain = next_EpsilonPl;
            iv_storage = next_iv_storage;

            // cout << "ASDPlasticMaterial::Runge_Kutta_45_Error_Control niter = " << niter << " maxStepError = " << maxStepError << " TolE = " << TolE << endl;;

            // ============================================================================================
            // Add the additional step: returning to the yield surface.
            // This algorithm is based on Crisfield(1996). Page 171. Section 6.6.3
            // After this step, the TrialStress(solution), TrialPlastic_Strain, and Stiffness will be updated to the yield surface.
            // ============================================================================================
            if (INT_OPT_return_to_yield_surface[this->getTag()])
            {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                // In the function below, depsilon_elpl is actually not used at all in xi_star_h_star
                // double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m_after_corrector,  TrialStress);
                double xi_star_h_star_after_corrector = yf.xi_star_h_star( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                double dLambda_after_corrector = yf_val_after_corrector / (
                                                     n_after_corrector.transpose() * Eelastic * m_after_corrector - xi_star_h_star_after_corrector
                                                 );
                TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial::Runge_Kutta_45_Error_Control  Numeric error!\n";
                printTensor("TrialStress = " , TrialStress);
                printTensor("CommitStress = " , CommitStress);
                printTensor("depsilon = " , depsilon);
                printTensor("dsigma   = " , dsigma);
                printTensor("intersection_stress = " , intersection_stress);
                printTensor4("Eelastic = " , Eelastic);
                printTensor4("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor("n = " , yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage) );
                printTensor("m = " , pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage) );
                cout << "xi_star_h_star  = " << yf.xi_star_h_star( depsilon_elpl, pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage),  TrialStress, iv_storage, parameters_storage) << endl;

                errorcode = -1;

                exit(-1);
            }

            ComputeTangentStiffness();
        }

        return 0;
    }

    // int Backward_Euler(const VoigtVector &strain_incr, bool debugrun = false, int NSteps = 1)
    // {
    //     using namespace ASDPlasticMaterialGlobals;
    //     int errorcode = 0;
    //     // // =====Sub-increments =====
    //     // // Problem is in the stiffness
    //     // double max_comp = Max_abs_Component(strain_incr)/ (double)NSteps;
    //     // double allowed_increment_magnitude = 1E-5;
    //     // if(max_comp > allowed_increment_magnitude){
    //     //     int required_Nstep = (int)(max_comp/ allowed_increment_magnitude)+1 ;
    //     //     Backward_Euler(strain_incr, false, required_Nstep);
    //     //     return 0;
    //     // }
    //     // // =====Sub-increments =====END
    //     internal_variables.revert();
    //     internal_variables.commit_tmp();
    //     static VoigtVector depsilon(3, 3, 0);
    //     static VoigtVector PredictorStress(3, 3, 0);
    //     const VoigtVector& sigma = CommitStress;
    //     depsilon(i, j) = strain_incr(i, j) / NSteps;
    //     TrialStress(i, j) = sigma(i, j);
    //     TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
    //     TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

    //     // Uncomment next line to enable the subincrement in Backward_Euler.
    //     // for (int step = 0; step < NSteps; step++)
    //     {
    //         VoigtMatrix& Eelastic = et(sigma);
    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

    //         dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

    //         TrialStress(i, j) +=  dsigma(i, j);
    //         PredictorStress(i, j) = TrialStress(i, j);

    //         double yf_val_start = yf(sigma);
    //         double yf_val_end = yf(PredictorStress);

    //         // printTensor (" CommitStress       " , CommitStress);
    //         // printTensor (" depsilon           " , depsilon);
    //         // printTensor (" dsigma             " , dsigma);
    //         // printTensor (" PredictorStress    " , PredictorStress);
    //         // cout <<      " yf_val_start        = " << yf_val_start << endl;
    //         // cout <<      " yf_val_end          = " << yf_val_end << endl;
    //         // printTensor (" TrialStress        " , TrialStress);

    //         if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //         {
    //             VoigtMatrix& Eelastic = et(TrialStress);
    //             Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
    //         }
    //         else  //Plasticity
    //         {
    //             static VoigtVector ResidualStress(3, 3, 0);
    //             static VoigtVector TrialStress_prev(3, 3, 0);
    //             ResidualStress *= 0;
    //             TrialStress_prev *= 0;
    //             double normResidualStress = -1;
    //             double stress_relative_error = this->stress_relative_tol * 10; //
    //             double dLambda = 0;
    //             double denominator = 0;
    //             int iteration_count = 0;

    //             double yf_PredictorStress = yf(PredictorStress);
    //             double yf_TrialStress = yf(TrialStress);

    //             // vonMises does NOT enter this part.
    //             // DruckerPrager requries this part.
    //             if (yf.hasCorner() && yf.in_Apex(PredictorStress))
    //             {
    //                 static VoigtVector small_stress(3, 3, 0.0);
    //                 small_stress *= 0;
    //                 // The small value 50*Pa refers to the lowest confinement test:
    //                 // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
    //                 double DP_k = yf.get_k();
    //                 double DP_p = 50 ;
    //                 // To make it on the yield surface, the q is equal to k*p.
    //                 double DP_q = DP_k * DP_p ;
    //                 // Assume the triaxial conditions sigma_2 = sigma_3.
    //                 small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
    //                 small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
    //                 small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

    //                 // (1) Update the trial stress
    //                 TrialStress(i, j) = small_stress(i, j);
    //                 // (2) Update the trial plastic strain
    //                 const VoigtVector& n = yf.df_dsigma_ij(small_stress);
    //                 const VoigtVector& m = pf(depsilon, small_stress);
    //                 const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
    //                 double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //                 double dLambda = yf_PredictorStress / denominator;
    //                 TrialPlastic_Strain(i, j) += dLambda * m(i, j);
    //                 // (3) Update the internal variables
    //                 internal_variables.evolve(dLambda, depsilon, m, TrialStress);
    //                 internal_variables.commit_tmp();
    //                 // (4) Update the stiffness.
    //                 // Backward_Euler use the inconsistent stiffness.
    //                 // Full_Backward_Euler use the consistent stiffness.
    //                 Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //                 // // ===================================================================
    //                 // // Working on consistent stiffness
    //                 // // ===================================================================
    //                 // // Construct the Tensor T
    //                 static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //                 IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //                 static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //                 dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //                 static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //                 Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //                 // // ===================================================================
    //                 // // Construct the Tensor H
    //                 static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //                 dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
    //                 static VoigtVector H(3, 3, 0.0);
    //                 H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //                 static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //                 // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
    //                 if ( !inverse4thTensor(Ts, invT) )
    //                 {
    //                     return 0;
    //                 }
    //                 static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //                 R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //                 // // ===================================================================
    //                 // // Construct the consistent stiffness
    //                 double denomin{0.0};
    //                 denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
    //                 Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;
    //                 // // ===================================================================
    //                 // // consistent stiffness  END
    //                 // // ===================================================================

    //                 // Skip the following iteration procedures.
    //                 return 0;
    //             }

    //             double f_relative_error = this->f_relative_tol * 10;

    //             while ((stress_relative_error > this-> stress_relative_tol ||
    //                     f_relative_error > this->f_relative_tol) &&
    //                     iteration_count < this->n_max_iterations
    //                   )
    //             {
    //                 iteration_count++;
    //                 // In this approach, we always start from
    //                 // PredictorStress and the initial committed internal variables.
    //                 internal_variables.revert_tmp();

    //                 const VoigtVector& n = yf.df_dsigma_ij(TrialStress);
    //                 const VoigtVector& m = pf(depsilon, TrialStress);
    //                 double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);
    //                 denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

    //                 //Compute the plastic multiplier
    //                 // if (abs(denominator) < MACHINE_EPSILON)
    //                 // {
    //                 //     cout << "CEP - denominator of plastic multiplier n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star = 0\n";
    //                 //     printTensor("m", m);
    //                 //     printTensor("n", n);
    //                 //     cout << "xi_star_h_star" << xi_star_h_star << endl;
    //                 //     cout << "denominator" << denominator << endl;
    //                 //     printTensor("depsilon", depsilon);
    //                 //     return -1;
    //                 // }
    //                 dLambda =  yf_PredictorStress / denominator;

    //                 TrialStress_prev(i, j) = TrialStress(i, j);

    //                 TrialStress(i, j) = PredictorStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
    //                 internal_variables.evolve(dLambda, depsilon, m, TrialStress);
    //                 // // ============================================================
    //                 // // Line search in the constitutive level
    //                 // // ============================================================
    //                 double subdLambda = dLambda;
    //                 while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
    //                 {
    //                     internal_variables.revert_tmp();
    //                     subdLambda *= 0.5;
    //                     TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
    //                     // Internal variable should be updated with subdLambda.
    //                     // The condition in while-loop will use the internal variables.
    //                     internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
    //                 }
    //                 // ============================================================
    //                 dLambda = subdLambda ;
    //                 yf_TrialStress = yf(TrialStress);

    //                 ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

    //                 normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

    //                 // update the stress and f relative error.
    //                 stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
    //                 f_relative_error = abs(yf_TrialStress / yf_PredictorStress);

    //                 double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
    //                 if (norm_trial_stress != norm_trial_stress || debugrun) //check for nan
    //                 {
    //                     cout << "=============================================================" << endl;
    //                     cout << "\nIteration # " << iteration_count <<  endl;
    //                     printTensor (" CommitStress       " , CommitStress);
    //                     double p, q, theta;
    //                     std::tie(p, q, theta) = getpqtheta(CommitStress);
    //                     fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
    //                     fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
    //                     fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
    //                     printTensor (" depsilon           " , depsilon);
    //                     printTensor (" dsigma             " , dsigma);
    //                     printTensor (" PredictorStress    " , PredictorStress);
    //                     std::tie(p, q, theta) = getpqtheta(PredictorStress);
    //                     fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
    //                     fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
    //                     fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
    //                     cout <<      " yf_val_start        = " << yf_val_start << endl;
    //                     cout <<      " yf_val_end          = " << yf_val_end << endl;
    //                     printTensor (" n                  " , n );
    //                     printTensor (" m                  " , m );
    //                     cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
    //                     cout <<      " denominator         = " << denominator << endl;
    //                     cout <<      " dLambda             = " << dLambda << endl;
    //                     printTensor (" TrialStress        " , TrialStress);
    //                     std::tie(p, q, theta) = getpqtheta(TrialStress);
    //                     fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
    //                     fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
    //                     fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

    //                     printTensor (" TrialStress_prev   " , TrialStress_prev);
    //                     std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
    //                     fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
    //                     fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
    //                     fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

    //                     printTensor (" ResidualStress     " , ResidualStress );
    //                     cout <<      " normResidualStress  = " << normResidualStress << endl;
    //                     cout <<      " stress_relative_error      = " << stress_relative_error << endl;
    //                     cout <<      " f_relative_error      = " << f_relative_error << endl;
    //                     cout <<      " iteration_count     = " << iteration_count << endl;
    //                     cout <<      " .........................................................." << endl;
    //                     cout <<      "  Internal variables:" << endl;
    //                     printTensor (" back_stress (alpha)      " , yf.get_alpha() );
    //                     cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
    //                     cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
    //                     cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
    //                     cout << "===============================================END===========" << endl;
    //                     cout << endl << endl;
    //                 }

    //             } // while for el-pl this step

    //             //Update the trial plastic strain.
    //             const VoigtVector& n = yf.df_dsigma_ij(TrialStress);
    //             const VoigtVector& m = pf(depsilon, TrialStress);
    //             const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);

    //             denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //             TrialPlastic_Strain(i, j) += dLambda * m(i, j);
    //             internal_variables.commit_tmp();

    //             Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //             // ===================================================================
    //             // Working on consistent stiffness
    //             // ===================================================================
    //             // Update the m and n:
    //             static VoigtVector nf(3, 3, 0.0);
    //             nf = yf.df_dsigma_ij(TrialStress);
    //             static VoigtVector mf(3, 3, 0.0);
    //             mf = pf(depsilon, TrialStress);
    //             // Construct the Tensor T
    //             static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //             IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //             static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //             dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //             static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //             Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //             // ===================================================================
    //             // Construct the Tensor H
    //             static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //             dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
    //             static VoigtVector H(3, 3, 0.0);
    //             H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //             static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //             // If cannot inverse T, return directly.
    //             // So "Stiffness" is the inconsistent stiffness tensor.
    //             if ( !inverse4thTensor(Ts, invT) )
    //             {
    //                 return 0;
    //             }
    //             static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //             R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //             // ===================================================================
    //             // Construct the consistent stiffness
    //             double denomin{0.0};
    //             double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
    //             denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

    //             Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;

    //             // ===================================================================
    //             // Working on consistent stiffness  END
    //             // ===================================================================

    //         } //else  //Plasticity
    //         internal_variables.commit_tmp();
    //     } //End the loop for (int step = 0; step < NSteps; step++)
    //     return errorcode;
    // }

    // int Backward_Euler_ddlambda(const VoigtVector &strain_incr, bool debugrun = false)
    // {
    //     using namespace ASDPlasticMaterialGlobals;
    //     int errorcode = 0;

    //     internal_variables.revert();
    //     internal_variables.commit_tmp();
    //     static VoigtVector depsilon(3, 3, 0);
    //     static VoigtVector PredictorStress(3, 3, 0);
    //     const VoigtVector& sigma = CommitStress;
    //     depsilon(i, j) = strain_incr(i, j) ;
    //     TrialStress(i, j) = sigma(i, j);
    //     TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
    //     TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

    //     VoigtMatrix& Eelastic = et(sigma);
    //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

    //     dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

    //     TrialStress(i, j) +=  dsigma(i, j);
    //     PredictorStress(i, j) = TrialStress(i, j);

    //     double yf_val_start = yf(sigma);
    //     double yf_val_end = yf(PredictorStress);

    //     // printTensor (" CommitStress       " , CommitStress);
    //     // printTensor (" depsilon           " , depsilon);
    //     // printTensor (" dsigma             " , dsigma);
    //     // printTensor (" PredictorStress    " , PredictorStress);
    //     // cout <<      " yf_val_start        = " << yf_val_start << endl;
    //     // cout <<      " yf_val_end          = " << yf_val_end << endl;
    //     // printTensor (" TrialStress        " , TrialStress);

    //     if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //     {
    //         VoigtMatrix& Eelastic = et(TrialStress);
    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
    //     }
    //     else  //Plasticity
    //     {
    //         static VoigtVector ResidualStress(3, 3, 0);
    //         static VoigtVector TrialStress_prev(3, 3, 0);
    //         ResidualStress *= 0;
    //         TrialStress_prev *= 0;
    //         double normResidualStress = -1;
    //         double stress_relative_error = this->stress_relative_tol * 10; //
    //         double dLambda = 0;
    //         double denominator = 0;
    //         int iteration_count = 0;

    //         double yf_PredictorStress = yf(PredictorStress);
    //         double yf_TrialStress = yf(TrialStress);
    //         // vonMises does NOT enter this part.
    //         // DruckerPrager requries this part.
    //         if (yf.hasCorner() && yf.in_Apex(PredictorStress))
    //         {
    //             static VoigtVector small_stress(3, 3, 0.0);
    //             small_stress *= 0;
    //             // The small value 50*Pa refers to the lowest confinement test:
    //             // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
    //             double DP_k = yf.get_k();
    //             double DP_p = 50 ;
    //             // To make it on the yield surface, the q is equal to k*p.
    //             double DP_q = DP_k * DP_p ;
    //             // Assume the triaxial conditions sigma_2 = sigma_3.
    //             small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
    //             small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
    //             small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

    //             // (1) Update the trial stress
    //             TrialStress(i, j) = small_stress(i, j);
    //             // (2) Update the trial plastic strain
    //             const VoigtVector& n = yf.df_dsigma_ij(small_stress);
    //             const VoigtVector& m = pf(depsilon, small_stress);
    //             const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
    //             double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //             double dLambda = yf_PredictorStress / denominator;
    //             TrialPlastic_Strain(i, j) += dLambda * m(i, j);
    //             // (3) Update the internal variables
    //             internal_variables.evolve(dLambda, depsilon, m, TrialStress);
    //             internal_variables.commit_tmp();
    //             // (4) Update the stiffness.
    //             Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //             // // =========================================================
    //             // // Working on consistent stiffness
    //             // // ========================================================
    //             // // Construct the Tensor T
    //             static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //             IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //             static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //             dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //             static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //             Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //             // // Construct the Tensor H
    //             static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //             dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
    //             static VoigtVector H(3, 3, 0.0);
    //             H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //             static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //             // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
    //             if ( !inverse4thTensor(Ts, invT) )
    //             {
    //                 return 0;
    //             }
    //             static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //             R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //             // // =======================================================
    //             // // Construct the consistent stiffness
    //             double denomin{0.0};
    //             denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
    //             Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;

    //             // // =====================================================
    //             // // Working on consistent stiffness  END
    //             // // =====================================================

    //             // skip the following iteration
    //             return 0;
    //         }//END yield surface APEX check

    //         double f_relative_error = this->f_relative_tol * 10;

    //         while ((stress_relative_error > this-> stress_relative_tol ||
    //                 f_relative_error > this->f_relative_tol) &&
    //                 iteration_count < this->n_max_iterations
    //               )
    //         {
    //             iteration_count++;
    //             internal_variables.revert_tmp();

    //             const VoigtVector& n = yf.df_dsigma_ij(TrialStress);
    //             const VoigtVector& m = pf(depsilon, TrialStress);
    //             double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);
    //             denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //             if (abs(denominator) < MACHINE_EPSILON)
    //             {
    //                 cout << "CEP - denominator of plastic multiplier  \
    //                            n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star = 0\n";
    //                 printTensor("m", m);
    //                 printTensor("n", n);
    //                 cout << "xi_star_h_star" << xi_star_h_star << endl;
    //                 cout << "denominator" << denominator << endl;
    //                 printTensor("depsilon", depsilon);
    //                 return -1;
    //             }
    //             // The first dLambda is based on the elastic predictor. The next dLambda's are updated by ddlambda
    //             TrialStress_prev(i, j) = TrialStress(i, j);

    //             dLambda =  yf_PredictorStress / denominator;
    //             TrialStress(i, j) = TrialStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
    //             internal_variables.evolve(dLambda, depsilon, m, TrialStress);

    //             // In the first step, TrialStress is the predictor stress.
    //             // =======Line Search ====
    //             // // Zeta is the cooefficent for line search.
    //             double subdLambda = dLambda;
    //             while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
    //             {
    //                 internal_variables.revert_tmp();
    //                 subdLambda *= 0.5;
    //                 TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
    //                 // update the internal variables with subdLambda and commit to update yield surface.
    //                 internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
    //             }
    //             dLambda = subdLambda;
    //             ResidualStress(i, j) =  TrialStress(i, j) - TrialStress_prev(i, j);
    //             // ==========================================
    //             // Works to calculate the ddlambda
    //             // ==========================================
    //             // Update the m and n:
    //             static VoigtVector n_new(3, 3, 0.0);
    //             n_new = yf.df_dsigma_ij(TrialStress);
    //             static VoigtVector m_new(3, 3, 0.0);
    //             m_new = pf(depsilon, TrialStress);
    //             // Construct the Tensor T
    //             static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0.0);
    //             IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //             static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //             dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //             static VoigtMatrix T_new(3, 3, 3, 3, 0.0);
    //             T_new(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //             static VoigtMatrix invTnew(3, 3, 3, 3, 0.0);
    //             if ( !inverse4thTensor(T_new, invTnew) )
    //             {
    //                 // cout<<"Singularity. Cannot inverse tensor T! "<<endl;
    //                 return 0;
    //             }

    //             static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //             dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m_new, TrialStress);
    //             static VoigtVector H_new(3, 3, 0.0);
    //             H_new(k, l) = m_new(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //             double const& xi_star_h_star_new = yf.xi_star_h_star( depsilon, m_new,  TrialStress);

    //             double denomina_ddlambda = 0.0;
    //             static VoigtVector E_times_H (3, 3, 0.0);
    //             Eelastic = et(TrialStress);
    //             E_times_H(i, j) = Eelastic(i, j, k, l) * H_new(k, l);
    //             denomina_ddlambda = n_new(r, s) * (E_times_H(p, q) * invTnew(p, q, r, s)) - xi_star_h_star_new ;
    //             double numerator_ddlambda = 0.0;
    //             numerator_ddlambda = yf(TrialStress) - n_new(r, s) * (ResidualStress(i, j) * invTnew(i, j, r, s)) ;
    //             double ddlambda =   numerator_ddlambda / denomina_ddlambda ;
    //             static VoigtVector dSigma_(3, 3, 0.0);
    //             dSigma_(r, s) = -
    //                             (
    //                                 ResidualStress(i, j) + ddlambda * Eelastic(i, j, p, q) * H_new(p, q)
    //                             )
    //                             * invTnew(i, j, r, s) ;
    //             TrialStress_prev(i, j) = TrialStress(i, j) ;
    //             TrialStress(i, j) =  TrialStress(i, j) + dSigma_(i, j) ;
    //             double yf_TrialStress_prev = yf(TrialStress_prev) ;
    //             internal_variables.evolve(ddlambda, depsilon, m_new, TrialStress);

    //             // =======Line Search ddlambda (double delta ) =======
    //             double subddlambda = ddlambda;
    //             double zeta = 1.0;
    //             while (yf(TrialStress) > yf_TrialStress_prev && subddlambda > MACHINE_EPSILON * 10)
    //             {
    //                 internal_variables.revert_tmp();
    //                 subddlambda *= 0.5 ;
    //                 zeta *= 0.5 ;
    //                 TrialStress(i, j) = TrialStress(i, j) + zeta * dSigma_(i, j) ;
    //                 // Internal variable should be updated with subddlambda.
    //                 // The condition in while-loop will use the internal variables.
    //                 internal_variables.evolve(subddlambda, depsilon, m, TrialStress);
    //             }

    //             dLambda = dLambda + subddlambda ;

    //             yf_TrialStress = yf(TrialStress);

    //             ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

    //             normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

    //             stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
    //             f_relative_error = abs(yf_TrialStress / yf_PredictorStress);
    //             debugrun = true;
    //             if (normResidualStress != normResidualStress || debugrun)  // check for NAN
    //             {
    //                 cout << "=============================================================" << endl;
    //                 cout << "\nIteration # " << iteration_count <<  endl;
    //                 printTensor (" CommitStress       " , CommitStress);
    //                 double p, q, theta;
    //                 std::tie(p, q, theta) = getpqtheta(CommitStress);
    //                 fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
    //                 printTensor (" depsilon           " , depsilon);
    //                 printTensor (" dsigma             " , dsigma);
    //                 printTensor (" PredictorStress    " , PredictorStress);
    //                 std::tie(p, q, theta) = getpqtheta(PredictorStress);
    //                 fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
    //                 cout <<      " yf_val_start        = " << yf_val_start << endl;
    //                 cout <<      " yf_val_end          = " << yf_val_end << endl;
    //                 printTensor (" n                  " , n );
    //                 printTensor (" m                  " , m );
    //                 cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
    //                 cout <<      " denominator         = " << denominator << endl;
    //                 cout <<      " dLambda             = " << dLambda << endl;
    //                 printTensor (" TrialStress        " , TrialStress);
    //                 std::tie(p, q, theta) = getpqtheta(TrialStress);
    //                 fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

    //                 printTensor (" TrialStress_prev   " , TrialStress_prev);
    //                 std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
    //                 fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

    //                 printTensor (" ResidualStress     " , ResidualStress );
    //                 cout <<      " normResidualStress  = " << normResidualStress << endl;
    //                 cout <<      " stress_relative_error      = " << stress_relative_error << endl;
    //                 cout <<      " f_relative_error      = " << f_relative_error << endl;
    //                 cout <<      " iteration_count     = " << iteration_count << endl;
    //                 cout <<      " .........................................................." << endl;
    //                 cout <<      "  Internal variables:" << endl;
    //                 printTensor (" back_stress (alpha)      " , yf.get_alpha() );
    //                 cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
    //                 cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
    //                 cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
    //                 cout << "===============================================END===========" << endl;
    //                 cout << endl << endl;
    //             }

    //         } // while for el-pl this step


    //         // ===============================================
    //         // Start working on the stiffness
    //         // ===============================================
    //         // Inconsistent stiffness
    //         //Update the trial plastic strain. and internal variables
    //         const VoigtVector& n = yf.df_dsigma_ij(TrialStress);
    //         const VoigtVector& m = pf(depsilon, TrialStress);
    //         const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  TrialStress);

    //         denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //         TrialPlastic_Strain(i, j) += dLambda * m(i, j);

    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //         // ===================================================================
    //         // Working on consistent stiffness
    //         // ===================================================================
    //         // Update the m and n:
    //         static VoigtVector nf(3, 3, 0.0);
    //         nf = yf.df_dsigma_ij(TrialStress);
    //         static VoigtVector mf(3, 3, 0.0);
    //         mf = pf(depsilon, TrialStress);
    //         // Construct the Tensor T
    //         static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //         IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //         static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //         dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //         static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //         Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //         // Construct the Tensor H
    //         static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //         dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
    //         static VoigtVector H(3, 3, 0.0);
    //         H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //         static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //         // If cannot inverse T, return directly.
    //         // So "Stiffness" is the inconsistent stiffness tensor.
    //         if ( !inverse4thTensor(Ts, invT) )
    //         {
    //             return 0;
    //         }
    //         static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //         R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //         // ===================================================================
    //         // Construct the consistent stiffness
    //         double denomin{0.0};
    //         double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
    //         denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

    //         Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;
    //         // ===================================================================
    //         // Working on consistent stiffness  END
    //         // ===================================================================

    //     } //else  //Plasticity
    //     internal_variables.commit();
    //     return errorcode;
    // }


    // int Backward_Euler_ddlambda_Subincrement(const VoigtVector &strain_incr, bool debugrun = false)
    // {
    //     using namespace ASDPlasticMaterialGlobals;
    //     int errorcode = 0;

    //     internal_variables.revert();
    //     internal_variables.commit_tmp();
    //     static VoigtVector depsilon(3, 3, 0);
    //     static VoigtVector PredictorStress(3, 3, 0);
    //     const VoigtVector& sigma = CommitStress;
    //     depsilon(i, j) = strain_incr(i, j) ;
    //     TrialStress(i, j) = sigma(i, j);
    //     TrialStrain(i, j) = CommitStrain(i, j) + depsilon(i, j);
    //     TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

    //     VoigtMatrix& Eelastic = et(sigma);
    //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

    //     dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

    //     TrialStress(i, j) +=  dsigma(i, j);
    //     PredictorStress(i, j) = TrialStress(i, j);

    //     double yf_val_start = yf(sigma);
    //     double yf_val_end = yf(PredictorStress);

    //     // printTensor (" CommitStress       " , CommitStress);
    //     // printTensor (" depsilon           " , depsilon);
    //     // printTensor (" dsigma             " , dsigma);
    //     // printTensor (" PredictorStress    " , PredictorStress);
    //     // cout <<      " yf_val_start        = " << yf_val_start << endl;
    //     // cout <<      " yf_val_end          = " << yf_val_end << endl;
    //     // printTensor (" TrialStress        " , TrialStress);

    //     if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //     {
    //         VoigtMatrix& Eelastic = et(TrialStress);
    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
    //     }
    //     else  //Plasticity
    //     {
    //         static VoigtVector ResidualStress(3, 3, 0);
    //         static VoigtVector TrialStress_prev(3, 3, 0);
    //         ResidualStress *= 0;
    //         TrialStress_prev *= 0;
    //         double normResidualStress = -1;
    //         double stress_relative_error = this->stress_relative_tol * 10; //
    //         double dLambda = 0;
    //         double denominator = 0;
    //         int iteration_count = 0;

    //         double yf_PredictorStress = yf(PredictorStress);
    //         double yf_TrialStress = yf(TrialStress);
    //         // vonMises does NOT enter this part.
    //         // DruckerPrager requries this part.
    //         if (yf.hasCorner() && yf.in_Apex(PredictorStress))
    //         {
    //             static VoigtVector small_stress(3, 3, 0.0);
    //             small_stress *= 0;
    //             // The small value 50*Pa refers to the lowest confinement test:
    //             // http://science.nasa.gov/science-news/science-at-nasa/1998/msad27may98_2/
    //             double DP_k = yf.get_k();
    //             double DP_p = 50 ;
    //             // To make it on the yield surface, the q is equal to k*p.
    //             double DP_q = DP_k * DP_p ;
    //             // Assume the triaxial conditions sigma_2 = sigma_3.
    //             small_stress(0, 0) = DP_p + 2. / 3.0 * DP_q;
    //             small_stress(1, 1) = DP_p - 1. / 3.0 * DP_q;
    //             small_stress(2, 2) = DP_p - 1. / 3.0 * DP_q;

    //             // (1) Update the trial stress
    //             TrialStress(i, j) = small_stress(i, j);
    //             // (2) Update the trial plastic strain
    //             const VoigtVector& n = yf.df_dsigma_ij(small_stress);
    //             const VoigtVector& m = pf(depsilon, small_stress);
    //             const double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  small_stress);
    //             double denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;
    //             double dLambda = yf_PredictorStress / denominator;
    //             TrialPlastic_Strain(i, j) += dLambda * m(i, j);
    //             // (3) Update the internal variables
    //             internal_variables.evolve(dLambda, depsilon, m, TrialStress);
    //             internal_variables.commit_tmp();
    //             // (4) Update the stiffness.
    //             Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m(p, q)) * (n(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //             // // =========================================================
    //             // // Working on consistent stiffness
    //             // // ========================================================
    //             // // Construct the Tensor T
    //             static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //             IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //             static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //             dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //             static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //             Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //             // // Construct the Tensor H
    //             static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //             dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m, TrialStress);
    //             static VoigtVector H(3, 3, 0.0);
    //             H(k, l) = m(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //             static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //             // // If cannot inverse T, return directly with the inconsistent stiffness tensor.
    //             if ( !inverse4thTensor(Ts, invT) )
    //             {
    //                 return 0;
    //             }
    //             static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //             R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //             // // =======================================================
    //             // // Construct the consistent stiffness
    //             double denomin{0.0};
    //             denomin = n(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star ;
    //             Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (n(r, s) * R(r, s, k, l))) / denomin;

    //             // // =====================================================
    //             // // Working on consistent stiffness  END
    //             // // =====================================================

    //             // skip the following iterations
    //             return 0;
    //         }//END yield surface APEX check

    //         // ==================================================================
    //         // From the start point(elastic predictor) to the first plastic corrector.
    //         // ==================================================================
    //         yf_PredictorStress = yf(PredictorStress);
    //         const VoigtVector& n = yf.df_dsigma_ij(PredictorStress);
    //         const VoigtVector& m = pf(depsilon, PredictorStress);
    //         double xi_star_h_star = yf.xi_star_h_star( depsilon, m,  PredictorStress);
    //         denominator = n(p, q) * Eelastic(p, q, r, s) * m(r, s) - xi_star_h_star;

    //         dLambda =  yf_PredictorStress / denominator;

    //         TrialStress_prev(i, j) = TrialStress(i, j);

    //         TrialStress(i, j) = PredictorStress(i, j) - dLambda * Eelastic(i, j, k, l) * m(k, l);
    //         internal_variables.evolve(dLambda, depsilon, m, TrialStress);
    //         // // ============================================================
    //         // // Line search in dlambda
    //         // // ============================================================
    //         double subdLambda = dLambda;
    //         while (yf(TrialStress) > yf_PredictorStress && subdLambda > MACHINE_EPSILON * 10)
    //         {
    //             internal_variables.revert_tmp();
    //             subdLambda *= 0.5;
    //             TrialStress(i, j) = PredictorStress(i, j) - subdLambda * Eelastic(i, j, k, l) * m(k, l);
    //             // Internal variable should be updated with subdLambda.
    //             // The condition in while-loop will use the internal variables.
    //             internal_variables.evolve(subdLambda, depsilon, m, TrialStress);
    //         }
    //         // ============================================================
    //         dLambda = subdLambda;
    //         ResidualStress(i, j) =  TrialStress(i, j) - TrialStress_prev(i, j);

    //         double f_relative_error = this->f_relative_tol * 10;

    //         // Iteration over the ddlambda
    //         while ((stress_relative_error > this-> stress_relative_tol ||
    //                 f_relative_error > this->f_relative_tol) &&
    //                 iteration_count < this->n_max_iterations
    //               )
    //         {
    //             iteration_count++;

    //             // ==========================================
    //             // Works to calculate the ddlambda
    //             // ==========================================
    //             // Update the m and n:
    //             static VoigtVector n_new(3, 3, 0.0);
    //             n_new = yf.df_dsigma_ij(TrialStress);
    //             static VoigtVector m_new(3, 3, 0.0);
    //             m_new = pf(depsilon, TrialStress);
    //             // Construct the Tensor T
    //             static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0.0);
    //             IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //             static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //             dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //             static VoigtMatrix T_new(3, 3, 3, 3, 0.0);
    //             T_new(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //             static VoigtMatrix invTnew(3, 3, 3, 3, 0.0);
    //             if ( !inverse4thTensor(T_new, invTnew) )
    //             {
    //                 // cout<<"Singularity. Cannot inverse tensor T! "<<endl;
    //                 return 0;
    //             }

    //             static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //             dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, m_new, TrialStress);
    //             static VoigtVector H_new(3, 3, 0.0);
    //             H_new(k, l) = m_new(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //             double const& xi_star_h_star_new = yf.xi_star_h_star( depsilon, m_new,  TrialStress);

    //             double denomina_ddlambda = 0.0;
    //             static VoigtVector E_times_H (3, 3, 0.0);
    //             Eelastic = et(TrialStress);
    //             E_times_H(i, j) = Eelastic(i, j, k, l) * H_new(k, l);
    //             denomina_ddlambda = n_new(r, s) * (E_times_H(p, q) * invTnew(p, q, r, s)) - xi_star_h_star_new ;
    //             double numerator_ddlambda = 0.0;
    //             numerator_ddlambda = yf(TrialStress) - n_new(r, s) * (ResidualStress(i, j) * invTnew(i, j, r, s)) ;
    //             double ddlambda =   numerator_ddlambda / denomina_ddlambda ;

    //             static VoigtVector dSigma_(3, 3, 0.0);
    //             dSigma_(r, s) = -
    //                             (
    //                                 ResidualStress(i, j) + ddlambda * Eelastic(i, j, p, q) * H_new(p, q)
    //                             )
    //                             * invTnew(i, j, r, s) ;

    //             TrialStress_prev(i, j) = TrialStress(i, j) ;
    //             TrialStress(i, j) =  TrialStress(i, j) + dSigma_(i, j) ;
    //             double yf_TrialStress_prev = yf(TrialStress_prev) ;
    //             internal_variables.evolve(ddlambda, depsilon, m_new, TrialStress);

    //             // =======Line Search ddlambda (double delta ) =======
    //             double subddlambda = ddlambda;
    //             double zeta = 1.0;
    //             while (yf(TrialStress) > yf_TrialStress_prev && subddlambda > MACHINE_EPSILON * 10)
    //             {
    //                 internal_variables.revert_tmp();
    //                 subddlambda *= 0.5 ;
    //                 zeta *= 0.5 ;
    //                 TrialStress(i, j) = TrialStress(i, j) + zeta * dSigma_(i, j) ;
    //                 // Internal variable should be updated with subddlambda.
    //                 // The condition in while-loop will use the internal variables.
    //                 internal_variables.evolve(subddlambda, depsilon, m, TrialStress);
    //             }

    //             dLambda = dLambda + subddlambda ;
    //             yf_TrialStress = yf(TrialStress);

    //             ResidualStress(i, j) = TrialStress(i, j) - TrialStress_prev(i, j);

    //             normResidualStress = sqrt(ResidualStress(i, j) * ResidualStress(i, j));

    //             stress_relative_error = normResidualStress / sqrt(TrialStress(i, j) * TrialStress(i, j));
    //             f_relative_error = abs(yf_TrialStress / yf_PredictorStress);
    //             // debugrun = true;
    //             if (normResidualStress != normResidualStress || debugrun)  // check for NAN
    //             {
    //                 cout << "=============================================================" << endl;
    //                 cout << "\nIteration # " << iteration_count <<  endl;
    //                 printTensor (" CommitStress       " , CommitStress);
    //                 double p, q, theta;
    //                 std::tie(p, q, theta) = getpqtheta(CommitStress);
    //                 fprintf(stderr, "CommitStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "CommitStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "CommitStress_theta= %16.8f\n"  , theta );
    //                 printTensor (" depsilon           " , depsilon);
    //                 printTensor (" dsigma             " , dsigma);
    //                 printTensor (" PredictorStress    " , PredictorStress);
    //                 std::tie(p, q, theta) = getpqtheta(PredictorStress);
    //                 fprintf(stderr, "PredictorStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "PredictorStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "PredictorStress_theta= %16.8f\n"  , theta );
    //                 cout <<      " yf_val_start        = " << yf_val_start << endl;
    //                 cout <<      " yf_val_end          = " << yf_val_end << endl;
    //                 printTensor (" n                  " , n );
    //                 printTensor (" m                  " , m );
    //                 cout <<      " xi_star_h_star      = " << xi_star_h_star << endl;
    //                 cout <<      " denominator         = " << denominator << endl;
    //                 cout <<      " dLambda             = " << dLambda << endl;
    //                 printTensor (" TrialStress        " , TrialStress);
    //                 std::tie(p, q, theta) = getpqtheta(TrialStress);
    //                 fprintf(stderr, "TrialStress_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "TrialStress_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "TrialStress_theta= %16.8f\n"  , theta );

    //                 printTensor (" TrialStress_prev   " , TrialStress_prev);
    //                 std::tie(p, q, theta) = getpqtheta(TrialStress_prev);
    //                 fprintf(stderr, "TrialStress_prev_p    = %16.8f\n"  , p );
    //                 fprintf(stderr, "TrialStress_prev_q    = %16.8f\n"  , q );
    //                 fprintf(stderr, "TrialStress_prev_theta= %16.8f\n"  , theta );

    //                 printTensor (" ResidualStress     " , ResidualStress );
    //                 cout <<      " normResidualStress  = " << normResidualStress << endl;
    //                 cout <<      " stress_relative_error      = " << stress_relative_error << endl;
    //                 cout <<      " f_relative_error      = " << f_relative_error << endl;
    //                 cout <<      " iteration_count     = " << iteration_count << endl;
    //                 cout <<      " .........................................................." << endl;
    //                 cout <<      "  Internal variables:" << endl;
    //                 printTensor (" back_stress (alpha)      " , yf.get_alpha() );
    //                 cout <<      " DP_k or VM_radius      " << yf.get_k() << endl;
    //                 cout << "back_stress(0,1) = " << (yf.get_alpha())(0, 1) << endl ;
    //                 cout << "PredictorStress(0,1) = " << PredictorStress(0, 1) << endl ;
    //                 cout << "===============================================END===========" << endl;
    //                 cout << endl << endl;
    //             }

    //         } // while for el-pl this step


    //         // ===============================================
    //         // Start working on the stiffness
    //         // ===============================================
    //         // Inconsistent stiffness
    //         //Update the trial plastic strain. and internal variables
    //         const VoigtVector& n_ = yf.df_dsigma_ij(TrialStress);
    //         const VoigtVector& m_ = pf(depsilon, TrialStress);
    //         const double xi_star_h_star_ = yf.xi_star_h_star( depsilon, m_,  TrialStress);

    //         denominator = n_(p, q) * Eelastic(p, q, r, s) * m_(r, s) - xi_star_h_star_;

    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * m_(p, q)) * (n_(r, s) * Eelastic(r, s, k, l) ) / denominator;

    //         // ===================================================================
    //         // Working on consistent stiffness
    //         // ===================================================================
    //         // Update the m and n:
    //         static VoigtVector nf(3, 3, 0.0);
    //         nf = yf.df_dsigma_ij(TrialStress);
    //         static VoigtVector mf(3, 3, 0.0);
    //         mf = pf(depsilon, TrialStress);
    //         // Construct the Tensor T
    //         static VoigtMatrix IdentityTensor4(3, 3, 3, 3, 0);
    //         IdentityTensor4(i, j, k, l) = kronecker_delta(i, j) * kronecker_delta(k, l);

    //         static VoigtMatrix dm_dsigma(3, 3, 3, 3, 0.0);
    //         dm_dsigma = pf.dm_over_dsigma(TrialStress);
    //         static VoigtMatrix Ts(3, 3, 3, 3, 0.0);
    //         Ts(i, j, k, l) = IdentityTensor4(i, k, j, l) + dLambda * Eelastic(i, j, p, q) * dm_dsigma(p, q, k, l);

    //         // Construct the Tensor H
    //         static VoigtVector dm_dq_star_h_star(3, 3, 0.0);
    //         dm_dq_star_h_star = pf.dm_over_dq_start_h_star(depsilon, mf, TrialStress);
    //         static VoigtVector H(3, 3, 0.0);
    //         H(k, l) = mf(k, l) + dLambda * dm_dq_star_h_star(k, l);

    //         static VoigtMatrix invT(3, 3, 3, 3, 0.0);
    //         // If cannot inverse T, return directly.
    //         // So "Stiffness" is the inconsistent stiffness tensor.
    //         if ( !inverse4thTensor(Ts, invT) )
    //         {
    //             return 0;
    //         }
    //         static VoigtMatrix R(3, 3, 3, 3, 0.0);
    //         R(p, q, r, s) = invT(k, l, p, q) * Eelastic(k, l, r, s);

    //         // ===================================================================
    //         // Construct the consistent stiffness
    //         double denomin{0.0};
    //         double const& xi_star_h_star_f = yf.xi_star_h_star( depsilon, mf,  TrialStress);
    //         denomin = nf(p, q) * R(p, q, r, s) * H(r, s) - xi_star_h_star_f ;

    //         Stiffness(i, j, k, l) = R(i, j, k, l) - ((R(i, j, p, q) * H(p, q)) * (nf(r, s) * R(r, s, k, l))) / denomin;
    //         // ===================================================================
    //         // Working on consistent stiffness  END
    //         // ===================================================================

    //     } //else  //Plasticity
    //     internal_variables.commit();
    //     return errorcode;
    // }


    // // The algorithm below Modified_Euler_Error_Control was implemented by Jose. 18Sep2016
    // int Modified_Euler_Error_Control(const VoigtVector &strain_incr)
    // {
    //     using namespace ASDPlasticMaterialGlobals;  // Brings indexes i,j,k,l,m,n,p,q into current scope

    //     int errorcode = 0;

    //     static VoigtVector depsilon(3, 3, 0);
    //     static VoigtVector start_stress(3, 3, 0.0);
    //     static VoigtVector end_stress(3, 3, 0.0);
    //     // const VoigtVector& sigma = CommitStress;

    //     //Zero out static variables to ensure nothing is there.
    //     depsilon *= 0;
    //     dsigma *= 0;
    //     start_stress *= 0;
    //     end_stress *= 0;
    //     intersection_stress *= 0;
    //     intersection_strain *= 0;

    //     depsilon(i, j) = strain_incr(i, j);         //This is needed because some function don't take const references

    //     //This makes sure we start with last committed values for internal_variables
    //     internal_variables.revert();
    //     TrialStress(i, j) = CommitStress(i, j);
    //     TrialStrain(i, j) = CommitStrain(i, j);
    //     TrialPlastic_Strain(i, j) = CommitPlastic_Strain(i, j);

    //     //Compute elasticity for current (committed) stress level and assign stiffness to be this.
    //     VoigtMatrix& Eelastic = et(TrialStress);
    //     Stiffness(i, j, k, l) = Eelastic(i, j, k, l);

    //     //Compute the elastic stress increment
    //     dsigma(i, j) = Eelastic(i, j, k, l) * depsilon(k, l);

    //     //Increment the trial stresses and strains
    //     TrialStress(i, j) += dsigma(i, j);
    //     TrialStrain(i, j) += depsilon(i, j);

    //     //Compute the values of yield function before and after elastic increment
    //     double yf_val_start = yf(CommitStress);
    //     double yf_val_end = yf(TrialStress);

    //     //These are the start and end stresses for the brent algorithm which
    //     //finds the yield surface intersection stress (taken from old NewTemplate3dEP classes)
    //     start_stress(i, j) = CommitStress(i, j);
    //     end_stress(i, j) = TrialStress(i, j);

    //     //Set the intersection stress at the start
    //     intersection_stress(i, j) = start_stress(i, j);

    //     //This pre_integration_callback_(...) provides an opportunity to exit before integration
    //     //in case its needed. For example when some materials go in tension.
    //     bool returns = false;
    //     int retval = pre_integration_callback_(depsilon, dsigma, TrialStress, Stiffness, yf_val_start, yf_val_end,  returns);
    //     if (returns)
    //     {
    //         return retval;
    //     }


    //     if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //     {
    //         //For an elastic step, compute the elasticity at the incremented stress,
    //         //for use in next step.
    //         VoigtMatrix& Eelastic = et(TrialStress);
    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l);
    //     }
    //     else  //Plasticity
    //     {
    //         //Some temporary variables used herein
    //         static VoigtVector sub_depsilon_elpl(3, 3, 0);
    //         static VoigtVector dsigma1(3, 3, 0);
    //         static VoigtVector dsigma2(3, 3, 0);
    //         static VoigtVector sigma_error(3, 3, 0);
    //         static VoigtVector m(3, 3, 0);

    //         sub_depsilon_elpl *= 0;
    //         dsigma1 *= 0;
    //         dsigma2 *= 0;
    //         sigma_error *= 0;
    //         m *= 0;

    //         //If the material was elastic but now is plastic, compute the intersection point
    //         //ie. the point at which the trial stress path intersects the yield surface.
    //         //Also compute the intersection strain increment and compute the plastic strain increment
    //         //from that.
    //         depsilon_elpl(i, j) = depsilon(i, j);
    //         if (yf_val_start < 0)
    //         {
    //             double intersection_factor = zbrentstress( start_stress, end_stress, 0.0, 1.0, TOLERANCE1 );
    //             intersection_stress(i, j) = start_stress(i, j) * (1 - intersection_factor) + end_stress(i, j) * intersection_factor;
    //             depsilon_elpl(i, j) = (1 - intersection_factor) * depsilon(i, j);
    //         }

    //         //Error controlled iteration after Sloan
    //         int count = 0;      //Number of iterations done
    //         double Tau = 0;     //Integration pseudo-time
    //         double dT = 1;      //Current integration pseudo-time step
    //         double q = 0;       //A factor tu multiply dT with for next step

    //         //Set both the start stress and the trial stress at the intersection stress
    //         //if trial fails, will return to start_stress.
    //         start_stress(i, j) = intersection_stress(i, j);
    //         TrialStress(i, j) = intersection_stress(i, j);

    //         internal_variables.commit_tmp(); // Stores temporary internal variables in static storage space

    //         while (Tau < 1)
    //         {
    //             //Compute current plastic strain sub-increment.
    //             sub_depsilon_elpl(i, j) = dT * depsilon_elpl(i, j);

    //             //Compute elasticity at current trial stress and determine the first elastic stress
    //             //increment based on this.
    //             Eelastic = et(TrialStress);
    //             dsigma1(i, j)  = Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);

    //             printTensor("  dsigma_elastic", dsigma1);

    //             //Compute normal to YF (n) and Plastic Flow direction (m)
    //             const VoigtVector& n1 = yf.df_dsigma_ij(TrialStress);
    //             const VoigtVector& m1 = pf(sub_depsilon_elpl, TrialStress);
    //             printTensor("  n1", n1);
    //             printTensor("  m1", m1);

    //             //Compute hardening parameter
    //             double xi_star_h_star1 = yf.xi_star_h_star( sub_depsilon_elpl, m1,  TrialStress);

    //             //Compute denominator of the plastic multiplier
    //             double den1 = n1(p, q) * Eelastic(p, q, r, s) * m1(r, s) - xi_star_h_star1;

    //             //Compute the plastic multiplier
    //             if (den1 == 0)
    //             {
    //                 cerr << "CEP - den = 0\n";
    //                 return -1;
    //             }
    //             double dLambda1 =  n1(i, j) * Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
    //             dLambda1 /= den1;

    //             //Correct the trial stress increment
    //             dsigma1(i, j) +=  - dLambda1 * Eelastic(i, j, k, l) * m1(k, l);

    //             //Now evolve material state
    //             internal_variables.evolve(dLambda1, sub_depsilon_elpl, m1, TrialStress);
    //             TrialStress(i, j) += dsigma1(i, j);

    //             //Now evaluate everything at the new trial stress
    //             Eelastic = et(TrialStress);
    //             dsigma2(i, j)  = Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
    //             const VoigtVector& n2 = yf.df_dsigma_ij(TrialStress);
    //             const VoigtVector& m2 = pf(sub_depsilon_elpl, TrialStress);
    //             double xi_star_h_star2 = yf.xi_star_h_star( sub_depsilon_elpl, m2,  TrialStress);
    //             double den2 = n2(p, q) * Eelastic(p, q, r, s) * m2(r, s) - xi_star_h_star2;
    //             if (den2 == 0)
    //             {
    //                 cerr << "CEP - den = 0\n";
    //                 return -1;
    //             }
    //             double dLambda2 =  n2(i, j) * Eelastic(i, j, k, l) * sub_depsilon_elpl(k, l);
    //             dLambda2 /= den2;
    //             dsigma2(i, j) += - dLambda2 * Eelastic(i, j, k, l) * m2(k, l);


    //             //Compute a better estimate of TrialStress by averaging the two increments
    //             TrialStress(i, j) = start_stress(i, j);
    //             sigma_error(i, j) = (dsigma2(i, j) - dsigma1(i, j)) / 2;
    //             TrialStress(i, j) += (dsigma1(i, j) + dsigma2(i, j)) / 2;

    //             //Compute the relative error based on the stress error and the trial stress
    //             double Relative_Error = sqrt(sigma_error(i, j) * sigma_error(i, j)) / sqrt(TrialStress(i, j) * TrialStress(i, j));

    //             double dLambda = 0;

    //             cout << "  den1 = " << den1 << ", dLambda1 = " << dLambda1 << endl;
    //             cout << "  den2 = " << den2 << ", dLambda2 = " << dLambda2 << endl;
    //             printTensor("  depsilon_elpl", depsilon_elpl);
    //             printTensor("  sub_depsilon_elpl", sub_depsilon_elpl);
    //             printTensor("  dsigma1", dsigma1);
    //             printTensor("  dsigma2", dsigma2);
    //             cout << "  dT = " << dT << " (T = " << Tau << " q = " << q << ")\n";
    //             cout << "  RE = " << Relative_Error << " (tol=" << this-> stress_relative_tol << ") ->";

    //             //Revert the variables to the start (use the temporary storage)
    //             internal_variables.revert_tmp();
    //             TrialStress(i, j) = start_stress(i, j);
    //             if (Relative_Error < this-> stress_relative_tol) // Accept this step
    //             {
    //                 cout << "Accept!\n";
    //                 Tau += dT;
    //                 q = fmin(0.8 * sqrt(this-> stress_relative_tol / Relative_Error), 2.0);
    //                 dT = q * dT;
    //                 dT = fmin(dT, 1.0 - Tau);

    //                 dLambda = dLambda1;//(dLambda1 + dLambda2) / 2;
    //                 internal_variables.evolve(dLambda, sub_depsilon_elpl, m1, TrialStress);
    //                 internal_variables.commit_tmp();
    //                 TrialStress(i, j) += dsigma1(i, j);//(dsigma1(i, j) + dsigma2(i, j)) / 2;


    //                 start_stress(i, j) = TrialStress(i, j);

    //                 //Update the trial plastic strain.
    //                 TrialPlastic_Strain(i, j) += dLambda * m1(i, j);
    //             }
    //             else //Reject step, take a smaller one!
    //             {
    //                 cout << "Reject!\n";
    //                 // cout << "Reject! R = " <<  Relative_Error << endl;
    //                 q = fmax(0.8 * sqrt(this-> stress_relative_tol / Relative_Error), 0.1);
    //                 dT = q * dT;
    //             }

    //             count += 1;
    //             if (count > this-> n_max_iterations)
    //             {
    //                 cerr << "ASDPlasticMaterial -- Modified euler failed to converge after " << count << " iterations with : \n"
    //                      << "  Relative_Error = " << Relative_Error  << endl
    //                      << "  Relative Tolerance = " << this-> stress_relative_tol  << endl
    //                      << "  q = " << q  << endl
    //                      << "  dT = " << dT  << endl
    //                      << "  dLambda = " << dLambda  << endl
    //                      << "  T = " << Tau  << endl;
    //                 printTensor("TrialStress = " , TrialStress);
    //                 printTensor("CommitStress = " , CommitStress);
    //                 printTensor("depsilon = " , depsilon);
    //                 printTensor("depsilon_elpl = " , depsilon_elpl);
    //                 printTensor("sub_depsilon_elpl = " , sub_depsilon_elpl);
    //                 printTensor("dsigma1 = " , dsigma1);
    //                 printTensor("dsigma2 = " , dsigma2);
    //                 printTensor("sigma_error = " , sigma_error);
    //                 printTensor("start_stress = " , start_stress);
    //                 printTensor("intersection_stress = " , intersection_stress);
    //                 return -1;
    //             }
    //         }

    //         //Stress correction back to the yield surface
    //         double yc = yf(TrialStress);
    //         int ys_correction_count = 0;
    //         printTensor("  TrialStress before correction", TrialStress);
    //         while (yc > this-> f_relative_tol && ys_correction_count < this-> n_max_iterations)
    //         {
    //             cout << "    +++ yf() = " << yc << endl;
    //             const VoigtVector& n = yf.df_dsigma_ij(TrialStress);
    //             const VoigtVector& m = pf(depsilon_elpl, TrialStress);
    //             double den = n(p, q) * Eelastic(p, q, r, s) * m(r, s);
    //             double dLambda_correction = yc / den;
    //             TrialStress(i, j) = TrialStress(i, j) - dLambda_correction * Eelastic(i, j, k, l) * m(k, l);
    //             printTensor("    +++ TrialStress correction", TrialStress);
    //             yc = yf(TrialStress);
    //             ys_correction_count++;
    //         }
    //         printTensor("  TrialStress after correction", TrialStress);


    //         // Once done compute the stiffness for the next step
    //         Eelastic = et(TrialStress);
    //         const VoigtVector& nf = yf.df_dsigma_ij(TrialStress);
    //         const VoigtVector& mf = pf(depsilon_elpl, TrialStress);
    //         double xi_star_h_star_f = yf.xi_star_h_star( depsilon_elpl, mf,  TrialStress);
    //         double denf = nf(p, q) * Eelastic(p, q, r, s) * mf(r, s) - xi_star_h_star_f;
    //         Stiffness(i, j, k, l) = Eelastic(i, j, k, l) - (Eelastic(i, j, p, q) * mf(p, q)) * (nf(r, s) * Eelastic(r, s, k, l) ) / denf;



    //         //Detect a NAN result and print diagnostics if one is found
    //         // double norm_trial_stress = TrialStress(i, j) * TrialStress(i, j);
    //         if (true)//norm_trial_stress != norm_trial_stress) //check for nan
    //         {
    //             // cout << "Nan Detected in Modified_Euler !\n";
    //             printTensor("TrialStress = " , TrialStress);
    //             printTensor("CommitStress = " , CommitStress);
    //             printTensor("depsilon = " , depsilon);
    //             printTensor("intersection_stress = " , intersection_stress);
    //             printTensor4("Eelastic = " , Eelastic);
    //             printTensor4("Stiffness = " , Stiffness);
    //             cout << "yf_val_start = " << yf_val_start << endl;
    //             cout << "yf_val_end = " << yf_val_end << endl;
    //             printTensor("n = " , nf );
    //             printTensor("m = " , mf );
    //             cout << "xi_star_h_star_f  = " << xi_star_h_star_f << endl;
    //             cout << "denf = " << denf << endl;
    //             // cout << "dLambda = " << dLambda << endl;
    //             // cout << "count = " << count << endl;
    //             // cout << "ys_correction_count = " << ys_correction_count << endl;

    //             // errorcode = -1;
    //         }

    //     }

    //     return errorcode;
    // }

//     template <typename U = T>
//     typename std::enable_if < !supports_pre_integration_callback<U>::value, int >::type
// // typename std::enable_if < !std::is_base_of<defines_pre_integration_callback, U>::value, int >::type
//     pre_integration_callback_(const VoigtVector &depsilon, const VoigtVector &dsigma,  const VoigtVector &TrialStress, const VoigtMatrix &Stiffness, double yf1, double yf2, bool & returns)
//     {
//         // cout << "pre_integration_callback_ disabled\n";
//         returns = false;
//         return 0;
//     }

//     template <typename U = T>
//     typename std::enable_if<supports_pre_integration_callback<U>::value, int>::type
// // typename std::enable_if<std::is_base_of<defines_pre_integration_callback, U>::value, int>::type
//     pre_integration_callback_(const VoigtVector &depsilon, const VoigtVector &dsigma, const VoigtVector &TrialStress, const VoigtMatrix &Stiffness, double yf1, double yf2, bool & returns)
//     {
//         // cout << "pre_integration_callback_ enabled\n";
//         return static_cast<U*>(this)->pre_integration_callback(depsilon, dsigma, TrialStress, Stiffness,  yf1,  yf2, returns);
//     }

private:
    // Routine used by yield_surface_cross to find the stresstensor at cross point
    //================================================================================
    // double zbrentstress(const VoigtVector& start_stress, const VoigtVector& end_stress, double x1, double x2, double tol) const {
    //     using namespace ASDPlasticMaterialGlobals;

    //     // Constants
    //     const double EPS_MULTIPLE = 2.0;
    //     const double TOL_MULTIPLE = 0.5;
    //     const double EPS = numeric_limits<double>::epsilon();

    //     // Lambda functions
    //     auto calculateSigma = [](const VoigtVector & start_stress, const VoigtVector & end_stress, double multiplier) {
    //         return start_stress * (1 - multiplier) + end_stress * multiplier;
    //     };

    //     auto calculateYf = [this](const VoigtVector & sigma) {
    //         return yf(sigma, iv_storage, parameters_storage);
    //     };

    //     // Rest of the function
    //     double a = x1;
    //     double b = x2;
    //     double c = 0.0;
    //     double d = 0.0;
    //     double e = 0.0;
    //     double fc = 0.0;
    //     double p = 0.0;
    //     double q = 0.0;
    //     double r = 0.0;
    //     double s = 0.0;
    //     double tol1 = 0.0;
    //     double xm = 0.0;

    //     VoigtVector sigma_a = calculateSigma(start_stress, end_stress, a);
    //     VoigtVector sigma_b = calculateSigma(start_stress, end_stress, b);

    //     double fa = calculateYf(sigma_a);
    //     double fb = calculateYf(sigma_b);

    //     if ( (fb * fa) > 0.0) {
    //         throw std::runtime_error("Root must be bracketed in ZBRENTstress");
    //     }

    //     fc = fb;

    //     for ( int iter = 1; iter <= ASDPlasticMaterial_MAXITER_BRENT; iter++ ) {
    //         if ( (fb * fc) > 0.0) {
    //             c = a;
    //             fc = fa;
    //             e = d = b - a;
    //         }

    //         if ( fabs(fc) < fabs(fb) ) {
    //             std::swap(a, b);
    //             std::swap(c, a);
    //             std::swap(fa, fb);
    //             std::swap(fc, fa);
    //         }

    //         tol1 = EPS_MULTIPLE * EPS * fabs(b) + TOL_MULTIPLE * tol;
    //         xm = 0.5 * (c - b);

    //         if ( fabs(xm) <= tol1 || fb == 0.0 ) {
    //             return b;
    //         }

    //         if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ) {
    //             s = fb / fa;

    //             if (a == c) {
    //                 p = 2.0 * xm * s;
    //                 q = 1.0 - s;
    //             } else {
    //                 q = fa / fc;
    //                 r = fb / fc;
    //                 p = s * ( 2.0 * xm * q * (q - r) - (b - a) * (r - 1.0) );
    //                 q = (q - 1.0) * (r - 1.0) * (s - 1.0);
    //             }

    //             if (p > 0.0) {
    //                 q = -q;
    //             }

    //             p = fabs(p);
    //             double min1 = 3.0 * xm * q - fabs(tol1 * q);
    //             double min2 = fabs(e * q);

    //             if (2.0 * p < std::min(min1, min2)) {
    //                 e = d;
    //                 d = p / q;
    //             } else {
    //                 d = xm;
    //                 e = d;
    //             }
    //         } else {
    //             d = xm;
    //             e = d;
    //         }

    //         a = b;
    //         fa = fb;

    //         if (fabs(d) > tol1) {
    //             b += d;
    //         } else {
    //             b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    //         }

    //         sigma_b = calculateSigma(start_stress, end_stress, b);
    //         fb = calculateYf(sigma_b);
    //     }

    //     throw std::runtime_error("Maximum iterations reached without convergence in ZBRENTstress");
    // }

    //More robust brent
    double zbrentstress(const VoigtVector & start_stress, const VoigtVector & end_stress, double x1, double x2, double tol) const
    {
        using namespace ASDPlasticMaterialGlobals;

        // Constants
        const double MAX_BRACKET_EXPANSION = 100.0; // Maximum factor to expand the bracket
        const double BRACKET_EXPANSION_FACTOR = 1.6; // Factor to expand the bracket each iteration
        const double EPS_MULTIPLE = 2.0;
        const double TOL_MULTIPLE = 0.5;
        const double EPS = std::numeric_limits<double>::epsilon();

        // Lambda functions
        auto calculateSigma = [](const VoigtVector & start, const VoigtVector & end, double multiplier) {
            return start * (1 - multiplier) + end * multiplier;
        };

        auto calculateYf = [this](const VoigtVector & sigma) {
            return yf(sigma, iv_storage, parameters_storage);
        };

        // Bracket Adjustment: Expand the bracket if necessary
        double a = x1;
        double b = x2;
        double c = x2;
        double d = 0;
        double e = 0.0;
        double fa = calculateYf(calculateSigma(start_stress, end_stress, a));
        double fb = calculateYf(calculateSigma(start_stress, end_stress, b));
        double fc = fb;


        for (double factor = 1.0; factor <= MAX_BRACKET_EXPANSION; factor *= BRACKET_EXPANSION_FACTOR) {
            if ((fb * fa) <= 0.0) {
                break; // Valid bracket found
            }
            if (std::abs(fa) < std::abs(fb)) {
                a -= (b - a) * factor;
                fa = calculateYf(calculateSigma(start_stress, end_stress, a));
            } else {
                b += (b - a) * factor;
                fb = calculateYf(calculateSigma(start_stress, end_stress, b));
            }
        }

        // Brent's Method Implementation
        if ((fb * fa) > 0.0) {
            throw std::runtime_error("Unable to find a valid bracket in zbrentstress");
        }

        for (int iter = 1; iter <= ASDPlasticMaterial_MAXITER_BRENT; iter++) {
            if ((fb * fc) > 0.0) {
                c = a;   // Rename a, b, c and adjust bounding interval d
                fc = fa;
                e = d = b - a;
            }

            if (std::abs(fc) < std::abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            double tol1 = EPS_MULTIPLE * EPS * std::abs(b) + TOL_MULTIPLE * tol;
            double xm = 0.5 * (c - b);

            if (std::abs(xm) <= tol1 || fb == 0.0) {
                return b;  // Convergence
            }

            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
                // Attempt inverse quadratic interpolation
                double s = fb / fa;
                double p, q;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) q = -q;
                p = std::abs(p);
                double min1 = 3.0 * xm * q - std::abs(tol1 * q);
                double min2 = std::abs(e * q);

                if (2.0 * p < std::min(min1, min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;  // Bounds decreasing too slowly, use bisection
                e = d;
            }

            a = b;
            fa = fb;

            if (std::abs(d) > tol1) {
                b += d;
            } else {
                b += (xm > 0.0 ? std::abs(tol1) : -std::abs(tol1));
            }

            fb = calculateYf(calculateSigma(start_stress, end_stress, b));

            // Check for NaN (Not a Number) values
            if (std::isnan(fb)) {
                throw std::runtime_error("zbrentstress: NaN encountered in function evaluation.");
            }
        }

        throw std::runtime_error("Maximum iterations reached without convergence in zbrentstress");
    }


protected:

    VoigtVector TrialStrain;
    VoigtVector TrialStress;
    VoigtVector TrialPlastic_Strain;

    VoigtVector CommitStress;
    VoigtVector CommitStrain;
    VoigtVector CommitPlastic_Strain;


    YieldFunctionType yf;
    ElasticityType    et;
    PlasticFlowType   pf;

    iv_storage_t iv_storage;
    parameters_storage_t parameters_storage;


protected:

    static std::map<int, ASDPlasticMaterial_Constitutive_Integration_Method> INT_OPT_constitutive_integration_method;     //
    static std::map<int, ASDPlasticMaterial_Tangent_Operator_Type> INT_OPT_tangent_operator_type;     //
    static std::map<int, double> INT_OPT_f_relative_tol;
    static std::map<int, double> INT_OPT_stress_relative_tol;
    static std::map<int, int> INT_OPT_n_max_iterations;
    static std::map<int, int> INT_OPT_return_to_yield_surface;

    bool first_step;

    static VoigtVector dsigma;
    static VoigtVector depsilon_elpl;    //Elastoplastic strain increment : For a strain increment that causes first yield, the step is divided into an elastic one (until yield) and an elastoplastic one.
    static VoigtVector intersection_stress;
    static VoigtVector intersection_strain;
    static VoigtMatrix Stiffness;


};

template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial_Constitutive_Integration_Method> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_constitutive_integration_method;
template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial_Tangent_Operator_Type> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_tangent_operator_type;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_f_relative_tol;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_stress_relative_tol;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_n_max_iterations;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial< E,  Y,  P,  tag>::INT_OPT_return_to_yield_surface;


template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial< E,  Y,  P,  tag>::dsigma;

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial< E,  Y,  P,  tag>::depsilon_elpl;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial< E,  Y,  P,  tag >::intersection_stress;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial< E,  Y,  P,  tag>::intersection_strain;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtMatrix ASDPlasticMaterial< E,  Y,  P,  tag>::Stiffness;  //Used to compute the yield surface intersection.


#endif
