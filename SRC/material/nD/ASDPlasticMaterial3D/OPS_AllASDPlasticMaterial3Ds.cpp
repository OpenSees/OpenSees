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

#ifdef _EIGEN3

#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <list>
#include "AllASDPlasticMaterial3Ds.h"

template <typename T>
void populate_ASDPlasticMaterial3D(T* instance);

typedef std::tuple<std::string, std::string, std::string, std::string> model_spec_t;


NDMaterial*  ASDPlasticMaterial3DFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models);


template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial3D(int instance_tag, const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, std::list<NDMaterial*> &instance_pointers, std::list<model_spec_t> &available_models);


void print_usage(void)
{
    opserr <<
       "nDMaterial ASDPlasticMaterial3D Error: Few arguments \n"
       "\n"
       "SYNTAX:\n"
       "nDMaterial ASDPlasticMaterial3D $tag \\ \n"
       "    YF_type \\ \n"
       "    PF_type \\ \n"
       "    EL_type \\ \n"
       "    IV_type \\ \n"
       "Begin_Internal_Variables \\ \n"
       "    IV1 (initial value(s)) \\ \n"
       "    IV2 (initial value(s)) \\ \n"
       "    ....\n"
       "End_Internal_Variables \\ \n"
       "Begin_Model_Parameters \\ \n"
       "    PARAM1 (value) \\\n"
       "    PARAM2 (value) \\\n"
       "    ....\n"
       "End_Model_Parameters \\ \n"
       "Begin_Integration_Options \\ \n"
       "    f_relative_tol (double value)\\ \n"
       "    stress_relative_tol (double value)\\ \n"
       "    n_max_iterations (int value)\\ \n"
       "    return_to_yield_surface (0 or 1)\\ \n"
       "    method (string) : Forward_Euler | Runge_Kutta_45_Error_Control\\ \n"
       "End_Integration_Options \\ \n"
       "\n";
}

void *OPS_AllASDPlasticMaterial3Ds(void)
{
    // some kudos
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDPlasticMaterial3D - Developed by: Jose Abell (UANDES), Massimo Petracca and Guido Camata (ASDEA Software Technology)\n";
        first_done = true;
    }

    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
        opserr << "numArgs = " << numArgs << endln;
    if (numArgs < 2) {
        print_usage();
        return nullptr;
    }

    int numData;

    int tag = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0)  {
        opserr << "nDMaterial ASDPlasticMaterial3D Error: invalid 'tag'.\n";
        return nullptr;
    }


    const char *yf_type = nullptr;
    const char *pf_type = nullptr;
    const char *el_type = nullptr;
    const char *iv_type = nullptr;

    // Now use conditional checks based on numArgs to assign values
    yf_type = numArgs >= 2 ? OPS_GetString() : " X ";
    pf_type = numArgs >= 3 ? OPS_GetString() : " X ";
    el_type = numArgs >= 4 ? OPS_GetString() : " X ";
    iv_type = numArgs >= 5 ? OPS_GetString() : " X ";


    cout << "Searching for instance with:\n";


    cout << "yf_type = " << yf_type << "\n";
    cout << "pf_type = " << pf_type << "\n";
    cout << "el_type = " << el_type << "\n";
    cout << "iv_type = " << iv_type << "\n";

    std::list<model_spec_t> available_models;

    NDMaterial* instance = ASDPlasticMaterial3DFactory(tag, yf_type, pf_type, el_type, iv_type, available_models);

    if(instance==nullptr)
    {

        bool matches_yf = false;
        bool matches_pf = false;
        bool matches_el = false;


        for(model_spec_t &model : available_models)
        {
            std::string model_yf_type = std::get<0>(model);
            std::string model_pf_type = std::get<1>(model);
            std::string model_el_type = std::get<2>(model);
            // std::string model_iv_type = std::get<3>(model);
            if (std::strcmp(yf_type, model_yf_type.c_str())==0)
                matches_yf = true;
            if (std::strcmp(pf_type, model_pf_type.c_str())==0)
                matches_pf = true;
            if (std::strcmp(el_type, model_el_type.c_str())==0)
                matches_el = true;
        }

        cout << endl;
        print_usage();
        cout << endl;

        cout << "ASDPlasticMaterial3D -- Material not found for input specification:\n";
        cout << "yf_type = " << yf_type << (matches_yf ? " :) " : " ") <<"\n";
        cout << "pf_type = " << pf_type << (matches_pf ? " :) " : " ") <<"\n";
        cout << "el_type = " << el_type << (matches_el ? " :) " : " ") <<"\n";
        cout << "iv_type = " << iv_type  <<"\n";

    }


    return instance;

}



NDMaterial*  ASDPlasticMaterial3DFactory(int instance_tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models)
{

    std::list<NDMaterial*> instance_pointers;


	#include "ASD_material_definitions.cpp"    

    //Search for the valid pointer and return that one
    for(auto instance : instance_pointers)
    {
        if(instance != nullptr)
        {
            return instance;
        }
    }

    return nullptr;
}






template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial3D(int instance_tag, 
        const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, std::list<NDMaterial*> &instance_pointers, std::list<model_spec_t> &available_models) {
    auto instance = new ASDPlasticMaterial3D<EL, YF, PF, ND_TAG_ASDPlasticMaterial3D>(instance_tag);


    available_models.push_back(std::make_tuple(
        instance->getYFName(),
        instance->getPFName(),
        instance->getELName(),
        instance->getIVName())
    );

    if(
        std::strcmp(yf_type,instance->getYFName().c_str())==0 &&
        std::strcmp(pf_type,instance->getPFName().c_str())==0 &&
        std::strcmp(el_type,instance->getELName().c_str())==0 &&
        std::strcmp(iv_type,instance->getIVName().c_str())==0 
        )
    {
        populate_ASDPlasticMaterial3D(instance);
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial3D \n\n\n";
        instance_pointers.push_back(static_cast<NDMaterial*>(instance));
        return static_cast<NDMaterial*>(instance);
    }  
    else
    {
        delete instance;
        instance_pointers.push_back(nullptr);        
        return nullptr;
    }
}








template <typename T>
void populate_ASDPlasticMaterial3D(T* instance)
{

    int get_one_value = 1;

    cout << "\n\nDefined internal variables: \n";
    auto iv_names = instance->getInternalVariablesNames();
    for_each_in_tuple(iv_names, [instance](auto & name)
    {
        std::cout << "   "  << name << " size = " << instance->getInternalVariableSizeByName(name) << std::endl;
    });

    cout << "\n\nDefined model parameters\n";
    auto parameter_names = instance->getParameterNames();
    for_each_in_tuple(parameter_names, [](auto & name)
    {
        std::cout << "   "  <<  name << std::endl;
    });

    // Default integration options
    int method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
    int tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
    double f_relative_tol = 1e-6; 
    double stress_relative_tol = 1e-6; 
    int n_max_iterations = 100;
    int return_to_yield_surface = 1;

    // Loop over input arguments
    while (OPS_GetNumRemainingInputArgs() > 0) {

        //Current command
        const char *cmd = OPS_GetString();


        // Specifying internal variables
        if (std::strcmp(cmd, "Begin_Internal_Variables") == 0)
        {
            cout << "\n\nReading internal variables from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double iv_values[6];
                const char *iv_name = OPS_GetString();

                if (std::strcmp(iv_name, "End_Internal_Variables") == 0)
                {
                    cout << "\n\n Done reading internal variables from input\n";
                    break;
                }

                int iv_size = instance->getInternalVariableSizeByName(iv_name);
                OPS_GetDouble(&iv_size, iv_values);
                cout << iv_name << " = ";
                for (int i = 0; i < iv_size; ++i)
                {
                    cout << iv_values[i] << " ";
                }
                cout << endl;
                instance->setInternalVariableByName(iv_name, iv_size, &iv_values[0]);
            }

        }


        // Specifying model parameters
        if (std::strcmp(cmd, "Begin_Model_Parameters") == 0)
        {
            cout << "\n\nReading parameters from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double param_value;
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Model_Parameters") == 0)
                {
                    cout << "\n\n Done reading parameters from input\n";
                    break;
                }

                OPS_GetDouble(&get_one_value, &param_value);
                cout << param_name << " = " << param_value << endl;
                instance->setParameterByName(param_name, param_value);
            }
        }


        // set_constitutive_integration_method(int method, double f_relative_tol, double stress_relative_tol, int n_max_iterations)
        if (std::strcmp(cmd, "Begin_Integration_Options") == 0)
        {
            cout << "\n\nReading Integration Options\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Integration_Options") == 0)
                {
                    cout << "\n\nDone reading Integration Options\n";
                    break;
                }

                if (std::strcmp(param_name, "f_relative_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &f_relative_tol);
                    cout << "   Setting f_relative_tol = " << f_relative_tol << endl;
                }

                if (std::strcmp(param_name, "stress_relative_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &stress_relative_tol);
                    cout << "   Setting stress_relative_tol = " << stress_relative_tol << endl;
                }

                if (std::strcmp(param_name, "n_max_iterations") == 0)
                {
                    OPS_GetInt(&get_one_value, &n_max_iterations);
                    cout << "   Setting n_max_iterations = " << n_max_iterations << endl;
                }
                
                if (std::strcmp(param_name, "return_to_yield_surface") == 0)
                {
                    OPS_GetInt(&get_one_value, &return_to_yield_surface);
                    cout << "   Setting return_to_yield_surface = " << return_to_yield_surface << endl;
                }

                if (std::strcmp(param_name, "integration_method") == 0)
                {
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Forward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
                    else
                    {
                        cout << "WARNING! Unrecognised ASDPlasticMaterial3D_Constitutive_Integration_Method name " << method_name << endl;
                        cout << "Defaulting to Runge_Kutta_45_Error_Control" << endl;
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
                    }
                    cout << "   Setting integration method = " << method_name << " method_int = " << method << endl;
                    
                }

                if (std::strcmp(param_name, "tangent_type") == 0)
                {
                    const char *tangent_type_name = OPS_GetString();
                    if (std::strcmp(tangent_type_name, "Elastic") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
                    else if (std::strcmp(tangent_type_name, "Continuum") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Continuum;
                    else if (std::strcmp(tangent_type_name, "Secant") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Secant;
                    else
                    {
                        cout << "WARNING! Unrecognised ASDPlasticMaterial3D_Tangent_Operator_Type name " << tangent_type_name << endl;
                        cout << "Defaulting to Elastic" << endl;
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
                    }
                    cout << "   Setting tangent type = " << tangent_type_name << " tangent int = " << tangent << endl;
                    
                }

            }
        }
    }

    instance->set_constitutive_integration_method(method, tangent, f_relative_tol, stress_relative_tol, n_max_iterations, return_to_yield_surface);
}


#endif // _EIGEN3
