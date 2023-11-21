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

#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <list>
#include "AllASDPlasticMaterials.h"

template <typename T>
void populate_ASDPlasticMaterial(T* instance);

NDMaterial*  ASDPlasticMaterialFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type);


template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial(int instance_tag, const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, std::list<NDMaterial*> &instance_pointers);

void *OPS_AllASDPlasticMaterials(void)
{
    // some kudos
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDPlasticMaterial - Developed by: Jose Abell (UANDES), Massimo Petracca and Guido Camata (ASDEA Software Technology)\n";
        first_done = true;
    }

    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr <<
               "nDMaterial ASDPlasticMaterial Error: Few arguments (< 3).\n"
               "nDMaterial ASDPlasticMaterial $tag "
               "-yf $YF_type <YF params> "
               "-pf $PF_type <PF params> "
               "-el $EL_type <EL params> "
               "-ev $EV_NUM $EV_type <EV params> "
               " .... repeat for all required EV's"
               "<-rho $rho> "
               "<-integrator $INTEGRATOR_TYPE> "
               // "<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha>"
               "\n";
        return nullptr;
    }

    int numData;

    int tag = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0)  {
        opserr << "nDMaterial ASDPlasticMaterial Error: invalid 'tag'.\n";
        return nullptr;
    }


    const char *yf_type = OPS_GetString();
    const char *pf_type = OPS_GetString();
    const char *el_type = OPS_GetString();
    const char *iv_type = OPS_GetString();


    // auto instance = TheAvailableProtos[std::make_tuple(yf_type,pf_type,el_type,iv_type)];


    cout << "Searching for instance with:\n";


    cout << "yf_type = " << yf_type << "\n";
    cout << "pf_type = " << pf_type << "\n";
    cout << "el_type = " << el_type << "\n";
    cout << "iv_type = " << iv_type << "\n";

    NDMaterial* instance = ASDPlasticMaterialFactory(tag, yf_type,pf_type,el_type,iv_type);

    if(instance==nullptr)
    {
        cout << "ASDPlasticMaterial -- Material not found for input specification:\n";
        cout << "yf_type = " << yf_type << "\n";
        cout << "pf_type = " << pf_type << "\n";
        cout << "el_type = " << el_type << "\n";
        cout << "iv_type = " << iv_type << "\n";
    }

    return instance;

}



NDMaterial*  ASDPlasticMaterialFactory(int instance_tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type)
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
NDMaterial* createASDPlasticMaterial(int instance_tag, 
        const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, std::list<NDMaterial*> &instance_pointers) {
    auto instance = new ASDPlasticMaterial<EL, YF, PF, ND_TAG_ASDPlasticMaterial>(instance_tag);
    if(
        std::strcmp(yf_type,instance->getYFName().c_str())==0 &&
        std::strcmp(pf_type,instance->getPFName().c_str())==0 &&
        std::strcmp(el_type,instance->getELName().c_str())==0 &&
        std::strcmp(iv_type,instance->getIVName().c_str())==0 
        )
    {
        populate_ASDPlasticMaterial(instance);
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial \n\n\n";
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
void populate_ASDPlasticMaterial(T* instance)
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
    int method = (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler;
    double f_relative_tol = 1e-6; 
    double stress_relative_tol = 1e-6; 
    int n_max_iterations = 100;

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

                if (std::strcmp(param_name, "method") == 0)
                {
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Forward_Euler") == 0)
                        method = (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
                    else
                    {
                        cout << "WARNING! Unrecognised ASDPlasticMaterial_Constitutive_Integration_Method name " << method_name << endl;
                        cout << "Defaulting to Forward_Euler" << endl;
                        method = (int) ASDPlasticMaterial_Constitutive_Integration_Method::Forward_Euler;
                    }
                    cout << "   Setting integration method = " << method_name << " method_int = " << method << endl;
                    
                }

            }
        }
    }

    instance->set_constitutive_integration_method(method, f_relative_tol, stress_relative_tol, n_max_iterations);
}
