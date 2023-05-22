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
#include "AllASDPlasticMaterials.h"

template <typename T>
void populate_ASDPlasticMaterial(T* instance);

NDMaterial*  ASDPlasticMaterialFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type);



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



    NDMaterial* instance = ASDPlasticMaterialFactory(tag, yf_type,pf_type,el_type,iv_type);

    if(instance==0)
    {
        cout << "ASDPlasticMaterial -- Material not found for input specification:\n";
        cout << "yf_type = " << yf_type << "\n";
        cout << "pf_type = " << pf_type << "\n";
        cout << "el_type = " << el_type << "\n";
        cout << "iv_type = " << iv_type << "\n";
    }

    return instance;

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
                    break;
                }

                OPS_GetDouble(&get_one_value, &param_value);
                cout << param_name << " = " << param_value << endl;
                instance->setParameterByName(param_name, param_value);
            }
        }
    }
}

NDMaterial*  ASDPlasticMaterialFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type)
{

    cout << "Searching for instance with:\n";


    cout << "yf_type = " << yf_type << "\n";
    cout << "pf_type = " << pf_type << "\n";
    cout << "el_type = " << el_type << "\n";
    cout << "iv_type = " << iv_type << "\n";

    if(
        std::strcmp(yf_type,"VonMises_YF")==0 &&
        std::strcmp(pf_type,"VonMises_PF")==0 &&
        std::strcmp(el_type,"LinearIsotropic3D_EL")==0 &&
        std::strcmp(iv_type,"BackStress(TensorLinearHardeningFunction):VonMisesRadius(ScalarLinearHardeningFunction):")==0
        )
    {
        auto instance = new 
            ASDPlasticMaterial <LinearIsotropic3D_EL,
                VonMises_YF<
                    BackStress<TensorLinearHardeningFunction>,
                    VonMisesRadius<ScalarLinearHardeningFunction>
                >,
                VonMises_PF<
                    BackStress<TensorLinearHardeningFunction>,
                    VonMisesRadius<ScalarLinearHardeningFunction>
                >,
            ND_TAG_ASDPlasticMaterial
            >(tag);
        populate_ASDPlasticMaterial(instance);
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial \n\n\n";
        return static_cast<NDMaterial*>(instance);
    }   
    else if(
    std::strcmp(yf_type,"DruckerPrager_YF")==0 &&
    std::strcmp(pf_type,"VonMises_PF")==0 &&
    std::strcmp(el_type,"LinearIsotropic3D_EL")==0 &&
    std::strcmp(iv_type,"BackStress(TensorLinearHardeningFunction):VonMisesRadius(ScalarLinearHardeningFunction):")==0
    )
    {
        auto instance = new 
            ASDPlasticMaterial <LinearIsotropic3D_EL,
                DruckerPrager_YF<
                    BackStress<TensorLinearHardeningFunction>,
                    VonMisesRadius<ScalarLinearHardeningFunction>
                >,
                VonMises_PF<
                    BackStress<TensorLinearHardeningFunction>,
                    VonMisesRadius<ScalarLinearHardeningFunction>
                >,
            ND_TAG_ASDPlasticMaterial
            >(tag);
        populate_ASDPlasticMaterial(instance);
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial \n\n\n";
        return static_cast<NDMaterial*>(instance);
    } 



    return 0;
}

























//Idea to automatically generate tensor product of types...

#include <tuple>
#include <utility>

template<typename T, typename S, typename U>
class TypeHolder
{
public:
    TypeHolder(T t_, S s_, U u_) : t(t_), s(s_), u(u_) {}

private:
    T t;
    S s;
    U u;
};

using the_types_T = std::tuple<double, int, char *>;
using the_types_S = std::tuple<int, double>;
using the_types_U = std::tuple<float, long>; // example types for U

template <typename Sequence1, typename Sequence2, typename Sequence3>
class TypeHolderGenerator;

template <std::size_t... I1, std::size_t... I2, std::size_t... I3>
class TypeHolderGenerator<std::index_sequence<I1...>, std::index_sequence<I2...>, std::index_sequence<I3...>>
{
public:
    std::tuple<TypeHolder<std::tuple_element_t<I1, the_types_T>, std::tuple_element_t<I2, the_types_S>, std::tuple_element_t<I3, the_types_U>>...> holders;
};

using HolderGenerator = TypeHolderGenerator<
    std::make_index_sequence<std::tuple_size<the_types_T>::value>,
    std::make_index_sequence<std::tuple_size<the_types_S>::value>,
    std::make_index_sequence<std::tuple_size<the_types_U>::value>>;