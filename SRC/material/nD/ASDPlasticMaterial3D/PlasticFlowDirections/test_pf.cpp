#include "VonMises_PF.h"
#include "../utuple_storage.h"
#include "../AllASDInternalVariableTypes.h"
#include "../AllASDHardeningFunctions.h"

#include <vector>


int main(void)
{

    // double H_scalar = 1;
    ScalarLinearHardeningParameter H_scalar(1.0);
    double initial_scalar = 1.0;
    using VMSL = VonMisesRadiusIV<ScalarLinearHardeningFunction>;
    VMSL kH(initial_scalar);

    TensorLinearHardeningParameter H_tensor(1.0);
    VoigtVector initial_tensor(0,0,0,0,0,0);
    using BSTL = BackStressIV<TensorLinearHardeningFunction>;
    BSTL alphaH(initial_tensor);

    // double H_tensor = 0;
    // LinearHardeningForTensor H_tensor(0);

    // LinearHardeningTensor_EV alphaH(initial_tensor);

    // using BS_t = BackStress<LinearHardeningTensor_EV>;
    // using m_t = VonMisesRadius<LinearHardeningScalar_EV>;

    // BS_t alpha(alphaH);
    // m_t m(kH);


    using VM = VonMises_PF<BSTL, VMSL>;
    VM pf;
 
    using concat_param_types = utuple_concat_type<VM::parameters_t>;
    using param_storage_t = utuple_storage<concat_param_types>;
    param_storage_t param_storage;

    using concat_iv_types = utuple_concat_type<VM::internal_variables_t>;
    using iv_storage_t = utuple_storage<concat_iv_types>;

    iv_storage_t iv_storage;
    iv_storage.set(kH);
    iv_storage.set(alphaH);


    std::cout << "[0] : " << std::get<0>(iv_storage.data) << "\n";
    std::cout << "[1] : " << std::get<1>(iv_storage.data) << "\n";

    std::cout << "alpha : " << iv_storage.get<BSTL>()   << "\n";
    std::cout << "alpha.trial : " << iv_storage.get<BSTL>().trial_value   << "\n";
    std::cout << "alpha.commit : " << iv_storage.get<BSTL>().committed_value   << "\n";
    std::vector<VoigtVector> stresses;

    VoigtVector depsilon;
    
    stresses.push_back({0., 0., 0., 0., 0., 0.});
    stresses.push_back({1., 1., 1., 0., 0., 0.});
    stresses.push_back({0., 0., 0., 1., 1., 1.});
    stresses.push_back({2., 1., 1., 0., 0., 0.});
    stresses.push_back({1., 1., 1., 1., 1., 1.});

    int index = 0;
    for (auto i = stresses.begin(); i != stresses.end(); ++i, index++)
    {
        VoigtVector sigma = *i;
        std::cout << "\nat sigma (" << index << ") = " << sigma.transpose() << std::endl;
        VoigtVector m = pf(depsilon, sigma, iv_storage, param_storage);
        std::cout << "   pf = " << m.transpose() << std::endl;
        std::cout << "tr(pf) = " << m.trace() << std::endl;
    }


    return 0;
}
