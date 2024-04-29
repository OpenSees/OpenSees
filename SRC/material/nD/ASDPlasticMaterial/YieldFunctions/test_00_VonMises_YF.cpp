#include "VonMises_YF.h"
#include "../utuple_storage.h"
#include "../AllASDInternalVariableTypes.h"
#include "../AllASDHardeningFunctions.h"
#include <vector>

int main(void)
{
    
    // LinearHardeningForScalar
    ScalarLinearHardeningParameter H_scalar(1.0);
    double initial_scalar = 1.0;
    using VMSL = VonMisesRadiusIV<ScalarLinearHardeningFunction>;
    VMSL kH(initial_scalar);
    
    // LinearHardeningForTensor
    TensorLinearHardeningParameter H_tensor(1.0);
    VoigtVector initial_tensor(0,0,0,0,0,0);
    using BSTL = BackStressIV<TensorLinearHardeningFunction>;
    BSTL alphaH(initial_tensor);


    using VM = VonMises_YF<BSTL, VMSL>;
    VM yf;
    
    //Setup the storage for the model parameters
    using parameter_storage_t = utuple_storage<tuple<ScalarLinearHardeningParameter, TensorLinearHardeningParameter>>;
    parameter_storage_t parameter_storage;
    parameter_storage.set(H_scalar);
    parameter_storage.set(H_tensor);

    //Setup the storage for the internal variables
    using concat_types = utuple_concat_type<VM::internal_variables_t>;
    using iv_storage_t = utuple_storage<concat_types>;
    iv_storage_t iv_storage;
    iv_storage.set(kH);
    iv_storage.set(alphaH);


    std::vector<VoigtVector> stresses;

    VoigtVector depsilon;
    // VoigtVector m;
    // VoigtVector sigma = {1., 1., 1., 1., 1., 1.0};

    double p = 10;
    double qmax = 2;
    double dq = 0.25;
    for (double q = 0; q < qmax; q += dq)
    {
        double sigma_xx = p;
        double sigma_yy = p;
        double sigma_zz = p + q;
        double sigma_xy = 0;
        double sigma_yz = 0;
        double sigma_xz = 0;
        stresses.push_back({sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz});
    }
    


    int index = 0;
    for (auto i = stresses.begin(); i != stresses.end(); ++i, index++)
    {

        VoigtVector sigma = *i;
        cout << "\nat sigma (" << index << ") = " << sigma.transpose() << endl;

        double yf_val = yf(sigma, iv_storage, parameter_storage);
        cout << "   yf_val = " << yf_val << endl;

        auto yf_der = yf.df_dsigma_ij(sigma, iv_storage, parameter_storage);
        cout << "   yf_der = " << yf_der.transpose() << endl;

        auto xi_star_h_star = yf.xi_star_h_star(depsilon, yf_der, sigma, iv_storage, parameter_storage);
        cout << "   xi_star_h_star = " << xi_star_h_star << endl;


    }

    // // VoigtVector alpha(1,1,1,1,1,1);
    // // VoigtScalar k(1);


    return 0;
}


//     VonMises_PF pf2(alphaH2, kH);

//     std::vector<VoigtVector> stresses;

//     VoigtVector depsilon;
//     // VoigtVector sigma = {1., 1., 1., 1., 1., 1.0};
    
//     stresses.push_back({0., 0., 0., 0., 0., 0.});
//     stresses.push_back({1., 1., 1., 0., 0., 0.});
//     stresses.push_back({0., 0., 0., 1., 1., 1.});
//     stresses.push_back({2., 1., 1., 0., 0., 0.});
//     stresses.push_back({1., 1., 1., 1., 1., 1.});

//     int index = 0;
//     for (auto i = stresses.begin(); i != stresses.end(); ++i, index++)
//     {

//         VoigtVector sigma = *i;
//         cout << "\nat sigma (" << index << ") = " << sigma.transpose() << endl;
//         VoigtVector m1 = pf1(depsilon, sigma);
//         VoigtVector m2 = pf2(depsilon, sigma);
//         cout << "   PF1 = " << m1.transpose() << endl;
//         cout << "tr(PF1) = " << m1.trace() << endl;
//         cout << "   PF2 = " << m2.transpose() << endl;
//         cout << "tr(PF2) = " << m2.trace() << endl;
//     }


//     return 0;
// }
