#include "VonMises_YF.h"
#include "../EvolvingVariables/LinearHardeningScalar_EV.h"
#include "../EvolvingVariables/LinearHardeningTensor_EV.h"

#include <vector>

int main(void)
{

    double H_scalar = 1;
    double initial_scalar = 1.0;
    LinearHardeningScalar_EV kH(H_scalar, initial_scalar);

    double H_tensor = 0;
    VoigtVector alpha01(0,0,0,0,0,0);
    VoigtVector alpha02(0,0,0,1,1,1);

    LinearHardeningTensor_EV alphaH1(H_tensor, alpha01);
    LinearHardeningTensor_EV alphaH2(H_tensor, alpha02);

    VonMises_YF yf1(alphaH1, kH);


    std::vector<VoigtVector> stresses;

    VoigtVector depsilon;
    VoigtVector m;
    // VoigtVector sigma = {1., 1., 1., 1., 1., 1.0};
    
    stresses.push_back({0., 0., 0., 0., 0., 0.});
    stresses.push_back({1., 1., 1., 0., 0., 0.});
    stresses.push_back({0., 0., 0., 1., 1., 1.});
    stresses.push_back({2., 1., 1., 0., 0., 0.});
    stresses.push_back({1., 1., 1., 1., 1., 1.});

    int index = 0;
    for (auto i = stresses.begin(); i != stresses.end(); ++i, index++)
    {

        VoigtVector sigma = *i;
        cout << "\nat sigma (" << index << ") = " << sigma.transpose() << endl;

        double yf1_val = yf1(sigma);
        // VoigtVector m2 = pf2(depsilon, sigma);
        cout << "   yf1_val = " << yf1_val << endl;

        auto yf1_der = yf1.df_dsigma_ij(sigma);
        cout << "   yf1_der = " << yf1_der.transpose() << endl;

        auto xi_star_h_star = yf1.xi_star_h_star(depsilon, m, sigma);
        cout << "   xi_star_h_star = " << xi_star_h_star << endl;


    }

    // VoigtVector alpha(1,1,1,1,1,1);
    // VoigtScalar k(1);


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
