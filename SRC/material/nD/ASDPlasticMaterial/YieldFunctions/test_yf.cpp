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
    VoigtVector initial_tensor1(0,0,0,0,0,0);
    VoigtVector initial_tensor2(0,0,0,1,1,1);

    LinearHardeningTensor_EV alphaH1(H_tensor, initial_tensor1);
    LinearHardeningTensor_EV alphaH2(H_tensor, initial_tensor2);

    VonMises_YF pf1(alphaH1, kH);

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
