#include "VonMisesLinearHardening.h"
#include <vector>

#include <StandardStream.h>
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main(void)
{
	// VonMisesLinearHardening(int tag_in, double k0_in, double H_alpha, double H_k, double E, double nu, double rho_)
	int tag_in = 1;
	double k0_in = 1;
	double H_alpha = 1;
	double H_k = 1;
	double E = 1;
	double nu = 0.0;
	double rho_ = 1000;
	VonMisesLinearHardening vm(tag_in, k0_in, H_alpha, H_k, E, nu, rho_);
	// double H_scalar = 0;
	// double k = 1.0;
	// LinearHardeningScalar_EV kH(H_scalar, k);

	// double H_tensor = 1;
	// VoigtVector alpha01(0,0,0,0,0,0);
	// VoigtVector alpha02(0.,0,0.5,0,0,0);

	// LinearHardeningTensor_EV alphaH1(H_tensor, alpha01.deviator()); //alpha should be deviatoric
	// LinearHardeningTensor_EV alphaH2(H_tensor, alpha02.deviator());

	// VonMises_YF yf1(alphaH1, kH);
	// VonMises_YF yf2(alphaH2, kH);


	// std::vector<VoigtVector> stresses;

	// VoigtVector depsilon;
	// VoigtVector m;
	// // VoigtVector sigma = {1., 1., 1., 1., 1., 1.0};

	// double p = 10;
	// double qmax = 2;
	// double dq = 0.25;
	// for (double q = 0; q < qmax; q += dq)
	// {
	//     double sigma_xx = p;
	//     double sigma_yy = p;
	//     double sigma_zz = p + q;
	//     double sigma_xy = 0;
	//     double sigma_yz = 0;
	//     double sigma_xz = 0;
	//     stresses.push_back({sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz});
	// }



	// int index = 0;
	// for (auto i = stresses.begin(); i != stresses.end(); ++i, index++)
	// {

	//     VoigtVector sigma = *i;
	//     cout << "\nat sigma (" << index << ") = " << sigma.transpose() << endl;

	//     double yf1_val = yf1(sigma);
	//     cout << "   yf1_val = " << yf1_val << endl;
	//     double yf2_val = yf2(sigma);
	//     cout << "   yf2_val = " << yf2_val << endl;

	//     auto yf1_der = yf1.df_dsigma_ij(sigma);
	//     cout << "   yf1_der = " << yf1_der.transpose() << endl;
	//     auto yf2_der = yf2.df_dsigma_ij(sigma);
	//     cout << "   yf2_der = " << yf2_der.transpose() << endl;

	//     auto xi_star_h_star1 = yf1.xi_star_h_star(depsilon, yf1_der, sigma);
	//     cout << "   xi_star_h_star1 = " << xi_star_h_star1 << endl;
	//     auto xi_star_h_star2 = yf2.xi_star_h_star(depsilon, yf1_der, sigma);
	//     cout << "   xi_star_h_star2 = " << xi_star_h_star2 << endl;


	// }

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
