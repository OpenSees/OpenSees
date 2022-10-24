#include "LinearIsotropic3D_EL.h"


int main(void)
{


    double E = 1000;
    double nu = 0.25;

    LinearIsotropic3D_EL elasticity(E, nu);

    VoigtVector depsilon = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    VoigtVector sigma = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    VoigtMatrix Et = elasticity(sigma);

    VoigtVector dsigma = Et * depsilon;

    cout << "\ndepsilon = \n" << depsilon.transpose() << endl;
    cout << "\nsigma = \n" << sigma.transpose() << endl;
    cout << "\ndsigma = \n" << dsigma.transpose() << endl;
    cout << "\nEt = \n" << Et.transpose() << endl;

    return 0;
}
