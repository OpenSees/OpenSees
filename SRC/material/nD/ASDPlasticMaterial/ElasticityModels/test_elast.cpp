#include "LinearIsotropic3D_EL.h"
#include "../utuple_storage.h"
// #include "../AllASDModelParameterTypes.h"

int main(void)
{


    YoungsModulus E(1000);
    PoissonsRatio nu(0.25);

    LinearIsotropic3D_EL elasticity;

    using concat_types = utuple_concat_type<LinearIsotropic3D_EL::parameters_t>;
    using parameter_storage_t = utuple_storage<concat_types>;

    parameter_storage_t parameter_storage;

    parameter_storage.set(E);
    parameter_storage.set(nu);


    VoigtVector depsilon = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    VoigtVector sigma = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    VoigtMatrix Et = elasticity(sigma, parameter_storage);

    VoigtVector dsigma = Et * depsilon;

    cout << "E = " << E << endl;
    cout << "nu = " << nu << endl;
    cout << "\ndepsilon = \n" << depsilon.transpose() << endl;
    cout << "\nsigma = \n" << sigma.transpose() << endl;
    cout << "\ndsigma = \n" << dsigma.transpose() << endl;
    cout << "\nEt = \n" << Et.transpose() << endl;

    return 0;
}
