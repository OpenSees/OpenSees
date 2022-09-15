#include "LinearIsotropic3D_EL.h"
#include "../../../ltensor/LTensor.h"


int main(void)
{

    Index<'i'> i;
    Index<'j'> j;
    Index<'k'> k;
    Index<'l'> l;

    double E = 1000;
    double nu = 0.25;

    LinearIsotropic3D_EL elasticity(E, nu);

    DTensor2 eps(3, 3, 0.0);
    DTensor2 sig(3, 3, 0.0);
    DTensor4 Et(3, 3, 3, 3, 0.0);

    Et = elasticity(sig);

    eps(0, 2) = 1;
    eps(2, 0) = 1;

    sig(i, j) = Et(i, j, k, l) * eps(k, l);

    cout << "eps = " << eps << endl;
    cout << "sig = " << sig << endl;
    // cout << "Et = " << Et << endl;

    return 0;
}
