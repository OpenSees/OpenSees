#include <iostream>

#include "../utuple_storage.h"
#include "../AllASDInternalVariableTypes.h"
#include "../AllASDModelParameterTypes.h"
#include "../AllASDHardeningFunctions.h"


#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

using namespace std;

int main(int argc, char const *argv[])
{
    ScalarLinearHardening H_scalar;
    double initial_scalar = 1.0;
    using VMSL = VonMisesRadiusIV<ScalarLinearHardening>;
    VMSL k(initial_scalar);

	cout << "k = " << k << endl;


    TensorLinearHardening H_tensor;
	VoigtVector initial_alpha;
    using BSTL = BackStressIV<TensorLinearHardening>;
    BSTL alpha(initial_alpha);

	cout << "alpha = " << alpha << endl;
	
	// Setup the storage fro the model parameters
	using parameter_storage_t = utuple_storage<tuple<ScalarLinearHardening, TensorLinearHardening>>;
    parameter_storage_t parameter_storage;
    parameter_storage.set(H_scalar);
    parameter_storage.set(H_tensor);

	// Setup the storage for internal variables
    using storage_t = utuple_storage<tuple<VMSL, BSTL>>;
    storage_t internal_variables_storage;
    internal_variables_storage.set(k);
    internal_variables_storage.set(alpha);

	VoigtVector depsilon 
		= VoigtVector::Zero();//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	VoigtVector m 
		= {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	VoigtVector sigma 
		= VoigtVector::Zero();//{0.0,, 0.0,, 0.0,, 0.0,, 0.0,, 0.0};


	// VoigtScalar h_k = k.getDerivative(depsilon, m, sigma, parameters);
	// VoigtVector h_alpha = alpha.getDerivative(depsilon, m, sigma, parameters);

	// cout << "depsilon = \n" << depsilon << endl;
	// cout << "m = \n" << m << endl;

	// cout << "LinearHardeningScalar_EV" << endl;
	// cout << "h_k = " << h_k << endl;
	// cout << "h_alpha = " << h_alpha << endl;





	// cout << "LinearHardeningTensor_EV alpha 1 H = " << H << " alpha0 = " << alpha0 << endl;
	// cout << "m = \n" << m << endl;
	// cout << "alpha0 = \n" << alpha0 << endl;
	// VoigtVector dalpha1 = alpha1.getDerivative(depsilon, m, sigma);
	// cout << "dalpha1 = " << dalpha1 << endl;

	// alpha0 = 2*kronecker_delta();
	// H = 10;
	// cout << "LinearHardeningTensor_EV alpha2 H = " << H << " alpha0 = " << alpha0 << endl;
	// cout << "m = \n" << m << endl;
	// LinearHardeningTensor_EV alpha2(H, alpha0);
	// VoigtVector dalpha2 = alpha2.getDerivative(depsilon, m, sigma);
	// cout << "dalpha2 = " << dalpha2 << endl;

	// alpha0.setZero();
	// H = 10;
	// m = {3.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
	// cout << "LinearHardeningTensor_EV alpha3 H = " << H << " alpha0 = " << alpha0 << endl;
	// cout << "m = \n" << m << endl;
	// LinearHardeningTensor_EV alpha3(H, alpha0);
	// VoigtVector dalpha3 = alpha3.getDerivative(depsilon, m, sigma);
	// cout << "dalpha3 = " << dalpha3 << endl;



	// // Test container of internal variables...;

	// double H_scalar = 2;
	// double initial_scalar = 1.0;
	// LinearHardeningScalar_EV aScalar(H_scalar, initial_scalar);

	// double H_tensor = 1000;
	// VoigtVector initial_tensor(1,2,3,4,5,6);

	// LinearHardeningTensor_EV aTensor(H_tensor, initial_tensor);

	// InternalVariablesList state(aScalar, aTensor);

	// cout << "Print state" << endl;

	// state.print(cout);

	return 0;
}