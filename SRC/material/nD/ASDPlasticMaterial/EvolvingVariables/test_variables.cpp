#include <iostream>
#include "LinearHardeningScalar_EV.h"
#include "LinearHardeningTensor_EV.h"
#include "../MaterialInternalVariables.h"

#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

typedef MaterialInternalVariables < LinearHardeningScalar_EV, LinearHardeningTensor_EV> InternalVariablesList;

int main(int argc, char const *argv[])
{
	
	LinearHardeningScalar_EV k(1.0, 1.0);

	cout << "k = " << k << endl;

	VoigtVector depsilon 
		= VoigtVector::Zero();//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	VoigtVector m 
		= {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	VoigtVector sigma 
		= VoigtVector::Zero();//{0.0,, 0.0,, 0.0,, 0.0,, 0.0,, 0.0};


	double h = k.getDerivative(depsilon, m, sigma);

	cout << "depsilon = \n" << depsilon << endl;
	cout << "m = \n" << m << endl;

	cout << "LinearHardeningScalar_EV" << endl;
	cout << "h = " << h << endl;

	VoigtVector kd = kronecker_delta();

	cout << "kd = " << kd << endl;


	VoigtVector alpha0;
	double H = 10;

	cout << "LinearHardeningTensor_EV alpha 1 H = " << H << " alpha0 = " << alpha0 << endl;
	cout << "m = \n" << m << endl;
	cout << "alpha0 = \n" << alpha0 << endl;
	LinearHardeningTensor_EV alpha1(H, alpha0);
	VoigtVector dalpha1 = alpha1.getDerivative(depsilon, m, sigma);
	cout << "dalpha1 = " << dalpha1 << endl;

	alpha0 = 2*kronecker_delta();
	H = 10;
	cout << "LinearHardeningTensor_EV alpha2 H = " << H << " alpha0 = " << alpha0 << endl;
	cout << "m = \n" << m << endl;
	LinearHardeningTensor_EV alpha2(H, alpha0);
	VoigtVector dalpha2 = alpha2.getDerivative(depsilon, m, sigma);
	cout << "dalpha2 = " << dalpha2 << endl;

	alpha0.setZero();
	H = 10;
	m = {3.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
	cout << "LinearHardeningTensor_EV alpha3 H = " << H << " alpha0 = " << alpha0 << endl;
	cout << "m = \n" << m << endl;
	LinearHardeningTensor_EV alpha3(H, alpha0);
	VoigtVector dalpha3 = alpha3.getDerivative(depsilon, m, sigma);
	cout << "dalpha3 = " << dalpha3 << endl;



	// Test container of internal variables...;

	double H_scalar = 2;
	double initial_scalar = 1.0;
	LinearHardeningScalar_EV aScalar(H_scalar, initial_scalar);

	double H_tensor = 1000;
	VoigtVector initial_tensor(1,2,3,4,5,6);

	LinearHardeningTensor_EV aTensor(H_tensor, initial_tensor);

	InternalVariablesList state(aScalar, aTensor);

	cout << "Print state" << endl;

	state.print(cout);

	return 0;
}