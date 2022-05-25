/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for a bar-slip model with local axial and bond-slip models
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#include <BarSlipMaterial2.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void*
OPS_BarSlipMaterial2(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 7) {
		opserr << "WARNING! Insufficient args in uniaxialMaterial BarSlip2" << endln;
		opserr << "want: uniaxialMaterial BarSlip2 tag? axialTag? bondTag? Ab? Cb? L? N? <-iter maxIterNo? maxTol?>" << endln;
		return 0;
	}

	int tags[3];
	double dData[3];
	int N;
	int maxIterNo = 50;
	double maxTol = 1.0e-6;
	bool dbl = false;

	int numData = 3;
	if (OPS_GetIntInput(&numData, tags) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial BarSlip2" << endln;
		return 0;
	}

	numData = 3;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid Ab, Cb, or L for uniaxialMaterial BarSlip2 " << tags[0] << endln;
		return 0;
	}

	numData = 1;
	if (OPS_GetIntInput(&numData, &N) != 0) {
		opserr << "WARNING! Invalid N for uniaxialMaterial BarSlip2 " << tags[0] << endln;
		return 0;
	}

	while (OPS_GetNumRemainingInputArgs() > 1) {
		const char* opt = OPS_GetString();

		if (strcmp(opt, "-iter") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, &maxIterNo) != 0) {
				opserr << "WARNING! Invalid maxIterNo for uniaxialMaterial BarSlip2 " << tags[0] << endln;
				return 0;
			}

			numData = 1;
			if (OPS_GetDoubleInput(&numData, &maxTol) != 0) {
				opserr << "WARNING! Invalid maxTol for uniaxialMaterial BarSlip2 " << tags[0] << endln;
				return 0;
			}
		}
	}

	// create pointers to reference material models
	UniaxialMaterial* axialMat = OPS_getUniaxialMaterial(tags[1]);

	if (axialMat == 0) {
		opserr << "WARNING! Reference material model " << tags[1] << " does not exist" << endln;
		opserr << "uniaxialMaterial BarSlip2: " << tags[0] << endln;

		delete axialMat;
	}

	UniaxialMaterial* bondMat = OPS_getUniaxialMaterial(tags[2]);

	if (axialMat == 0) {
		opserr << "WARNING! Reference material model " << tags[2] << " does not exist" << endln;
		opserr << "uniaxialMaterial BarSlip2: " << tags[0] << endln;

		delete bondMat;
	}

	// copy material models
	int aMatsNo, bMatsNo;

	if (!free) {
		aMatsNo = N;
		bMatsNo = N - 1;
	}
	else {
		aMatsNo = N - 1;
		bMatsNo = N;
	}

	UniaxialMaterial** aMats = new UniaxialMaterial * [aMatsNo];
	UniaxialMaterial** bMats = new UniaxialMaterial * [bMatsNo];

	for (int i = 0; i < aMatsNo; i++)
		aMats[i] = axialMat;

	for (int i = 0; i < bMatsNo; i++)
		bMats[i] = bondMat;

	// allocate material model
	theMaterial = new BarSlipMaterial2(tags[0], aMats, bMats, dData[0], dData[1], dData[2], N, 0.0, maxIterNo, maxTol, false);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type BarSlip2\n";
		return 0;
	}

	delete[] aMats;
	delete[] bMats;

	return theMaterial;
}


BarSlipMaterial2::BarSlipMaterial2(int tag, UniaxialMaterial** aMats, UniaxialMaterial** bMats, double ab, double cb, double l, int n, double Lc, int iter, double tol, bool dbl)
	:UniaxialMaterial(tag, MAT_TAG_BarSlip2),
	Ab(ab), Cb(cb), L(l), N(n), lc(Lc), maxIterNo(iter), maxTol(tol), doubleSided(dbl),
	lambda((cb / ab) * l / (2.0 * (n - 1.0))), dx(l / (n - 1.0)),
	s_L_T(0.0), s_L_C(0.0), tau_L(0.0), sig_L_T(0.0), sig_L_C(0.0),
	E_tan_T(0.0), E_tan_C(0.0),
	s_u_T(0), s_u_C(0), s_ub(0), eps_T(0), eps_C(0), tau_u(0), sig(0),
	J_o(2 * n - 2, 2 * n -2), J_T(2 * n - 2, 2 * n - 2), J_C(2 * n - 2, 2 * n - 2), H_inv(n, n),
	axialMats(0), bondMats(0)
{
	// create pointers for state variables
	s_u_T = new double[N - 2];
	s_u_C = new double[N - 2];
	s_ub = new double[N - 2];
	tau_u = new double[N - 2];
	eps_T = new double[N];
	eps_C = new double[N];
	sig = new double[N];

	if (s_u_T == 0 || s_u_C == 0 || s_ub == 0 || eps_T == 0 || eps_C == 0 || tau_u == 0 || sig == 0) {
		opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
		opserr << "- ran out of memory for state variables\n";
		exit(-1);
	}

	for (int i = 0; i < N - 2; i++)		// zero the pointer values
		s_u_T[i] = s_u_C[i] = s_ub[i] = tau_u[i] = 0.0;

	for (int i = 0; i < N; i++)		// zero the pointer values
		eps_T[i] = eps_C[i] = sig[i] = 0.0;

	// create arrays to store copies of material models
	axialMats = new UniaxialMaterial* [N];

	if (axialMats == 0) {
		opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
		opserr << "- ran out of memory for axial material model array of size: " << N << endln;
		exit(-1);
	}

	bondMats = new UniaxialMaterial* [N - 1];	// excluding the bar's embedded end, assuming its slip is zero

	if (bondMats == 0) {
		opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
		opserr << "- ran out of memory for bond matetrial model array of size: " << N-1 << endln;
		exit(-1);
	}

	// store copies of material models into their arrays
	for (int i = 0; i < N; i++) {
		axialMats[i] = aMats[i]->getCopy();

		if (axialMats[i] == 0) {
			opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
			opserr << "- failed to get copy of axial material model\n";
			exit(-1);
		}
	}
	
	for (int i = 0; i < N - 1; i++) {
		bondMats[i] = bMats[i]->getCopy();

		if (bondMats[i] == 0) {
			opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
			opserr << "- failed to get copy of bond material model\n";
			exit(-1);
		}
	}

	// create [H_inv] matrix
	Matrix H(N, N); H.Zero();

	double A = 1.0 + pow(lc / dx, 2.0);
	double B = -0.5 * pow(lc / dx, 2.0);

	H(0, 0) = H(N - 1, N - 1) = 1.0;
	for (int i = 1; i < N - 1; i++) {
		H(i, i) = A;
		H(i, i - 1) = H(i, i + 1) = B;
	}

	if (H.Invert(H_inv) != 0) {
		opserr << "FATAL BarSlipMaterial2::BarSlipMaterial2() ";
		opserr << "- failed to invert [H]\n";
		exit(-1);
	}

	// create [F_su]
	Matrix F_su(N - 2, N); F_su.Zero();
	Matrix F_suT(N, N - 2); F_suT.Zero();	// transpose of [F_su]

	for (int i = 0; i < N - 2; i++)
		F_su(i, i + 1) = F_suT(i + 1, i) = 1.0;

	//// create the initial Jacobian matrix
	J_o.Zero();

	// submatrix: [F_su][H_inv][F_suT]
	Matrix FHF(N - 2, N - 2);
	FHF = F_su * (H_inv * F_suT);

	for (int i = 0; i < N - 2; i++)
		for (int j = 0; j < N - 2; j++)
			J_o(i, j) = FHF(i, j);

	// submatrix: -[W_u]
	for (int i = 0; i < N - 2; i++) {
		J_o(i, N - 2) = J_o(i, N - 1 + i) = -0.5 * dx;

		for (int j = 1; j <= i; j++)
			J_o(i, N - 2 + j) = -dx;
	}

	// submatrix: -[W_L]
	J_o(N - 2, N - 2) = J_o(N - 2, 2 * N - 3) = -0.5 * dx;
	for (int j = 1; j < N - 1; j++)
		J_o(N - 2, N - 2 + j) = -dx;

	// submatrix: -[F_tu][K_t]
	for (int j = 0; j < N - 2; j++)
		J_o(N - 1 + j, j) = J_o(N + j, j) = -lambda * (bondMats[j]->getInitialTangent());

	// submatrix: [F_s][K_s]
	J_o(N - 1, N - 2) = -(axialMats[0]->getInitialTangent());
	J_o(2 * N - 3, 2 * N - 3) = axialMats[N - 1]->getInitialTangent();

	for (int j = 0; j < N - 2; j++) {
		J_o(N - 1 + j, N - 1 + j) = axialMats[j + 1]->getInitialTangent();
		J_o(N + j, N - 1 + j) = -(axialMats[j + 1]->getInitialTangent());
	}

	J_T = J_C = J_o;

	// print [J_o]
	/*for (int i = 0; i < 2 * N - 2; i++) {
		for (int j = 0; j < 2 * N - 2; j++)
			opserr << J_o(i, j) << '\t';
		opserr << endln;
	}
	opserr << endln;*/

	// initialize E_tan
	E_tan_T = E_tan_C = this->getInitialTangent();
}


BarSlipMaterial2::BarSlipMaterial2()
	:UniaxialMaterial(0, MAT_TAG_BarSlip2),
	L(0.0), N(0), lc(0.0), maxIterNo(0), maxTol(0.0), doubleSided(0),
	lambda(0.0), dx(0.0),
	s_L_T(0.0), s_L_C(0.0), tau_L(0.0), sig_L_T(0.0), sig_L_C(0.0), 
	E_tan_T(0.0), E_tan_C(0.0),
	s_u_T(0), s_u_C(0), s_ub(0), eps_T(0), eps_C(0), tau_u(0), sig(0),
	J_o(0,0), J_T(0, 0), J_C(0, 0), H_inv(0, 0),
	axialMats(0), bondMats(0)
{

}


BarSlipMaterial2::~BarSlipMaterial2()
{
	// delete the state variable arrays
	if (s_u_T != 0)
		delete[] s_u_T;
	if (s_u_C != 0)
		delete[] s_u_C;
	if (s_ub != 0)
		delete[] s_ub;
	if (tau_u != 0)
		delete[] tau_u;
	if (eps_T != 0)
		delete[] eps_T;
	if (eps_C != 0)
		delete[] eps_C;
	if (sig != 0)
		delete[] sig;
	
	// delete the material arrays
	for (int i = 0; i < N; i++)
		if (axialMats[i] != 0)
			delete axialMats[i];

	for (int i = 0; i < N - 1; i++)
		if (bondMats[i] != 0)
			delete bondMats[i];

	if (axialMats != 0)
		delete[] axialMats;

	if (bondMats != 0)
		delete[] bondMats;
}


int
BarSlipMaterial2::setTrialStrain(double strain, double strainRate)
{
	if (!doubleSided)
		s_L_T = strain;	// slip at bar's free end
	else
		s_L_T = strain / 2.;

	if (fabs(s_L_T - s_L_C) <= DBL_EPSILON) {
		sig_L_T = sig_L_C;
		E_tan_T = E_tan_C;

		return 0;
	}

	// compute bond stress at end IP
	bondMats[N - 2]->setTrialStrain(s_L_T, 0.0);		// set slip at bar's free end's IP
	tau_L = bondMats[N - 2]->getStress();	// bond stress at bar's free end's IP

	// iteratibe solution
	Vector y(2 * N - 2);		// unknowns vector
	Vector y_pre(2 * N - 2);	// previous iteration's unknowns vector
	Vector dy(2 * N - 2);		// unknowns vector correction (via Newton-Raphson)
	Vector R(2 * N - 2);		// error vector

	double R_ele;	// used in calculation of {R} elements

	double ds;	// slip calculation error norm
	double df;	// force calculation error norm

	int converged = -1;		// equals 0 when converged for the whole step

	// for subdivision
	double s_L;
	double r = 1.0;
	double r_pre = 0.0;
	double dr = 1.0;

	const double subdiv = 10.0;
	const int maxDivNo = 4;

	bool convergence;

	Matrix J_pre = J_C;

	double *s_u_pre = new double[N - 2];
	double* eps_pre = new double[N];

	for (int i = 0; i < N - 2; i++)
		s_u_pre[i] = s_u_C[i];

	for (int i = 0; i < N; i++)
		eps_pre[i] = eps_C[i];

	for (int divNo = 0; divNo < maxDivNo; divNo++) {

		// current s_L considering subdivision
		s_L = s_L_T - (1.0 - r) * (s_L_T - s_L_C);

		convergence = false;

		for (int m = 0; m < 2; m++) {	// three methods with combinations of initial and tangent Jacobians

			// initial guesses for {y} elements
			for (int i = 0; i < N - 2; i++)
				s_u_T[i] = s_u_pre[i];

			for (int i = 0; i < N; i++)
				eps_T[i] = eps_pre[i];

			// iterations
			for (int iter = 0; iter < maxIterNo; iter++) {

				//// update material models with unknown inputs and get their output stresses
				for (int i = 0; i < N - 2; i++) {
					bondMats[i]->setTrialStrain(s_u_T[i], 0.0);
					tau_u[i] = bondMats[i]->getStress();
				}

				for (int i = 0; i < N; i++) {
					axialMats[i]->setTrialStrain(eps_T[i], 0.0);
					sig[i] = axialMats[i]->getStress();
				}

				//// compute nonlocal slip trials
				for (int i = 0; i < N - 2; i++) {
					if (lc == 0.0)
						s_ub[i] = s_u_T[i];
					else {
						s_ub[i] = H_inv(i + 1, N - 1) * s_L;

						for (int j = 0; j < N - 2; j++)
							s_ub[i] += H_inv(i + 1, j + 1) * s_u_T[j];
					}
				}

				//// compute {R} elements
				// {s_ub} - [W_u]{eps}
				for (int i = 0; i < N - 2; i++) {
					R_ele = s_ub[i] - 0.5 * dx * (eps_T[0] + eps_T[i + 1]);

					for (int j = 1; j <= i; j++)
						R_ele -= (dx * eps_T[j]);

					R(i) = R_ele;
				}

				// s_L - [W_L]{eps}
				R_ele = s_L - 0.5 * dx * (eps_T[0] + eps_T[N - 1]);

				for (int j = 1; j <= N - 2; j++)
					R_ele -= (dx * eps_T[j]);

				R(N - 2) = R_ele;

				// [F_s]{sig} - [F_tu]{tau_u} - tau_L {F_tL}
				R(N - 1) = -sig[0] + sig[1] - lambda * tau_u[0];
				R(2 * N - 3) = -sig[N - 2] + sig[N - 1] - lambda * (tau_u[N - 3] + tau_L);

				for (int i = 0; i < N - 3; i++)
					R(N + i) = -sig[i + 1] + sig[i + 2] - lambda * (tau_u[i] + tau_u[i + 1]);

				//// check convergence and update solution if needed
				// compute error norms
				ds = df = 0.0;

				for (int i = 0; i < N - 1; i++) {
					ds += (R(i) * R(i));
					df += (R(N - 1 + i) * R(N - 1 + i));
				}

				if (sqrt(ds) <= fmax(0.01 * maxTol * s_L, maxTol) && sqrt(df) <= fmax(0.01 * maxTol * tau_L, maxTol)) {	// converged
					convergence = true;
					break;
				}
				else if (iter < maxIterNo - 1) {	// not converged
					//// form {y}
					for (int i = 0; i < N - 2; i++)
						y(i) = s_u_T[i];

					for (int i = 0; i < N; i++)
						y(N - 2 + i) = eps_T[i];

					//// update [J]
					if (m == 0 || (m == 2 && iter > 0)) {
						// submatrix: -[F_tu][K_t]
						for (int j = 0; j < N - 2; j++) {
							J_T(N - 1 + j, j) = J_T(N + j, j) = -lambda * (bondMats[j]->getTangent());
						}

						// submatrix: [F_s][K_s]
						J_T(N - 1, N - 2) = -(axialMats[0]->getTangent());
						J_T(2 * N - 3, 2 * N - 3) = axialMats[N - 1]->getTangent();

						for (int j = 0; j < N - 2; j++) {
							J_T(N - 1 + j, N - 1 + j) = axialMats[j + 1]->getTangent();
							J_T(N + j, N - 1 + j) = -(axialMats[j + 1]->getTangent());
						}
					}
					else if (m == 1)
						J_T = J_pre;

					// print [J]
					/*if (iter > maxIterNo - 4) {
						for (int i = 0; i < 2 * N - 2; i++) {
							for (int j = 0; j < 2 * N - 2; j++)
								opserr << J(i, j) << '\t';
							opserr << endln;
						}
						opserr << endln;
					}*/

					//// Newton-Raphson correction
					if (J_T.Solve(R, dy) < 0) {
						opserr << "WARNING! BarSlipMaterial2::setTrialStrain() - tag: " << this->getTag()
							<< "\ncould not invert Jacobian\n";
						break;
					}

					// maximum correction size control
					double gamma = 1.0;
					/*double dy_max = 1.0e-4;
					if (iter > maxIterNo / 2)
						for (int i = 0; i < 2 * N - 2; i++)
							if (fabs(dy(i)) > dy_max)
								gamma = fmin(gamma, dy_max / fabs(dy(i)));*/

					y -= gamma * dy;

					/*if (iter > maxIterNo / 2) {
						y += y_pre;
						y /= 2.0;
					}
					
					// save {y}
					y_pre = y;*/

					/*if (iter > maxIterNo - 6) {
						opserr << '\n' << iter << endln;
						for (int i = 0; i < 2 * N - 2; i++)
							opserr << y(i) << '\t';
						opserr << endln;

						for (int i = 0; i < 2 * N - 2; i++)
							opserr << R(i) << '\t';
						opserr << endln;
					}*/

					//// extract corrected unknowns
					for (int i = 0; i < N - 2; i++)
						s_u_T[i] = y(i);

					for (int i = 0; i < N; i++)
						eps_T[i] = y(N - 2 + i);
				}
			}	// end iterations loop

			// break if converged for current s_L from subdivision
			if (convergence)
				break;
		}	// end iteration method loop

		if (convergence) {
			if (r == 1.0) {
				converged = 0;
				break;
			}
			else {
				// increase subdivision ratio
				divNo--;
				r_pre = r;
				r = fmax(r + dr, 1.0);

				//opserr << r << endln;

				// save previous converged solution
				for (int i = 0; i < N - 2; i++)
					s_u_pre[i] = s_u_T[i];

				for (int i = 0; i < N; i++)
					eps_pre[i] = eps_T[i];

				J_pre = J_T;
			}
		}
		else {
			dr /= subdiv;
			r = r_pre + dr;
			//opserr << r << endln;
		}
	}	// end while loop

	if (converged != 0)
		opserr << "\nWARNING! BarSlipMaterial2::setTrialStrain() - tag: " << this->getTag()
			<< "\nNo convergence after " << maxIterNo << " iterations - Error norm: " << R.Norm()
			<< "\nTrial slip: " << s_L_T << " - Slip increment: " << s_L_T - s_L_C << endln;
	else {
		sig_L_T = sig[N - 1];
		E_tan_T = (fabs(s_L_T - s_L_C) < 1.0e-16) ? E_tan_C : ((sig_L_T - sig_L_C) / (s_L_T - s_L_C));
	}

	// delete dynamic pointer arrays
	if (s_u_pre != 0)
		delete[] s_u_pre;

	if (eps_pre != 0)
		delete[] eps_pre;

	return converged;
}


double
BarSlipMaterial2::getStrain(void)
{
	return s_L_T;
}


double
BarSlipMaterial2::getStress(void)
{
	return sig_L_T;
}


double
BarSlipMaterial2::getTangent(void)
{
	return E_tan_T;
}


Vector
BarSlipMaterial2::getAxialStrains(void)
{
	Vector vecN(N);

	for (int i = 0; i < N; i++)
		vecN(i) = eps_T[i];

	return vecN;
}


Vector
BarSlipMaterial2::getAxialStresses(void)
{
	Vector vecN(N);

	for (int i = 0; i < N; i++)
		vecN(i) = sig[i];

	return vecN;
}


Vector
BarSlipMaterial2::getSlips(void)
{
	Vector vecN(N);
	
	vecN(0) = 0.0;
	vecN(N - 1) = s_L_T;
	for (int i = 0; i < N - 2; i++)
		vecN(i + 1) = s_u_T[i];

	return vecN;
}


Vector
BarSlipMaterial2::getBondStresses(void)
{
	Vector vecN(N);
	
	vecN(0) = 0.0;
	vecN(N - 1) = tau_L;
	for (int i = 0; i < N - 2; i++)
		vecN(i + 1) = tau_u[i];

	return vecN;
}


double
BarSlipMaterial2::getInitialTangent(void)
{
	return (axialMats[N - 1]->getInitialTangent());
}


int
BarSlipMaterial2::commitState(void)
{
	int err0, err = 0;

	// invoke commitState() on each MaterialModel object
	for (int i = 0; i < N; i++) {

		err0 = axialMats[i]->commitState();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::commitState() ";
			opserr << "- MaterialModel failed to commitState():";
			axialMats[i]->Print(opserr);

			err += err0;
		}
	}

	for (int i = 0; i < N - 1; i++) {

		err0 = bondMats[i]->commitState();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::commitState() ";
			opserr << "- MaterialModel failed to commitState():";
			bondMats[i]->Print(opserr);

			err += err0;
		}
	}

	// commit state variables
	J_C = J_T;
	s_L_C = s_L_T;
	sig_L_C = sig_L_T;
	E_tan_C = E_tan_T;

	for (int i = 0; i < N - 2; i++)
		s_u_C[i] = s_u_T[i];

	for (int i = 0; i < N; i++)
		eps_C[i] = eps_T[i];
	
	return err;
}


int
BarSlipMaterial2::revertToLastCommit(void)
{
	int err0, err = 0;

	// invoke revertToLastCommit() on each MaterialModel object
	for (int i = 0; i < N; i++) {

		err0 = axialMats[i]->revertToLastCommit();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::revertToLastCommit() ";
			opserr << "- MaterialModel failed to revertToLastCommit():";
			axialMats[i]->Print(opserr);

			err += err0;
		}
	}

	for (int i = 0; i < N - 1; i++) {

		err0 = bondMats[i]->revertToLastCommit();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::revertToLastCommit() ";
			opserr << "- MaterialModel failed to revertToLastCommit():";
			bondMats[i]->Print(opserr);

			err += err0;
		}
	}

	// state variables
	s_L_T = s_L_C;
	sig_L_T = sig_L_C;
	E_tan_T = E_tan_C;

	return err;
}


int
BarSlipMaterial2::revertToStart(void)
{
	int err0, err = 0;

	// zero state variables
	s_L_C = s_L_T = 0.0;
	sig_L_C = sig_L_T = 0.0;
	E_tan_C = E_tan_T = this->getInitialTangent();
	tau_L = 0.0;

	for (int i = 0; i < N - 2; i++)		// zero the pointer values
		s_u_T[i] = s_u_C[i] = tau_u[i] = 0.0;

	for (int i = 0; i < N; i++)		// zero the pointer values
		eps_T[i] = eps_C[i] = sig[i] = 0.0;

	// invoke revertToStart() on each MaterialModel object
	for (int i = 0; i < N; i++) {

		err0 = axialMats[i]->revertToStart();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::revertToStart() ";
			opserr << "- MaterialModel failed to revertToStart():";
			axialMats[i]->Print(opserr);

			err += err0;
		}
	}

	for (int i = 0; i < N - 1; i++) {

		err0 = bondMats[i]->revertToStart();

		if (err0 != 0) {
			opserr << "WARNING BarSlipMaterial2::revertToStart() ";
			opserr << "- MaterialModel failed to revertToStart():";
			bondMats[i]->Print(opserr);

			err += err0;
		}
	}

	return err;
}


UniaxialMaterial *
BarSlipMaterial2::getCopy(void)
{
	BarSlipMaterial2 *theCopy = new BarSlipMaterial2(this->getTag(), axialMats, bondMats, Ab, Cb, L, N, lc, maxIterNo, maxTol, doubleSided);

	theCopy->s_L_T = s_L_T;
	theCopy->s_L_C = s_L_C;
	theCopy->tau_L = tau_L;
	theCopy->sig_L_T = sig_L_T;
	theCopy->sig_L_C = sig_L_C;
	theCopy->E_tan_T = E_tan_T;
	theCopy->E_tan_C = E_tan_C;
	theCopy->H_inv = H_inv;

	for (int i = 0; i < N - 2; i++) {
		theCopy->s_u_T[i] = s_u_T[i];
		theCopy->s_u_C[i] = s_u_C[i];
		theCopy->tau_u[i] = tau_u[i];
	}

	for (int i = 0; i < N; i++) {
		theCopy->eps_T[i] = eps_T[i];
		theCopy->eps_C[i] = eps_C[i];
		theCopy->sig[i] = sig[i];
	}

	return theCopy;
}


int
BarSlipMaterial2::sendSelf(int cTag, Channel &theChannel)
{
	opserr << "BarSlipMaterial2::sendSelf() - not implemented!\n";
	return 0;
}


int
BarSlipMaterial2::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	opserr << "BarSlipMaterial2::recvSelf() - not implemented!\n";
	return 0;
}


Response*
BarSlipMaterial2::setResponse(const char** argv, int argc,
	OPS_Stream& theOutput)
{
	Response* theResponse = 0;

	if ((strcmp(argv[0], "endAxialStress") == 0) ||
		(strcmp(argv[0], "stress") == 0) ||
		(strcmp(argv[0], "tangent") == 0) ||
		(strcmp(argv[0], "endSlip") == 0) ||
		(strcmp(argv[0], "strain") == 0) ||
		(strcmp(argv[0], "endAxialStressSlip") == 0) ||
		(strcmp(argv[0], "stressStrain") == 0) ||
		(strcmp(argv[0], "endAxialStressSlipTangent") == 0) ||
		(strcmp(argv[0], "stressStrainTangent") == 0) ||
		(strstr(argv[0], "axialStrains") != 0) ||
		(strstr(argv[0], "axialStresses") != 0) ||
		(strstr(argv[0], "slips") != 0) ||
		(strstr(argv[0], "bondStresses") != 0)) {

		theOutput.tag("UniaxialMaterialOutput");
		theOutput.attr("matType", this->getClassType());
		theOutput.attr("matTag", this->getTag());

		// end stress
		if (strcmp(argv[0], "endAxialStress") == 0 ||
			strcmp(argv[0], "stress") == 0) {
			theOutput.tag("ResponseType", "sig_L");
			theResponse = new MaterialResponse(this, 1, this->getStress());
		}
		// tangent
		else if (strcmp(argv[0], "tangent") == 0) {
			theOutput.tag("ResponseType", "k");
			theResponse = new MaterialResponse(this, 2, this->getTangent());
		}

		// end slip
		else if (strcmp(argv[0], "endSlip") == 0 ||
			strcmp(argv[0], "strain") == 0) {
			theOutput.tag("ResponseType", "s_L");
			theResponse = new MaterialResponse(this, 3, this->getStrain());
		}

		// end stress-slip
		else if (strcmp(argv[0], "endAxialStressSlip") == 0 ||
			strcmp(argv[0], "stressStrain") == 0) {
			theOutput.tag("ResponseType", "sig_L");
			theOutput.tag("ResponseType", "s_L");
			theResponse = new MaterialResponse(this, 4, Vector(2));
		}

		// end stress-slip-tangent
		else if (strcmp(argv[0], "endAxialStressSlipTangent") == 0 ||
			strcmp(argv[0], "stressStrainTangent") == 0) {
			theOutput.tag("ResponseType", "sig_L");
			theOutput.tag("ResponseType", "s_L");
			theOutput.tag("ResponseType", "k");
			theResponse = new MaterialResponse(this, 5, Vector(3));
		}

		// axial strains
		else if (strcmp(argv[0], "axialStrains") == 0) {
			theOutput.tag("ResponseType", "eps_x");
			theResponse = new MaterialResponse(this, 6, this->getAxialStrains());
		}

		// axial stresses
		else if (strcmp(argv[0], "axialStresses") == 0) {
			theOutput.tag("ResponseType", "sig_x");
			theResponse = new MaterialResponse(this, 7, this->getAxialStresses());
		}

		// slips
		else if (strcmp(argv[0], "slips") == 0) {
			theOutput.tag("ResponseType", "s_x");
			theResponse = new MaterialResponse(this, 8, this->getSlips());
		}

		// bond stresses
		else if (strcmp(argv[0], "bondStresses") == 0) {
			theOutput.tag("ResponseType", "tau_x");
			theResponse = new MaterialResponse(this, 9, this->getBondStresses());
		}

		theOutput.endTag();
	}

	return theResponse;

}

int
BarSlipMaterial2::getResponse(int responseID, Information& matInfo)
{
	Vector vec2(2);
	Vector vec3(3);

	switch (responseID) {
	case 1:
		matInfo.setDouble(this->getStress());
		return 0;

	case 2:
		matInfo.setDouble(this->getTangent());
		return 0;

	case 3:
		matInfo.setDouble(this->getStrain());
		return 0;

	case 4:
		vec2(0) = this->getStress();
		vec2(1) = this->getStrain();

		matInfo.setVector(vec2);
		return 0;

	case 5:
		vec3(0) = this->getStress();
		vec3(1) = this->getStrain();
		vec3(2) = this->getTangent();

		matInfo.setVector(vec3);
		return 0;

	case 6:
		matInfo.setVector(this->getAxialStrains());
		return 0;

	case 7:
		matInfo.setVector(this->getAxialStresses());
		return 0;

	case 8:
		matInfo.setVector(this->getSlips());
		return 0;

	case 9:
		matInfo.setVector(this->getBondStresses());
		return 0;

	default:
		return -1;
	}
}

void
BarSlipMaterial2::Print(OPS_Stream &s, int flag)
{
	s << "BarSlip2 material tag: " << this->getTag() << endln;
	s << "  axialMaterial tag: " << axialMats[0]->getTag() << endln;
	s << "  bondMaterial tag:  " << bondMats[0]->getTag() << endln;
	s << "  bar area:          " << Ab << endln;
	s << "  bar circumference: " << Cb << endln;
	s << "  bar length:        " << L << endln;
	s << "  number of IPs:     " << N << endln;
}
