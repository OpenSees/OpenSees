/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for a mechanical bar-buckling model
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#include <BarBucklingMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void*
OPS_BarBucklingMaterial(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 7) {
		opserr << "WARNING! Insufficient args in uniaxialMaterial BarBuckling" << endln;
		opserr << "want: uniaxialMaterial BarBuckling tag? refMatTag? D? sl? <-imp r?> <-min minEps?> <-iter maxIterNo? maxTol?>" << endln;
		return 0;
	}

	int tags[2];
	double dData[2];
	double imp = 0.001, minEps = 0.0;
	int maxIterNo = 50;
	double maxTol = 1.0e-6;
	bool dbl = false;

	int numData = 2;
	if (OPS_GetIntInput(&numData, tags) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial BarBuckling" << endln;
		return 0;
	}

	numData = 2;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid D or sl for uniaxialMaterial BarBuckling " << tags[0] << endln;
		return 0;
	}

	while (OPS_GetNumRemainingInputArgs() > 1) {
		const char* opt = OPS_GetString();

		if (strcmp(opt, "-imp") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &imp) != 0) {
				opserr << "WARNING! Invalid imp for uniaxialMaterial BarBuckling " << tags[0] << endln;
				return 0;
			}
		}
		else if (strcmp(opt, "-min") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &minEps) != 0) {
				opserr << "WARNING! Invalid minEps for uniaxialMaterial BarBuckling " << tags[0] << endln;
				return 0;
			}
		}
		else if (strcmp(opt, "-iter") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, &maxIterNo) != 0) {
				opserr << "WARNING! Invalid maxIterNo for uniaxialMaterial BarBuckling " << tags[0] << endln;
				return 0;
			}

			numData = 1;
			if (OPS_GetDoubleInput(&numData, &maxTol) != 0) {
				opserr << "WARNING! Invalid maxTol for uniaxialMaterial BarBuckling " << tags[0] << endln;
				return 0;
			}
		}
	}

	// create pointers to reference material models
	UniaxialMaterial* theMat = OPS_getUniaxialMaterial(tags[1]);

	if (theMat == 0) {
		opserr << "WARNING! Reference material model " << tags[1] << " does not exist" << endln;
		opserr << "uniaxialMaterial BarBuckling: " << tags[0] << endln;

		delete theMat;
	}

	// allocate material model
	theMaterial = new BarBucklingMaterial(tags[0], theMat, dData[0], dData[1], imp, minEps, maxIterNo, maxTol);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type BarBuckling\n";
		return 0;
	}

	return theMaterial;
}


BarBucklingMaterial::BarBucklingMaterial(int tag, UniaxialMaterial* mat, double d, double sl, double imp, double eps, int iter, double tol)
	:UniaxialMaterial(tag, MAT_TAG_BarBuckling),
	D(d), LD(sl), L(d * sl), IMP(imp), minBcklEps(eps), maxIterNo(iter), maxTol(tol),
	eps_T(0.0), eps_C(0.0), sig_T(0.0), sig_C(0.0), eps_min(0.0),
	eps_eT(0.0), eps_eC(0.0), eps_qT(0.0), eps_qC(0.0), phi_eT(0.0), phi_eC(0.0),
	N_e(0.0), N_q(0.0), M_e(0.0),
	J_o(3, 3), w_m(0.0), E_tanT(0.0), E_tanC(0.0), mats(0)
{
	if (minBcklEps > 0.0) {
		opserr << "WARNING BarBucklingMaterial::BarBucklingMaterial() ";
		opserr << "- minimum buckling strain is positive - continued using 0.0" << endln;
		minBcklEps = 0.0;
	}

	// compute bar area
	barA = PI * D * D / 4.0;

	//// assign values for IP locations and weights (10-point Gauss-Lobatto quadrature)
	secX[0] = 0.0;
	secX[1] = 0.020116522950;
	secX[2] = 0.065306533725;
	secX[3] = 0.130518762550;
	secX[4] = 0.208680260575;
	secX[5] = 0.291319739425;
	secX[6] = 0.369481237450;
	secX[7] = 0.434693466275;
	secX[8] = 0.479883477050;
	secX[9] = 0.5;

	secL[0] = 0.005555555555 * L;
	secL[1] = 0.033326497700 * L;
	secL[2] = 0.056222335525 * L;
	secL[3] = 0.073010670900 * L;
	secL[4] = 0.081884940275 * L;
	secL[5] = 0.081884940275 * L;
	secL[6] = 0.073010670900 * L;
	secL[7] = 0.056222335525 * L;
	secL[8] = 0.033326497700 * L;
	secL[9] = 0.005555555555 * L;

	secY[0] = -0.5 * D;
	secY[1] = -0.4597669541 * D;
	secY[2] = -0.3693869326 * D;
	secY[3] = -0.2389624749 * D;
	secY[4] = -0.0826394788 * D;
	secY[5] = 0.0826394788 * D;
	secY[6] = 0.2389624749 * D;
	secY[7] = 0.3693869326 * D;
	secY[8] = 0.4597669541 * D;
	secY[9] = 0.5 * D;

	secA[0] = 0.0015564028 * D * D;
	secA[1] = 0.0266733508 * D * D;
	secA[2] = 0.0758342135 * D * D;
	secA[3] = 0.1278425307 * D * D;
	secA[4] = 0.1607925839 * D * D;
	secA[5] = 0.1607925839 * D * D;
	secA[6] = 0.1278425307 * D * D;
	secA[7] = 0.0758342135 * D * D;
	secA[8] = 0.0266733508 * D * D;
	secA[9] = 0.0015564028 * D * D;

	// create arrays to store copies of material models
	mats = new UniaxialMaterial* [N + 1];

	if (mats == 0) {
		opserr << "FATAL BarBucklingMaterial::BarBucklingMaterial() ";
		opserr << "- ran out of memory for material model array of size: " << N + 1 << endln;
		exit(-1);
	}

	// store copies of material models into their arrays
	for (int i = 0; i < N + 1; i++) {
		mats[i] = mat->getCopy();

		if (mats[i] == 0) {
			opserr << "FATAL BarBucklingMaterial::BarBucklingMaterial() ";
			opserr << "- failed to get copy of material model\n";
			exit(-1);
		}
	}

	// initialize material models
	if (this->updateSectionStrains() != 0) {
		opserr << "FATAL BarBucklingMaterial::BarBucklingMaterial() ";
		opserr << "- failed to initialize material models\n";
		exit(-1);
	}

	E_tanT = E_tanC = mats[0]->getInitialTangent();

	//// create the initial Jacobian matrix
	J_o.Zero();

	double N2eps_e, N2phi_e, N2eps_q, M2eps_e, M2phi_e;

	this->getSectionTangents(N2eps_e, N2phi_e, N2eps_q, M2eps_e, M2phi_e);

	// row 1
	J_o(0, 0) = ((eps_min <= minBcklEps ? IMP : 0.0) * L) * N2eps_e + 2.0 * M2eps_e;
	J_o(0, 1) = 0.0;
	J_o(0, 2) = ((eps_min <= minBcklEps ? IMP : 0.0) * L) * N2phi_e + 2.0 * M2phi_e;

	// row 2
	J_o(1, 0) = -N2eps_e * cos((eps_min <= minBcklEps ? IMP : 0.0) * PI);
	J_o(1, 1) = N2eps_q;
	J_o(1, 2) = -N2phi_e * cos((eps_min <= minBcklEps ? IMP : 0.0) * PI);

	// row 3
	J_o(2, 0) = -L / 2.0;
	J_o(2, 1) = -L / 2.0;
	J_o(2, 2) = 0.0;
}


BarBucklingMaterial::BarBucklingMaterial()
	:UniaxialMaterial(0, MAT_TAG_BarBuckling),
	LD(0.0), D(0.0), L(0.0), IMP(0.0), minBcklEps(0.0), maxIterNo(0), maxTol(0.0),
	eps_T(0.0), eps_C(0.0), sig_T(0.0), sig_C(0.0), eps_min(0.0),
	eps_eT(0.0), eps_eC(0.0), eps_qT(0.0), eps_qC(0.0), phi_eT(0.0), phi_eC(0.0),
	N_e(0.0), N_q(0.0), M_e(0.0),
	J_o(3, 3), w_m(0.0), E_tanT(0.0), E_tanC(0.0), barA(0.0), mats(0)
{
	//// assign values for IP locations and weights (10-point Gauss-Lobatto quadrature)
	secX[0] = 0.0;
	secX[1] = 0.020116522950;
	secX[2] = 0.065306533725;
	secX[3] = 0.130518762550;
	secX[4] = 0.208680260575;
	secX[5] = 0.291319739425;
	secX[6] = 0.369481237450;
	secX[7] = 0.434693466275;
	secX[8] = 0.479883477050;
	secX[9] = 0.5;

	secL[0] = 0.005555555555 * L;
	secL[1] = 0.033326497700 * L;
	secL[2] = 0.056222335525 * L;
	secL[3] = 0.073010670900 * L;
	secL[4] = 0.081884940275 * L;
	secL[5] = 0.081884940275 * L;
	secL[6] = 0.073010670900 * L;
	secL[7] = 0.056222335525 * L;
	secL[8] = 0.033326497700 * L;
	secL[9] = 0.005555555555 * L;

	secY[0] = -0.5 * D;
	secY[1] = -0.4597669541 * D;
	secY[2] = -0.3693869326 * D;
	secY[3] = -0.2389624749 * D;
	secY[4] = -0.0826394788 * D;
	secY[5] = 0.0826394788 * D;
	secY[6] = 0.2389624749 * D;
	secY[7] = 0.3693869326 * D;
	secY[8] = 0.4597669541 * D;
	secY[9] = 0.5 * D;

	secA[0] = 0.0015564028 * D * D;
	secA[1] = 0.0266733508 * D * D;
	secA[2] = 0.0758342135 * D * D;
	secA[3] = 0.1278425307 * D * D;
	secA[4] = 0.1607925839 * D * D;
	secA[5] = 0.1607925839 * D * D;
	secA[6] = 0.1278425307 * D * D;
	secA[7] = 0.0758342135 * D * D;
	secA[8] = 0.0266733508 * D * D;
	secA[9] = 0.0015564028 * D * D;
}


BarBucklingMaterial::~BarBucklingMaterial()
{
	// delete the material array
	for (int i = 0; i < N + 1; i++)
		if (mats[i] != 0)
			delete mats[i];

	if (mats != 0)
		delete[] mats;
}


int
BarBucklingMaterial::setTrialStrain(double strain, double strainRate)
{
	eps_T = strain;
	double u_e = eps_T * L;		// trial end displacement

	if (fabs(eps_T - eps_C) <= DBL_EPSILON) {
		sig_T = sig_C;
		E_tanT = E_tanC;

		// update material models anyway
		eps_eT = eps_eC;
		eps_qT = eps_qC;
		phi_eT = phi_eC;

		if (this->updateSectionStrains() != 0) {
			opserr << "WARNING! BarBucklingMaterial::setTrialStrain() - tag: " << this->getTag()
				<< "\nfailed to invoke setSectionStrains()\n";
			return -1;
		}

		return 0;
	}

	// iteratibe solution variables
	static Vector z(3);		// unknowns vector
	static Vector dz(3);	// unknowns vector correction (via Newton-Raphson)
	static Vector R(3);		// residuals vector
	static Matrix J(3, 3);	// Jacobian matrix

	double imp, u_e_int, eps_sum, eps_dif, val1, val2;
	double N2eps_e, N2phi_e, N2eps_q, M2eps_e, M2phi_e;
	double w2eps_e, w2phi_e, w2eps_q, u2eps_e, u2phi_e, u2eps_q;

	double* dL = new double[N];
	double* CS = new double[N];
	double* SN = new double[N];

	double tol;
	int convergence = -1;

	for (int m = 0; m < 1; m++) {	// three methods with combinations of initial and tangent Jacobians

		// initial Jacobian
		J = J_o;

		// initial guesses for {z} elements
		eps_eT = eps_eC;
		eps_qT = eps_qC;
		phi_eT = phi_eC;

		// iterations
		for (int iter = 0; iter < maxIterNo; iter++) {

			// set section strains
			if (this->updateSectionStrains() != 0) {
				opserr << "WARNING! BarBucklingMaterial::setTrialStrain() - tag: " << this->getTag()
					<< "\nfailed to invoke setSectionStrains()\n";
				return -1;
			}

			// get section stresses
			this->updateSectionStresses();

			//// compute {R}
			// integrate section strains over bar's half length to obtain w_m and u_e
			eps_sum = eps_eT + eps_qT;
			eps_dif = eps_eT - eps_qT;
			val1 = phi_eT * L / (2.0 * PI);

			w_m = 0.0;			// mid-length deflection (excluding initial imperfection)
			u_e_int = -L;		// end axial displacement obtained through integral

			for (int i = 0; i < N; i++) {
				dL[i] = (1.0 + 0.5 * (eps_sum + eps_dif * cos(4.0 * PI * secX[i]))) * secL[i];
				val2 = val1 * sin(2.0 * PI * secX[i]);

				CS[i] = cos(val2);
				SN[i] = sin(val2);

				w_m += dL[i] * SN[i];
				u_e_int += 2.0 * dL[i] * CS[i];
			}

			// determine imperfection
			if (minBcklEps == 0.0)
				imp = IMP;
			else
				imp = (eps_min < minBcklEps ? (fmin((eps_min - minBcklEps) / (0.1 * minBcklEps), 1.0) * IMP) : 0.0);	//(fmin((eps_min - minBcklEps) / (0.1 * minBcklEps), 1.0) * IMP)

			// determine {R} elements
			R(0) = (w_m + imp * L) * N_e + 2.0 * M_e;
			R(1) = N_q - N_e * cos(val1 + imp * PI);
			R(2) = u_e - u_e_int;

			//// check convergence and update solution if needed
			if (iter < maxIterNo / 2)
				tol = maxTol;
			else if (iter < 3 * maxIterNo / 4)
				tol = 10. * maxTol;
			else
				tol = 100. * maxTol;

			if (fabs(R(0)) <= fmax(0.01 * tol * fabs(M_e), tol)
				&& fabs(R(1)) <= fmax(0.01 * tol * fabs(N_e), tol)
				&& fabs(R(2)) <= fmax(0.01 * tol * fabs(u_e), tol)) {	// converged
				convergence = 0;
				break;
			}
			else if (iter < maxIterNo - 1) {	// not converged
				// form {z}
				z(0) = eps_eT;
				z(1) = eps_qT;
				z(2) = phi_eT;

				// update [J]
				if (m == 0 || (m == 1 && iter > 0)) {
					// compute section tangents
					this->getSectionTangents(N2eps_e, N2phi_e, N2eps_q, M2eps_e, M2phi_e);

					// compute other derivatives
					w2eps_e = w2eps_q = w2phi_e = 0.0;
					u2eps_e = u2eps_q = u2phi_e = 0.0;

					for (int i = 0; i < N; i++) {
						val1 = (1.0 + cos(4.0 * PI * secX[i])) * secL[i];
						val2 = (1.0 - cos(4.0 * PI * secX[i])) * secL[i];

						w2eps_e += 0.5 * val1 * SN[i];
						w2eps_q += 0.5 * val2 * SN[i];
						w2phi_e += (L / (2.0 * PI)) * dL[i] * sin(2.0 * PI * secX[i]) * CS[i];

						u2eps_e += val1 * CS[i];
						u2eps_q += val2 * CS[i];
						u2phi_e -= (L / PI) * dL[i] * sin(2.0 * PI * secX[i]) * SN[i];
					}

					// assemble [J]
					J(0, 0) = w2eps_e * N_e + (w_m + imp * L) * N2eps_e + 2.0 * M2eps_e;
					J(0, 1) = w2eps_q * N_e;
					J(0, 2) = w2phi_e * N_e + (w_m + imp * L) * N2phi_e + 2.0 * M2phi_e;

					val1 = phi_eT * L / (2.0 * PI) + imp * PI;

					J(1, 0) = -N2eps_e * cos(val1);
					J(1, 1) = N2eps_q;
					J(1, 2) = -N2phi_e * cos(val1) + N_e * (L / (2.0 * PI)) * sin(val1);

					J(2, 0) = -u2eps_e;
					J(2, 1) = -u2eps_q;
					J(2, 2) = -u2phi_e;
				}

				// Newton-Raphson correction
				if (J.Solve(R, dz) < 0) {
					opserr << "WARNING! BarBucklingMaterial::setTrialStrain() - tag: " << this->getTag()
						<< "\ncould not invert Jacobian\n";
					return -1;
				}

				double gamma = 1.0;

				for (int i = 0; i < 3; i++)
					if (fabs(dz(i)) > 0.0001)
						gamma = fmin(gamma, 0.0001 / fabs(dz(i)));

				z -= gamma * dz;

				/*if (iter > maxIterNo - 6) {
					opserr << "iter " << iter << endln;
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++)
							opserr << J(i, j) << '\t';
						opserr << endln;
					}
					opserr << "\nR: " << R(0) << '\t' << R(1) << '\t' << R(2) << endln;
					opserr << "\nz: " << z(0) << '\t' << z(1) << '\t' << z(2) << endln;
					opserr << endln;
				}*/

				// extract corrected unknowns
				eps_eT = z(0);
				eps_qT = z(1);
				phi_eT = z(2);
			}
		}

		if (convergence == 0)
			break;
	}

	if (convergence != 0)
		opserr << "WARNING! BarBucklingMaterial::setTrialStrain() - tag: " << this->getTag()
			<< "\nNo convergence after " << maxIterNo << " iterations - error norm: " << R.Norm() << endln;
	else {
		// trial average stress
		sig_T = N_e / barA;
		
		// trial tangent
		/*double F_tan = 0.0;

		static Matrix K_sec(3, 3);

		K_sec(0, 0) = N2eps_e;
		K_sec(0, 1) = N2phi_e;
		K_sec(0, 2) = 0.0;
		K_sec(1, 0) = M2eps_e;
		K_sec(1, 1) = M2phi_e;
		K_sec(1, 2) = 0.0;
		K_sec(2, 0) = 0.0;
		K_sec(2, 1) = 0.0;
		K_sec(2, 2) = N2eps_q;

		static Vector B(3);
		B(0) = 1.0;
		B(1) = -(w_m + imp * L) / 2.0;
		B(2) = cos(phi_eT * L / (2.0 * PI) + imp * PI);

		static Vector FB(3);

		if (K_sec.Solve(B, FB) < 0) {
			opserr << "WARNING! BarBucklingMaterial::setTrialStrain() - tag: " << this->getTag()
				<< "\ncould not invert [K_sec]\n";
			return -1;
		}

		F_tan = u2eps_e * FB(0) + u2phi_e * FB(1) + u2eps_q * FB(2);
		E_tanT = (L / barA) / F_tan;*/

		E_tanT = (sig_T - sig_C) / (eps_T - eps_C);

		if (isnan(E_tanT))
			E_tanT = E_tanC;
	}

	// delete dynamic memory
	if (dL)
		delete[] dL;

	if (CS)
		delete[] CS;

	if (SN)
		delete[] SN;

	return convergence;
}


double
BarBucklingMaterial::getStrain(void)
{
	return eps_T;
}


double
BarBucklingMaterial::getStress(void)
{
	return sig_T;
}


double
BarBucklingMaterial::getTangent(void)
{
	return E_tanT;
}


double
BarBucklingMaterial::getInitialTangent(void)
{
	return (mats[0]->getInitialTangent());
}


int
BarBucklingMaterial::commitState(void)
{
	int err0, err = 0;

	// invoke commitState() on each MaterialModel object
	for (int i = 0; i < N + 1; i++) {

		err0 = mats[i]->commitState();

		if (err0 != 0) {
			opserr << "WARNING BarBucklingMaterial::commitState() ";
			opserr << "- MaterialModel failed to commitState():";
			mats[i]->Print(opserr);

			err += err0;
		}
	}

	// save minimum reached strain
	eps_min = fmin(eps_T,eps_min);

	// commit state variables
	eps_C = eps_T;
	sig_C = sig_T;
	eps_eC = eps_eT;
	eps_qC = eps_qT;
	phi_eC = phi_eT;
	E_tanC = E_tanT;

	return err;
}


int
BarBucklingMaterial::revertToLastCommit(void)
{
	int err0, err = 0;

	// invoke revertToLastCommit() on each MaterialModel object
	for (int i = 0; i < N + 1; i++) {

		err0 = mats[i]->revertToLastCommit();

		if (err0 != 0) {
			opserr << "WARNING BarBucklingMaterial::revertToLastCommit() ";
			opserr << "- MaterialModel failed to revertToLastCommit():";
			mats[i]->Print(opserr);

			err += err0;
		}
	}

	// state variables
	eps_T = eps_C;
	sig_T = sig_C;
	eps_eT = eps_eC;
	eps_qT = eps_qC;
	phi_eT = phi_eC;
	E_tanT = E_tanC;

	return err;
}


int
BarBucklingMaterial::revertToStart(void)
{
	int err0, err = 0;

	// invoke revertToStart() on each MaterialModel object
	for (int i = 0; i < N + 1; i++) {

		err0 = mats[i]->revertToStart();

		if (err0 != 0) {
			opserr << "WARNING BarBucklingMaterial::revertToStart() ";
			opserr << "- MaterialModel failed to revertToStart():";
			mats[i]->Print(opserr);

			err += err0;
		}
	}

	// zero state variables
	eps_min = eps_min = 0.0;
	eps_T = eps_C = 0.0;
	sig_T = sig_C = 0.0;
	eps_eT = eps_eC = 0.0;
	eps_qT = eps_qC = 0.0;
	phi_eT = phi_eC = 0.0;
	E_tanC = E_tanT = this->getInitialTangent();

	return err;
}


UniaxialMaterial *
BarBucklingMaterial::getCopy(void)
{
	BarBucklingMaterial *theCopy = new BarBucklingMaterial(this->getTag(), mats[0], D, LD, IMP, minBcklEps, maxIterNo, maxTol);

	theCopy->eps_T = eps_T;
	theCopy->eps_C = eps_C;
	theCopy->sig_T = sig_T;
	theCopy->sig_C = sig_C;
	theCopy->eps_min = eps_min;
	theCopy->eps_eT = eps_eT;
	theCopy->eps_eC = eps_eC;
	theCopy->eps_qT = eps_qT;
	theCopy->eps_qC = eps_qC;
	theCopy->phi_eT = phi_eT;
	theCopy->phi_eC = phi_eC;
	theCopy->E_tanT = E_tanT;
	theCopy->E_tanC = E_tanC;

	return theCopy;
}


int
BarBucklingMaterial::sendSelf(int cTag, Channel &theChannel)
{
	opserr << "BarBucklingMaterial::sendSelf() - not implemented!\n";
	return 0;
}


int
BarBucklingMaterial::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	opserr << "BarBucklingMaterial::recvSelf() - not implemented!\n";
	return 0;
}


Response*
BarBucklingMaterial::setResponse(const char** argv, int argc,
	OPS_Stream& theOutput)
{
	Response* theResponse = 0;

	if ((strcmp(argv[0], "stress") == 0) ||
		(strcmp(argv[0], "tangent") == 0) ||
		(strcmp(argv[0], "strain") == 0) ||
		(strcmp(argv[0], "stressStrain") == 0) ||
		(strcmp(argv[0], "stressStrainTangent") == 0) ||
		(strstr(argv[0], "midDeflection") == 0) ||
		(strstr(argv[0], "sectionStrains") == 0) ||
		(strstr(argv[0], "sectionStresses") == 0)) {

		theOutput.tag("UniaxialMaterialOutput");
		theOutput.attr("matType", this->getClassType());
		theOutput.attr("matTag", this->getTag());

		// average stress
		if (strcmp(argv[0], "stress") == 0) {
			theOutput.tag("ResponseType", "sig");
			theResponse = new MaterialResponse(this, 1, this->getStress());
		}
		// tangent
		else if (strcmp(argv[0], "tangent") == 0) {
			theOutput.tag("ResponseType", "tan");
			theResponse = new MaterialResponse(this, 2, this->getTangent());
		}

		// average strain
		else if (strcmp(argv[0], "strain") == 0) {
			theOutput.tag("ResponseType", "eps");
			theResponse = new MaterialResponse(this, 3, this->getStrain());
		}

		// average stress-strain
		else if (strcmp(argv[0], "stressStrain") == 0) {
			theOutput.tag("ResponseType", "sig");
			theOutput.tag("ResponseType", "eps");
			theResponse = new MaterialResponse(this, 4, Vector(2));
		}

		// average stress-slip-tangent
		else if (strcmp(argv[0], "stressStrainTangent") == 0) {
			theOutput.tag("ResponseType", "sig");
			theOutput.tag("ResponseType", "eps");
			theOutput.tag("ResponseType", "tan");
			theResponse = new MaterialResponse(this, 5, Vector(3));
		}

		// transverse deflection at bar mid-length
		else if (strcmp(argv[0], "midDeflection") == 0) {
			theOutput.tag("ResponseType", "w_mid");
			theResponse = new MaterialResponse(this, 6, w_m);
		}

		// representative section strains
		else if (strcmp(argv[0], "sectionStrains") == 0) {
			theOutput.tag("ResponseType", "eps_e");
			theOutput.tag("ResponseType", "phi_e");
			theOutput.tag("ResponseType", "eps_q");
			theResponse = new MaterialResponse(this, 7, Vector(3));
		}

		// representative section strains
		else if (strcmp(argv[0], "sectionStresses") == 0) {
			theOutput.tag("ResponseType", "N_e");
			theOutput.tag("ResponseType", "M_e");
			theOutput.tag("ResponseType", "N_q");
			theResponse = new MaterialResponse(this, 8, Vector(3));
		}

		theOutput.endTag();
	}

	return theResponse;

}

int
BarBucklingMaterial::getResponse(int responseID, Information& matInfo)
{
	static Vector vec2(2);
	static Vector vec3(3);

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
		matInfo.setDouble(w_m);
		return 0;

	case 7:
		vec3(0) = eps_eT;
		vec3(1) = phi_eT;
		vec3(2) = eps_qT;

		matInfo.setVector(vec3);
		return 0;

	case 8:
		vec3(0) = N_e;
		vec3(1) = M_e;
		vec3(2) = N_q;

		matInfo.setVector(vec3);
		return 0;

	default:
		return -1;
	}
}

void
BarBucklingMaterial::Print(OPS_Stream &s, int flag)
{
	s << "BarBuckling material tag: " << this->getTag() << endln;
	s << "  reference material tag: " << mats[0]->getTag() << endln;
	s << "  bar diameter:           " << D << endln;
	s << "  bar slenderness:        " << LD << endln;
	s << "  middle imperfection:    " << IMP << endln;
	s << "  min. buckling strain:   " << minBcklEps << endln;
	s << "  number of IPs:          " << N << endln;
}


//// Private Methods
int 
BarBucklingMaterial::updateSectionStrains(void)
{
	// NOTE: first N material models are for end/middle fiber sections,
	// while the (N+1)th model is for quarter-length section (with zero curvature)
	int err0, err = 0;

	// invoke setTrialStrain() on each MaterialModel object
	for (int i = 0; i < N; i++) {		// fiber sections
		err0 = mats[i]->setTrialStrain(eps_eT - secY[i] * phi_eT, 0.0);

		if (err0 != 0) {
			opserr << "WARNING BarBucklingMaterial::setSectionStrains() ";
			opserr << "- MaterialModel failed to setTrialStrain():";
			mats[i]->Print(opserr);

			err += err0;
		}
	}

	err += mats[N]->setTrialStrain(eps_qT, 0.0);		// quarter-length section (with constant strain over area)

	return err;
}


void 
BarBucklingMaterial::updateSectionStresses(void)
{
	// NOTE: first N material models are for end/middle fiber sections,
	// while the (N+1)th model is for quarter-length section (with zero curvature)

	// fiber section
	N_e = M_e = 0.0;
	double df;

	for (int i = 0; i < N; i++) {		// fiber sections
		df = secA[i] * mats[i]->getStress();
		N_e += df;
		M_e -= secY[i] * df;
	}

	// quarter-length section (constant stress)
	N_q = barA * mats[N]->getStress();
}


void
BarBucklingMaterial::getSectionTangents(double& N2eps_e, double& N2phi_e, double& N2eps_q, double& M2eps_e, double& M2phi_e)
{
	// NOTE: first N material models are for end / middle fiber sections,
	// while the (N+1)th model is for quarter-length section (with zero curvature)

	// fiber section
	N2eps_e = N2phi_e = M2phi_e = 0.0;
	double dk;

	for (int i = 0; i < N; i++) {		// fiber sections
		dk = secA[i] * mats[i]->getTangent();
		N2eps_e += dk;
		N2phi_e -= secY[i] * dk;
		M2phi_e += secY[i] * secY[i] * dk;
	}

	M2eps_e = N2phi_e;

	// quarter-length section (constant stress)
	N2eps_q = barA * mats[N]->getTangent();
}