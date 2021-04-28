/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                                                                                                                                           
// Written: vavgen  (Evangelos Avgenakis)
// Some preexisting code from basic OpenSees elements was used
// Created: July 2020
// Revision: A
//
// Description: This element can be used to describe the rocking motion of 2d deformable bodies,
// 		either elastic or inelastic, under static or dynamic loading. Apart from the deformability
// 		along the length of the member, the element is able to account for deformability near the
// 		contact area, where nonlinear stress distributions develop and sections do not remain plane.
// 		Furthermore, the element is able to account for constraints along the length of the rocking
// 		member imposed by other structural members, as well as sliding and upthrow.
//
// References:
// 	1. Avgenakis E. and Psycharis I.N. (2017) “Modeling of Rocking Elastic Flexible Bodies under Static
//  		Loading Considering the Nonlinear Stress Distribution at Their Base.” Journal of Structural
//		Engineering 143(7): 04017051.
//	2. Avgenakis, E. and Psycharis, I. N. (2019) “Determination of the nonlinear displacement distribution
//		of the semi-infinite strip–Application to deformable rocking bodies.” International Journal
//		of Solids and Structures, 170, 22-37.
//	3. Avgenakis E. and Psycharis I.N. (2020) “Modeling of inelastic rocking bodies under cyclic loading.”
//		Journal of Engineering Mechanics 146(4): 04020020.
// 	4. Avgenakis E. and Psycharis I.N. (2020) “An integrated macroelement formulation for the dynamic
//		response of inelastic deformable rocking bodies.” Earthquake Engineering and Structural Dynamics,
//      49(11), 1072-1094.

#ifndef RockingBC_h
#define RockingBC_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <fstream>
#include <Renderer.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>

// #include "Eigen/Dense" // REMOVE EIGEN

using Vec = std::vector<double>;
using Vecint = std::vector<int>;
using VecVecOS = std::vector< Vector >;
using VecVecXd = std::vector<Vector>;
using VecVec = std::vector< std::vector<double> >;
using VecVecint = std::vector< std::vector<int> >;
using VecVecVec = std::vector< std::vector< std::vector<double> > >;
using VecMatOS = std::vector< Matrix >;

class Channel;
class Information;
class CrdTransf;
class Response;
class Renderer;

class RockingBC : public Element
{
public:
	RockingBC();
	RockingBC(int tag, int Nd1, int Nd2, int Nw,
		double e, double Nu, double Sy, double bb, double ww, double Mu, double Convlim, int Maxtries, double Af, double Aflim, double Convlimmult, int Usecomstiff, int Useshear, int Blevery,
		double NLimN, double NLimT, double DTlim, int ErrorifNexceeds, int UseUelNM);
	~RockingBC();

	const char *getClassType(void) const { return "RockingBC"; };

	int getNumExternalNodes(void) const;
	const ID &getExternalNodes(void);
	Node **getNodePtrs(void);

	int getNumDOF(void);
	void setDomain(Domain *theDomain);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int initialize(Node *nodeIPointer, Node *nodeJPointer);
	int update(void);
	double getInitialLength(void);
	int state_determination(void);
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getDamp(void);
	const Matrix &getMass(void);

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);

	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &info);

	const Matrix & inverse3x3matrix(Matrix &A) const;
	int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes = 0, int numModes = 0);
private:

	Node *theNodes[2];

	ID  connectedExternalNodes;

	//CrdTransf *theCoordTransf;

	int compElemtLengthAndOrient(void);
	void compTransfMatrixLocalGlobal(Matrix &Tlg);
	Node *nodeIPtr, *nodeJPtr;

	const Vector &getLocalTrialDisp(void);
	const Vector &getLocalIncrDisp(void);
	const Vector &getLocalIncrDeltaDisp(void);

	const Vector &getGlobalResistingForce(const Vector &localForce);
	const Matrix &getGlobalStiffMatrix(const Matrix &localStiff);
	const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);

	const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

	double E{}, sy{}, G{}, ey{};
	double B{}, w{}, A{}, I{}, b{};
	double alpha{}, nu{}, mu{};
	double EI{}, EA{}, GA{};
	int Nw{};

	double committedtime{0};
	double curtime{};
	double Dt{};
	double Dtprev{};
	bool is_analysis_dynamic(void);
	double getDt(void);
	double beta_Dt{};

	Matrix K = Matrix(6, 6);
	Vector P = Vector(6);
	Matrix Tlg = Matrix(6, 6);
	Vector pg = Vector(6);
	Matrix kg = Matrix(6, 6);

	Matrix fe = Matrix(6, 6);
	Matrix ke = Matrix(6, 6);
	Matrix kecommit = Matrix(6, 6);
	Matrix kepr = Matrix(6, 6);
	Vector Fe = Vector(6);
	Vector Fecommit = Vector(6);
	Vector Fepr = Vector(6);

	double cosTheta, sinTheta; // direction cosines of undeformed element wrt to global system 
	double L;                  // undeformed element length

	Vector ue = Vector(6);                 // local displacements
	Vector uecommit;           // commited local displacements
	Vector uepr;               // previous local displacements
	Vector due = Vector(6);
	Vector Due = Vector(6);

	Vector W;
	Vector Winit;
	Vector Wpr;
	Vector Wcommit;
	Vector DW;

	Vector Yw, Yup, Up, Yup_com, Up_com, Ys_com, S_com;
	Matrix dW_due,dW_due_com,dW_due_pr;
	VecVec Ysi, Si, Yupi, Upi, Ysi_com, Si_com, Yupi_com, Upi_com, Ua_pos;
	VecMatOS dYsi_dW{}; VecMatOS dSi_dW{};

	Vec ysi_new, si_new, yupi_new, upi_new;

	Vector Ys, S, Ud;
	Matrix dYs_dW, dS_dW;
	Matrix dF_due;
	Matrix dUd_dW, dUd_due;

	Matrix DFe_Due = Matrix(6, 6);

	double N{ 0. }, M{ 0. }, Q{ 0. }, t{ 0. };

	Vector Nints{}, Mints{};
	Matrix dNints_dW{}, dMints_dW{};

	double sL{ 0. };
	double sLpr{ 0. };
	double sLcommit{ 0. };

	bool noshear{ true };

	Vector dueV = Vector(6);
	Vector ueV = Vector(6);
	Vector DueV = Vector(6);

	Matrix UN{}, UM{};

	std::ofstream Yup_file;
	std::ofstream Up_file;
	std::ofstream Ys_file;
	std::ofstream S_file;

	int tagtag = 0;

	int maxtries{};
	double af{};
	double convlim{};
	double aflim{};
	int usecomstiff{};
	int useshear{};
	double convlimmult{};
	int blevery{};
	double NlimN{};
	double NlimT{};
	double Dtlim{};
	int errorifNexceeds{};
	double Dtmax{}, DtmaxN{}, DtmaxT{};
	int useUelNM{};

	int comcount = 0;
	int dyncount = 0;
	bool hasreverted = false;
	bool isdynamic = false;

	Vector DWcommit;
	double Dtcommit;

	Matrix fr = Matrix(3,3);

	std::vector< std::vector<int> > Ys_cats{};
	Vecint Ys_cats_dist;
	Vector Upl, Ua, Ec, El;
	Matrix dUa_dW;
	Vector Youter = Vector(2);
	Matrix dYouter_dW;
	Vector dN_dW, dM_dW, dQ_dW, dt_dW;
	Vector dt_due = Vector(6);
	Vector utn = Vector(2);
	Vector ut = Vector(2);
	Vector urf = Vector(2);
	Matrix dutn_dYouter = Matrix(2, 2);
	Matrix dut_due = Matrix(2, 6);
	Matrix durf_due = Matrix(2, 6);
	Matrix dut_dW, durf_dW;
	Vector Uel;
	Matrix dUel_dW;

	Vector un = Vector(3);
	Matrix dun_dW;
	Matrix dun_due = Matrix(3, 6);
	Vector ues = Vector(6);
	Matrix dues_dW;
	Matrix dues_due = Matrix(6, 6);
	Vector dsL_dW;
	Vector dsL_due = Vector(6);
	Vector dQ_due = Vector(6);

	Vector dw1_due;
	Vector dw2_due;
	Vector dr_due;

	double k1{}, k2{};
	double w1{}, w2{}, r{}, lim1{}, lim2{}, cval{}, gQ{};

	Vector dlim1_dW, dlim2_dW, dgQ_dW;
	Vector dlim1_due = Vector(6);
	Vector dlim2_due = Vector(6);

	Matrix Tn = Matrix(3, 6);
	Matrix dun_dues = Matrix(3, 6);
	Matrix dutn_dW;

	double th2{};
	Vector dth2_dW;
	Vector dth2_due = Vector(6);

	Matrix frr = Matrix(2, 2);
	Vector urth = Vector(2);

	Matrix durth_due, durth_dW;
	Vector Fn2 = Vector(2);
	Matrix dFn2_dW;

	Matrix CC;
	Matrix BB = Matrix(2, 2);
	Matrix CB;

	Vector utar = Vector(2);
	Matrix dutar_due;
	Matrix dutar_dW;

	Vector Ut;
	Matrix dUt_dW;
	Matrix dUt_due;

	Vector Urf;
	Matrix dUrf_dW;
	Matrix dUrf_due;

	Vector Utar;
	Matrix dUtar_dW;
	Matrix dUtar_due;

	Vector Fn = Vector(3);
	Vector Fn_com = Vector(3);
	Vector Fntot = Vector(3);
	Matrix dFn_dW;
	Matrix dFn_due = Matrix(3, 6);
	Matrix dFntot_dW;
	Matrix dFntot_due = Matrix(3, 6);
	Matrix TF1 = Matrix(3, 3);
	Matrix TF;
	Vector FnNN = Vector(3);
	Vector FnVec = Vector(3);
	Vector FnVec_com = Vector(3);

	Vector FnD = Vector(3);
	Vector FnD_com = Vector(3);
	Matrix dFnD_dW;
	Matrix dFnD_due = Matrix(3, 6);

	Vector Fnntot;
	Matrix dFnntot_dW;
	Matrix dFnntot_due;
	Vector Fes = Vector(6);
	Matrix dFes_dW;
	Matrix dFes_due = Matrix(6, 6);

	Vector FeV = Vector(6);
	Matrix dFeV_dW;
	Matrix dFeV_due = Matrix(6, 6);

	Vector due5_due = Vector(6);

	Vector dNtot_dW;
	Vector dPA_dW;
	Vector dPB_dW;
	Vector dPA_due;

	Vector Ks, Ks_com, Ydks, Kup, Kup_com, DS, Dks, DDKs;
	Matrix dKs_dW, dYdks_dW, dDks_dW, dDDKs_dW, dDS_dW;
	Vec rnotfound{};
	std::vector<int> rfoundi{};
	std::vector<int> ifound{};
	std::vector<int> inotfound{};
	int rnfi{ 0 }; int ifi{ 0 };
	Vector rnotfoundvec;
	Matrix Unf;
	Matrix dUnf_dR;
	Matrix UB{};
	Matrix dUB_dR{};
	Vec UB_R{};
	Matrix UBnew{};
	Matrix dUBnew_dR{};
	Vec UBnew_R{};
	Vector Im1, Jm1;
	Vector Uel_com;

	double Ntot{}, N_com{}, Q_com{}, PA{}, PB{};
	double ND_com{}, QD_com{};

	bool usespecslidmode{ true };
	int slidmode_com{ 0 }, slidmode{ 0 }, newslidmode{ 0 }, slidmode_init{ 0 };
	Vecint slidingmodes{};
	Vecint slidingmodes_try{};

	double Fst{ 0 };
	double forceratioN{ 0 };
	double forceratioT{ 0 };
	double forceratioNmax{ 0 };
	double forceratioTmax{ 0 };
	int triesfromcommitstate{ 0 };

	Vector find_in_dist(const Vector& X, const Vector& Y, const Vector& Xf);
	void simplify_dist_up(const Vector& X, const Vector& Y, const Vector& Xw, Vector& Xnew, Vector& Ynew);

	void W_to_ua_upl();
	void W_to_ua_upl_K();
	void Youter_calc();

	void NM_calc();
	void NM_calc_YS();
	void NM_calc_Fncom();

	void fr_calc();
	void sL_Q_t_calc();
	void un_calc();
	void ut_calc();
	void urf_calc();

	void Uel_NM_calc();
	void Uel_K_calc();

	void disp_comb();
	void forces();

	void WZ_solve();

	int NL_solve_dyn();
	void writedbgfile();
	
	// SIS_funcs
	
	const double SISfunclim{ 1.0e-15 };
	
	double YMXLOGYMX(double y, double p);
	double OMXYLOGOMXYOXY(double xy);
	double J2(double yp);
	double OMXATANYMOOXMO(double y, double p);
	double OMYLOGSQ(double y, double p);

	double I_FA(double y, double p);
	double J_FA(double y, double p);
	double I_FB(double y, double p);
	double J_FB(double y, double p);
	double I_FP(double y, double p);
	double J_FP(double y, double p);
	double I_FP_alt(double y, double p);

	double I_calc(double y, double r);
	double J_calc(double y, double r);

	double FAa(double y, double p);
	double dFAa_dp(double y, double p);

	double FA(double y, double p);
	double FB(double y, double p);
	double FP(double y, double p);

	double D_FA(double y, double p);
	double D_FB(double y, double p);
	double D_FP(double y, double p);

	double I_FAb(double y, double p);
	double J_FAb(double y, double p);

	double Ib_calc(double y, double r);
	double Jb_calc(double y, double r);

	double pImJ_FA(double y, double p);
	double pImJ_FB(double y, double p);
	double pImJ_FP(double y, double p);

	double pImJ_calc(double y, double r);

	double pImJ_FA_nochecks(double y, double p);
	double pImJ_FB_nochecks(double y, double p);

	void Imat_calc(const Vector& Y, const Vector& R, Matrix& Imat);
	void Jmat_calc(const Vector& Y, const Vector& R, Matrix& Jmat);
	void Im1_calc(const Vector& Y, Vector& Im1);
	void Jm1_calc(const Vector& Y, Vector& Jm1);

	void Imatb_calc(const Vector& Y, const Vector& R, Matrix& Imat);
	void Jmatb_calc(const Vector& Y, const Vector& R, Matrix& Jmat);
	void Im1b_calc(const Vector& Y, Vector& Im1);
	void Jm1b_calc(const Vector& Y, Vector& Jm1);

	void pImJmat_calc(const Vector& Y, const Vector& R, Matrix& pImJmat);

	void Usgm_trapz(const Vector& Yw, Matrix& Usgm);
	void triangle_dispslope_disps(const Vector& R, const Vector& Y, Matrix& U, Matrix& dU_dR);
	void triangle_dispslope_disps_givenMat1(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR);
	void triangle_dispslope_disps_2(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR);

	void UNM_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U);
	void UNM_rect(const Vector& R, const Vector& Yw, Matrix& U);
	void UNM_calc(const Vector& Yw, Matrix& UN, Matrix& UM);

	void UNMb_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U);
	void UNMb_rect(const Vector& R, const Vector& Yw, Matrix& U);
	void UNMb_calc(const Vector& Yw, Matrix& UN, Matrix& UM);
	
	// SCfuncs
	
	double inline SC_A(double x) { return 2.436222252877402 - 2.3818059387327604*x + 0.7078998718614156*x*x; };
	double inline SC_DA(double x) { return -2.3818059387327604 + 2 * 0.7078998718614156*x; };

	double inline SC_B(double x) { return 0.6982001887951753 - 1.098308073905204*x + 1.9266756798514126*x*x + -1.1270666845181774*x*x*x + 0.688867046041808*x*x*x*x; };
	double inline SC_DB(double x) { return -1.098308073905204 + 2 * 1.9266756798514126*x + 3 * (-1.1270666845181774)*x*x + 4 * 0.688867046041808*x*x*x; };

	double inline SC_C(double x) { return 0.8134604447686402*pow(1. - x, 3.770057533864266) + 1.; };
	double inline SC_DC(double x) { return -0.8134604447686402*3.770057533864266*pow(1. - x, 2.770057533864266); };

	double inline SC_D(double x) { return (1. - x)*(2.340417693163326 - 1.9592356132890616*x + 0.8914260492531663*x*x); };
	double inline SC_DD(double x) { return -2.340417693163326 - (x - 1.)*(-1.9592356132890616 + 2 * 0.8914260492531663*x) + 1.9592356132890616*x - 0.8914260492531663*x*x; };

	double inline SC_E(double x) { return 1.4043226196463283 + 0.1302424508017461*pow(1. - x, 3.6564163357661053) - 0.0549296131209048*x; };
	double inline SC_DE(double x) { return -0.0549296131209048 - 0.1302424508017461*3.6564163357661053*pow(1. - x, 2.6564163357661053); };

	double inline SC_F(double x) { return 0.4343458286281541*x + 3.107476490749382*x*x + (-6.967836976078876)*x*x*x + 6.501720103798543*x*x*x*x + (-2.284276614857206)*x*x*x*x*x; };
	double inline SC_DF(double x) { return 0.4343458286281541 + 2 * 3.107476490749382*x + 3 * (-6.967836976078876)*x*x + 4 * 6.501720103798543*x*x*x + 5 * (-2.284276614857206)*x*x*x*x; };

	void Dt_calc(const Vector& P, double& d, Vector& dddP);
	void Rt_calc(const Vector& P, double& th, Vector& dthdP);

	void se_shear_1der(const Vector& Youter, Vector& Ut, Matrix& dUt_dYouter);
	
	// region
	
	void Up_interval_split(const Vector& Yup, const Vector& Up, const Vector& Yw,
		VecVec& Yup_ints, VecVec& Up_ints);
		
	Vector interval_join(const VecVec& X_ints);
	Matrix interval_join(const VecMatOS& X_ints);
	Vector array_join(const VecVec& X_ints);
	Matrix array_join(const VecMatOS& X_ints);

	void commony(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB);

	void interval_interior(double wl, double wr, double ey, double dy, const Vec& up_com, const Vec& yup_com,
		const Vec& ys_com, const Vec& s_com, double beta_Dt,
		Vec& ys_new, Vec& s_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new,
		Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr,
		Vec& ua_pos);

	void NM_calc_int(const Vec& Ys, const Matrix& dYs_dW, const Vec& S, const Matrix& dS_dW, double& N, double& M, Vector& dN_dW, Vector& dM_dW);

	void interval_dists(const Vector& Yw, const Vector& W, const VecVec& Yupi_com, const VecVec& Upi_com, const VecVec& Ysi_com, const VecVec& Si_com, double ey, double beta_Dt,
		VecVec& Ysi, VecVec& Si, VecVec& Yupi_new, VecVec& Upi_new,
		VecVecint& Ys_cats, Vector& Nints, Vector& Mints, Matrix& dNints_dW, Matrix& dMints_dW, VecVec& Ua_pos, VecMatOS& dYsi_dW, VecMatOS& dSi_dW);

	void critpoints(const Vec& y, const Vec& s, int rinit, int rend, Vecint& cp);

	void int_bilin(const Vecint& ys_cats, const Vec& ys, const Vec& s, const Vec& yup, const Vec& up, const Vec& ua_pos, double ey,
		Vec& ys_new, Vec& s_new, Vec& yup_new, Vec& up_new);

	void Up_interval_split_K(const Vector& Yup, const Vector& Up, const Vector& Kup, const Vector& Yw,
		VecVecOS& Yup_ints, VecVecOS& Up_ints, VecVecOS& Kup_ints);

	void commony_K(const Vector& ya, const Vector& fa, const Vector& ka, const Vector& yb, const Vector& fb, const Vector& kb, Vec& Y, Vec& FA, Vec& FB, Vec& KA, Vec& KB);

	void interval_interior_K(double wl, double wr, double ey, double dy, const Vector& up_com, const Vector& yup_com, const Vector& kup_com,
		const Vector& ys_com, const Vector& s_com, const Vector& ks_com, double beta_Dt,
		Vec& ys_new, Vec& s_new, Vec& ks_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new, Vec& kup_new,
		Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr, Vec& dks_new_dwl, Vec& dks_new_dwr,
		Vec& ydks, Vec& dks, Vec& dydks_dwl, Vec& dydks_dwr, Vec& ddks_dwl, Vec& ddks_dwr,
		Vec& ds, Vec& dds_dwl, Vec& dds_dwr);

	void interval_dists_K(const Vector& Yw, const Vector& W, const Vector& Yup_com, const Vector& Up_com, const Vector& Kup_com, const Vector& Ys_com, const Vector& S_com, const Vector& Ks_com, double ey, double beta_Dt,
		Vector& Ys, Vector& S, Vector& Ks, Vector& Yup_new, Vector& Up_new, Vector& Kup_new,
		Matrix& dYs_dW, Matrix& dS_dW, Matrix& dKs_dW, Vecint& Ys_cats, Vector& Ydks, Vector& Dks, Matrix& dYdks_dW, Matrix& dDks_dW, Vector& DS, Matrix& dDS_dW);

	void Ys_cats_dist_calc(const VecVecint& Ys_cats, Vecint& Ys_cats_dist);

	// bilin

	void commony_BL(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB);
	bool distintersec(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q);
	bool twobilinintersec(double y1, double y2, double p1, double p2, double q1, double q2, double yp, double p0, double yq, double q0);
	void NM_BL(const Vec& Y, const Vec& S, double& N, double& M, double& Nd, double& Md);
	bool bilinable(double Nd, double Md, double y1, double y2, double BILINLIM = 1.0e-18);
	void bilindist(const Vec& Y, const Vec& S, double Nd, double Md, Vec& Ybl, Vec& Sbl, double BILINLIM = 1.0e-18);
	bool bilin_two(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q, Vec& YPn, Vec& Pn, Vec& YQn, Vec& Qn);
	bool bilin_one(const Vec& YP, const Vec& P, Vec& YPn, Vec& Pn);



};

#endif
