/* ****************************************************************** **
**    Rocking 2D Beam - Column Element    							  **
** ****************************************************************** */

// Written by: Evangelos Avgenakis
// Some preexisting code from basic OpenSees elements was used

#ifndef RockingBC_h
#define RockingBC_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <vector>
#include <fstream>

#include "RockingBC_SISfuncs.h"
#include "RockingBC_SCfuncs.h"
#include "RockingBC_region.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

// #include "Eigen/Dense" // REMOVE EIGEN

using Vec = std::vector<double>;
using VecVecXd = std::vector<Vector>;
using VecVecOS = std::vector< Vector >;
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

private:

	Node *theNodes[2];

	ID  connectedExternalNodes;

	CrdTransf *theCoordTransf;

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

};

#endif
