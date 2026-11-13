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

// Code Developed by: Shiye Wang (shiyewang@tongji.edu.cn, shiye.wang@epfl.ch, 465317046@qq.com)

#ifndef CFSTsteel_h
#define CFSTsteel_h

#include <UniaxialMaterial.h>
#include <vector>
#include <cmath>
#include <OPS_Globals.h>
#include <elementAPI.h>

class CFSTsteel : public UniaxialMaterial
{
public:
    CFSTsteel(int tag,
        double E0, double Fy,
        double b2, double siglb, double a, double sigr, double b3,
        double Lambda_y, double Lambda_lb, double Lambda_a, double Lambda_yneg, double Lambda_re,
        double C_y, double C_lb, double C_a, double C_yneg, double C_re,
        double qInf, double b,
        std::vector<double> cK, std::vector<double> gammaK);

    CFSTsteel(void);
    virtual ~CFSTsteel();

public:
    const char* getClassType(void) const { return "CFSTsteel"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial* getCopy(void);
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);
    void Print(OPS_Stream& s, int flag = 0);
    double fxfunction1(double x);
    double fxfunction2(double x);
    double fxfunction3(double x);
    double fdxfunction1(double x);
    double fdxfunction2(double x);
    double fdxfunction3(double x);

private:
    std::vector<double*> stateVariables();
    void returnMapping(double strain_inc);
    void calculateStiffness();
    template <typename T> int sgn(T val);
    double Fy = 0.0;
    double E0 = 0.0;
    double b2 = 0.0;
    double siglb = 0.0;
    double a = 0.0;
    double sigr = 0.0;
    double b3 = 0.0;
    double Lambda_y = 0.0;
    double Lambda_lb = 0.0;
    double Lambda_a = 0.0;
    double Lambda_yneg = 0.0;
    double Lambda_re = 0.0;
    double C_y = 0.0;
    double C_lb = 0.0;
    double C_a = 0.0;
    double C_yneg = 0.0;
    double C_re = 0.0;
    double Eb_pos = 0.0, Eb_neg = 0.0;
    double epsy_pos = 0.0, epsy_neg = 0.0;
    double sigy_project_pos = 0.0, sigy_project_neg = 0.0;
    double epslb = 0.0;
    double slope_lb = 0.0;
    double slope_r = 0.0;
    double siglb_project = 0.0;
    double sigr_project = 0.0;
    double Eb_pos_j_1 = 0.0, cEb_pos_j_1 = 0.0;
    double sigy_project_pos_j_1 = 0.0, csigy_project_pos_j_1 = 0.0;
    double Eb_neg_j_1 = 0.0, cEb_neg_j_1 = 0.0;
    double sigy_project_neg_j_1 = 0.0, csigy_project_neg_j_1 = 0.0;
    double epslb_j_1 = 0.0, cepslb_j_1 = 0.0;
    double siglb_j_1 = 0.0, csiglb_j_1 = 0.0;
    double a_j_1 = 0.0, ca_j_1 = 0.0;
    double sigr_j_1 = 0.0, csigr_j_1 = 0.0;
    double slope_lb_j_1 = 0.0, cslope_lb_j_1 = 0.0;
    double slope_r_j_1 = 0.0, cslope_r_j_1 = 0.0;
    double siglb_project_j_1 = 0.0, csiglb_project_j_1 = 0.0;
    double siglb_project_j = 0.0, csiglb_project_j = 0.0;
    double sigr_project_j_1 = 0.0, csigr_project_j_1 = 0.0;
    double Di = 0.0, cDi = 0.0;
    double Di_1 = 0.0, cDi_1 = 0.0;
    double Di_2 = 0.0, cDi_2 = 0.0;
    double epsreversal = 0.0, cepsreversal = 0.0;
    double sigreversal = 0.0, csigreversal = 0.0;
    double epsreversal_1 = 0.0, cepsreversal_1 = 0.0;
    double sigreversal_1 = 0.0, csigreversal_1 = 0.0;
    double epsreversal_2 = 0.0, cepsreversal_2 = 0.0;
    double sigreversal_2 = 0.0, csigreversal_2 = 0.0;
    double sig = 0.0, csig = 0.0;
    double sig_1 = 0.0, csig_1 = 0.0;
    double sig_2 = 0.0, csig_2 = 0.0;
    double eps = 0.0, ceps = 0.0;
    double eps_1 = 0.0, ceps_1 = 0.0;
    double eps_2 = 0.0, ceps_2 = 0.0;
    double kon = 0.0, ckon = 0.0;
    double kon_1 = 0.0, ckon_1 = 0.0;
    double kon_2 = 0.0, ckon_2 = 0.0;
    double lbstage = 0.0, clbstage = 0.0;
    double beta_y_j_1 = 0.0, cbeta_y_j_1 = 0.0;
    double beta_lb_j_1 = 0.0, cbeta_lb_j_1 = 0.0;
    double beta_a_j_1 = 0.0, cbeta_a_j_1 = 0.0;
    double beta_yneg_j_1 = 0.0, cbeta_yneg_j_1 = 0.0;
    double beta_re_j_1 = 0.0, cbeta_re_j_1 = 0.0;
    double curve_epslb = 0.0, ccurve_epslb = 0.0;
    double curve_epsr = 0.0, ccurve_epsr = 0.0;
    double beta_y_j = 0.0, beta_lb_j = 0.0, beta_a_j = 0.0, beta_yneg_j = 0.0, beta_re_j = 0.0;
    double cbeta_y_j = 0.0, cbeta_lb_j = 0.0, cbeta_a_j = 0.0, cbeta_yneg_j = 0.0, cbeta_re_j = 0.0;
    double Excursion_Flag = 0.0, cExcursion_Flag = 0.0;
    double Reversal_Flag = 0.0, cReversal_Flag = 0.0;
    double Energy_Excrsni_1 = 0.0, cEnergy_Excrsni_1 = 0.0;
    double Energy_Flag = 0.0, cEnergy_Flag = 0.0;
    double Energy_Excrsn = 0.0, cEnergy_Excrsn = 0.0;
    double Energy_total = 0.0, cEnergy_total = 0.0;
    double Ref_Energy_y = 0.0;
    double Ref_Energy_lb = 0.0;
    double Ref_Energy_a = 0.0;
    double Ref_Energy_yneg = 0.0;
    double Ref_Energy_re = 0.0;
    double e = 0.0, ce = 0.0;
    double e_1 = 0.0, ce_1 = 0.0;
    double e_2 = 0.0, ce_2 = 0.0;
    double Yield_Flag = 0.0, cYield_Flag = 0.0;
    double Bucklingdegree = 0.0, cBucklingdegree = 0.0;
    double Buckling_flag = 0.0, cBuckling_flag = 0.0;
    double Stiffness_neg_Flag = 0.0, cStiffness_neg_Flag = 0.0;
    double eps_inflection = 0.0, ceps_inflection = 0.0;
    double sig_inflection = 0.0, csig_inflection = 0.0;
    double e_inflection = 0.0, ce_inflection = 0.0;
    double Minus = 0.0, cMinus = 0.0;
    double Minus_1 = 0.0, cMinus_1 = 0.0;
    double Minus_Flag = 0.0, cMinus_Flag = 0.0;
    double Minus_Flag_1 = 0.0, cMinus_Flag_1 = 0.0;
    double Minus_Flag_2 = 0.0, cMinus_Flag_2 = 0.0;
    double epsdiatance = 0.0, cepsdiatance = 0.0;
    double E_r_j_1 = 0.0, cE_r_j_1 = 0.0;
    double Stiffness_neg = 0.0, cStiffness_neg = 0.0;
    double sig_iso = 0.0, csig_iso = 0.0;
    double exp_point = 0.0, cexp_point = 0.0;
    double sig_lb_re = 0.0, csig_lb_re = 0.0;
    double kon_re = 0.0, ckon_re = 0.0;
    double eps_re = 0.0, ceps_re = 0.0;
    double sig_re = 0.0, csig_re = 0.0;
    double e_re = 0.0, ce_re = 0.0;
    double deps_re = 0.0, cdeps_re = 0.0;
    double epsreversal_re = 0.0, cepsreversal_re = 0.0;
    double sigreversal_re = 0.0, csigreversal_re = 0.0;
    double reloading_Flag = 0.0, creloading_Flag = 0.0;
    double sig_Trial = 0.0, csig_Trial = 0.0;
    double eps_max = 0.0, ceps_max = 0.0;
    double eps_max_Flag = 0.0, ceps_max_Flag = 0.0;
    double lbstage_re = 0.0, clbstage_re = 0.0;
    double Minus_Flag_re = 0.0, cMinus_Flag_re = 0.0;
    double sig_boundary = 0.0, csig_boundary = 0.0;
    double alpha_boundary = 0.0, calpha_boundary = 0.0;
    double Damage_Flag = 0.0, cDamage_Flag = 0.0;
    double lbstrength = 0.0, clbstrength = 0.0;
    const int N_BASIC_PARAMS = 4;
    const int N_PARAM_PER_BACK = 2;
    const double RETURN_MAP_TOL = 10.0e-10;
    const int MAXIMUM_ITERATIONS = 1000;
    double qInf = 0.0;
    double b = 0.0;
    double stiffnessInitial = 0.0;
    std::vector<double> cK;
    std::vector<double> gammaK;
    int nBackstresses = 0;
    double strainConverged = 0.0;
    double strainTrial = 0.0;
    double strainPEqConverged = 0.0;
    double strainPEqTrial = 0.0;
    double stressConverged = 0.0;
    double stressTrial = 0.0;
    std::vector<double> alphaKConverged;
    std::vector<double> alphaKTrial;
    double stiffnessConverged = 0.0;
    double stiffnessTrial = 0.0;
    double flowDirection = 0.0;
    bool plasticLoading = false;

    double Newtoniteration(double (CFSTsteel::* fxfunction)(double), double (CFSTsteel::* fdxfunction)(double), double x0) {
        double fx = (this->*fxfunction)(x0);
        double fdx = (this->*fdxfunction)(x0);
        double h = x0 - (fx / fdx);
        int i = 0;
        while (std::abs(h - x0) > 10.0e-10) {
            x0 = h;
            i++;
            fdx = (this->*fdxfunction)(x0);
            h = x0 - ((this->*fxfunction)(x0) / fdx);
            if (i > 1000) {
                break;
            }
        }
        return h;
    }
};

#endif
