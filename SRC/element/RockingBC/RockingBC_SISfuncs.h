/* ****************************************************************** **
**    RockingBC Semi-infinite strip functions						  **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#ifndef RockingBC_SISfuncs_h
#define RockingBC_SISfuncs_h

#include <cmath>
#include <iostream>

#include <Matrix.h>
#include <Vector.h>

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

#endif