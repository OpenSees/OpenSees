/* ****************************************************************** **
**    RockingBC region functions											      **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#ifndef RockingBC_region_h
#define RockingBC_region_h

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include <Matrix.h>
#include <Vector.h>

using Vec = std::vector<double>;
using Vecint = std::vector<int>;
using VecVecOS = std::vector< Vector >;
using VecVec = std::vector< std::vector<double> >;
using VecVecint = std::vector< std::vector<int> >;
using VecVecVec = std::vector< std::vector< std::vector<double> > >;
using VecMatOS = std::vector< Matrix >;

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

#endif