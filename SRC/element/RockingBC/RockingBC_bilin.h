/* ****************************************************************** **
**    RockingBC bilin functions											      **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#ifndef RockingBC_bilin_h
#define RockingBC_bilin_h

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using Vec = std::vector<double>;

void commony_BL(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB);
bool distintersec(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q);
bool twobilinintersec(double y1, double y2, double p1, double p2, double q1, double q2, double yp, double p0, double yq, double q0);
void NM_BL(const Vec& Y, const Vec& S, double& N, double& M, double& Nd, double& Md);
bool bilinable(double Nd, double Md, double y1, double y2, double BILINLIM = 1.0e-18);
void bilindist(const Vec& Y, const Vec& S, double Nd, double Md, Vec& Ybl, Vec& Sbl, double BILINLIM = 1.0e-18);
bool bilin_two(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q, Vec& YPn, Vec& Pn, Vec& YQn, Vec& Qn);
bool bilin_one(const Vec& YP, const Vec& P, Vec& YPn, Vec& Pn);

#endif