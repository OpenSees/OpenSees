/* ****************************************************************** **
**    RockingBC bilin functions											      **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#include "RockingBC_bilin.h"

void commony_BL(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB)
{
	Y.clear();
	FA.clear();
	FB.clear();

	int ia = 0;
	int ib = 0;
	while ((ia < ya.size() - 1) || (ib < yb.size() - 1))
	{
		if (ya[ia] == yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib]);
			ia += 1;
			ib += 1;
		}
		else if (ya[ia] < yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib - 1] + (ya[ia] - yb[ib - 1]) / (yb[ib] - yb[ib - 1])*(fb[ib] - fb[ib - 1]));
			ia += 1;
		}
		else {
			Y.push_back(yb[ib]);
			FB.push_back(fb[ib]);
			FA.push_back(fa[ia - 1] + (yb[ib] - ya[ia - 1]) / (ya[ia] - ya[ia - 1])*(fa[ia] - fa[ia - 1]));
			ib += 1;
		}
	}
	Y.push_back(ya[ya.size() - 1]);
	FA.push_back(fa[fa.size() - 1]);
	FB.push_back(fb[fb.size() - 1]);

	return;
}

bool distintersec(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q)
{
	static Vec Y{};
	static Vec PT{};
	static Vec QT{};
	commony_BL(YP, P, YQ, Q, Y, PT, QT);

	int sgn = 0;
	int s = 0;
	for (size_t i = 0; i != PT.size(); i++) {
		if (PT[i] < QT[i]) {
			s = -1;
		} else if(PT[i] == QT[i]) {
			s = 0;
		} else {
			s = 1;
		}
		if (s!=0 && sgn == 0) {
			sgn = s;
		}
		else if (s!=0 && s!=sgn) {
			return true;
		}
	}
	return false;
}

bool twobilinintersec(double y1, double y2, double p1, double p2, double q1, double q2, double yp, double p0, double yq, double q0)
{
	double pe{};
	double qe{};
	if (yp <= yq)
	{
		pe = p0 + (yq - yp) / (y2 - yp)*(p2 - p0);
		qe = q1 + (yp - y1) / (yq - y1)*(q0 - q1);
	} else {
		qe = q0 + (yp - yq) / (y2 - yq)*(q2 - q0);
		pe = p1 + (yq - y1) / (yp - y1)*(p0 - p1);
	}
	if (p1 <= q1 && p2 <= q2 && p0 <= qe && pe <= q0) {
		return false;
	}
	else if (p1 >= q1 && p2 >= q2 && p0 >= qe && pe >= q0) {
		return false;
	}
	else {
		return true;
	}
}

void NM_BL(const Vec& Y, const Vec& S, double& N, double& M, double& Nd, double& Md)
{
	N = 0.;
	M = 0.;

	for (size_t i = 0; i != Y.size()-1; i++) {
		N += (Y[i + 1] - Y[i])*(S[i + 1] + S[i]) / 2.;
		M += (Y[i + 1] - Y[i])*(2 * S[i] * Y[i] + S[i] * Y[i + 1] + S[i + 1] * Y[i] + 2 * S[i + 1] * Y[i + 1]) / 6.;
	}

	double Nlin = (Y[Y.size()-1] - Y[0])*(S[S.size() -1] + S[0]) / 2.;
	double Mlin = (Y[Y.size() -1] - Y[0])*(2 * S[0] * Y[0] + S[0] * Y[Y.size() -1] + S[S.size() -1] * Y[0] + 2 * S[S.size() -1] * Y[Y.size() -1]) / 6.;

	Nd = N - Nlin;
	Md = M - Mlin;
}

bool bilinable(double Nd, double Md, double y1, double y2, double BILINLIM) {
	if (fabs(Nd) < BILINLIM && fabs(Md) > BILINLIM) {
		return false;
	} else if ((fabs(Nd) < BILINLIM && fabs(Md) < BILINLIM) || (2. * y1 + y2 < 3. * Md / Nd && 3. * Md / Nd < y1 + 2. * y2)) {
		return true;
	} else {
		return false;
	}
}

void bilindist(const Vec& Y, const Vec& S, double Nd, double Md, Vec& Ybl, Vec& Sbl, double BILINLIM)
{
	Ybl.clear();
	Sbl.clear();
	if (fabs(Nd) < BILINLIM && fabs(Md) < BILINLIM) {
		Ybl = { Y[0], Y[Y.size() - 1] };
		Sbl = { S[0], S[S.size() - 1] };
		return;
	}

	double s = 2. * Nd / (Y[Y.size()-1] - Y[0]);
	double y0 = 3. * Md / Nd - Y[0] - Y[Y.size() -1];
	double k = (S[S.size() -1] - S[0]) / (Y[Y.size() -1] - Y[0]);
	double snew = S[0] + (y0 - Y[0])*k + s;

	Ybl = { Y[0], y0, Y[Y.size() - 1] };
	Sbl = { S[0], snew, S[S.size() - 1] };
	return;

}

bool bilin_two(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q, Vec& YPn, Vec& Pn, Vec& YQn, Vec& Qn)
{
	double NP{}, MP{}, NPd{}, MPd{}, NQ{}, MQ{}, NQd{}, MQd{};
	NM_BL(YP, P, NP, MP, NPd, MPd);
	NM_BL(YQ, Q, NQ, MQ, NQd, MQd);

	if (!bilinable(NPd, MPd, YP[0], YP[YP.size() - 1]) || !bilinable(NQd, MQd, YQ[0], YQ[YQ.size() - 1])) {
		return false;
	}

	bilindist(YP, P, NPd, MPd, YPn, Pn);
	bilindist(YQ, Q, NQd, MQd, YQn, Qn);

	double yp{}, p0{}, yq{}, q0{};
	if (YPn.size() == 3) {
		yp = YPn[1];
		p0 = Pn[1];
	}
	else {
		yp = 0.5*(YPn[0] + YPn[1]);
		p0 = 0.5*(Pn[0] + Pn[1]);
	}
	if (YQn.size() == 3) {
		yq = YQn[1];
		q0 = Qn[1];
	}
	else {
		yq = 0.5*(YQn[0] + YQn[1]);
		q0 = 0.5*(Qn[0] + Qn[1]);
	}

	bool intersec = twobilinintersec(YPn[0], YPn[YPn.size() - 1], Pn[0], Pn[Pn.size() - 1], Qn[0], Qn[Qn.size() - 1], yp, p0, yq, q0);
	if (!intersec) {
		return true;
	}
	else {
		return false;
	}

}

bool bilin_one(const Vec& YP, const Vec& P, Vec& YPn, Vec& Pn) {
	double NP{}, MP{}, NPd{}, MPd{};
	NM_BL(YP, P, NP, MP, NPd, MPd);

	if (!bilinable(NPd, MPd, YP[0], YP[YP.size() - 1])) {
		return false;
	}

	bilindist(YP, P, NPd, MPd, YPn, Pn);
	return true;

}

