/* ****************************************************************** **
**    RockingBC region functions											      **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#include "RockingBC_region.h"
#include "RockingBC_bilin.h"

void Up_interval_split(const Vector& Yup, const Vector& Up, const Vector& Yw,
	VecVec& Yup_ints, VecVec& Up_ints)
{
	static std::vector<int> Yind{};
	Yind.clear();
	int iy = 0;
	for (size_t iw = 0; iw != Yw.Size(); iw++) {
		while (true) {
			if (Yup[iy] == Yw[iw]) {
			//if (std::fabs(Yup[iy]-Yw[iw])<1.0e-12) {
				Yind.push_back(iy);
				iy += 1;
				break;
			}
			iy += 1;
		}
	}
	
	Yup_ints.clear();
	Up_ints.clear();
	for (size_t i = 0; i != Yind.size() - 1; i++) {
		Vec X1{};
		for (size_t j = Yind[i]; j != Yind[i + 1] +1; j++) {
			X1.push_back(Up[j]);
		}
		Up_ints.push_back(X1);
		Vec X2{};
		for (size_t j = Yind[i]; j != Yind[i + 1] +1; j++) {
			X2.push_back(Yup[j]);
		}
		Yup_ints.push_back(X2);
	}
	return;
}

Vector interval_join(const VecVec& X_ints) {
	static Vec X{};
	X.clear();

	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t j = 0; j != X_ints.at(i).size() - 1; j++) {
			X.push_back(X_ints[i][j]);
		}
	}
	X.push_back(X_ints[X_ints.size()-1][X_ints.at(X_ints.size() - 1).size() - 1]);

	static Vector XX;
	XX = Vector(X.size());

	for (size_t i = 0; i != X.size(); i++) {
		XX[i] = X[i];
	}
	return XX;
}

Matrix interval_join(const VecMatOS& X_ints) {
	static std::vector<int> vecints;
	vecints.clear();
	vecints.push_back(0);

	for (size_t i = 0; i != X_ints.size(); i++) {
		vecints.push_back(vecints[vecints.size() - 1] + X_ints[i].noRows()-1);
	}
	
	static Matrix res;
	res = Matrix(vecints[vecints.size() - 1] + 1, X_ints.at(0).noCols());
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t k = 0; k != X_ints.at(i).noRows()-1; k++) {
			for (size_t l = 0; l != X_ints.at(i).noCols(); l++) {
				res(vecints[i] + k, l) = X_ints[i](k, l);
			}
		}
	}
	
	const Matrix& mm = X_ints[X_ints.size() - 1];
	for (size_t l = 0; l != mm.noCols(); l++) {
		res(res.noRows() - 1, l) = mm(mm.noRows() - 1, l);
	}
	return res;
}

Vector array_join(const VecVec& X_ints) {
	Vec X{};
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t j = 0; j != X_ints.at(i).size(); j++) {
			X.push_back(X_ints[i][j]);
		}
	}

	Vector XX(X.size());
	for (size_t i = 0; i != X.size(); i++) {
		XX[i] = X[i];
	}
	return XX;
}

Matrix array_join(const VecMatOS& X_ints) {
	std::vector<int> vecints{0};
	for (size_t i = 0; i != X_ints.size(); i++) {
		vecints.push_back(vecints[vecints.size() - 1] + X_ints[i].noRows());
	}
	Matrix res = Matrix(vecints[vecints.size() - 1], X_ints.at(0).noCols());
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t k = 0; k != X_ints.at(i).noRows(); k++) {
			for (size_t l = 0; l != X_ints.at(i).noCols(); l++) {
				res(vecints[i] + k, l) = X_ints[i](k, l);
			}
		}
	}

	return res;
}

void commony(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB)
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
		} else if (ya[ia] < yb[ib]) {
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

void interval_interior(double wl, double wr, double ey, double dy, const Vec& up_com, const Vec& yup_com,
	const Vec& ys_com, const Vec& s_com, double beta_Dt,
	Vec& ys_new, Vec& s_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new, 
	Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr, Vec& ua_pos)
{
	
	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = dy / pi;
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}
	double eyn = ey*DU_DS;

	static Vec Y{};
	static Vec Up{};
	static Vec S{};

	commony(yup_com,up_com,ys_com,s_com,Y,Up,S);

	// Plastic displacements differences
	static Vec Upd; Upd.clear();
	double yline{};
	double kyline = (Up[Up.size()-1]-Up[0])/(Y[Y.size()-1]-Y[0]);
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		yline = Up[0]+(Y[iy]-Y[0])*kyline;
		Upd.push_back(yline - Up[iy]);
	}
	
	// Limits
	static Vec Slim; Slim.clear();
	static Vec Slimn; Slimn.clear();
	for (size_t i = 0; i != Y.size(); i++)
	{
		Slim.push_back(S[i]*DAMPC);
		Slimn.push_back(S[i]*DAMPC*DU_DS);
	}

	// Edge stress normalization
	double wln{};
	double dwln_dwl{};
	double wrn{};
	double dwrn_dwr{};
	if (wl >= Slim[0]) {
		wln = Slimn[0]+(wl-Slim[0]);
		dwln_dwl = 1.0;
	} 
	else if (wl>=ey){
		wln = wl*DU_DS;
		dwln_dwl = DU_DS;
	}
	else {
		wln = eyn + (wl - ey);
		dwln_dwl = 1.0;
	}

	if (wr >= Slim[Slim.size()-1]) {
		wrn = Slimn[Slimn.size()-1]+(wr-Slim[Slimn.size()-1]);
		dwrn_dwr = 1.0;
	}
	else if (wr >= ey) {
		wrn = wr*DU_DS;
		dwrn_dwr = DU_DS;
	}
	else {
		wrn = eyn + (wr - ey);
		dwrn_dwr = 1.0;
	}

	double kwnline=(wrn-wln)/dy;
    double dkwnline_dwl=-dwln_dwl/dy;
    double dkwnline_dwr=dwrn_dwr/dy;

	// Plastic displacements into stresses insertion
	static Vec Wn; Wn.clear();
	static Vec dWn_dwl; dWn_dwl.clear();
	static Vec dWn_dwr; dWn_dwr.clear();
	double wline{};
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		Wn.push_back(wln + (Y[iy] - Y[0])*kwnline + Upd[iy]);
		dWn_dwl.push_back(dwln_dwl + (Y[iy] - Y[0])*dkwnline_dwl);
		dWn_dwr.push_back((Y[iy] - Y[0])*dkwnline_dwr);
	}

	// Crossings
	static Vec Yf{}; Yf.clear();
	static Vec Wnf{}; Wnf.clear();
	static Vec Upf{}; Upf.clear();
	static Vec Slimnf{}; Slimnf.clear();
	static Vec dYf_dwl{}; dYf_dwl.clear();
	static Vec dWnf_dwl{}; dWnf_dwl.clear();
	static Vec dSlimnf_dwl{}; dSlimnf_dwl.clear();
	static Vec dYf_dwr{}; dYf_dwr.clear();
	static Vec dWnf_dwr{}; dWnf_dwr.clear();
	static Vec dSlimnf_dwr{}; dSlimnf_dwr.clear();
	
	double wnf1{}, wnf2{};
	bool wnf1found = false;
	bool wnf2found = false;
	double yf1{}, upf1{}, slimnf1{}, dyf1_dwl{}, dyf1_dwr{}, dwnf1_dwl{}, dwnf1_dwr{}, dslimnf1_dwl{}, dslimnf1_dwr{};
	double yf2{}, upf2{}, slimnf2{}, dyf2_dwl{}, dyf2_dwr{}, dwnf2_dwl{}, dwnf2_dwr{}, dslimnf2_dwl{}, dslimnf2_dwr{};
	for (size_t i = 0; i != Wn.size()-1; i++)
	{
        Wnf.push_back(Wn[i]);
        Yf.push_back(Y[i]);
        Upf.push_back(Up[i]);
        Slimnf.push_back(Slimn[i]);
        dWnf_dwl.push_back(dWn_dwl[i]);
        dWnf_dwr.push_back(dWn_dwr[i]);
        dYf_dwl.push_back(0.);
        dYf_dwr.push_back(0.);   
        dSlimnf_dwl.push_back(0.);
        dSlimnf_dwr.push_back(0.);
		wnf1found = false;
		wnf2found = false;

		if ((Wn[i]<Slimn[i] && Wn[i+1]>Slimn[i+1]) || (Wn[i]>Slimn[i] && Wn[i+1]<Slimn[i+1])) {
			wnf1found = true;
            yf1=Y[i]-(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]);
            wnf1=Wn[i]-(Wn[i]-Slimn[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])*(Wn[i+1]-Wn[i]);
            upf1=Up[i]+(yf1-Y[i])/(Y[i+1]-Y[i])*(Up[i+1]-Up[i]);
            slimnf1=Slimn[i]+(yf1-Y[i])/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dyf1_dwl=-(Y[i+1]-Y[i])*(dWn_dwl[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])+(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/pow((Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]),2)*(dWn_dwl[i+1]-dWn_dwl[i]);
            dyf1_dwr=-(Y[i+1]-Y[i])*(dWn_dwr[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])+(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/pow((Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]),2)*(dWn_dwr[i+1]-dWn_dwr[i]);
            dwnf1_dwl=dWn_dwl[i]+(dyf1_dwl)/(Y[i+1]-Y[i])*(Wn[i+1]-Wn[i])+(yf1-Y[i])/(Y[i+1]-Y[i])*(dWn_dwl[i+1]-dWn_dwl[i]);
            dwnf1_dwr=dWn_dwr[i]+(dyf1_dwr)/(Y[i+1]-Y[i])*(Wn[i+1]-Wn[i])+(yf1-Y[i])/(Y[i+1]-Y[i])*(dWn_dwr[i+1]-dWn_dwr[i]);
            dslimnf1_dwl=(dyf1_dwl)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dslimnf1_dwr=(dyf1_dwr)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
		}
		if ((Wn[i]<eyn && Wn[i+1]>eyn) || (Wn[i]>eyn && Wn[i+1]<eyn)) {
			wnf2found = true;
            yf2=Y[i]-(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i]);
            wnf2=eyn;
            upf2=Up[i]+(yf2-Y[i])/(Y[i+1]-Y[i])*(Up[i+1]-Up[i]);
            slimnf2=Slimn[i]+(yf2-Y[i])/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dyf2_dwl=-(Y[i+1]-Y[i])*(dWn_dwl[i])/(Wn[i+1]-Wn[i])+(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i])/(Wn[i+1]-Wn[i])*(dWn_dwl[i+1]-dWn_dwl[i]);
            dyf2_dwr=-(Y[i+1]-Y[i])*(dWn_dwr[i])/(Wn[i+1]-Wn[i])+(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i])/(Wn[i+1]-Wn[i])*(dWn_dwr[i+1]-dWn_dwr[i]);
            dwnf2_dwl=0.;
            dwnf2_dwr=0.;
            dslimnf2_dwl=(dyf2_dwl)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dslimnf2_dwr=(dyf2_dwr)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
		}
		
		if (wnf1found && !wnf2found) {
            Wnf.push_back(wnf1);
            Yf.push_back(yf1);
            Upf.push_back(upf1);
            Slimnf.push_back(slimnf1);
            dWnf_dwl.push_back(dwnf1_dwl);
            dYf_dwl.push_back(dyf1_dwl);
            dWnf_dwr.push_back(dwnf1_dwr);
            dYf_dwr.push_back(dyf1_dwr);
            dSlimnf_dwl.push_back(dslimnf1_dwl);
            dSlimnf_dwr.push_back(dslimnf1_dwr);
		}
		if (wnf2found && !wnf1found) {
            Wnf.push_back(wnf2);
            Yf.push_back(yf2);
            Upf.push_back(upf2);
            Slimnf.push_back(slimnf2);
            dWnf_dwl.push_back(dwnf2_dwl);
            dYf_dwl.push_back(dyf2_dwl);
            dWnf_dwr.push_back(dwnf2_dwr);
            dYf_dwr.push_back(dyf2_dwr);
            dSlimnf_dwl.push_back(dslimnf2_dwl);
            dSlimnf_dwr.push_back(dslimnf2_dwr);
		}
		if (wnf1found && wnf2found) {
			if (yf1 <= yf2) {
                Wnf.push_back(wnf1);
                Yf.push_back(yf1);
                Upf.push_back(upf1);
                Slimnf.push_back(slimnf1);
                Wnf.push_back(wnf2);
                Yf.push_back(yf2);
                Upf.push_back(upf2);
                Slimnf.push_back(slimnf2);
                dWnf_dwl.push_back(dwnf1_dwl);
                dYf_dwl.push_back(dyf1_dwl);
                dWnf_dwr.push_back(dwnf1_dwr);
                dYf_dwr.push_back(dyf1_dwr);
                dWnf_dwl.push_back(dwnf2_dwl);
                dYf_dwl.push_back(dyf2_dwl);
                dWnf_dwr.push_back(dwnf2_dwr);
                dYf_dwr.push_back(dyf2_dwr);
                dSlimnf_dwl.push_back(dslimnf1_dwl);
                dSlimnf_dwr.push_back(dslimnf1_dwr);
                dSlimnf_dwl.push_back(dslimnf2_dwl);
                dSlimnf_dwr.push_back(dslimnf2_dwr);
			}
			else {
                Wnf.push_back(wnf2);
                Yf.push_back(yf2);
                Upf.push_back(upf2);
                Slimnf.push_back(slimnf2);
                Wnf.push_back(wnf1);
                Yf.push_back(yf1);
                Upf.push_back(upf1);   
                Slimnf.push_back(slimnf1);
                dWnf_dwl.push_back(dwnf2_dwl);
                dYf_dwl.push_back(dyf2_dwl);
                dWnf_dwr.push_back(dwnf2_dwr);
                dYf_dwr.push_back(dyf2_dwr);
                dWnf_dwl.push_back(dwnf1_dwl);
                dYf_dwl.push_back(dyf1_dwl);
                dWnf_dwr.push_back(dwnf1_dwr);
                dYf_dwr.push_back(dyf1_dwr);
                dSlimnf_dwl.push_back(dslimnf2_dwl);
                dSlimnf_dwr.push_back(dslimnf2_dwr);
                dSlimnf_dwl.push_back(dslimnf1_dwl);
                dSlimnf_dwr.push_back(dslimnf1_dwr);
			}
		}
	}
	
	//std::cout << Eigen::Map<Vector>(&s_com[0], s_com.size()).transpose() << std::endl;
	
    Wnf.push_back(Wn[Wn.size()-1]);
    Yf.push_back(Y[Y.size()-1]);   
    Upf.push_back(Up[Up.size()-1]);
    Slimnf.push_back(Slimn[Slimn.size()-1]);
    dWnf_dwl.push_back(dWn_dwl[dWn_dwl.size()-1]);
    dWnf_dwr.push_back(dWn_dwr[dWn_dwr.size()-1]);
    dYf_dwl.push_back(0.);
    dYf_dwr.push_back(0.);
    dSlimnf_dwl.push_back(0.);
    dSlimnf_dwr.push_back(0.);
	
	//Categorization

	static std::vector<int> intcats{};
	intcats.clear();
	for (size_t i = 0; i != Wnf.size()-1; i++)
	{
		if (0.5*(Wnf[i] + Wnf[i + 1]) > 0.5*(Slimnf[i] + Slimnf[i + 1])) {
			intcats.push_back(0);
		}
		else if (0.5*(Wnf[i] + Wnf[i + 1]) > eyn) {
			intcats.push_back(1);
		}
		else {
			intcats.push_back(2);
		}
	}

	//Separation into stresses, plastic displacements
	
	static Vec Sf_new{}; Sf_new.clear();
	static Vec dSf_new_dwl{}; dSf_new_dwl.clear();
	static Vec dSf_new_dwr{}; dSf_new_dwr.clear();
	static Vec Upf_new{}; Upf_new.clear();
	static Vec Ua_pos{}; Ua_pos.clear();

	for (size_t i = 0; i != Wnf.size(); i++) {
		if (Wnf[i] > Slimnf[i]) {
            Sf_new.push_back(Slimnf[i]/DU_DS);
			Upf_new.push_back(Upf[i]);
			Ua_pos.push_back(Wnf[i] - Slimnf[i]);
            dSf_new_dwl.push_back(dSlimnf_dwl[i]/DU_DS);
            dSf_new_dwr.push_back(dSlimnf_dwr[i]/DU_DS);
		}
		else if (Wnf[i] > eyn) {
            Sf_new.push_back(Wnf[i]/DU_DS);
			Upf_new.push_back(Upf[i]);
			Ua_pos.push_back(0.);
            dSf_new_dwl.push_back(dWnf_dwl[i]/DU_DS);
            dSf_new_dwr.push_back(dWnf_dwr[i]/DU_DS);
		}
		else {
            Sf_new.push_back(ey);
			Upf_new.push_back(Upf[i]+Wnf[i]-eyn);
			Ua_pos.push_back(0.);
            dSf_new_dwl.push_back(0.);
            dSf_new_dwr.push_back(0.);
		}
	}
	
	// Simplification
	
	ys_new.clear();
	s_new.clear();
	ys_cats.clear();
	dys_new_dwl.clear();
	dys_new_dwr.clear();
	ds_new_dwl.clear();
	ds_new_dwr.clear();
	yup_new.clear();
	up_new.clear();
	ua_pos.clear();

	ys_new.push_back(Yf[0]);
	s_new.push_back(Sf_new[0]);	
	
	dys_new_dwl.push_back(dYf_dwl[0]);
	ds_new_dwl.push_back(dSf_new_dwl[0]);	
	dys_new_dwr.push_back(dYf_dwr[0]);
	ds_new_dwr.push_back(dSf_new_dwr[0]);	
	ys_cats.push_back(intcats[0]);
	ua_pos.push_back(Ua_pos[0]);
	
	for (size_t i = 1; i != Yf.size()-1; i++) {
		if (intcats[i-1]==2 && intcats[i]==2) {
			continue;
		}
		ys_new.push_back(Yf[i]);
		s_new.push_back(Sf_new[i]);
		ys_cats.push_back(intcats[i]);
		dys_new_dwl.push_back(dYf_dwl[i]);
		ds_new_dwl.push_back(dSf_new_dwl[i]);	
		dys_new_dwr.push_back(dYf_dwr[i]);
		ds_new_dwr.push_back(dSf_new_dwr[i]);	
		ua_pos.push_back(Ua_pos[i]);
	}
	ys_new.push_back(Yf[Yf.size()-1]);
	s_new.push_back(Sf_new[Sf_new.size()-1]);
	dys_new_dwl.push_back(dYf_dwl[dYf_dwl.size()-1]);
	ds_new_dwl.push_back(dSf_new_dwl[dSf_new_dwl.size()-1]);
	dys_new_dwr.push_back(dYf_dwr[dYf_dwr.size()-1]);
	ds_new_dwr.push_back(dSf_new_dwr[dSf_new_dwr.size()-1]);	
	ua_pos.push_back(Ua_pos[Ua_pos.size()-1]);
	
	yup_new.push_back(Yf[0]);
	up_new.push_back(Upf_new[0]);	
	for (size_t i = 1; i != Yf.size()-1; i++) {
		if (intcats[i - 1] == 2 && intcats[i] == 2) {
			continue;
		}
		yup_new.push_back(Yf[i]);
		up_new.push_back(Upf_new[i]);
	}
	yup_new.push_back(Yf[Yf.size()-1]);
	up_new.push_back(Upf_new[Upf_new.size()-1]);

	return;
}

void interval_dists(const Vector& Yw, const Vector& W, const VecVec& Yupi_com, const VecVec& Upi_com, const VecVec& Ysi_com, const VecVec& Si_com, double ey, double beta_Dt,
	VecVec& Ysi, VecVec& Si, VecVec& Yupi_new, VecVec& Upi_new,
	VecVecint& Ys_cats, Vector& Nints, Vector& Mints, Matrix& dNints_dW, Matrix& dMints_dW, VecVec& Ua_pos, VecMatOS& dYsi_dW, VecMatOS& dSi_dW)
{
	
	VecVec dys_dwl_list( W.Size() - 1, std::vector<double>{} );
	VecVec ds_dwl_list( W.Size() - 1, std::vector<double>{} );
	VecVec dys_dwr_list( W.Size() - 1, std::vector<double>{} );
	VecVec ds_dwr_list( W.Size() - 1, std::vector<double>{} );
	
	for (size_t i = 0; i != W.Size()-1; i++) {

		interval_interior(W[i], W[i + 1], ey, Yw[i + 1] - Yw[i], Upi_com[i], Yupi_com[i],Ysi_com[i], Si_com[i], beta_Dt,
			Ysi[i], Si[i], Ys_cats[i], Yupi_new[i], Upi_new[i],
			dys_dwl_list[i], dys_dwr_list[i], ds_dwl_list[i], ds_dwr_list[i],
			Ua_pos[i]);
	}

	static Vector dNdW{}, dMdW{};
	
	for (size_t i = 0; i != W.Size() - 1; i++) {
		
		Vec dwl_dW(W.Size()); dwl_dW[i] = 1.0;
		Vec dwr_dW(W.Size()); dwr_dW[i+1] = 1.0;
		Matrix dys_dW = Matrix(dys_dwl_list[i].size(), W.Size());
		Matrix ds_dW = Matrix(ds_dwl_list[i].size(), W.Size());
		for (size_t l = 0; l != W.Size(); l++) {
			for (size_t k = 0; k != dys_dwl_list[i].size(); k++) {
				dys_dW(k, l) += dys_dwl_list[i][k] * dwl_dW[l];
				dys_dW(k, l) += dys_dwr_list[i][k] * dwr_dW[l];
				ds_dW(k, l) += ds_dwl_list[i][k] * dwl_dW[l];
				ds_dW(k, l) += ds_dwr_list[i][k] * dwr_dW[l];
			}
		}
		dYsi_dW[i] = dys_dW;
		dSi_dW[i] = ds_dW;
		NM_calc_int(Ysi[i], dys_dW, Si[i], ds_dW, Nints[i], Mints[i], dNdW, dMdW);
		for (size_t j = 0; j != W.Size(); j++) {
			dNints_dW(i,j) = dNdW(j);
			dMints_dW(i,j) = dMdW(j);
		}
	}
	return;

}

void NM_calc_int(const Vec& Ys, const Matrix& dYs_dW, const Vec& S, const Matrix& dS_dW, double& N, double& M, Vector& dN_dW, Vector& dM_dW)
{
	N = 0;
	M = 0;
	dN_dW = Vector(dYs_dW.noCols());
	dM_dW = Vector(dS_dW.noCols());

	for (size_t i = 0; i != Ys.size() - 1; i++)
	{

		N += (Ys[i + 1] - Ys[i])*(S[i] + S[i + 1]) / 2.;
		M += (Ys[i + 1] - Ys[i])*(2 * S[i] * Ys[i] + S[i] * Ys[i + 1] + S[i + 1] * Ys[i] + 2 * S[i + 1] * Ys[i + 1]) / 6.;

		for (size_t j= 0 ; j != dN_dW.Size(); j++) {
			dN_dW(j) += (-S[i] / 2. - S[i + 1] / 2.) * dYs_dW(i,j) + (S[i] / 2. + S[i + 1] / 2.) * dYs_dW(i + 1,j) + (Ys[i + 1] / 2. - Ys[i] / 2.) * dS_dW(i,j) + (Ys[i + 1] / 2. - Ys[i] / 2.) * dS_dW(i + 1,j);
			dM_dW(j) += (-(S[i] * Ys[i]) / 3. - (S[i] * Ys[i + 1]) / 6. - (S[i + 1] * Ys[i]) / 6. - (S[i + 1] * Ys[i + 1]) / 3. - ((2. * S[i] + S[i + 1]) * (Ys[i] - Ys[i + 1])) / 6.) * dYs_dW(i,j) +
				((S[i] * Ys[i]) / 3. + (S[i] * Ys[i + 1]) / 6. + (S[i + 1] * Ys[i]) / 6. + (S[i + 1] * Ys[i + 1]) / 3. - ((S[i] + 2. * S[i + 1]) * (Ys[i] - Ys[i + 1])) / 6.) * dYs_dW(i + 1,j) +
				(-((Ys[i] - Ys[i + 1]) * (2. * Ys[i] + Ys[i + 1])) / 6.) * dS_dW(i,j) + (-((Ys[i] - Ys[i + 1]) * (Ys[i] + 2. * Ys[i + 1])) / 6.) * dS_dW(i + 1,j);
		}

	}
	return;
}

void critpoints(const Vec& y, const Vec& s, int rinit, int rend, Vecint& cp)
{
	cp.clear();

	for (size_t i = rinit + 1; i != rend; i++) {
		// Slope change
		if ((s[i] - s[i - 1])*(s[i + 1] - s[i]) <= 0 && (s[i] - s[i - 1] != 0 || s[i + 1] - s[i] != 0))
		{
			cp.push_back(i);
			continue;
		}
		// Curvature change
		//if (i != rend - 1) {
		//	if (((s[i + 1] - s[i])*(y[i] - y[i - 1])*(y[i + 2] - y[i + 1]) - (s[i] - s[i - 1])*(y[i + 1] - y[i])*(y[i + 2] - y[i + 1]))*(
		//		(s[i + 2] - s[i + 1])*(y[i + 1] - y[i])*(y[i] - y[i - 1]) - (s[i + 1] - s[i])*(y[i] - y[i - 1])*(y[i + 2] - y[i + 1])) < 0) {
		//		cp.push_back(i);
		//		cp.push_back(i + 1);
		//	}
		//}
	}
	return;
}

void int_bilin(const Vecint& ys_cats, const Vec& ys, const Vec& s, const Vec& yup, const Vec& up, const Vec& ua_pos, double ey,
	Vec& ys_new, Vec& s_new, Vec& yup_new, Vec& up_new)
{

	if ((ys.size() != yup.size()) || (ys_cats.size() != ys.size() - 1) || (ys.size() != ua_pos.size())) {
		ys_new = ys;
		s_new = s;
		yup_new = yup;
		up_new = up;
		return;
	}

	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = (ys[ys.size()-1] - ys[0]) / pi;

	static Vec w{}; w.clear();
	static Vec sm{}; sm.clear();
	for (size_t i = 0; i != s.size(); i++) {
		w.push_back(s[i] * DU_DS + ua_pos[i]);
		sm.push_back(s[i] * DU_DS);
	}

	static Vecint regi{}; regi.clear();
	regi.push_back(0);
	for (size_t i = 0; i != ys_cats.size()-1; i++) {
		if (ys_cats[i + 1] != ys_cats[i]) {
			regi.push_back(i + 1);
		}
	}
	regi.push_back(ys.size() - 1);

	static VecVecint regs{}; regs.clear();
	static Vecint regs_cats{}; regs_cats.clear();
	for (size_t i = 0; i != regi.size()-1; i++) {
		if (ys_cats[regi[i]] != 2) {
			Vecint v = { regi[i], regi[i + 1] };
			regs.push_back(v);
			regs_cats.push_back(ys_cats[regi[i]]);
		}
	}

	static VecVecint regs2{}; regs2.clear();
	static Vecint regs2_cats{}; regs2_cats.clear();
	static Vecint m{}; m.clear();
	static Vecint cp{}; cp.clear();
	for (size_t ir = 0; ir != regs.size(); ir++) {
		Vecint r = regs[ir];
		m.clear();
		m.push_back(r[0]);
		critpoints(ys, s, r[0], r[1], cp);
		m.insert(m.end(), cp.begin(), cp.end());
		critpoints(yup, up, r[0], r[1], cp);
		m.insert(m.end(), cp.begin(), cp.end());
		m.push_back(r[1]);
		sort(m.begin(), m.end());
		m.erase(unique(m.begin(), m.end()), m.end());
		for (size_t j = 0; j != m.size() - 1; j++) {
			Vecint v = { m[j],m[j + 1] };
			regs2.push_back(v);
			regs2_cats.push_back(regs_cats[ir]);
		}
	}
	static VecVecint i_s_bl{}; i_s_bl.clear();
	static VecVecint i_up_bl{}; i_up_bl.clear();
	static VecVec s_bl{}; s_bl.clear();
	static VecVec ys_bl{}; ys_bl.clear();
	static VecVec up_bl{}; up_bl.clear();
	static VecVec yup_bl{}; yup_bl.clear();

	static Vec ysp{};
	static Vec smp{};
	static Vec sp{};
	static Vec wp{};
	static Vec ss{};
	static Vec ys_try{}, sm_try{}, yw_try{}, w_try{}, s_try{};
	static Vecint v{};

	for (size_t ir = 0; ir != regs2.size(); ir++) {
		Vecint r = regs2[ir];
		if (regs2_cats[ir] == 0) {
			ysp = Vec(ys.begin() + r[0], ys.begin() + r[1] + 1);
			smp = Vec(sm.begin() + r[0], sm.begin() + r[1] + 1);
			wp = Vec(w.begin() + r[0], w.begin() + r[1] + 1);
			ys_try.clear(); sm_try.clear(); yw_try.clear(); w_try.clear();
			bool suc = bilin_two(ysp, smp, ysp, wp, ys_try, sm_try, yw_try, w_try);
			if (suc) {

				if (yw_try.size() == 2) {
					v = { r[0],r[1] };
					i_s_bl.push_back(v);
					i_up_bl.push_back(v);
					ys_bl.push_back(ys_try);
					ss.clear();
					for (size_t i = 0; i != sm_try.size(); i++) {
						ss.push_back(sm_try[i] / DU_DS);
					}
					s_bl.push_back(ss);

					ss = { yup[r[0]],yup[r[r.size() - 1]] };
					yup_bl.push_back(ss);
					ss = { up[r[0]],up[r[r.size() - 1]] };
					up_bl.push_back(ss);
				}
				else {
					double yupn0 = yw_try[1];
					double upn0 = (up[r[0]] + (yupn0 - yup[r[0]]) / (yup[r[1]] - yup[r[0]])*(up[r[1]] - up[r[0]])) + (w_try[0] + (yupn0 - yw_try[0]) / (yw_try[yw_try.size() - 1] - yw_try[0])*(w_try[w_try.size() - 1] - w_try[0])) - (w_try[1]);
					
					if (sm_try[1] <= 0 && sm_try[1] / DU_DS >= ey && upn0 <= 0) {
						v = { r[0],r[1] };
						i_s_bl.push_back(v);
						i_up_bl.push_back(v);
						ys_bl.push_back(ys_try);
						ss.clear();
						for (size_t i = 0; i != sm_try.size(); i++) {
							ss.push_back(sm_try[i] / DU_DS);
						}
						s_bl.push_back(ss);

						ss = { yup[r[0]],yupn0,yup[r[1]] };
						yup_bl.push_back(ss);
						ss = { up[r[0]],upn0,up[r[1]] };
						up_bl.push_back(ss);
					}
				}
			}
		}
		else {
			ysp = Vec(ys.begin() + r[0], ys.begin() + r[1] + 1);
			sp = Vec(s.begin() + r[0], s.begin() + r[1] + 1);
			ys_try.clear(); s_try.clear();
			bool suc = bilin_one(ysp, sp, ys_try, s_try);
			if (suc) {
				if (ys_try.size() == 2) {

					v= { r[0],r[1] };
					i_s_bl.push_back(v);
					i_up_bl.push_back(v);
					ys_bl.push_back(ys_try);
					s_bl.push_back(s_try);

					ss = { yup[r[0]],yup[r[1]] };
					yup_bl.push_back(ss);
					ss = { up[r[0]],up[r[1]] };
					up_bl.push_back(ss);
				}
				else {
					double yupn0 = ys_try[1];
					double upn0 = (up[r[0]] + (yupn0 - yup[r[0]]) / (yup[r[r.size() - 1]] - yup[r[0]])*(up[r[r.size() - 1]] - up[r[0]])) - (s_try[1] - (s_try[0] + (yupn0 - ys_try[0]) / (ys_try[ys_try.size() - 1] - ys_try[0])*(s_try[s_try.size() - 1] - s_try[0])))*DU_DS;
					
					if (s_try[1] <= 0 && s_try[1] >= ey && upn0 <= 0) {
						v = { r[0],r[1] };
						i_s_bl.push_back(v);
						i_up_bl.push_back(v);
						ys_bl.push_back(ys_try);
						s_bl.push_back(s_try);

						ss = { yup[r[0]],yupn0,yup[r[1]] };
						yup_bl.push_back(ss);
						ss = { up[r[0]],upn0,up[r[1]] };
						up_bl.push_back(ss);
					}
				}
			}

		}
	}

	ys_new.clear();
	s_new.clear();
	int ir = 0;
	int i = 0;
	while (i < ys.size()) {
		if (ir < i_s_bl.size() && (i == i_s_bl[ir][0] || i == i_s_bl[ir][0] + 1)) {
			if (i == i_s_bl[ir][0]) {
				for (size_t k = 0; k != ys_bl[ir].size(); k++) {
					ys_new.push_back(ys_bl[ir][k]);
					s_new.push_back(s_bl[ir][k]);
				}
			}
			else {
				for (size_t k = 1; k != ys_bl[ir].size(); k++) {
					ys_new.push_back(ys_bl[ir][k]);
					s_new.push_back(s_bl[ir][k]);
				}
			}
			i = i_s_bl[ir][1] + 1;
			ir += 1;
		}
		else {
			ys_new.push_back(ys[i]);
			s_new.push_back(s[i]);
			i += 1;
		}
	}

	yup_new.clear();
	up_new.clear();
	ir = 0;
	i = 0;
	while (i < yup.size()) {
		if (ir < i_up_bl.size() && (i == i_up_bl[ir][0] || i == i_up_bl[ir][0] + 1)) {
			if (i == i_up_bl[ir][0]) {
				for (size_t k = 0; k != yup_bl[ir].size(); k++) {
					yup_new.push_back(yup_bl[ir][k]);
					up_new.push_back(up_bl[ir][k]);
				}
			}
			else {
				for (size_t k = 1; k != yup_bl[ir].size(); k++) {
					yup_new.push_back(yup_bl[ir][k]);
					up_new.push_back(up_bl[ir][k]);
				}
			}
			i = i_up_bl[ir][1] + 1;
			ir += 1;
		}
		else {
			yup_new.push_back(yup[i]);
			up_new.push_back(up[i]);
			i += 1;
		}
	}
}

void Up_interval_split_K(const Vector& Yup, const Vector& Up, const Vector& Kup, const Vector& Yw,
	VecVecOS& Yup_ints, VecVecOS& Up_ints, VecVecOS& Kup_ints)
{
	static std::vector<int> Yind{};
	Yind.clear();

	int iy = 0;
	for (size_t iw = 0; iw != Yw.Size(); iw++) {
		while (true) {
			if (Yup[iy] == Yw[iw]) {
				//if (std::fabs(Yup[iy]-Yw[iw])<1.0e-15) {
				Yind.push_back(iy);
				iy += 1;
				break;
			}
			iy += 1;
		}
	}

	Yup_ints.clear();
	Up_ints.clear();
	Kup_ints.clear();
	for (size_t i = 0; i != Yind.size() - 1; i++) {
		//Up_ints.push_back(Up.segment(Yind[i], Yind[i + 1] - Yind[i] + 1));
		//Yup_ints.push_back(Yup.segment(Yind[i], Yind[i + 1] - Yind[i] + 1));
		//Kup_ints.push_back(Kup.segment(Yind[i], Yind[i + 1] - Yind[i]));

		Vector upint(Yind[i + 1] - Yind[i] + 1);
		Vector yupint(Yind[i + 1] - Yind[i] + 1);
		for (size_t j = 0; j != Yind[i + 1] - Yind[i] + 1; j++) {
			upint(j) = Up(Yind[i] + j);
			yupint(j) = Yup(Yind[i] + j);
		}
		Vector kupint(Yind[i + 1] - Yind[i]);
		for (size_t j = 0; j != Yind[i + 1] - Yind[i]; j++) {
			kupint(j) = Kup(Yind[i] + j);
		}
		Up_ints.push_back(upint);
		Yup_ints.push_back(yupint);
		Kup_ints.push_back(kupint);
	}

	return;
}

void commony_K(const Vector& ya, const Vector& fa, const Vector& ka, const Vector& yb, const Vector& fb, const Vector& kb, Vec& Y, Vec& FA, Vec& FB, Vec& KA, Vec& KB)
{
	Y.clear();
	FA.clear();
	FB.clear();
	KA.clear();
	KB.clear();

	int ia = 0;
	int ib = 0;
	while ((ia < ya.Size() - 1) || (ib < yb.Size() - 1))
	{
		if (ya[ia] == yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib]);
			KA.push_back(ka[ia]);
			KB.push_back(kb[ib]);
			ia += 1;
			ib += 1;
		}
		else if (ya[ia] < yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib - 1] + (ya[ia] - yb[ib - 1]) / (yb[ib] - yb[ib - 1]) * (fb[ib] - fb[ib - 1]));
			KA.push_back(ka[ia]);
			KB.push_back(kb[ib - 1]);
			ia += 1;
		}
		else {
			Y.push_back(yb[ib]);
			FB.push_back(fb[ib]);
			FA.push_back(fa[ia - 1] + (yb[ib] - ya[ia - 1]) / (ya[ia] - ya[ia - 1]) * (fa[ia] - fa[ia - 1]));
			KB.push_back(kb[ib]);
			KA.push_back(ka[ia - 1]);
			ib += 1;
		}
	}
	Y.push_back(ya[ya.Size() - 1]);
	FA.push_back(fa[fa.Size() - 1]);
	FB.push_back(fb[fb.Size() - 1]);

	return;
}

void interval_interior_K(double wl, double wr, double ey, double dy, const Vector& up_com, const Vector& yup_com, const Vector& kup_com,
	const Vector& ys_com, const Vector& s_com, const Vector& ks_com, double beta_Dt,
	Vec& ys_new, Vec& s_new, Vec& ks_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new, Vec& kup_new,
	Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr, Vec& dks_new_dwl, Vec& dks_new_dwr,
	Vec& ydks, Vec& dks, Vec& dydks_dwl, Vec& dydks_dwr, Vec& ddks_dwl, Vec& ddks_dwr,
	Vec& ds, Vec& dds_dwl, Vec& dds_dwr)
{
	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = dy / pi;
	double DAMPC = beta_Dt / (1.0 + beta_Dt);
	double eyn = ey * DU_DS;

	static Vec Y{};
	static Vec Up{};
	static Vec S{};
	static Vec KUp{};
	static Vec KS{};

	commony_K(yup_com, up_com, kup_com, ys_com, s_com, ks_com, Y, Up, S, KUp, KS);

	// Plastic displacements differences
	static Vec Upd; Upd.clear();
	static Vec KUpd; KUpd.clear();
	double yline{};
	double kyline = (Up[Up.size() - 1] - Up[0]) / (Y[Y.size() - 1] - Y[0]);
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		yline = Up[0] + (Y[iy] - Y[0]) * kyline;
		Upd.push_back(yline - Up[iy]);
	}
	for (size_t iy = 0; iy != Y.size() - 1; iy++)
	{
		KUpd.push_back(kyline - KUp[iy]);
	}

	// Limits
	static Vec Slim; Slim.clear();
	static Vec Slimn; Slimn.clear();
	static Vec KSlim; KSlim.clear();
	for (size_t i = 0; i != Y.size(); i++)
	{
		Slim.push_back(S[i] * DAMPC);
		Slimn.push_back(S[i] * DAMPC * DU_DS);
	}
	for (size_t i = 0; i != Y.size() - 1; i++)
	{
		KSlim.push_back(KS[i] * DAMPC);
	}

	// Edge stress normalization
	double wln{};
	double dwln_dwl{};
	double wrn{};
	double dwrn_dwr{};
	if (wl >= Slim[0]) {
		wln = Slimn[0] + (wl - Slim[0]);
		dwln_dwl = 1.0;
	}
	else if (wl >= ey) {
		wln = wl * DU_DS;
		dwln_dwl = DU_DS;
	}
	else {
		wln = eyn + (wl - ey);
		dwln_dwl = 1.0;
	}

	if (wr >= Slim[Slim.size() - 1]) {
		wrn = Slimn[Slimn.size() - 1] + (wr - Slim[Slimn.size() - 1]);
		dwrn_dwr = 1.0;
	}
	else if (wr >= ey) {
		wrn = wr * DU_DS;
		dwrn_dwr = DU_DS;
	}
	else {
		wrn = eyn + (wr - ey);
		dwrn_dwr = 1.0;
	}

	double kwnline = (wrn - wln) / dy;
	double dkwnline_dwl = -dwln_dwl / dy;
	double dkwnline_dwr = dwrn_dwr / dy;

	// Plastic displacements into stresses insertion
	static Vec Wn; Wn.clear();
	static Vec dWn_dwl; dWn_dwl.clear();
	static Vec dWn_dwr; dWn_dwr.clear();
	static Vec KWn; KWn.clear();
	static Vec dKWn_dwl; dKWn_dwl.clear();
	static Vec dKWn_dwr; dKWn_dwr.clear();
	double wline{};
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		Wn.push_back(wln + (Y[iy] - Y[0]) * kwnline + Upd[iy]);
		dWn_dwl.push_back(dwln_dwl + (Y[iy] - Y[0]) * dkwnline_dwl);
		dWn_dwr.push_back((Y[iy] - Y[0]) * dkwnline_dwr);
	}
	for (size_t iy = 0; iy != Y.size() - 1; iy++)
	{
		KWn.push_back(kwnline + KUpd[iy]);
		dKWn_dwl.push_back(dkwnline_dwl);
		dKWn_dwr.push_back(dkwnline_dwr);
	}

	// Crossings
	static Vec Yf{}; Yf.clear();
	static Vec Wnf{}; Wnf.clear();
	static Vec Upf{}; Upf.clear();
	static Vec Slimnf{}; Slimnf.clear();
	static Vec dYf_dwl{}; dYf_dwl.clear();
	static Vec dWnf_dwl{}; dWnf_dwl.clear();
	static Vec dSlimnf_dwl{}; dSlimnf_dwl.clear();
	static Vec dYf_dwr{}; dYf_dwr.clear();
	static Vec dWnf_dwr{}; dWnf_dwr.clear();
	static Vec dSlimnf_dwr{}; dSlimnf_dwr.clear();
	static Vec KWnf{}; KWnf.clear();
	static Vec KUpf{}; KUpf.clear();
	static Vec KSlimf{}; KSlimf.clear();
	static Vec dKWnf_dwl{}; dKWnf_dwl.clear();
	static Vec dKWnf_dwr{}; dKWnf_dwr.clear();

	double wnf1{}, wnf2{};
	bool wnf1found = false;
	bool wnf2found = false;
	double yf1{}, upf1{}, slimnf1{}, dyf1_dwl{}, dyf1_dwr{}, dwnf1_dwl{}, dwnf1_dwr{}, dslimnf1_dwl{}, dslimnf1_dwr{};
	double yf2{}, upf2{}, slimnf2{}, dyf2_dwl{}, dyf2_dwr{}, dwnf2_dwl{}, dwnf2_dwr{}, dslimnf2_dwl{}, dslimnf2_dwr{};
	for (size_t i = 0; i != Wn.size() - 1; i++)
	{
		Wnf.push_back(Wn[i]);
		Yf.push_back(Y[i]);
		Upf.push_back(Up[i]);
		Slimnf.push_back(Slimn[i]);
		dWnf_dwl.push_back(dWn_dwl[i]);
		dWnf_dwr.push_back(dWn_dwr[i]);
		dYf_dwl.push_back(0.);
		dYf_dwr.push_back(0.);
		dSlimnf_dwl.push_back(0.);
		dSlimnf_dwr.push_back(0.);
		KWnf.push_back(KWn[i]);
		dKWnf_dwl.push_back(dKWn_dwl[i]);
		dKWnf_dwr.push_back(dKWn_dwr[i]);
		KUpf.push_back(KUp[i]);
		KSlimf.push_back(KSlim[i]);
		wnf1found = false;
		wnf2found = false;

		if ((Wn[i]<Slimn[i] && Wn[i + 1]>Slimn[i + 1]) || (Wn[i] > Slimn[i] && Wn[i + 1] < Slimn[i + 1])) {
			wnf1found = true;
			yf1 = Y[i] - (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]);
			wnf1 = Wn[i] - (Wn[i] - Slimn[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) * (Wn[i + 1] - Wn[i]);
			upf1 = Up[i] + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (Up[i + 1] - Up[i]);
			slimnf1 = Slimn[i] + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dyf1_dwl = -(Y[i + 1] - Y[i]) * (dWn_dwl[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / pow((Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]), 2) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dyf1_dwr = -(Y[i + 1] - Y[i]) * (dWn_dwr[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / pow((Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]), 2) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dwnf1_dwl = dWn_dwl[i] + (dyf1_dwl) / (Y[i + 1] - Y[i]) * (Wn[i + 1] - Wn[i]) + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dwnf1_dwr = dWn_dwr[i] + (dyf1_dwr) / (Y[i + 1] - Y[i]) * (Wn[i + 1] - Wn[i]) + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dslimnf1_dwl = (dyf1_dwl) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dslimnf1_dwr = (dyf1_dwr) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
		}
		if ((Wn[i]<eyn && Wn[i + 1]>eyn) || (Wn[i] > eyn&& Wn[i + 1] < eyn)) {
			wnf2found = true;
			yf2 = Y[i] - (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]);
			wnf2 = eyn;
			upf2 = Up[i] + (yf2 - Y[i]) / (Y[i + 1] - Y[i]) * (Up[i + 1] - Up[i]);
			slimnf2 = Slimn[i] + (yf2 - Y[i]) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dyf2_dwl = -(Y[i + 1] - Y[i]) * (dWn_dwl[i]) / (Wn[i + 1] - Wn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]) / (Wn[i + 1] - Wn[i]) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dyf2_dwr = -(Y[i + 1] - Y[i]) * (dWn_dwr[i]) / (Wn[i + 1] - Wn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]) / (Wn[i + 1] - Wn[i]) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dwnf2_dwl = 0.;
			dwnf2_dwr = 0.;
			dslimnf2_dwl = (dyf2_dwl) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dslimnf2_dwr = (dyf2_dwr) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
		}

		if (wnf1found && !wnf2found) {
			Wnf.push_back(wnf1);
			Yf.push_back(yf1);
			Upf.push_back(upf1);
			Slimnf.push_back(slimnf1);
			dWnf_dwl.push_back(dwnf1_dwl);
			dYf_dwl.push_back(dyf1_dwl);
			dWnf_dwr.push_back(dwnf1_dwr);
			dYf_dwr.push_back(dyf1_dwr);
			dSlimnf_dwl.push_back(dslimnf1_dwl);
			dSlimnf_dwr.push_back(dslimnf1_dwr);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
		if (wnf2found && !wnf1found) {
			Wnf.push_back(wnf2);
			Yf.push_back(yf2);
			Upf.push_back(upf2);
			Slimnf.push_back(slimnf2);
			dWnf_dwl.push_back(dwnf2_dwl);
			dYf_dwl.push_back(dyf2_dwl);
			dWnf_dwr.push_back(dwnf2_dwr);
			dYf_dwr.push_back(dyf2_dwr);
			dSlimnf_dwl.push_back(dslimnf2_dwl);
			dSlimnf_dwr.push_back(dslimnf2_dwr);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
		if (wnf1found && wnf2found) {
			if (yf1 <= yf2) {
				Wnf.push_back(wnf1);
				Yf.push_back(yf1);
				Upf.push_back(upf1);
				Slimnf.push_back(slimnf1);
				Wnf.push_back(wnf2);
				Yf.push_back(yf2);
				Upf.push_back(upf2);
				Slimnf.push_back(slimnf2);
				dWnf_dwl.push_back(dwnf1_dwl);
				dYf_dwl.push_back(dyf1_dwl);
				dWnf_dwr.push_back(dwnf1_dwr);
				dYf_dwr.push_back(dyf1_dwr);
				dWnf_dwl.push_back(dwnf2_dwl);
				dYf_dwl.push_back(dyf2_dwl);
				dWnf_dwr.push_back(dwnf2_dwr);
				dYf_dwr.push_back(dyf2_dwr);
				dSlimnf_dwl.push_back(dslimnf1_dwl);
				dSlimnf_dwr.push_back(dslimnf1_dwr);
				dSlimnf_dwl.push_back(dslimnf2_dwl);
				dSlimnf_dwr.push_back(dslimnf2_dwr);
			}
			else {
				Wnf.push_back(wnf2);
				Yf.push_back(yf2);
				Upf.push_back(upf2);
				Slimnf.push_back(slimnf2);
				Wnf.push_back(wnf1);
				Yf.push_back(yf1);
				Upf.push_back(upf1);
				Slimnf.push_back(slimnf1);
				dWnf_dwl.push_back(dwnf2_dwl);
				dYf_dwl.push_back(dyf2_dwl);
				dWnf_dwr.push_back(dwnf2_dwr);
				dYf_dwr.push_back(dyf2_dwr);
				dWnf_dwl.push_back(dwnf1_dwl);
				dYf_dwl.push_back(dyf1_dwl);
				dWnf_dwr.push_back(dwnf1_dwr);
				dYf_dwr.push_back(dyf1_dwr);
				dSlimnf_dwl.push_back(dslimnf2_dwl);
				dSlimnf_dwr.push_back(dslimnf2_dwr);
				dSlimnf_dwl.push_back(dslimnf1_dwl);
				dSlimnf_dwr.push_back(dslimnf1_dwr);
			}
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
	}

	//std::cout << Eigen::Map<Vector>(&s_com[0], s_com.size()).transpose() << std::endl;

	Wnf.push_back(Wn[Wn.size() - 1]);
	Yf.push_back(Y[Y.size() - 1]);
	Upf.push_back(Up[Up.size() - 1]);
	Slimnf.push_back(Slimn[Slimn.size() - 1]);
	dWnf_dwl.push_back(dWn_dwl[dWn_dwl.size() - 1]);
	dWnf_dwr.push_back(dWn_dwr[dWn_dwr.size() - 1]);
	dYf_dwl.push_back(0.);
	dYf_dwr.push_back(0.);
	dSlimnf_dwl.push_back(0.);
	dSlimnf_dwr.push_back(0.);

	//Categorization

	static std::vector<int> intcats{};
	intcats.clear();
	for (size_t i = 0; i != Wnf.size() - 1; i++)
	{
		if (0.5 * (Wnf[i] + Wnf[i + 1]) > 0.5 * (Slimnf[i] + Slimnf[i + 1])) {
			intcats.push_back(0);
		}
		else if (0.5 * (Wnf[i] + Wnf[i + 1]) > eyn) {
			intcats.push_back(1);
		}
		else {
			intcats.push_back(2);
		}
	}

	//Separation into stresses, plastic displacements

	static Vec Sf_new{}; Sf_new.clear();
	static Vec dSf_new_dwl{}; dSf_new_dwl.clear();
	static Vec dSf_new_dwr{}; dSf_new_dwr.clear();
	static Vec Upf_new{}; Upf_new.clear();

	for (size_t i = 0; i != Wnf.size(); i++) {
		if (Wnf[i] > Slimnf[i]) {
			Sf_new.push_back(Slimnf[i] / DU_DS);
			Upf_new.push_back(Upf[i]);
			dSf_new_dwl.push_back(dSlimnf_dwl[i] / DU_DS);
			dSf_new_dwr.push_back(dSlimnf_dwr[i] / DU_DS);
		}
		else if (Wnf[i] > eyn) {
			Sf_new.push_back(Wnf[i] / DU_DS);
			Upf_new.push_back(Upf[i]);
			dSf_new_dwl.push_back(dWnf_dwl[i] / DU_DS);
			dSf_new_dwr.push_back(dWnf_dwr[i] / DU_DS);
		}
		else {
			Sf_new.push_back(ey);
			Upf_new.push_back(Upf[i] + Wnf[i] - eyn);
			dSf_new_dwl.push_back(0.);
			dSf_new_dwr.push_back(0.);
		}
	}

	static Vec DSf{}; DSf.clear();
	static Vec dDSf_dwl{}; dDSf_dwl.clear();
	static Vec dDSf_dwr{}; dDSf_dwr.clear();

	for (size_t i = 0; i != Sf_new.size(); i++) {
		DSf.push_back(Sf_new[i] - Slimnf[i] / DU_DS);
		dDSf_dwl.push_back(dSf_new_dwl[i] - dSlimnf_dwl[i] / DU_DS);
		dDSf_dwr.push_back(dSf_new_dwr[i] - dSlimnf_dwr[i] / DU_DS);
	}

	//Slopes

	static Vec KSf_new{}; KSf_new.clear();
	static Vec dKSf_new_dwl{}; dKSf_new_dwl.clear();
	static Vec dKSf_new_dwr{}; dKSf_new_dwr.clear();
	static Vec KUpf_new{}; KUpf_new.clear();
	static Vec DKSf{}; DKSf.clear();
	static Vec dDKSf_dwl{}; dDKSf_dwl.clear();
	static Vec dDKSf_dwr{}; dDKSf_dwr.clear();

	for (size_t i = 0; i != intcats.size(); i++) {
		if (intcats[i] == 0) {
			KSf_new.push_back(KSlimf[i]);
			dKSf_new_dwl.push_back(0.);
			dKSf_new_dwr.push_back(0.);
			KUpf_new.push_back(KUpf[i]);
			DKSf.push_back(0.);
			dDKSf_dwl.push_back(0.);
			dDKSf_dwr.push_back(0.);
		}
		else if (intcats[i] == 1) {
			KSf_new.push_back(KWnf[i] / DU_DS);
			dKSf_new_dwl.push_back(dKWnf_dwl[i] / DU_DS);
			dKSf_new_dwr.push_back(dKWnf_dwr[i] / DU_DS);
			KUpf_new.push_back(KUpf[i]);
			DKSf.push_back(KSf_new[i] - KSlimf[i]);
			dDKSf_dwl.push_back(dKSf_new_dwl[i]);
			dDKSf_dwr.push_back(dKSf_new_dwr[i]);
		}
		else {
			KSf_new.push_back(0.);
			dKSf_new_dwl.push_back(0.);
			dKSf_new_dwr.push_back(0.);
			KUpf_new.push_back(KUpf[i] + KWnf[i]);
			DKSf.push_back(KSf_new[i] - KSlimf[i]);
			dDKSf_dwl.push_back(dKSf_new_dwl[i]);
			dDKSf_dwr.push_back(dKSf_new_dwr[i]);
		}
	}

	// Simplification

	static double REMLIM = 1.0e-16;

	ys_new.clear();
	s_new.clear();
	ks_new.clear();
	ys_cats.clear();
	dys_new_dwl.clear();
	dys_new_dwr.clear();
	ds_new_dwl.clear();
	ds_new_dwr.clear();
	dks_new_dwl.clear();
	dks_new_dwr.clear();
	yup_new.clear();
	up_new.clear();
	kup_new.clear();
	ydks.clear();
	dydks_dwl.clear();
	dydks_dwr.clear();
	dks.clear();
	ddks_dwl.clear();
	ddks_dwr.clear();
	ds.clear();
	dds_dwl.clear();
	dds_dwr.clear();

	ys_new.push_back(Yf[0]);
	s_new.push_back(Sf_new[0]);
	ks_new.push_back(KSf_new[0]);
	dys_new_dwl.push_back(dYf_dwl[0]);
	ds_new_dwl.push_back(dSf_new_dwl[0]);
	dks_new_dwl.push_back(dKSf_new_dwl[0]);
	dys_new_dwr.push_back(dYf_dwr[0]);
	ds_new_dwr.push_back(dSf_new_dwr[0]);
	dks_new_dwr.push_back(dKSf_new_dwr[0]);
	ys_cats.push_back(intcats[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(KSf_new[i - 1] - KSf_new[i]) <= REMLIM) {
			continue;
		}
		ys_new.push_back(Yf[i]);
		s_new.push_back(Sf_new[i]);
		ks_new.push_back(KSf_new[i]);
		ys_cats.push_back(intcats[i]);
		dys_new_dwl.push_back(dYf_dwl[i]);
		ds_new_dwl.push_back(dSf_new_dwl[i]);
		dks_new_dwl.push_back(dKSf_new_dwl[i]);
		dys_new_dwr.push_back(dYf_dwr[i]);
		ds_new_dwr.push_back(dSf_new_dwr[i]);
		dks_new_dwr.push_back(dKSf_new_dwr[i]);
	}
	ys_new.push_back(Yf[Yf.size() - 1]);
	s_new.push_back(Sf_new[Sf_new.size() - 1]);
	dys_new_dwl.push_back(dYf_dwl[dYf_dwl.size() - 1]);
	ds_new_dwl.push_back(dSf_new_dwl[dSf_new_dwl.size() - 1]);
	dys_new_dwr.push_back(dYf_dwr[dYf_dwr.size() - 1]);
	ds_new_dwr.push_back(dSf_new_dwr[dSf_new_dwr.size() - 1]);

	yup_new.push_back(Yf[0]);
	up_new.push_back(Upf_new[0]);
	kup_new.push_back(KUpf_new[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(KUpf_new[i - 1] - KUpf_new[i]) <= REMLIM) {
			continue;
		}
		yup_new.push_back(Yf[i]);
		up_new.push_back(Upf_new[i]);
		kup_new.push_back(KUpf_new[i]);
	}
	yup_new.push_back(Yf[Yf.size() - 1]);
	up_new.push_back(Upf_new[Upf_new.size() - 1]);

	ydks.push_back(Yf[0]);
	dks.push_back(DKSf[0]);
	ds.push_back(DSf[0]);
	dydks_dwl.push_back(dYf_dwl[0]);
	ddks_dwl.push_back(dDKSf_dwl[0]);
	dydks_dwr.push_back(dYf_dwr[0]);
	ddks_dwr.push_back(dDKSf_dwr[0]);
	dds_dwl.push_back(dDSf_dwl[0]);
	dds_dwr.push_back(dDSf_dwr[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(DKSf[i - 1] - DKSf[i]) <= REMLIM) {
			continue;
		}
		ydks.push_back(Yf[i]);
		ds.push_back(DSf[i]);
		dks.push_back(DKSf[i]);
		dydks_dwl.push_back(dYf_dwl[i]);
		ddks_dwl.push_back(dDKSf_dwl[i]);
		dydks_dwr.push_back(dYf_dwr[i]);
		ddks_dwr.push_back(dDKSf_dwr[i]);
		dds_dwl.push_back(dDSf_dwl[i]);
		dds_dwr.push_back(dDSf_dwr[i]);
	}
	ydks.push_back(Yf[Yf.size() - 1]);
	ds.push_back(DSf[DSf.size() - 1]);
	dydks_dwl.push_back(dYf_dwl[dYf_dwl.size() - 1]);
	dydks_dwr.push_back(dYf_dwr[dYf_dwr.size() - 1]);
	dds_dwl.push_back(dDSf_dwl[dDSf_dwl.size() - 1]);
	dds_dwr.push_back(dDSf_dwr[dDSf_dwr.size() - 1]);

	return;
}

void interval_dists_K(const Vector& Yw, const Vector& W, const Vector& Yup_com, const Vector& Up_com, const Vector& Kup_com, const Vector& Ys_com, const Vector& S_com, const Vector& Ks_com, double ey, double beta_Dt,
	Vector& Ys, Vector& S, Vector& Ks, Vector& Yup_new, Vector& Up_new, Vector& Kup_new,
	Matrix& dYs_dW, Matrix& dS_dW, Matrix& dKs_dW, Vecint& Ys_cats, Vector& Ydks, Vector& Dks, Matrix& dYdks_dW, Matrix& dDks_dW, Vector& DS, Matrix& dDS_dW)
{
	static VecVecOS Yup_ints{}; Yup_ints.clear();
	static VecVecOS Up_ints{}; Up_ints.clear();
	static VecVecOS Kup_ints{}; Kup_ints.clear();
	static VecVecOS Ys_ints{}; Ys_ints.clear();
	static VecVecOS S_ints{}; S_ints.clear();
	static VecVecOS Ks_ints{}; Ks_ints.clear();

	Up_interval_split_K(Yup_com, Up_com, Kup_com, Yw, Yup_ints, Up_ints, Kup_ints);
	Up_interval_split_K(Ys_com, S_com, Ks_com, Yw, Ys_ints, S_ints, Ks_ints);

	//std::cout << Ys_ints[3].transpose() << std::endl;
	//std::cout << S_ints[3].transpose() << std::endl;
	//std::cout << Ks_ints[3].transpose() << std::endl;

	VecVec ys_list(W.Size() - 1, std::vector<double>{});
	VecVec s_list(W.Size() - 1, std::vector<double>{});
	VecVec ks_list(W.Size() - 1, std::vector<double>{});
	VecVec yup_new_list(W.Size() - 1, std::vector<double>{});
	VecVec up_new_list(W.Size() - 1, std::vector<double>{});
	VecVec kup_new_list(W.Size() - 1, std::vector<double>{});

	VecVec dys_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dys_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVecint ys_cats_list(W.Size() - 1, std::vector<int>{});

	VecVec ydks_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_list(W.Size() - 1, std::vector<double>{});
	VecVec dydks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dydks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec ddks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec ddks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec dds_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dds_dwr_list(W.Size() - 1, std::vector<double>{});

	for (size_t i = 0; i != W.Size() - 1; i++) {

		interval_interior_K(W[i], W[i + 1], ey, Yw[i + 1] - Yw[i], Up_ints[i], Yup_ints[i], Kup_ints[i],
			Ys_ints[i], S_ints[i], Ks_ints[i], beta_Dt,
			ys_list[i], s_list[i], ks_list[i], ys_cats_list[i], yup_new_list[i], up_new_list[i], kup_new_list[i],
			dys_dwl_list[i], dys_dwr_list[i], ds_dwl_list[i], ds_dwr_list[i], dks_dwl_list[i], dks_dwr_list[i],
			ydks_list[i], dks_list[i], dydks_dwl_list[i], dydks_dwr_list[i], ddks_dwl_list[i], ddks_dwr_list[i],
			ds_list[i], dds_dwl_list[i], dds_dwr_list[i]);
	}

	Ys = interval_join(ys_list);
	S = interval_join(s_list);
	Ks = array_join(ks_list);
	Yup_new = interval_join(yup_new_list);
	Up_new = interval_join(up_new_list);
	Kup_new = array_join(kup_new_list);

	Ydks = interval_join(ydks_list);
	Dks = array_join(dks_list);
	DS = interval_join(ds_list);

	static VecMatOS dYs_dW_list{}; dYs_dW_list.clear();
	static VecMatOS dS_dW_list{}; dS_dW_list.clear();
	static VecMatOS dKs_dW_list{}; dKs_dW_list.clear();
	static VecMatOS dYdks_dW_list{}; dYdks_dW_list.clear();
	static VecMatOS dDks_dW_list{}; dDks_dW_list.clear();
	static VecMatOS dDS_dW_list{}; dDS_dW_list.clear();

	for (size_t i = 0; i != W.Size() - 1; i++) {

		Vec dwl_dW(W.Size()); dwl_dW[i] = 1.0;
		Vec dwr_dW(W.Size()); dwr_dW[i + 1] = 1.0;
		Matrix dys_dW = Matrix(dys_dwl_list[i].size(), W.Size());
		Matrix ds_dW = Matrix(ds_dwl_list[i].size(), W.Size());
		Matrix dks_dW = Matrix(dks_dwl_list[i].size(), W.Size());
		Matrix dydks_dW = Matrix(dydks_dwl_list[i].size(), W.Size());
		Matrix ddks_dW = Matrix(ddks_dwl_list[i].size(), W.Size());
		Matrix dds_dW = Matrix(dds_dwl_list[i].size(), W.Size());
		for (size_t l = 0; l != W.Size(); l++) {
			for (size_t k = 0; k != dys_dwl_list[i].size(); k++) {
				dys_dW(k, l) += dys_dwl_list[i][k] * dwl_dW[l];
				dys_dW(k, l) += dys_dwr_list[i][k] * dwr_dW[l];
				ds_dW(k, l) += ds_dwl_list[i][k] * dwl_dW[l];
				ds_dW(k, l) += ds_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dks_dwl_list[i].size(); k++) {
				dks_dW(k, l) += dks_dwl_list[i][k] * dwl_dW[l];
				dks_dW(k, l) += dks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dydks_dwl_list[i].size(); k++) {
				dydks_dW(k, l) += dydks_dwl_list[i][k] * dwl_dW[l];
				dydks_dW(k, l) += dydks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != ddks_dwl_list[i].size(); k++) {
				ddks_dW(k, l) += ddks_dwl_list[i][k] * dwl_dW[l];
				ddks_dW(k, l) += ddks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dds_dwl_list[i].size(); k++) {
				dds_dW(k, l) += dds_dwl_list[i][k] * dwl_dW[l];
				dds_dW(k, l) += dds_dwr_list[i][k] * dwr_dW[l];
			}
		}

		dYs_dW_list.push_back(dys_dW);
		dS_dW_list.push_back(ds_dW);
		dKs_dW_list.push_back(dks_dW);
		dYdks_dW_list.push_back(dydks_dW);
		dDks_dW_list.push_back(ddks_dW);
		dDS_dW_list.push_back(dds_dW);
	}
	dYs_dW = interval_join(dYs_dW_list);
	dS_dW = interval_join(dS_dW_list);
	dKs_dW = array_join(dKs_dW_list);
	dYdks_dW = interval_join(dYdks_dW_list);
	dDks_dW = array_join(dDks_dW_list);
	dDS_dW = interval_join(dDS_dW_list);

	Ys_cats.clear();
	for (size_t i = 0; i != ys_cats_list.size(); i++) {
		for (size_t j = 0; j != ys_cats_list[i].size(); j++) {
			Ys_cats.push_back(ys_cats_list[i][j]);
		}
	}

	return;

}

void Ys_cats_dist_calc(const VecVecint& Ys_cats, Vecint& Ys_cats_dist) {
	Ys_cats_dist.clear();
	for (size_t i = 0; i != Ys_cats.size(); i++) {
		for (size_t j = 0; j != Ys_cats[i].size(); j++) {
			Ys_cats_dist.push_back(Ys_cats[i][j]);
		}
	}
}