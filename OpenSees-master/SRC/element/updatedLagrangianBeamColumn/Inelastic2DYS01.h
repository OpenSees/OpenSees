#ifndef Inelastic2DYS01_H
#define Inelastic2DYS01_H

#include "InelasticYS2DGNL.h"

class Inelastic2DYS01 : public InelasticYS2DGNL
{
public:
	Inelastic2DYS01(int tag, double A, double E, double Iz,
                    int Nd1, int Nd2,
                    YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
                    int rf_algo=-1, bool islinear=false, double rho=0.0);

	~Inelastic2DYS01();

protected:
	virtual void getLocalStiff(Matrix &K);

private:
double A;
double Iz;
double E;
double damageFactorEnd1;
double damageFactorEnd2;
double fpeak;
};

#endif

