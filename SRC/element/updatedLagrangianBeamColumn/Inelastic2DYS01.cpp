// Inelastic2DYS01.cpp
//////////////////////////////////////////////////////////////////////

#include "Inelastic2DYS01.h"
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_Inelastic2DYS01)
{
    //cerr << "Press key to continue...\n";
    //cin.get();

    if (OPS_GetNumRemainingInputArgs() < 9)
    {
	opserr << "WARNING insufficient arguments\n";
	opserr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

	return 0;
    }

    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid element2dYS int inputs" << endln;
    }
    
    int tag = idata[0];
    int ndI = idata[1];
    int ndJ = idata[2];

    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid element2dYS double inputs" << endln;
    }
    
    double A = data[0];
    double E = data[1];
    double I = data[2];

    
    numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid element2dYS int inputs" << endln;
    }
    
    int ysID1 = idata[0];
    int ysID2 = idata[1];
    int rf_algo = idata[2];

    YieldSurface_BC *theYS1 = OPS_getYieldSurface_BC(ysID1);
    if(theYS1 == 0)
    {
	opserr << "WARNING element2dYS: " << tag << "\n";
	opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
	return 0;
    }

    YieldSurface_BC *theYS2 = OPS_getYieldSurface_BC(ysID2);
    if(theYS2 == 0)
    {
	opserr << "WARNING element2dYS: " << tag << "\n";
	opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
	return 0;
    }

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated
    return new Inelastic2DYS01(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Inelastic2DYS01::Inelastic2DYS01(int tag, double a, double e, double iz,
                                 int Nd1, int Nd2,
				                 YieldSurface_BC *ysEnd1,
				                 YieldSurface_BC *ysEnd2,
				                 int rf_algo, bool islinear, double rho
				                 )
			    :InelasticYS2DGNL (tag, Nd1, Nd2, ysEnd1, ysEnd2,
			                   rf_algo, islinear, rho),
			     A(a), E(e), Iz(iz),
			     damageFactorEnd1(0.0), damageFactorEnd2(0.0), fpeak(0.0)
{
  massDof = A*L*rho;
  massDof = massDof/2;

}

Inelastic2DYS01::~Inelastic2DYS01()
{
	//does nothing
}

// very simple element
// provides the elastic stiffness based on the average degraded
// Iz at each ends
void Inelastic2DYS01::getLocalStiff(Matrix &K)
{

//	double i1Factor = (1 - damageFactorEnd1); // this->getDegradeFactor(1);
//	double i2Factor = (1 - damageFactorEnd2); // this->getDegradeFactor(2);
//
//	double iFactor = (i1Factor + i2Factor)/2;
//	iFactor = 1.0;
	double iz = Iz; //*iFactor;
	
    double	EIbyL = E*iz/L;

    K(0, 1) = K(0, 2) = K(0, 4) = K(0, 5)=0;
    K(1, 0) = K(1, 3) =0;
    K(2, 0) = K(2, 3) =0;
    K(3, 1) = K(3, 2) = K(3, 4) = K(3, 5)=0;
    K(4, 0) = K(4, 3) =0;
    K(5, 0) = K(5, 3) =0;

	K(0,0) = K(3,3) = (A/iz)*(EIbyL);
	K(0,3) = K(3,0) = (-A/iz)*(EIbyL);
	K(1,1) = K(4,4) = (12/(L*L))*(EIbyL);
	K(1,4) = K(4,1) = (-12/(L*L))*(EIbyL);
	K(1,2) = K(2,1) = K(1,5) = K(5,1) = (6/L)*(EIbyL);
	K(2,4) = K(4,2) = K(4,5) = K(5,4) = (-6/L)*(EIbyL);
	K(2,2) = K(5,5) = 4*(EIbyL);
	K(2,5) = K(5,2) = 2*(EIbyL);

}//getLocalStiff



