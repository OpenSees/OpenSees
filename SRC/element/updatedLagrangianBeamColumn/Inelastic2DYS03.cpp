// Inelastic2DYS03.cpp
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "Inelastic2DYS03.h"
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_Inelastic2DYS03)
{
    if (OPS_GetNumRemainingInputArgs() < 9)
    {
	opserr << "WARNING insufficient arguments\n";
	opserr << "element element2dYS03 tag? Nd1? Nd2? A_ten? A_com? E? IzPos? IzNeg? ysID1? ysID2? algo?";

	return 0;
    }

    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid element2dYS int inputs" << endln;
	return 0;
    }
    
    int tag = idata[0];
    int ndI = idata[1];
    int ndJ = idata[2];

    double data[5];
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid element2dYS double inputs" << endln;
	return 0;
    }

    double aTens = data[0];
    double aComp = data[1];
    double E = data[2];
    double Ipos = data[3];
    double Ineg = data[4];


    numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid element2dYS int inputs" << endln;
	return 0;
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

    return new Inelastic2DYS03(tag, aTens, aComp, E,
			       Ipos, Ineg, ndI, ndJ,
			       theYS1, theYS2, rf_algo);

}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Inelastic2DYS03::Inelastic2DYS03(int tag, double a_ten, double a_com,
                 double e, double iz_pos, double iz_neg,
                 int Nd1, int Nd2,
				 YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
				 int rf_algo, bool islinear, double rho)

			:InelasticYS2DGNL(tag, Nd1, Nd2,
			ysEnd1, ysEnd2, rf_algo, islinear, rho),
			Atens(a_ten), Acomp(a_com), E(e), IzPos(iz_pos), IzNeg(iz_neg), ndisp(6), ndisp_hist(6)
{
	ndisp_hist.Zero();
	ndisp.Zero();
}

Inelastic2DYS03::~Inelastic2DYS03()
{
	//does nothing
}


void Inelastic2DYS03::getLocalStiff(Matrix &K)
{
	double L1,L2,I1,I2,A;
	Vector ndisp_inc(6);

    getIncrNaturalDisp(ndisp_inc);
	ndisp = ndisp_hist + ndisp_inc;

	// getTrialnaturalDisp(ndisp);
	
    opserr << ndisp;
    // opserr << ndisp_hist;
    // opserr << ndisp_inc;
    opserr << "\a";

	if(ndisp(2)*ndisp(5) < 0  || fabs(ndisp(2)*ndisp(5)) < 1e-10) {	//if single curvature
	  L1 = L;
	  L2 = 0;
	  if(ndisp(2) > 0 || ndisp(5) < 0)
	    I1 = I2 = IzNeg;
	  else I1 = I2 = IzPos;
	} else {		//double curvature
		
		if((fabs(ndisp(2)) + fabs(ndisp(5)) < 1e-10))
			L1 = 0;
		else
			L1 = (fabs(ndisp(2))*L) / (fabs(ndisp(2)) + fabs(ndisp(5)));
			L2 = L - L1;
		if(ndisp(2) > 0) {
			I1 = IzNeg;
			I2 = IzPos;
		} else {
			I1 = IzPos;
			I2 = IzNeg;
		}
	}

	opserr << L1 << "  " << L2 << "\n";
	
	if(ndisp(3) < 0) //element is in compression
		A = Acomp;
	else					 //element is in tension
		A = Atens;

	//some factors used in stiffness matrix
	double X1 = 
I2*I2*L1*L1*L1*L1+4*I2*L1*L1*L1*L2*I1+6*I2*L1*L1*L2*L2*I1+4*I2*L1*L2*L2*L2*I1+L2*L2*L2*L2*I1*I1;
	double X2 = (I2*I1*(I2*L1*L1+2*L2*I1*L1+L2*L2*I1))/(X1);
	double X3 = ((L1*I2+L2*I1)*I2*I1)/(X1);
	double X4 = ((I2*L1*L1+2*I2*L1*L2+L2*L2*I1)*I2*I1)/(X1);

	//zeros in stiffness matrix
    K(0,1) = K(0,2) = K(0,4) = K(0,5)=0;
    K(1,0) = K(1,3) = 0;
    K(2,0) = K(2,3) = 0;
    K(3,1) = K(3,2) = K(3,4) = K(3,5)=0;
    K(4,0) = K(4,3) = 0;
    K(5,0) = K(5,3) = 0;

	//axial terms
  	K(0,0) = K(3,3) = (A*E)/(L); 
  	K(0,3) = K(3,0) = (-A*E)/(L);

	//shear and moment terms
  	K(1,1) = K(4,4) = 12*E*X3;
  	K(1,4) = K(4,1) = -12*E*X3;
  	K(1,2) = K(2,1) = 6*E*X2;
	K(1,5) = K(5,1) = 6*E*X4;
  	K(2,4) = K(4,2) = -6*E*X2;
	K(4,5) = K(5,4) = -6*E*X4;
  	K(2,2) = (4*E*I2*I1*(I2*L1*L1*L1 + 3*L2*I1*L1*L1 + 3*L2*L2*I1*L1 + 
              L2*L2*L2*I1))/X1;
	K(5,5) = (4*E*I2*I1*(I2*L1*L1*L1 + 3*I2*L1*L1*L2 + 3*I2*L1*L2*L2 + 
              L2*L2*L2*I1))/X1;
  	K(2,5) = K(5,2) = 2*E*I2*I1*(I2*L1*L1*L1 + 3*I2*L1*L1*L2 + 3*L2*L2*I1*L1 
             + L2*L2*L2*I1)/X1;

   opserr << "\nInelastic2DYS03::getLocalStiff(..) = \n" << K;
}//getLocalStiff


int Inelastic2DYS03::commitState()
{
	// first let the super classes do their stuff
    this->InelasticYS2DGNL::commitState();
	// now set the commit natural disps
	ndisp_hist = ndisp;
	return 0;

}

