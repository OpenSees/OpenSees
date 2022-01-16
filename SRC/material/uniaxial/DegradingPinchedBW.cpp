#include <elementAPI.h>
#include <DegradingPinchedBW.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>


void * OPS_ADD_RUNTIME_VPV(OPS_DegradingPinchedBW)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData1[1];
  double dData[18];
  int    iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DegradingPinchedBW tag" << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new DegradingPinchedBW(iData1[0], dData[0], dData[1], dData[2], 
  dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], 
  dData[10], dData[11], dData[12], dData[13], dData[14], dData[15],
  dData[16], dData[17], iData2[0]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type DegradingPinchedBW\n";
    return 0;
  }

  return theMaterial;
}



DegradingPinchedBW::DegradingPinchedBW(int tag, 
       double p_m,
	   double p_Fy,
	   double p_xu,
	   double p_alpha,
	   double p_ko,
	   double p_n,
	   double p_eta,
	   double p_beta,
	   double p_rhoeps,
	   double p_rhox,
	   double p_phi,
	   double p_deltak,
	   double p_deltaf,
	   double p_sigma,
	   double p_u,
	   double p_epsp,
	   double p_rhop,
	   double ptolerance,
	   int pMaxNumIter)
 :UniaxialMaterial(tag, 0),
  xmax(0.0), xmaxp(0.0), m(p_m), Fy(p_Fy), xu(p_xu),alpha(p_alpha), 
  ko(p_ko), n(p_n), eta(p_eta), beta(p_beta), rhoeps(p_rhoeps), rhox(p_rhox), 
  phi(p_phi), deltak(p_deltak), deltaf(p_deltaf), sigma(p_sigma), 
  u(p_u), epsp(p_epsp), rhop(p_rhop), tolerance(ptolerance),
  maxNumIter(pMaxNumIter)
{
  // Initialize variables
  this->revertToStart();
}


DegradingPinchedBW::~DegradingPinchedBW()
{
  // does nothing
}

double 
DegradingPinchedBW::signum(double value)
{
  if (value > 0.0) {
    return 1.0;
  }
  else {
    return -1.0;
  }
}


int 
DegradingPinchedBW::setTrialStrain (double strain, double strainRate)
{
	// Set trial strain and compute strain increment
	Tstrain = strain;
	double dStrain = Tstrain - Cstrain;



	// Initial declarations (make sure not to declare class variables here!)
	double TDamIndex, Tpk, TA, Tbetak, Tbetaf, TPow, Psi, Phi, f, Tznew, Tzold;
	double Te_z, TDamIndex_z, Tpk_z, TA_z, Tbetak_z, Tbetaf_z, TPow_z, Phi_z, f_z;
	double Te_x, TDamIndex_x, Tpk_x, TA_x, Tbetak_x, Tbetaf_x, Phi_x, f_x;
	double Te_z_x, TDamIndex_z_x, Tpk_z_x, TA_z_x, Tbetak_z_x, Tbetaf_z_x, Phi_z_x, f_z_x;
	double gp, feps, kp;


	// Newton-Raphson scheme to solve for z_{i+1} = z1
	int count = 0;
	double startPoint = 0.01;
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;  
	while ( ( fabs(Tzold-Tznew) > tolerance ) && count < maxNumIter) {

  
		//Evaluate function f
		gp = exp((-0.5)*pow(Cstrain/sigma, u));
		feps = 1.0 - exp((-0.5)*pow(Ce/epsp, 8));
		kp = ko*(1.0 - feps*gp*rhop);
		Te = Ce + (1.0-alpha)*kp/m*dStrain*Tz;
		if ((fabs(Tstrain)) >= (fabs(xmax))) {
			xmax = fabs(Tstrain);
		}
		else {
			xmax = xmaxp;
		}
		TDamIndex = rhoeps*Te*m/(Fy*xu) + rhox*fabs(xmax)/xu;
		Tpk = exp(-phi*TDamIndex);
		TA = exp(-deltak*TDamIndex*Tpk);
		Tbetak = beta*TA;
		Tbetaf = exp(n*deltaf*TDamIndex);
		Psi = (eta + signum(dStrain*Tz));
		TPow = pow(fabs(Tz),(n));
		Phi = TA - TPow*Tbetak*Tbetaf*Psi;
		f = Tz - Cz - Phi*dStrain;

		// Evaluate function derivative f_z ( _z = derivative with respect to z_n+1 )
		Te_z = (1.0-alpha)*kp*dStrain/m;
		TDamIndex_z = rhoeps*Te_z*m/(Fy*xu);
		Tpk_z = Tpk*(-phi*TDamIndex_z);
		TA_z = TA*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z);
		Tbetak_z = beta*TA_z;
		Tbetaf_z = Tbetaf*(n*deltaf*TDamIndex_z);
		if (Tz == 0.0) {
			TPow = 0.0;
			TPow_z = 0.0;
		}
		else {
			TPow = pow(fabs(Tz),(n));
			TPow_z = n*pow(fabs(Tz),(n-1))*signum(Tz);
		}
		Phi_z = TA_z - (TPow_z*Tbetak*Tbetaf + TPow*Tbetak_z*Tbetaf + TPow*Tbetak*Tbetaf_z)*Psi;
		f_z = 1.0 - Phi_z*dStrain;
		
		
		// Issue warning if derivative is zero
		if ( fabs(f_z)<1.0e-10 ) {
			opserr << "WARNING: DegradingPinchedBW::setTrialStrain() -- zero derivative " << endln
				<< " in Newton-Raphson scheme" << endln;
		}

		// Take a Newton step
		Tznew = Tz - f/f_z;


		// Update the root (but the keep the old for convergence check)
		Tzold = Tz; 
		Tz = Tznew;

		// Update counter
		count++;


		// Issue warning if we didn't converge
		if (count == maxNumIter) {

			opserr << "WARNING: DegradingPinchedBW::setTrialStrain() -- did not" << endln
				<< " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
				<< " and norm: " << fabs(Tzold-Tznew) << endln;
		}

		// Compute stress
		
		Tstress = alpha*kp*Tstrain + (1-alpha)*kp*Tz;


		// Compute deterioration parameters
	
		Te = Ce + (1-alpha)*kp*dStrain*Tz/m;
		TDamIndex = rhoeps*Te*m/(Fy*xu) + rhox*fabs(xmax)/xu;
		Tpk = exp(-phi*TDamIndex);
		TA = exp(-deltak*TDamIndex*Tpk);
		Tbetak = beta*TA;
		Tbetaf = exp(n*deltaf*TDamIndex);
		if (Tz == 0.0) {
			TPow = 0.0;
			TPow_z = 0.0;
		}
		else {
			TPow = pow(fabs(Tz),(n));
			TPow_z = n*pow(fabs(Tz),(n-1))*signum(Tz);
		}
        Psi = eta + signum(dStrain*Tz);
		Phi =  TA - TPow*Tbetak*Tbetaf*Psi;
		
		// Compute tangent ( _x = derivative with respect to x_n+1 )
		
		// Compute the derivative of f with respect to x_n+1
		if (Tz != 0.0) {
			Te_x = (1 - alpha)*kp*Tz/m;
			if (xmax == Tstrain) {
				TDamIndex_x = rhoeps*Te_x*m/(Fy*xu) + rhox/xu;
			}
			else {
				TDamIndex_x = rhoeps*Te_x*m/(Fy*xu);
			}
			Tpk_x = Tpk*(-phi*TDamIndex_x);
			TA_x = TA*(-deltak*TDamIndex_x*Tpk - deltak*TDamIndex*Tpk_x);
			Tbetak_x = beta*TA_x;
			Tbetaf_x = Tbetaf*(n*deltaf*TDamIndex_x);
			Phi_x = TA_x - (TPow*Tbetak_x*Tbetaf + TPow*Tbetak*Tbetaf_x)*Psi;
			f_x = -Phi - Phi_x*dStrain;
			
		// Compute the derivative of f_z with respect to x_n+1
			Te_z_x = (1 - alpha)*kp/m;
			TDamIndex_z_x = rhoeps*Te_z_x*m/(Fy*xu);
			Tpk_z_x = Tpk_x*(-phi*TDamIndex_z) + Tpk*(-phi*TDamIndex_z_x);
			TA_z_x = TA_x*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z) - TA*(deltak*TDamIndex_z_x*Tpk + deltak*TDamIndex_z*Tpk_x + deltak*TDamIndex_x*Tpk_z + deltak*TDamIndex*Tpk_z_x);
			Tbetak_z_x = beta*TA_z_x;
			Tbetaf_z_x = Tbetaf_x*(n*deltaf*TDamIndex_z) + Tbetaf*(n*deltaf*TDamIndex_z_x);
			Phi_z_x = TA_z_x - (TPow_z*Tbetak*Tbetaf_x + TPow_z*Tbetak_x*Tbetaf + TPow*Tbetaf_z_x*Tbetak + TPow*Tbetaf_z*Tbetak_x + TPow*Tbetaf_x*Tbetak_z + TPow*Tbetaf*Tbetak_z_x)*Psi;
			f_z_x = -Phi_z - Phi_z_x*dStrain ;
		
		
		// Compute tangent	
			double DzDeps = -(f_x*f_z - f*f_z_x)/(pow(f_z,2));
			Ttangent = alpha*kp + (1-alpha)*kp*DzDeps;
			//Ttangent = Tstress/Tstrain;
			
		}
		else {
			Ttangent = alpha*ko + (1-alpha)*ko;
		}
		
	

	}
	
	
    return 0;
}

double 
DegradingPinchedBW::getStress(void)
{
    return Tstress;
}

double 
DegradingPinchedBW::getInitialTangent(void)
{
    return ( alpha*ko + (1-alpha)*ko );
}


double 
DegradingPinchedBW::getTangent(void)
{
    return Ttangent;
}

double 
DegradingPinchedBW::getStrain(void)
{
    return Tstrain;
}

int 
DegradingPinchedBW::commitState(void)
{
    // Commit trial history variables
	xmaxp = xmax;
	Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
DegradingPinchedBW::revertToLastCommit(void)
{
	xmax = xmaxp;
    return 0;
}

int 
DegradingPinchedBW::revertToStart(void)
{
    xmaxp = 0.0;
	Tstrain = 0.0;
	Cstrain = 0.0;
	Tz = 0.0;
	Cz = 0.0;
	Te = 0.0;
	Ce = 0.0;
	Tstress = 0.0;
	Ttangent = alpha*ko + (1-alpha)*ko;

    return 0;
}

UniaxialMaterial *
DegradingPinchedBW::getCopy(void)
{
    DegradingPinchedBW *theCopy =
	new DegradingPinchedBW(this->getTag(), m, Fy, xu, alpha, ko, n, eta, beta, rhoeps, 
	rhox, phi, deltak, deltaf, sigma, u, epsp, rhop, tolerance, maxNumIter);
    theCopy->xmax = xmax;
	theCopy->xmaxp = xmaxp;
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->Te = Te;
    theCopy->Ce = Ce;
    theCopy->Tstress = Tstress;
	theCopy->Ttangent = Ttangent;

    return theCopy;
}

int 
DegradingPinchedBW::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
  static Vector data(2);
  data(0) = this->getTag();
  data(1) = xmaxp;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "DegradingPinchedBW::sendSelf() - failed to send data\n";

  return res;
}

int 
DegradingPinchedBW::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
  static Vector data(2);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  this->setTag(int(data(0)));
  xmaxp = data(1);
  this->revertToLastCommit();
    
  return res;
}

void 
DegradingPinchedBW::Print(OPS_Stream &s, int flag)
{
    s << "DegradingPinchedBW, tag: " << this->getTag() << endln;
	s << "  m: " << m << endln;
	s << "  Fy: " << Fy << endln;
	s << "  xu: " << xu << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
    s << "  n: " << n << endln;
    s << "  eta: " << eta << endln;
    s << "  beta: " << beta << endln;
	s << "  rhoeps: " << rhoeps << endln;
	s << "  rhox: " << rhox << endln;
	s << "  phi: " << phi << endln;
	s << "  deltak: " << deltak << endln;
	s << "  deltaf: " << deltaf << endln;
	s << "  sigma: " << sigma << endln;
	s << "  u: " << u << endln;
	s << "  epsp: " << epsp << endln;
	s << "  rhop: " << rhop << endln;
}


