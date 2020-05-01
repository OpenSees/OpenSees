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

// Written in Matlab: Thanh Do
// Created: 07/16
                                                                        
#include <PlasticDamageConcrete3d.h>           
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>

static Vector Iv6(6); 
static Matrix Ivp(6,6); 
static Matrix Idp(6,6); 
static Matrix I(6,6);
static Matrix Id(6,6); 


void *
OPS_NewPlasticDamageConcrete3d(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 5 || numArgs > 9) {
    opserr << "Want: nDMaterial PlasticDamageConcrete3d $tag $E $nu $ft $fc <$beta $Ap $An $Bn>\n";
    return 0;	
  }
  
  int iData[1];
  double dData[8];
  dData[4] = 0.6;
  dData[5] = 0.5;
  dData[6] = 2.0;
  dData[7] = 0.75;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
    return 0;
  }
  
  numData = numArgs - 1;;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new PlasticDamageConcrete3d(iData[0], 
					    dData[0], dData[1], dData[2],dData[3], 
					    dData[4], dData[5], dData[6],dData[7]);

  return theMaterial;
}


PlasticDamageConcrete3d::PlasticDamageConcrete3d(int tag,
						 double _e, 
						 double _nu, 
						 double _ft,
						 double _fc, 
						 double _beta, 
						 double _Ap, 
						 double _An, 
						 double _Bn)
  :NDMaterial(tag,ND_TAG_PlasticDamageConcrete3d),
 E(_e), nu(_nu), ft(_ft), fc(_fc), beta(_beta), Ap(_Ap), An(_An), Bn(_Bn),
 eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
 epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
 Ce(6,6), C(6,6), Ccommit(6,6)
{
  eps.Zero();
  sig.Zero();
  sige.Zero();
  eps_p.Zero();
  sigeP.Zero();
  Ce.Zero();

  // additional material parameters
  double G   = E/2/(1+nu);        //shear modulus
  double  K   = E/3/(1-2*nu);     // bulk  modulus

  Iv6.Zero(); Iv6(0)=1.;Iv6(1)=1.;Iv6(2)=1.;

  Ivp.Zero(); 
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++)
      Ivp(i,j)=1.;

  Idp.Zero(); 
  I.Zero(); 
  Id.Zero(); 
  for (int i=0; i<6; i++) {
    Idp(i,i) = 1.;
    if (i<3) {
      I(i,i) = 1.0;
      Id(i,i) = 1.0;
    } else {
      I(i,i) = 0.5;
      Id(i,i) = 0.5;
    }
  }
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      Id(i,j)=Idp(i,j)-1/3.;
      Idp(i,j) = Id(i,j);
    }

  Ce.addMatrix(0.0, Ivp, K);
  Ce.addMatrix(1.0,  Id, 2.*G);
  
  C = Ce;
  
  double f2c = 1.16*fc;
  double k = sqrt(2.0)*(f2c - fc)/(2.*f2c - fc);

  //      % initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k+sqrt(2.0))*fc/sqrt(3.0));
      
  rp = rp0;
  rn = rn0;
  dp = 0.;
  dn = 0.;

  this->commitState();
}

PlasticDamageConcrete3d::PlasticDamageConcrete3d()
  :NDMaterial (0, ND_TAG_PlasticDamageConcrete3d),
   eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
   epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
   Ce(6,6), C(6,6), Ccommit(6,6)
{

}

PlasticDamageConcrete3d::~PlasticDamageConcrete3d ()
{

}

/*
%% Stress invariant function: octahedral normal and shear stresses
function [sigoct, tauoct] = StrsInvar (sig)

% normal stress
sigoct = (sig(1) + sig(2) + sig(3))/3;

% shear stress
J2 = ((sig(1) - sig(2))^2 + (sig(1) - sig(3))^2 + (sig(2) - sig(3))^2)/6 +...
     (sig(4))^2 + (sig(5))^2 + (sig(6))^2;
tauoct = (2/3*J2)^0.5;
end

%% Stress decomposition function: algebraic approach
function [sigpos, signeg, Qpos, Qneg] = StrsDecA (sig)

% positive and negative stress tensors
sigpos = (sig + abs(sig))/2;
signeg = sig - sigpos;

% projection tensors
Qpos = diag((1+sign(sig))/2);
Qneg = eye(6) - Qpos;
end
*/

void
StrsInvar(const Vector &sig, double &sigoct, double &tauoct) {
  sigoct = (sig(0) + sig(1) + sig(2))/3.;
  double J2 = (pow((sig(0) - sig(1)),2) + pow((sig(0) - sig(2)),2) + pow((sig(1) - sig(2)),2))/6. + 
    pow(sig(3),2) + pow(sig(4),2) + pow(sig(5),2);
  tauoct = sqrt(2./3.*J2);  
}

void StrsDecA(const Vector &sig, Vector &sigpos, Vector &signeg, Matrix &Qpos, Matrix &Qneg) {

  Qneg.Zero();
  Qpos.Zero();
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      sigpos(i) = sig(i);
      signeg(i) = 0.;
      Qpos(i,i) = 1;
      Qneg(i,i) = 0;
    } else if (sig(i)  < -1.e-8){
      sigpos(i) = 0.;
      signeg(i) = sig(i);
      Qpos(i,i) = 0;
      Qneg(i,i) = 1;
    } else {
      sigpos(i) = sig(i)/2.;
      signeg(i) = sig(i)/2.;
      Qpos(i,i) = 0.5;
      Qneg(i,i) = 0.5;
    }
  }
}

int
PlasticDamageConcrete3d::setTrialStrain(Vector const&v1, Vector const&v2){
  return this->setTrialStrain(v1);
}

int
PlasticDamageConcrete3d::setTrialStrain (const Vector &strain)
{
  //  opserr << "PlasticDamageConcrete3d::setTrialStrain: " << strain << endln;

  // bunch of Vectors and Matrices used in the method
  static Vector Depse_tr(6);
  static Vector Deps(6);  
  static Vector sige_tr(6);
  static Vector sigpos(6);
  static Vector signeg(6);
  static Matrix Qpos(6,6);
  static Matrix Qneg(6,6);
  static Vector L_tr(6);
  static Vector L_tr_temp(6);  
  static Vector Dnrm_Dsig(6);
  static Vector Dlam_Dsig(6);
  static Vector Dnrm_Deps(6);
  static Matrix Dsigpos_Deps(6,6);
  static Matrix Dsigneg_Deps(6,6);
  static Vector Ddp_Deps(6);
  static Vector Ddn_Deps(6);
  static Matrix Cbar(6,6);

  double f2c = 1.16*fc;
  double k = sqrt(2.0)*(f2c - fc)/(2.*f2c - fc);
  // initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k+sqrt(2.0))*fc/sqrt(3.0));

  double tol = 1.0e-5;
  
  // retrieve history variables
  eps_p = eps_pCommit;
  sigeP = sigePCommit;
  sigeP = sigeCommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

  // current strain
  eps = strain;

  // incremental strain
  Depse_tr = eps - eps_p;
  Deps = eps - epsCommit;

  //  opserr << "Depse_tr: " << Depse_tr;
  //  opserr << "Deps: " << Deps;

  // PLASTIC part
  // elastic trial stress
  sige_tr = sigeP + Ce*Deps;
  //  opserr << "Ce*Deps: " << Ce*Deps;
  //  opserr << "sigeP: " << sigeP;
  //  opserr << "sige_tr: " << sige_tr;

  // decomposition of trial stress tensor
  StrsDecA(sige_tr, sigpos, signeg, Qpos, Qneg);

  // compute stress variants of the negative stress tensor
  double sigoct, tauoct;
  StrsInvar(signeg, sigoct, tauoct);

  /*
  opserr << "sige_tr: " << sige_tr << endln;
  opserr << "sige: " << sige << endln;
  opserr << "sigeP: " << sigeP << endln;
  opserr << "signeg: " << signeg << endln;
  */

  // check tauneg
  //  opserr << "K: " << k << " sigoct: " << sigoct << " tauoct: " << tauoct << endln;
  double taun = sqrt(sqrt(3.0)*(k*sigoct + tauoct)); // % negative equivalent stress

  // Correction
  //  opserr << "taun1: " << taun << endln;

  if ((taun - rn) <= (tol*rn0)) { // elastic state, accept trial response
    sige = sige_tr;                                                    
    Cbar = Ce;                                                          
  } else {
    //  norm of trial effective stress
    double nrm = sqrt(pow(sige_tr(0),2) + pow(sige_tr(1),2) + pow(sige_tr(2),2) + 
		      2*pow(sige_tr(3),2) + 2*pow(sige_tr(4),2) + 2*pow(sige_tr(5),2)); 
    L_tr = sige_tr; L_tr/=nrm;    // normalized trial effective stress
    
    //  Deps_p = beta*E*(L_tr'*Deps)*Depse_tr/nrm;     // plastic strain increment
    double L_trDotDeps = L_tr ^ Deps;
    static Vector Deps_p(6);
    Deps_p = Depse_tr;
    Deps_p *= beta*E*L_trDotDeps/nrm;

    double lam    = 1 - beta*E/nrm * L_trDotDeps;      //  scale factor
    
    //    opserr << "lam: " << lam << endln;

    sige   = sige_tr; sige *= lam;                   //  corrected effective stress
        
    // check damage
    StrsDecA(sige, sigpos, signeg, Qpos, Qneg);         //  decompose the effective stress        
    StrsInvar(signeg, sigoct, tauoct);                  //  find octahedral stresses
    //    opserr << "k:" << k << "sigoct: " << sigoct << " tauoct: " << tauoct << endln;
    taun = sqrt(sqrt(3.)*(k*sigoct + tauoct));          //  negative equivalent stress

    /*
    opserr << "taun: " << taun << endln;
    opserr << "rn: " << rn << endln;
    opserr << "L_trDotDeps: " << L_trDotDeps << endln;
    */

    if ((taun - rn <= tol*rn0) || (L_trDotDeps <= 0)) { //  no damage or sige and eps in different direction
	sige = sige_tr;
	Cbar = Ce;            
    } else {
      eps_p  = eps_p + Deps_p;                          //  update plastic strain
        
      // tangent in effective space, Cbar
      //  L_tr_temp = [L_tr(1:3); 2*L_tr(4:6)];
      for (int i=0; i<3; i++) L_tr_temp(i) = L_tr(i);
      for (int i=3; i<6; i++) L_tr_temp(i) = 2*L_tr(i);

      double Dlam_Dnrm = 2*beta*E/pow(nrm,3)*(sige_tr^Deps);

      Dnrm_Dsig = L_tr_temp;
      Dlam_Dsig = Deps; Dlam_Dsig *= -beta*E/(nrm*nrm); // Dlam_Dsig = -beta*E/(nrm*nrm)*Deps;
      static Vector Dlam_Deps(6);
      Dlam_Deps = L_tr; Dlam_Deps *= -beta*E/nrm;       // Dlam_Deps = -beta*E/nrm*L_tr;
      // Dlam_Deps = Dlam_Dnrm * Ce * Dnrm_Dsig + Ce*Dlam_Dsig + Dlam_Deps; 
      Dlam_Deps = Dlam_Dnrm * (Ce * Dnrm_Dsig) + Ce*Dlam_Dsig + Dlam_Deps;

      /*
      opserr << "lam: " << lam << endln;
      opserr << "Ce: " << Ce << endln;
      opserr << "sige_tr: " << sige_tr << endln;
      opserr << "Dlam_Deps: " << Dlam_Deps << endln;
      */

      Cbar = lam*Ce + sige_tr % Dlam_Deps;   
    }
  }


  // DAMAGE part
  // decompose into positive and negative effective stress tensor
  StrsDecA(sige, sigpos, signeg, Qpos, Qneg);    // decompose the effective stress  
    
  // calculate equivalent stresses
  static Vector tmp(6);
  Ce.Solve(sigpos, tmp);
  double taup = sqrt(sigpos^tmp);                // positive equivalent stress



  StrsInvar(signeg, sigoct, tauoct);             // find octahedral stresses
  taun = sqrt((sqrt(3.)*(k*sigoct + tauoct)));   // negative equivalent stress

  double Ddp_Drp = 0.;

  // positive damage
  if ((taup - rp) <= (tol*rp0)) {                // no positive damage
    Ddp_Drp = 0;
  } else {                                       // positive damage evolves
    rp = taup;                                   // update rp = max(taup, rp)
    dp = 1 - rp0/rp * exp(Ap*(1 - rp/rp0));
    //    opserr << "dp: " << dp << " rpo: " << rp0 << " rp: " << rp << "Ap: " << Ap << endln;

    Ddp_Drp =  (Ap*rp + rp0)/(rp*rp) * exp(Ap*(1 - rp/rp0));               
    dp = dp*(1-tol);                             // cap the damage variable 
    Ddp_Drp = Ddp_Drp*(1-tol);        
    if (dp > 1-tol) {
      dp = 1- tol; Ddp_Drp = 0;
    }
  }

  //  opserr << "db: " << dp << endln;
  //  opserr << "Ddp_Drp: " << Ddp_Drp << endln;

  // negative damage
  double Ddn_Drn;
  if (taun - rn <= tol*rn0) {                    // no negative damage
    Ddn_Drn = 0;
  }  else {                                      // negative damage evolves
    rn = taun;                                   // update rn
    dn = 1 - rn0/rn*(1-An) - An*exp(Bn*(1 - rn/rn0));
    Ddn_Drn = rn0/(rn*rn)*(1-An) + An*Bn/rn0*exp(Bn*(1 - rn/rn0));
    dn = dn*(1-tol);                             // cap the damage variable
    Ddn_Drn = Ddn_Drn*(1-tol);
    if (dn > 1-tol) {
      dn = 1- tol; Ddn_Drn = 0;
    }
  }

  // stress update
  sig = (1-dp)*sigpos + (1-dn)*signeg;

  // TANGENT
  Dsigpos_Deps = Qpos*Cbar;
  Dsigneg_Deps = Qneg*Cbar;

  /*  
  opserr << "Qpos: " << Qpos;
  opserr << "Qneg: " << Qneg;
  opserr << "Cbar: " << Cbar;
  */
  
  static Vector s(6);
  s = Idp*signeg;                    // deviatoric stress

  //  opserr << "Idp: " << Idp;
  //  opserr << "signeg: " << signeg;

  // norm of deviatoric stress
  double nrms = sqrt( pow(s(0),2) + pow(s(1),2) + pow(s(2),2) +  
		      2*pow(s(3),2) + 2*pow(s(4),2) + 2*pow(s(5),2));

  static Vector n(6);
  if (nrms <= tol) 
    n.Zero(); 
  else {
    n = s; n/=nrms;
  }

  static Vector Dtaup_Dsigpos(6);
  static Vector Dtaun_Dsigneg(6);

  if (taup <= tol) {
    Dtaup_Dsigpos.Zero();   //  Dtaup_Dsigpos = zeros(6,1); 
  } else  {
    Ce.Solve(sigpos, tmp);  //  Dtaup_Dsigpos = (Ce\sigpos)/taup; end
    Dtaup_Dsigpos = tmp; Dtaup_Dsigpos/=taup;
  }

  if (taun <= tol) {
    Dtaun_Dsigneg.Zero();
  } else {
    double Dtaun_Dsigoct = pow(3,0.25) * k/2/sqrt(k*sigoct + tauoct);
    double Dtaun_Dtauoct = pow(3,0.25) /2/sqrt(k*sigoct + tauoct);
    static Vector Dsigoct_Dsigneg(6);
    Dsigoct_Dsigneg = Iv6; Dsigoct_Dsigneg/=3.;
    static Vector Dtauoct_Dsigneg(6);
    
    Dtauoct_Dsigneg = n; Dtauoct_Dsigneg/=sqrt(3.0);
    /*
    opserr << "Dtaun_Dsigoct: " << Dtaun_Dsigoct << endln;
    opserr << "Dsigoct_Dsigneg: " << Dsigoct_Dsigneg << endln;
    opserr << "Dtaun_Dtauoct: " << Dtaun_Dtauoct << endln;
    opserr << "Dtauoct_Dsigneg: " << Dtauoct_Dsigneg << endln;
    opserr << "n: " << n;
    */
    Dtaun_Dsigneg = Dtaun_Dsigoct * Dsigoct_Dsigneg + Dtaun_Dtauoct * Dtauoct_Dsigneg;
  }

  // Ddp_Deps = Ddp_Drp * Dsigpos_Deps' * Dtaup_Dsigpos;
  Ddp_Deps = Dsigpos_Deps ^ Dtaup_Dsigpos;
  Ddp_Deps *= Ddp_Drp; 

  Ddn_Deps = Dsigneg_Deps ^ Dtaun_Dsigneg;
  Ddn_Deps *= Ddn_Drn;
  
  C = (1-dp)*Dsigpos_Deps + (1-dn)*Dsigneg_Deps - 
    sigpos % Ddp_Deps - signeg % Ddn_Deps;  

  /*
  opserr << "dp: " << dp << endln;
  opserr << "Dsigpos_Deps: " << Dsigpos_Deps << endln;
  opserr << "dn: " << dn << endln;
  opserr << "Dsigneg_Deps: " << Dsigneg_Deps << endln;
  opserr << "sigpos: " << sigpos;
  opserr << "signeg: " << signeg;

  opserr << "Ddp_Deps: " << Ddp_Deps;
  opserr << "Ddn_Deps: " << Ddn_Deps;
  opserr << "Ddp_Drp: " << Ddp_Drp;
  opserr << "Ddn_Drn: " << Ddn_Drn;
  opserr << "Dtaup_Dsigpos: " << Dtaup_Dsigpos;
  opserr << "Dtaun_Dsigneg: " << Dtaun_Dsigneg;

  opserr << "taup: " << taup << endln;
  opserr << "taun: " << taun << endln;
  opserr << "rp: " << rp << endln;		
  opserr << "rn: " << rn << endln;		
  opserr << "rp0: " << rp0 << endln;		
  opserr << "rn0: " << rn0 << endln;		
  opserr << "dn: " << dn << endln;		
  opserr << "dp: " << dp << endln;		

  opserr << "TANGENT: " << C << endln;
  opserr << "Stress: " << sig << endln;
  */

  return 0;
}

int
PlasticDamageConcrete3d::setTrialStrainIncr (const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
PlasticDamageConcrete3d::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
PlasticDamageConcrete3d::getTangent (void)
{
  return C;
}

const Matrix&
PlasticDamageConcrete3d::getInitialTangent (void)
{
  return Ce;
}

const Vector&
PlasticDamageConcrete3d::getStress (void)
{
  return sig;
}

const Vector&
PlasticDamageConcrete3d::getStrain (void)
{
  return eps;
}

int
PlasticDamageConcrete3d::commitState (void)
{
  rpCommit = rp;
  rnCommit = rn;
  dpCommit = dp;
  dnCommit = dn;
  epsCommit = eps;
  sigCommit = sig;
  sigeCommit = sige;
  eps_pCommit = eps_p;
  sigePCommit = sige;
  return 0;
}

int
PlasticDamageConcrete3d::revertToLastCommit (void)
{
  C = Ccommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;
  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  sigeP = sigePCommit;
  return 0;
}

int
PlasticDamageConcrete3d::revertToStart (void)
{
  eps.Zero();
  sig.Zero();
  sige.Zero();
  eps_p.Zero();
  sigeP.Zero();
  Ce.Zero();

  return 0;
}

NDMaterial*
PlasticDamageConcrete3d::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
    PlasticDamageConcrete3d *theCopy =
      new PlasticDamageConcrete3d (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
    return theCopy;  
  } else {
    return 0;
  }
}

NDMaterial*
PlasticDamageConcrete3d::getCopy (void)
{
  PlasticDamageConcrete3d *theCopy =
    new PlasticDamageConcrete3d (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
  return theCopy;
}

const char*
PlasticDamageConcrete3d::getType (void) const
{
  return "ThreeDimensional";
}

int
PlasticDamageConcrete3d::getOrder (void) const
{
  return 6;
}


int 
PlasticDamageConcrete3d::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3d::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PlasticDamageConcrete3d::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3d::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
PlasticDamageConcrete3d::Print(OPS_Stream &s, int flag) {
  opserr << "PlasticDamageConcrete3d: " << this->getTag();
  opserr << "strain: " << eps;
  opserr << "strain: " << sig;
  opserr << "tangent: " << C;
}       
